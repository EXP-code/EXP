#include <math.h>
#include <pthread.h>

#include "expand.H"
#include <localmpi.H>

#include <gaussQ.H>
#include <massmodel.H>

#include <UserDiffuse.H>

static pthread_mutex_t randlock = PTHREAD_MUTEX_INITIALIZER;

UserDiffuse::UserDiffuse(const YAML::Node& conf) : ExternalForce(conf)
{
  id = "Two-body relaxation";

  name = "";			// Default component name
  nfreq = 10000000;		// An integer bigger than any reasonable
				// number of steps

  logr = true;			// Use logarithmic radial intervals
  rmin = 1.0e-3;		// Minimum radius for histogram
  rmax = 2.0;			// Maximum radius for histogram
  logL = 5.7;			// LogLambda
  pmass = 1.0e-6;		// Particle mass for relaxation
  clip = 100.0;			// Number of sigma for clipping

  numr = 200;			// Number of radial points in model
  numv = 64;			// Number of velocity points in 
				// coef grid
  nume = 128;			// Number of energy integration knots

  use_file = false;		// Use
  modfile = "SLGridSph.model";	// model file rather than histogram
  diverge = 0;			// Use power divergence for model
  diverge_rfac = 1.0;		// with this exponent

  check_ev = false;		// Deep debugging statement (don't use this)

  initialize();

  if (name.size()>0) {
				// Look for the fiducial component
    bool found = false;
    for (auto c : comp->components) {
      if ( !name.compare(c->name) ) {
	c0 = c;
	found = true;
      break;
      }
    }

    if (!found) {
      cerr << "Process " << myid << ": can't find desired component <"
	   << name << ">" << endl;
      MPI_Abort(MPI_COMM_WORLD, 35);
    }

  }
  else
    c0 = NULL;

  Gamma = 4.0*M_PI*pmass*logL;
  
  userinfo();
}

UserDiffuse::~UserDiffuse()
{
  delete model;
  delete jq;
}

void UserDiffuse::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  if (c0)
    cout << "** User routine 2-BODY DIFFUSION initialized for component: <" << name << ">";
  else
    cout << "** User routine 2-BODY DIFFUSION disabled: no component specified";
  
  cout << ", mass = " << pmass;
  cout << ", logL = " << logL;
  cout << ", seed = " << seed;
  if (use_file)
    cout << ", using model file <" << modfile << ">";
  cout << ", (Rmin, Rmax) = (" << rmin << ", " << rmax << ")";
  if (logr)
    cout << ", using log radial grid";
  if (check_ev)
    cout << ", direction debug stats are *ON*";
  cout << endl;

  print_divider();
}

void UserDiffuse::initialize()
{
  try {
    if (conf["name"])           name               = conf["name"].as<string>();
    if (conf["pmass"])          pmass              = conf["pmass"].as<double>();
    if (conf["clip"])           clip               = conf["clip"].as<double>();
    if (conf["logL"])           logL               = conf["logL"].as<double>();
    if (conf["seed"])           seed               = conf["seed"].as<int>();
    if (conf["nfreq"])          nfreq              = conf["nfreq"].as<int>();
    if (conf["rmin"])           rmin               = conf["rmin"].as<double>();
    if (conf["rmax"])           rmax               = conf["rmax"].as<double>();
    if (conf["logr"])           logr               = conf["logr"].as<bool>();
    if (conf["numr"])           numr               = conf["numr"].as<int>();
    if (conf["numv"])           numv               = conf["numv"].as<int>();
    
    if (conf["use_file"])       use_file           = conf["use_file"].as<bool>();
    if (conf["modfile"])        modfile            = conf["modfile"].as<string>();
    if (conf["diverge"])        diverge            = conf["diverge"].as<int>();
    if (conf["diverge_rfac"])   diverge_rfac       = conf["diverge_rfac"].as<double>();
    if (conf["check_ev"])       check_ev           = conf["check_ev"].as<bool>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in UserDiffuse: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }

  if (check_ev) {
    ev_mean_th.resize(nthrds);
    ev_disp_th.resize(nthrds);
    ev_numb_th.resize(nthrds);

    for (int i=0; i<nthrds; i++) {
      ev_mean_th[i].resize(6);
      ev_disp_th[i].resize(6);
    }

    if (myid==0) {
				// Open output stream for writing
      check_file = outdir + runtag + ".check_ev";
      std::ofstream out(check_file.c_str(), ios::out | ios::app);
      if (out.good()) {
	out.setf(ios::left);
	out << "# " << setw(14) << "Time"
	    << "| " << setw(14) << "Mean (1*V)"
	    << "| " << setw(14) << "Disp (1*V)"
	    << "| " << setw(14) << "Mean (2*V)"
	    << "| " << setw(14) << "Disp (2*V)"
	    << "| " << setw(14) << "Mean (3*V)"
	    << "| " << setw(14) << "Disp (3*V)"
	    << "| " << setw(14) << "Mean (1*2)"
	    << "| " << setw(14) << "Disp (1*2)"
	    << "| " << setw(14) << "Mean (1*3)"
	    << "| " << setw(14) << "Disp (1*3)"
	    << "| " << setw(14) << "Mean (2*3)"
	    << "| " << setw(14) << "Disp (2*3)"
	    << "| " << setw(14) << "Number"
	    << endl;
	out << "#-1-" << setw(12) << setfill('-') << "-"
	    << "|-2-" << setw(12) << setfill('-') << "-"
	    << "|-3-" << setw(12) << setfill('-') << "-"
	    << "|-4-" << setw(12) << setfill('-') << "-"
	    << "|-5-" << setw(12) << setfill('-') << "-"
	    << "|-6-" << setw(12) << setfill('-') << "-"
	    << "|-7-" << setw(12) << setfill('-') << "-"
	    << "|-8-" << setw(12) << setfill('-') << "-"
	    << "|-9-" << setw(12) << setfill('-') << "-"
	    << "|-10" << setw(12) << setfill('-') << "-"
	    << "|-11" << setw(12) << setfill('-') << "-"
	    << "|-12" << setw(12) << setfill('-') << "-"
	    << "|-13" << setw(12) << setfill('-') << "-"
	    << "|-14" << setw(12) << setfill('-') << "-"
	    << endl;
      }
    }
  }

}


void UserDiffuse::determine_acceleration_and_potential(void)
{
  if (!c0) return;

  if (multistep && mlevel>0) return;

#if HAVE_LIBCUDA==1		// Cuda compatibility
  getParticlesCuda(c0);
#endif

  if (!(this_step % nfreq)) {
    compute_model();
    compute_diffuse();
  }

  if (check_ev) {
    for (int i=0; i<nthrds; i++) {
      for (int k=0; k<6; k++) ev_mean_th[i][k] = ev_disp_th[i][k] = 0.0;
      ev_numb_th[i] = 0;
    }
  }

  exp_thread_fork(false);

  if (check_ev) {

    vector<double> ev_mean(6, 0.0), ev_disp(6, 0.0);
    vector<double> tt_mean(6, 0.0), tt_disp(6, 0.0);
    int ev_numb = 0;
    int tt_numb = 0;

    for (int i=0; i<nthrds; i++) {
      for (int k=0; k<6; k++) {
	ev_mean[k] += ev_mean_th[i][k];
	ev_disp[k] += ev_disp_th[i][k];
      }
      ev_numb += ev_numb_th[i];
    }

				// Get contribution from all processes
    MPI_Reduce(&ev_mean[0], &tt_mean[0], 6, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ev_disp[0], &tt_disp[0], 6, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ev_numb, &tt_numb, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (myid==0) {
				// Open output stream for writing
      std::ofstream out(check_file.c_str(), ios::out | ios::app);
      if (out.good()) {
	out.setf(ios::left);
	out << setw(16) << tnow;
	for (int k=0; k<6; k++) {
	  if (ev_numb>1)
	    out << setw(16) << ev_mean[k]/ev_numb
		<< setw(16) << sqrt((ev_disp[k]-ev_mean[k]*ev_mean[k]/ev_numb)/(ev_numb-1));
	  else
	    out << setw(16) << 0.0
		<< setw(16) << 0.0;
	}
	out << setw(16) << ev_numb << endl;
      } else
	cerr << "UserDiffuse: error opening <" << check_file << "> for append\n";
    }
  }

}


void * UserDiffuse::determine_acceleration_and_potential_thread(void * arg) 
{
  int nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  double rr, vv, vdotr;
  double dvpara1, dvpara2, dvperp2;
  vector<double> er(3), ev(3), e1(3), e2(3);

  for (int i=nbeg; i<nend; i++) {

				// Compute radius and velocity scalars
    rr = vv = vdotr = 0.0;
    for (int k=0; k<3; k++) {
      er[k] = cC->Pos(i, k) - c0->center[k];
      ev[k] = cC->Vel(i, k);
      rr += er[k]*er[k];
      vv += ev[k]*ev[k];
      vdotr += er[k]*ev[k];
    }

    rr = sqrt(rr+1.0e-12);
    vv = sqrt(vv+1.0e-12);

				// Compute unit vectors
    for (int k=0; k<3; k++) {
      er[k] /= rr;
      ev[k] /= vv;
    }

				// Subtract velocity direction
				// radial direction
    double norm = 0.0;
    for (int k=0; k<3; k++) {
      e1[k] = er[k] - vdotr*ev[k]/(vv*rr);
      norm += e1[k]*e1[k];
    }
    for (int k=0; k<3; k++) e1[k] /= sqrt(norm+1.0e-12);


				// Unit vec perp to ${\hat r}$ and ${\hat v}$
    /*
    e2[0] = er[1]*e1[2] - er[2]*e1[1];
    e2[1] = er[2]*e1[0] - er[0]*e1[2];
    e2[2] = er[0]*e1[1] - er[1]*e1[0];
    */
				// Unit vec perp to ${\hat e1}$ and ${\hat v}$
    e2[0] = ev[1]*e1[2] - ev[2]*e1[1];
    e2[1] = ev[2]*e1[0] - ev[0]*e1[2];
    e2[2] = ev[0]*e1[1] - ev[1]*e1[0];

				// Get diffusion coefficients
    get_coefs(rr, vv, dvpara1, dvpara2, dvperp2);
    

				// Normal random variates
    pthread_mutex_lock(&randlock);
    double r1 = nrand(random_gen);
    double r2 = nrand(random_gen);
    double r3 = urand(random_gen);
    pthread_mutex_unlock(&randlock);
				// Clip

    if (fabs(r1)>clip) r1 = copysign(clip, r1);
    if (fabs(r2)>clip) r2 = copysign(clip, r2);

				// Velocity changes
    double deltaVpara = dvpara1*dtime + r1*sqrt(dvpara2*dtime);
    double deltaVperp = r2*sqrt(dvperp2*dtime);

    double angle = 2.0*M_PI*r3;	// Angle in tangent velocity plane
    double cosA = cos(angle);
    double sinA = sin(angle);
				// Update particle velocity
    for (int k=0; k<3; k++)
      cC->AddVel(i, k, 
	deltaVpara*ev[k]      +	// Parallel direction
	deltaVperp*e1[k]*cosA +	// Perpendicular direction #1
	deltaVperp*e2[k]*sinA 	// Perpendicular direction #2
		      );
				// Debug directions
    if (check_ev) {

      double tmp;
				// 1*V
      tmp = 0.0;
      for (int k=0; k<3; k++) tmp += cC->Vel(i, k) * ev[k];
      tmp /= vv;
      ev_mean_th[id][0] += tmp;
      ev_disp_th[id][0] += tmp*tmp;
				// 2*V
      tmp = 0.0;
      for (int k=0; k<3; k++) tmp += cC->Vel(i, k) * e1[k];
      tmp /= vv;
      ev_mean_th[id][1] += tmp;
      ev_disp_th[id][1] += tmp*tmp;
				// 3*V
      tmp = 0.0;
      for (int k=0; k<3; k++) tmp += cC->Vel(i, k) * e2[k];
      tmp /= vv;
      ev_mean_th[id][2] += tmp;
      ev_disp_th[id][2] += tmp*tmp;
				// 1*2
      tmp = 0.0;
      for (int k=0; k<3; k++) tmp += ev[k] * e1[k];
      ev_mean_th[id][3] += tmp;
      ev_disp_th[id][3] += tmp*tmp;
				// 1*3
      tmp = 0.0;
      for (int k=0; k<3; k++) tmp += ev[k] * e2[k];
      ev_mean_th[id][4] += tmp;
      ev_disp_th[id][4] += tmp*tmp;
				// 2*3
      tmp = 0.0;
      for (int k=0; k<3; k++) tmp += e1[k] * e2[k];
      ev_mean_th[id][5] += tmp;
      ev_disp_th[id][5] += tmp*tmp;
				// Increment counter
      ev_numb_th[id]++;
    }

  }

  return (NULL);
}


void UserDiffuse::compute_model()
{
  if (use_file) {
    model = new SphericalModelTable(modfile, diverge, diverge_rfac);
    model->setup_df(800);
    return;
  }


  // ===================================================================
  // Accumulate density histogram
  // ===================================================================

  int indx;
  double r, dr;

  if (logr)
    dr = (log(rmax) - log(rmin))/(numr - 1);
  else
    dr = (rmax - rmin)/(numr - 1);

  
  vector<double> histo(numr, 0.0);


  for (unsigned n=0; n<c0->Number(); n++)
    {
      r = 0.0;
      for (int i=0; i<3; i++) 
	r += (c0->Pos(n, i)-c0->center[i])*(c0->Pos(n, i)-c0->center[i]);
      r = sqrt(r);

      if (r>rmin && r<rmax) {
	if (logr)
	  indx = (int)( (log(r)-log(rmin))/dr );
	else
	  indx = (int)( (r-rmin)/dr );

	histo[indx] += c0->Mass(n);
      }
      
    }


  // ===================================================================
  // Compute model
  // ===================================================================

  const double SMALL = 1.0e-12;
  vector<double> R(numr), D(numr), M(numr), P(numr), P1(numr);
  vector<double> delr(numr), r1(numr), r2(numr), ld(numr);

  for (int i=0; i<numr; i++) {

    if (logr) {
      r1[i] = rmin*exp(dr*i);
      r2[i] = rmin*exp(dr*(i+1));
      R[i] = rmin*exp(dr*(0.5+i));
      delr[i] = r2[i] - r1[i];
    } else {
      r1[i] = rmin + dr*i;
      r2[i] = rmin + dr*(i+1);
      R[i] = rmin + dr*(0.5+i);
      delr[i] = dr;
    }

    D[i] = (double)histo[i] / 
      (4.0*M_PI/3.0 * (pow(r2[i], 3.0) - pow(r1[i], 3.0)));
    ld[i] = log(D[i]+SMALL);
  }


  for (int i=0; i<numr; i++) {

    if (i==0) M[i] = 0.0;
    else M[i] = M[i-1] + 2.0*M_PI*delr[i]*(D[i]*r2[i]*r2[i] + D[i-1]*r1[i]*r1[i]);
    
    if (i==0) P1[i] = 4.0*M_PI*delr[i]*D[i]*r2[i];
    else P1[i] = P1[i-1] + 2.0*M_PI*delr[i]*(D[i]*r2[i] + D[i-1]*r1[i]);
  }
  
  for (int i=numr-1; i>0; i--) 
    P[i] = -M[i]/R[i] - (P1[numr-1] - 0.5*(P1[i] + P1[i-1]));
  P[0] = -M[0]/R[0]  - (P1[numr-1] - 0.5*P1[0]);


  model = new SphericalModelTable(numr, &R[0]-1, &D[0]-1, &M[0]-1, &P[0]-1,
				  0, 0.0, 0, "model from histogram");
  model->setup_df(800);
  
				// Debug: print out model
  if (myid==0) {
    model->print_model("test_model.dat");
    model->print_df("test_df.dat");
  }

}


void UserDiffuse::compute_diffuse()
{
  // ===================================================================
  // Create grids (if necessary)
  // =================================================================== 
  
  if (R.size() == 0) {

    R = vector<double>(numr);
    Vmax = vector<double>(numr);
    dVpara1 = vector<dvector>(numr);
    dVpara2 = vector<dvector>(numr);
    dVperp2 = vector<dvector>(numr);

    for (int i=0; i<numr; i++) {
      dVpara1[i] = vector<double>(numv);
      dVpara2[i] = vector<double>(numv);
      dVperp2[i] = vector<double>(numv);
    }

    jq = new LegeQuad(nume);
  }

  Rmin = model->get_min_radius();
  Rmax = model->get_max_radius();
  Emin = model->get_pot(Rmin);
  Emax = model->get_pot(Rmax);

  if (logr) {
    delR = (log(Rmax) - log(Rmin))/numr;
    for (int i=0; i<numr; i++) R[i] = Rmin*exp(delR*i);
  }
  else {
    delR = (Rmax - Rmin)/numr;
    for (int i=0; i<numr; i++) R[i] = Rmin + delR*i;
  }


  // ===================================================================
  // Dump grid to file
  // =================================================================== 

  std::ofstream out;
  if (myid==0) {
    string outfile = outdir + "diffusion.grid" + runtag;
    out.open(outfile);
    if (out.good()) {
      out.setf(ios::left);

      out << setw(15) << "# Radius" 
	  << setw(15) << "Velocity"
	  << setw(15) << "D(para)"
	  << setw(15) << "D(para^2)"
	  << setw(15) << "D(perp^2)"
	  << endl;

      char c = out.fill('-');
      out << setw(15) << "#-1" 
	  << setw(15) << "|-2"
	  << setw(15) << "|-3"
	  << setw(15) << "|-4"
	  << setw(15) << "|-5"
	  << endl;
      out.fill(c);
    }
  }

  // ===================================================================
  // Compute diffusion coefficients
  // =================================================================== 
  
  double F0, F1, F2;
  double E, EE, phi, dV, V, df, f1, vrat;

				// Expressions from Cohn 1979
  delV = 1.0/(numv-1);

  for (int i=0; i<numr; i++) {
    
    phi = model->get_pot(R[i]);
    Vmax[i] = sqrt(2.0*(Emax - phi));
    dV = Vmax[i]/numv;

    for (int j=0; j<numv; j++) {

      V = dV*(j+0.5);
      EE = 0.5*V*V + phi;

      F0 = F1 = F2 = 0.0;

      for (int k=0; k<jq->get_n(); k++) {

	F0 += 2.0*model->distf(EE + (Emax - EE)*jq->knot(k), 0.5) * 
	  jq->weight(k);
      
	E = phi + (EE - phi)*jq->knot(k);
	df = 2.0*model->distf(E, 0.5);
	vrat = (E - phi)/(EE - phi);
	f1 = df*sqrt(vrat);
      
	F1 += f1 * jq->weight(k);
	F2 += f1*vrat * jq->weight(k);
      }
      
      F0 *= 4.0*M_PI*Gamma*(Emax - EE);
      F1 *= 4.0*M_PI*Gamma*(EE - phi);
      F2 *= 4.0*M_PI*Gamma*(EE - phi);
	
      dVpara1[i][j] = -2.0*F1/V;
      dVpara2[i][j] = 2.0*(F0 + F2)/3.0;
      dVperp2[i][j] = 2.0*(2.0*F0 + 3.0*F1 - F2)/3.0;
    
      if (myid==0 and out.good())
	out << setw(15) << R[i] 
	    << setw(15) << V
	    << setw(15) << dVpara1[i][j]
	    << setw(15) << dVpara2[i][j]
	    << setw(15) << dVperp2[i][j]
	    << endl;

    } // End V loop

    if (myid==0 and out.good()) out << endl;

  } // End R loop

  if (myid==0) {
    try {
      out.close();
    }
    catch (const ofstream::failure& e) {
      std::cout << "UserDiffuse: exception closing file <"
		<< outdir + "diffusion.grid" + runtag
		<< ": " << e.what() << std::endl;
    }
  }

}


void UserDiffuse::get_coefs(double r, double v,
			    double& dvpara1, double& dvpara2, double& dvperp2)
{
  
  // ===================================================================
  // Linear interpolation coefficients
  // =================================================================== 

  int indxR, indxV;
  double a0, a1, b0, b1;
  
				// Enforce upper and lower bounds
#ifdef DEBUG
  if (r>1.1*Rmax)
    cout << "Process " << myid << ": R=" << r
	 << " V=" << v
	 << " Rmin=" << Rmin << " Rmax=" << Rmax
	 << ": radius out of bounds" << endl;
#endif      
  r = max<double>(r, Rmin);
  r = min<double>(r, Rmax);

  if (logr)
    indxR = (int)( (log(r) - log(Rmin))/delR );
  else
    indxR = (int)( (r - Rmin)/delR );

  indxR = min<int>(indxR, numr-2);

  if (logr) {
    a0 = (log(R[indxR+1]) - log(r))/delR;
    a1 = (log(r) - log(R[indxR  ]))/delR;
  } else {
    a0 = (R[indxR+1] - r)/delR;
    a1 = (r - R[indxR  ])/delR;
  }

				// Get Vmax
  double vmax = a0*Vmax[indxR] + a1*Vmax[indxR+1];

				// Enforce upper and lower bounds
  double vrel = v/vmax;
#ifdef DEBUG
  if (vrel>1.4)
    cout << "Process " << myid << ": R=" << r
	 << " V=" << v << " Vmax=" << vmax
	 << ": velocity out of bounds" << endl;
#endif      
  vrel = max<double>(vrel, 0.0);
  vrel = min<double>(vrel, 1.0);

  indxV = (int)( vrel/delV-0.5 );
  indxV = min<int>(indxV, numv-2);

  b0 = (delV*(1.5+indxV) - vrel)/delV;
  b1 = (vrel - delV*(0.5+indxV))/delV;


				// Do the interpolation on all quatities
  dvpara1 = 
    b0*(a0*dVpara1[indxR][indxV]   + a1*dVpara1[indxR+1][indxV]  ) +
    b1*(a0*dVpara1[indxR][indxV+1] + a1*dVpara1[indxR+1][indxV+1]) ;

  dvpara2 = 
    b0*(a0*dVpara2[indxR][indxV]   + a1*dVpara2[indxR+1][indxV]  ) +
    b1*(a0*dVpara2[indxR][indxV+1] + a1*dVpara2[indxR+1][indxV+1]) ;

  dvperp2 = 
    b0*(a0*dVperp2[indxR][indxV]   + a1*dVperp2[indxR+1][indxV]  ) +
    b1*(a0*dVperp2[indxR][indxV+1] + a1*dVperp2[indxR+1][indxV+1]) ;

				// Interpolation to unphysical values?
  if (dvpara2 < 0.0 || dvperp2 < 0.0) 
    {
#ifdef DEBUG
      cout << "Process " << myid << ": R=" << r
	   << " V=" << v << " Vmax=" << vmax
	   << " dvpara1, dvpara2, dvperp2="
	   << dvpara1 << ", " << dvpara2 << ", " << dvperp2 << endl;
#endif      
      dvpara1 = dvpara2 = dvperp2 = 0.0;
    }

}


extern "C" {
  ExternalForce *makerDiffuse(const YAML::Node& conf)
  {
    return new UserDiffuse(conf);
  }
}

class proxydiffuse { 
public:
  proxydiffuse()
  {
    factory["userdiffuse"] = makerDiffuse;
  }
};

proxydiffuse p;
