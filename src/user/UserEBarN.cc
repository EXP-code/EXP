#include <sstream>
#include <cstdlib>
#include <cmath>

#include "expand.H"
#include <localmpi.H>
#include <gaussQ.H>
#include <UserEBarN.H>
#include <Timer.H>

static Timer timer_tot, timer_thrd;
static bool timing = false;

UserEBarN::UserEBarN(const YAML::Node& conf) : ExternalForce(conf)
{
  id = "RotatingBarWithMonopole";

  length = 1.0;			// Bar length
  bratio = 0.5;			// Ratio of b to a
  cratio = 0.1;			// Ratio of c to b
  amplitude = 0.3;		// Bar quadrupole amplitude
  angmomfac = 1.0;		// Artifically change the total bar ang mom
  barmass = 1.0;		// Total bar mass
  monoamp = 1.0;		// Monopole amplitude
  Ton = -20.0;			// Turn on start time
  Toff = 200.0;			// Turn off start time
  TmonoOn = -20.0;		// Turn on start time for monopole
  TmonoOff = 200.0;		// Turn off start time monopole
  DeltaT = 1.0;			// Turn on duration
  DeltaMonoT = 1.0;		// Turn on duration for monopole
  DOmega = 0.0;			// Change in pattern speed
  dtom = -1.0;			// Width of forced bar slow down
  T0 = 10.0;			// Center of pattern speed change
  Fcorot  = 1.0;		// Corotation factor
  fixed = false;		// Constant pattern speed
  alpha = 5.0;			// Variable sharpness bar potential
  monopole = true;		// Use the monopole part of the potential
  monopole_follow = true;	// Follow monopole center
  monopole_onoff = false;	// To apply turn-on and turn-off to monopole
  monopole_frac = 1.0;		// Fraction of monopole to turn off
  quadrupole_frac = 1.0;	// Fraction of quadrupole to turn off

  bartype = powerlaw;		// Use the the powerlaw density profile
  modelp = 0.0;			// Either the power law exponent or the exponential scale length
  rsmin  = 0.0001;		// The minimum radius for the quadrupole table in major axis units
  rsmax  = 2;	                // The maximum radius for the quadrupole table in major axis units
  numt   = 1000;                // Number of table entries for quadrupole table
  N      = 100;			// Number of integration knots

				// Output file name
  filename = outdir + "BarRot." + runtag;

  firstime = true;
  omega0 = -1.0;

  ctr_name = "";		// Default component for com
  angm_name = "";		// Default component for angular momentum
  
  initialize();

  if (ctr_name.size()>0) {
				// Look for the fiducial component for
				// centering
    bool found = false;
    for (auto c : comp->components) {
      if ( !ctr_name.compare(c->name) ) {
	c0 = c;
	found = true;
      break;
      }
    }

    if (!found) {
      cerr << "Process " << myid << ": can't find desired component <"
	   << ctr_name << ">" << endl;
      MPI_Abort(MPI_COMM_WORLD, 35);
    }

  }
  else
    c0 = NULL;

  if (angm_name.size()>0) {
				// Look for the fiducial component
				// for angular momentum
    bool found = false;
    for (auto c : comp->components) {
      if ( !angm_name.compare(c->name) ) {
	c1 = c;
	found = true;
      break;
      }
    }

    if (!found) {
      cerr << "Process " << myid << ": can't find desired component <"
	   << angm_name << ">" << endl;
      MPI_Abort(MPI_COMM_WORLD, 35);
    }

  }
  else
    c1 = NULL;

  // Zero monopole variables
  teval = vector<double>(multistep+1, tnow);
  for (int k=0; k<3; k++) bps[k] = vel[k] = acc[k] = 0.0;

  // Assign working vectors for each thread
  tacc.resize(nthrds);
  for (int n=0; n<nthrds; n++) tacc[n].resize(3);

  userinfo();

  // Only turn on bar timing for extreme debugging levels
  if (VERBOSE>49) timing = true;

  gq = std::make_shared<LegeQuad>(N);
}

UserEBarN::~UserEBarN()
{
  // Nothing
}

void UserEBarN::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine ROTATING BAR with MONOPOLE initialized, " ;
  if (fixed) {
    cout << "prescribed pattern speed with"
	 << " domega=" <<  DOmega << " and t0=" << T0;
    if (dtom>0) cout << ", dT_om=" << dtom;
    cout << ", ";
  }
  else
    cout << "initial corotation fraction, ";
  cout << " amplitude=" << amplitude << ", ";

  if (fabs(angmomfac-1.0)>1.0e-10)
    cout << " ang mom factor=" << angmomfac << ", ";
    
  if (omega0<0.0)
    cout << "initial pattern speed to be computed, ";
  else
    cout << "initial pattern speed " << omega0 << ", ";

  cout << "quadrupole fraction=" << quadrupole_frac
       << ", Ton=" << Ton << ", Toff=" << Toff << ", DeltaT=" << DeltaT 
       << ", ";
  if (monopole) {
    if (monopole_onoff)
      cout << "using monopole (amplitude=" << monoamp 
	   << ") with turn-on/off with fraction=" 
	   << monopole_frac << ", TmonoOn=" << TmonoOn
	   << ", TmonoOff=" << TmonoOff << ", DeltaMonoT=" << DeltaMonoT
	   << ", ";
    else
      cout << "using monopole (amplitude=" << monoamp << "), ";
    if (monopole_follow)
      cout << "self-consistent monopole centering, ";
    else
      cout << "monopole center fixed, ";
  }
  else
    cout << "without monopole, ";

  cout << "bar softness (alpha)=" << alpha << ", ";

  switch (bartype) {
  case powerlaw:
    cout << "powerlaw density profile with exponent=" << modelp << ", ";
    break;
  case ferrers:
    cout << "Ferrers density profile with exponent=" << modelp << ", ";
    break;
  case expon:
    cout << "Exponential density profile with scale length=" << modelp << ", ";
    break;
  }

  if (c0) 
    cout << "center on component <" << ctr_name << ">, ";
  else
    cout << "using inertial center, ";
  if (c1) 
    cout << "angular momentum from <" << angm_name << ">" << endl;
  else
    cout << "no initial angular momentum (besides bar)" << endl;

  print_divider();
}

void UserEBarN::initialize()
{
  try {
    if (conf["ctrname"])        ctr_name           = conf["ctrname"].as<string>();
    if (conf["angmname"])       angm_name          = conf["angmname"].as<string>();
    if (conf["length"])         length             = conf["length"].as<double>();
    if (conf["bratio"])         bratio             = conf["bratio"].as<double>();
    if (conf["cratio"])         cratio             = conf["cratio"].as<double>();
    if (conf["amp"])            amplitude          = conf["amp"].as<double>();
    if (conf["angmomfac"])      angmomfac          = conf["angmomfac"].as<double>();
    if (conf["barmass"])        barmass            = conf["barmass"].as<double>();
    if (conf["Ton"])            Ton                = conf["Ton"].as<double>();
    if (conf["Toff"])           Toff               = conf["Toff"].as<double>();
    if (conf["TmonoOn"])        TmonoOn            = conf["TmonoOn"].as<double>();
    if (conf["TmonoOff"])       TmonoOff           = conf["TmonoOff"].as<double>();
    if (conf["DeltaT"])         DeltaT             = conf["DeltaT"].as<double>();
    if (conf["DeltaMonoT"])     DeltaMonoT         = conf["DeltaMonoT"].as<double>();
    if (conf["DOmega"])         DOmega             = conf["DOmega"].as<double>();
    if (conf["dtom"])           dtom               = conf["dtom"].as<double>();
    if (conf["T0"])             T0                 = conf["T0"].as<double>();
    if (conf["Fcorot"])         Fcorot             = conf["Fcorot"].as<double>();
    if (conf["omega"])          omega0             = conf["omega"].as<double>();
    if (conf["fixed"])          fixed              = conf["fixed"].as<bool>();
    if (conf["self"])           fixed              = conf["self"].as<bool>();
    if (conf["alpha"])          alpha              = conf["alpha"].as<double>();
    if (conf["monopole"])       monopole           = conf["monopole"].as<bool>();
    if (conf["follow"])         monopole_follow    = conf["follow"].as<bool>();
    if (conf["onoff"])          monopole_onoff     = conf["onoff"].as<bool>();
    if (conf["monofrac"])       monopole_frac      = conf["monofrac"].as<double>();
    if (conf["quadfrac"])       quadrupole_frac    = conf["quadfrac"].as<double>();
    if (conf["monoamp"])        monoamp            = conf["monoamp"].as<double>();
    if (conf["filename"])       filename           = conf["filename"].as<string>();
    
    if (conf["modelp"])         modelp             = conf["modelp"].as<double>();
    if (conf["rmin"])           rsmin              = conf["rmin"].as<double>();
    if (conf["rmax"])           rsmax              = conf["rmax"].as<double>();
    if (conf["numt"])           numt               = conf["numt"].as<int>();
    
    if (conf["bartype"]) {
      switch(conf["bartype"].as<int>()) {
      case powerlaw:
	bartype = powerlaw;
	break;
      case ferrers:
	bartype = ferrers;
	break;
      case expon:
	bartype = expon;
	break;
      default:
	if (myid==0)
	  std::cerr << "UnderEBarN: no such bar profile="
		    << conf["bartype"].as<int>() << std::endl;
	MPI_Abort(MPI_COMM_WORLD, 36);
      }
    }
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in UserEBarN: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }
}


vector<double> a(3);

double find_fct(double u, vector<double>& z)
{
  double ans = -1.0;
  for (int k=0; k<3; k++) ans += z[k]*z[k]/(a[k]*a[k] + u);
  return ans;
}

// Solve for u or lambda
double solve(vector<double> x, double m2)
{
  ZBrent< vector<double> > zbrent;

  double r2 = 0.0;
  for (int k=0; k<3; k++) {
    x[k] /= m2;
    r2 += x[k]*x[k];
  }
  double umax = r2-a[2]*a[2];

  if (umax<0.0) {
    cerr << "out of bounds!!!" << endl;
    return -1.0;
  }
  
  double ans = 0.0;
  ZBrentReturn ret = zbrent.find(find_fct, x, 0.0, umax, 1.0e-10, ans);
  
  switch (ret) {
  case Bracket:
    cout << "Root is not bracketed" << endl;
    break;

  case Iteration:
    cout << "Number of iterations exceeded" << endl;
    break;

  case Good:
    break;
  }

  return ans;
}

double UserEBarN::Density(vector<double> x)
{
  double m2=0.0, rho=0.0;
  for (int k=0; k<3; k++) m2 += x[k]*x[k]/(a[k]*a[k]);

  if (m2>1.0) return 0.0;

  switch (bartype){
  case powerlaw:
    rho = rho0*pow(m2, modelp);
    break;
  case ferrers:
    rho = rho0*pow(1.0 - m2, modelp);
    break;
  case expon:
    rho =  rho0*exp(-a[0]*sqrt(m2)/modelp)/sqrt(m2);
    break;
  }

  return rho/(M_PI*a[0]*a[1]*a[2]);
}

double UserEBarN::RhoBar(double r) 
{
  vector<double> z(3);
  const int nphi = 200, ntheta = 100;
  const double dphi = 0.5*M_PI/nphi;

  double cosx, sinx, phi, ans=0.0;

  if (gt==0) gt = std::make_shared<LegeQuad>(ntheta);

  for (int i=0; i<nphi; i++) {
    phi = dphi*i;
    for (int j=0; j<ntheta; j++) {
      cosx = gt->knot(j);
      sinx = sqrt(1.0 - cosx*cosx);

      z[0] = r*sinx*cos(phi);
      z[1] = r*sinx*sin(phi);
      z[2] = r*cosx;

      ans += gt->weight(j)*dphi * Density(z);
    }
  }

  return ans*8.0/(4.0*M_PI);
}


double UserEBarN::Potential(vector<double> x)
{
  double mshell=0.0, ans=0.0, u, d, t, denom, m2;
  double tmin, tmax=0.5*M_PI;

  double ellip = 0.0;
  for (int k=0; k<3; k++) ellip += x[k]*x[k]/(a[k]*a[k]);
  ellip -= 1.0;

  if (ellip<0.0) {		// Inside
    tmin = 0.0;
  } else {			// Outside
    tmin = atan(solve(x, 1.0));
  }

  for (int i=0; i<N; i++) {
    t = tmin + (tmax - tmin)*gq->knot(i);
    u = tan(t);
    d = cos(t);
    d = 1.0/(d*d);
    
    m2 = 0.0;
    denom = 1.0;
    for (int k=0; k<3; k++) {
      m2 += x[k]*x[k]/(a[k]*a[k] + u);
      denom *= a[k]*a[k]+u;
    }

    switch (bartype) {
    case powerlaw:
      mshell = rho0/(modelp+1.0)*(1.0 - pow(m2, modelp+1.0));
      break;
    case ferrers:
      mshell = rho0/(modelp+1.0)*(1.0 - pow(1.0-m2, modelp+1.0));
      break;
    case expon:
      mshell = 2.0*rho0*modelp/a[0]*(exp(-a[0]/modelp) - exp(-a[0]*sqrt(m2)/modelp));
      break;
    }

    ans += d*gq->weight(i) * mshell/sqrt(denom);
  }

  return ans*(tmax - tmin);
}


double UserEBarN::U22(double r) 
{
  vector<double> z(3);
  const int nphi = 200, ntheta = 100;
  const double dphi = 2.0*M_PI/nphi;

  double cosx, sinx, phi, ans=0.0;

  if (gt==0) gt = std::make_shared<LegeQuad>(ntheta);

  for (int i=0; i<nphi; i++) {
    phi = dphi*i;
    for (int j=0; j<ntheta; j++) {
      cosx = 2.0*(gt->knot(j) - 0.5);
      sinx = sqrt(1.0 - cosx*cosx);

      z[0] = r*sinx*cos(phi);
      z[1] = r*sinx*sin(phi);
      z[2] = r*cosx;

      ans += 2.0*gt->weight(j)*dphi * Potential(z) * sinx*sinx*cos(2.0*phi);
    }
  }

  return ans*numfac;
}


void UserEBarN::Inertia(vector<double>& I)
{
  I = vector<double>(3, 0.0);
  vector<double> z(3);
  double fac;

  for (int i=0; i<N; i++) {
    z[0] = a[0]*gq->knot(i);

    for (int j=0; j<N; j++) {
      z[1] = a[1]*gq->knot(j);

      for (int k=0; k<N; k++) {
	z[2] = a[2]*gq->knot(k);

	fac = gq->weight(i)*gq->weight(j)*gq->weight(k) * Density(z); 
	
	I[0] += fac * (z[1]*z[1] + z[2]*z[2]);
	I[1] += fac * (z[0]*z[0] + z[2]*z[2]);
	I[2] += fac * (z[0]*z[0] + z[1]*z[1]);
      }
    }
  }

  for (int k=0; k<3; k++) I[k] *= 8.0*a[0]*a[1]*a[2];
}

bool UserEBarN::quadpot(double r, double& p, double& f, double& fr, double& M)
{
  if (r<rr.front()) {
    p = pp.front();
    f = fr = 0.0;
    M = mm.front();
    return false;
  }

  if (r>=rr.back()) {
    p = pp.back()*rr.back()/r;
    f = uu.back()*pow(rr.back()/r, 3.0);
    fr = -3.0*f/r;
    M = mm.back();
    return false;
  }

  double lr = log(r);
  int indx = max<int>(0, min<int>((int)floor((lr - lrmin)/ldr), numt-2));
  int indx2 = max<int>(1, indx);
  double dlr = (lr - (lrmin + ldr*indx2))/ldr;

  double a = (lrmin + ldr*(indx+1) - lr)/ldr;
  double b = (lr - lrmin - ldr*indx)/ldr;

  p = a*pp[indx] + b*pp[indx+1];
  f = a*uu[indx] + b*uu[indx+1];
  fr = (uu[indx2+1]*(dlr+0.5) - 2.0*uu[indx]*dlr + uu[indx2-1]*(dlr-0.5))/(r*ldr);
  M = a*mm[indx] + b*mm[indx+1];

  return true;
}


void UserEBarN::determine_acceleration_and_potential(void)
{
  if (timing) timer_tot.start();
				// Write to bar state file, if true
  bool update = false;

  if (c1) c1->get_angmom();	// Tell component to compute angular momentum
  // cout << "Process " << myid << ": Lz=" << c1->angmom[2] << endl; // debug  

  if (firstime) {
    
    if (omega0 < 0.0) {

      double R=length*Fcorot;
      double phi, theta=0.5*M_PI;
      double dens0, potl0, dens, potl, potr, pott, potp;
      double avg=0.0;
      
      for (int n=0; n<8; n++) {
	phi = 2.0*M_PI/8.0 * n;

	for (auto c : comp->components) {
	
	  if (c->force->geometry == PotAccel::sphere || 
	      c->force->geometry == PotAccel::cylinder) {
	  
	    ((Basis*)c->force)->
	      determine_fields_at_point_sph(R, theta, phi,
					    &dens0, &potl0,
					    &dens, &potl, &potr, &pott, &potp);
	  
	    avg += potr/8.0;
	  }
	}
      }

      omega0 = sqrt(avg/R);
    }

    if (dtom>0.0)
      omega = omega0*(1.0 + DOmega*0.5*(1.0 + erf( (tnow - T0)/dtom )));
    else
      omega = omega0*(1.0 + DOmega*(tnow - T0*0.5));

    a[0] = length;
    a[1] = bratio*a[0];
    a[2] = cratio*a[1];

    // *M_PI*a[0]*a[1]*a[2] sucked into the leading factor for the 
    // gravitational potential
    //
    switch (bartype) {
    case expon:
      rho0 = a[0]*a[0]*barmass/(4.0*modelp*modelp) / 
	( 1.0 - (1.0 + a[0]/modelp)*exp(-a[0]/modelp) ); 
      break;
    case powerlaw:
      rho0 = (2.0*modelp + 3.0)*barmass/4.0; 
      break;
    case ferrers:
      rho0 = 2.0*exp(lgamma(2.5+modelp) - lgamma(1.5) - lgamma(1.0+modelp))*barmass/4.0;
      break;
    }

    lrmin = log(rsmin);
    lrmax = log(rsmax);
    ldr = (lrmax - lrmin)/(numt-1);
    double r, msum=0.0, psum=0.0;

    vector<double> rr1(numt, 0.0);
    vector<double> rb1(numt, 0.0), rb(numt, 0.0);
    vector<double> ub1(numt, 0.0);

    rr = vector<double>(numt, 0.0);
    mm = vector<double>(numt, 0.0);
    pp = vector<double>(numt, 0.0);
    uu = vector<double>(numt, 0.0);

				// Evaluate the table in parallel
				// 
    int jbeg = (myid+0)*numt/numprocs;
    int jend = (myid+1)*numt/numprocs;

    for (int j=jbeg; j<jend; j++) {
      rr1[j] = r = exp(lrmin + ldr*j);
      rb1[j] = RhoBar(r);
      ub1[j] = U22(r);
    }

				// Combine evaluations from all processes
				// 
    MPI_Allreduce(&rr1[0], &rr[0], numt, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&rb1[0], &rb[0], numt, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&ub1[0], &uu[0], numt, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				// Compute the mass and potential
				// 
    for (int j=1; j<numt; j++) {
      msum += 4.0*M_PI/2.0 * ldr * 
	(rb[j-1]*pow(rr[j-1], 3) + rb[j]*pow(rr[j], 3));

      psum += 4.0*M_PI/2.0 * ldr * 
	(rb[j-1]*pow(rr[j-1], 2) + rb[j]*pow(rr[j], 2));

      mm[j] = msum;
      pp[j] = psum;
    }

    for (int j=0; j<numt; j++)
      pp[j] = -mm[j]/rr[j] - (pp.back() - pp[j]);

    if (myid==0) {

      if (VERBOSE>3) {
	string name = outdir + "UserEBarN." + runtag + ".debug";
	ofstream out(name.c_str());

	out << "# a1, a2, a3=" << a[0] << ", " << a[1] << ", " << a[2] << endl;
	out << "# Expon=" << expon << endl;
	out << "# Param=" << modelp << endl;
	out << "# Mass="  << barmass << endl;
	out << "# Ampl="  << amplitude << endl;
	out << "# I_3="   << 0.2*barmass*(a[0]*a[0] + a[1]*a[1]) << endl;
	out << "#" << endl;
	for (int j=0; j<numt; j++)
	  out << setw(18) << rr[j]
	      << setw(18) << mm[j] 
	      << setw(18) << pp[j] 
	      << setw(18) << uu[j] << endl;
      }

      cout << "====================================================\n";
      cout << "Computed quadrupole fit for an ellipsoid using\n";
      switch (bartype) {
      case powerlaw:
	cout << "powerlaw density profile with exponent=" << modelp << endl;
	break;
      case ferrers:
	cout << "Ferrers density profile with exponent=" << modelp << endl;
	break;
      case expon:
	cout << "Exponential density profile with scale length=" << modelp << endl;
	break;
      }
      cout << "Rho0=" << rho0 << endl;
      cout << "Mass=" << barmass << endl;
      cout << "Ampl=" << amplitude << endl;
      cout << "I_3=" << 0.2*barmass*(a[0]*a[0] + a[1]*a[1]) << endl;
      cout << "Omega(0)=" << omega0 << endl;
      cout << "====================================================\n";
      
      name = outdir + filename;
      name += ".barstat";

      if (!restart) {
	ofstream out(name.c_str(), ios::out | ios::app);

	out << setw(15) << "# Time"
	    << setw(15) << "Phi"
	    << setw(15) << "Omega"
	    << setw(15) << "L_z(Bar)"
	    << setw(15) << "L_z(PS)"
	    << setw(15) << "Amp"
	    << setw(15) << "x"
	    << setw(15) << "y"
	    << setw(15) << "z"
	    << setw(15) << "u"
	    << setw(15) << "v"
	    << setw(15) << "w"
	    << setw(15) << "ax"
	    << setw(15) << "ay"
	    << setw(15) << "az"
	    << endl;
      }
    }

    vector<double> I;
    Inertia(I);
    Iz = I[2] * angmomfac;
    Lz = Iz * omega;

    if (c1)
      Lz0 = c1->angmom[2];
    else
      Lz0 = 0.0;

    posang = 0.0;
    lastomega = omega;
    lasttime = tnow;
    
    if (restart) {

      if (myid == 0) {
	
	// Backup up old file
	string backupfile = outdir + name + ".bak";
	string command("cp ");
	command += outdir + name + " " + backupfile;
	if (system(command.c_str()) == -1) {
	  std::cerr << "UserEBarN: error in executing <"
		    << command << ">" << endl;
	}

	// Open new output stream for writing
	std::ofstream out(string(outdir+name).c_str());
	if (!out) {
	  throw FileCreateError(outdir+name, "UserEbarN: error opening new log file",
				__FILE__, __LINE__);
	}
	
	// Open old file for reading
	ifstream in(backupfile.c_str());
	if (!in) {
	  throw FileCreateError(backupfile, "UserEbarN: error opening original log file",
				__FILE__, __LINE__);
	}

	const int linesize = 1024;
	char line[linesize];
	
	in.getline(line, linesize); // Discard header
	in.getline(line, linesize); // Next line

	double Lzp,  am1;
	bool firstime1 = true;
	while (in) {
	  istringstream ins(line);

	  lastomega = omega;

	  ins >> lasttime;
	  ins >> posang;
	  ins >> omega;
	  ins >> Lz;
	  ins >> Lzp;
	  ins >> am1;
	  ins >> bps[0];
	  ins >> bps[1];
	  ins >> bps[2];
	  ins >> vel[0];
	  ins >> vel[1];
	  ins >> vel[2];
	  ins >> acc[0];
	  ins >> acc[1];
	  ins >> acc[2];

	  if (firstime1) {
	    Lz0 = Lzp;
	    firstime1 = false;
	  }

	  if (lasttime >= tnow) break;

	  out << line << "\n";

	  in.getline(line, linesize); // Next line
	}

	cout << "UserEBarN: restart at T=" << lasttime 
	     << " with PosAng=" << posang
	     << ", Omega=" << omega
	     << ", Lz=" << Lz
	     << ", Lz0=" << Lz0
	     << endl;

      }

      MPI_Bcast(&lasttime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&posang, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&omega, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&lastomega, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&Lz, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&Lz0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&bps[0], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&vel[0], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&acc[0], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      // Recompute Lz from log output
      if (c1) Lz = Lz - Lz0 + c1->angmom[2];
    }

    firstime = false;
    update = true;

  } else if (mlevel==0) {

    if (!fixed) {
      if (c1)
	omega = (Lz + Lz0 - c1->angmom[2])/Iz;
      else
	omega = Lz/Iz;
    }
    else {
      if (dtom>0.0)
	omega = omega0*(1.0 + DOmega*0.5*(1.0 + erf( (tnow - T0)/dtom )));
      else
	omega = omega0*(1.0 + DOmega*(tnow - T0*0.5));
    }
    
    if ( fabs(tnow-lasttime) > 2.0*DBL_EPSILON) {
      posang += 0.5*(omega + lastomega)*(tnow - lasttime);
      lastomega = omega;
      lasttime = tnow;
      update = true;
    }
  }

				// Zero thread variables
  for (int n=0; n<nthrds; n++) {
    for (int k=0; k<3; k++) tacc[n][k] = 0.0;
  }

  if (timing) timer_thrd.start();
  exp_thread_fork(false);
  if (timing) timer_thrd.stop();

				// Get full contribution from all threads
  for (int k=0; k<3; k++) acc[k] = acc1[k] = 0.0;
  for (int n=0; n<nthrds; n++) {
    for (int k=0; k<3; k++) acc1[k] += tacc[n][k]/barmass;
  }

				// Get contribution from all processes
  MPI_Allreduce(&acc1[0], &acc[0], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				// Backward Euler
  if (mlevel==0 && monopole && monopole_follow) {
    for (int k=0; k<3; k++) {
      bps[k] += vel[k] * (tnow - teval[mlevel]);
      vel[k] += acc[k] * (tnow - teval[mlevel]);
    }
    for (unsigned m=mlevel; m<=multistep; m++) teval[m] = tnow;
  }

  if (myid==0 && update) 
    {
      ofstream out(string(outdir+name).c_str(), ios::out | ios::app);
      out.setf(ios::scientific);

      out << setw(15) << tnow
	  << setw(15) << posang
	  << setw(15) << omega;

      if (c1)
	out << setw(15) << Lz + Lz0 - c1->angmom[2]
	    << setw(15) << c1->angmom[2];
      else
	out << setw(15) << Lz
	    << setw(15) << 0.0;

      if (amplitude==0.0)
	out << setw(15) <<  0.0;
      else
	out << setw(15) << amplitude/fabs(amplitude) *  
	  0.5*(1.0 + erf( (tnow - Ton )/DeltaT )) *
	  0.5*(1.0 - erf( (tnow - Toff)/DeltaT ));

      for (int k=0; k<3; k++) out << setw(15) << bps[k];
      for (int k=0; k<3; k++) out << setw(15) << vel[k];
      for (int k=0; k<3; k++) out << setw(15) << acc[k];
      
      out << endl;
    }

  if (timing) {
    timer_tot.stop();
    cout << setw(20) << "Bar total: "
	 << setw(18) << timer_tot.getTime()  << endl
	 << setw(20) << "Bar threads: "
	 << setw(18) << timer_thrd.getTime() << endl;
    timer_tot.reset();
    timer_thrd.reset();
  }
}


void * UserEBarN::determine_acceleration_and_potential_thread(void * arg) 
{
  int id = *((int*)arg), nbodies, nbeg, nend, indx;
  double fac=0.0, ffac=0.0, dfac=0.0, amp=0.0, pot=0.0, dpot=0.0;
  double xx, yy, zz, rr, pp=0.0, extpot, M0=0.0, P0=0.0;
  vector<double> pos(3), acct(3); 
  double cos2p = cos(2.0*posang);
  double sin2p = sin(2.0*posang);

  double fraction_on =   0.5*(1.0 + erf( (tnow - Ton )/DeltaT )) ;
  double fraction_off =  0.5*(1.0 - erf( (tnow - Toff)/DeltaT )) ;

  double quad_onoff = 
    fraction_on*( (1.0 - quadrupole_frac) + quadrupole_frac * fraction_off );

  double mono_fraction = 
    0.5*(1.0 + erf( (tnow - TmonoOn )/DeltaMonoT )) *
    0.5*(1.0 - erf( (tnow - TmonoOff)/DeltaMonoT )) ;

  double mono_onoff = 
    (1.0 - monopole_frac) + monopole_frac*mono_fraction;

  if (amplitude==0.0) 
    amp = 0.0;
  else
    amp = amplitude * quad_onoff;


#if HAVE_LIBCUDA==1		// Cuda compatibility
  getParticlesCuda(cC);
#endif

  for (unsigned lev=mlevel; lev<=multistep; lev++) {

    nbodies = cC->levlist[lev].size();
    nbeg = nbodies*(id  )/nthrds;
    nend = nbodies*(id+1)/nthrds;

    for (int i=nbeg; i<nend; i++) {

      indx = cC->levlist[lev][i];
      
      for (int k=0; k<3; k++) pos[k] = cC->Pos(indx, k);

      if (c0)
	for (int k=0; k<3; k++) pos[k] -= c0->center[k];
      else if (monopole)
	for (int k=0; k<3; k++) pos[k] -= bps[k];
    
      xx = pos[0];
      yy = pos[1];
      zz = pos[2];
      rr = sqrt( xx*xx + yy*yy + zz*zz );

				// Variable sharpness potential
      if (quadpot(rr, P0, pot, dpot, M0)) {
	fac = pot/(rr*rr);
	dfac = (dpot - 2.0*pot/rr)/(rr*rr*rr);
	ffac = -amp*numfac;
	pp = (xx*xx - yy*yy)*cos2p + 2.0*xx*yy*sin2p;
      } else {
	ffac = 0.0;
      }

				// Quadrupole acceleration
      acct[0] = ffac*
	( 2.0*( xx*cos2p + yy*sin2p)*fac + dfac*pp*xx );
    
      acct[1] = ffac*
	( 2.0*(-yy*cos2p + xx*sin2p)*fac + dfac*pp*yy );

      acct[2] = ffac*
	( dfac*pp*zz );
    
				// Quadrupole potential
      extpot = -ffac*pp*fac;
    

				// Monopole contribution
      if (monopole) {

	M0 *= monoamp;
	
	if (monopole_onoff) M0 *= mono_onoff;
	
	for (int k=0; k<3; k++) {
				// Add monopole acceleration
	  acct[k] += -M0*pos[k]/(rr*rr*rr);

				// Force on bar (via Newton's 3rd law)
	  tacc[id][k] += -cC->Mass(indx) * acct[k];
	}

				// Monopole potential
	extpot += P0;
      }

				// Add bar acceleration to particle
      cC->AddAcc(indx, acct);
    
				// Add external potential
      cC->AddPotExt(indx, extpot);

    }
  }

  return (NULL);
}


extern "C" {
  ExternalForce *makerEBarN(const YAML::Node& conf)
  {
    return new UserEBarN(conf);
  }
}

class proxyebarn { 
public:
  proxyebarn()
  {
    factory["userebarn"] = makerEBarN;
  }
};

static proxyebarn p;
