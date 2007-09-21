#include <math.h>
#include <sstream>

#include "expand.h"
#include <localmpi.h>
#include <UserEBarN.H>
#include <Timer.h>
static Timer timer_tot(true), timer_thrd(true);
static bool timing = false;

UserEBarN::UserEBarN(string &line) : ExternalForce(line)
{
  id = "RotatingBarWithMonopole";

  length = 1.0;			// Bar length
  bratio = 0.5;			// Ratio of b to a
  cratio = 0.1;			// Ratio of c to b
  amplitude = 0.3;		// Bar quadrupole amplitude
  angmomfac = 1.0;		// Artifically change the total bar ang mom
  barmass = 1.0;		// Total bar mass
  Ton = -20.0;			// Turn on start time
  Toff = 200.0;			// Turn off start time
  TmonoOn = -20.0;		// Turn on start time for monopole
  TmonoOff = 200.0;		// Turn off start time monopole
  DeltaT = 1.0;			// Turn on duration
  DeltaMonoT = 1.0;		// Turn on duration for monopole
  DOmega = 0.0;			// Change in pattern speed
  tom0 = 1.0;			// Midpoint of forced bar slow down
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

  expon = false;		// Use the exponential disk density rather than power law
  modelp = 0.0;			// Either the power law exponent or the exponential scale length
  rsmin  = 0.001;		// The minimum radius for the quadrupole table in major axis units
  rsmax  = 100;	                // The maximum radius for the quadrupole table in major axis units
  numt   = 400;                 // Number of table entries for quadrupole table
  N      = 100;			// Number of integration knots

				// Output file name
  filename = outdir + "BarRot." + runtag;

  firstime = true;
  omega0 = -1.0;

  ctr_name = "";		// Default component for com
  angm_name = "";		// Default component for angular momentum
  
  ellip = 0;

  initialize();

  if (ctr_name.size()>0) {
				// Look for the fiducial component for
				// centering
    bool found = false;
    list<Component*>::iterator cc;
    Component *c;
    for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
      c = *cc;
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
    list<Component*>::iterator cc;
    Component *c;
    for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
      c = *cc;
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
  tacc = new double* [nthrds];
  for (int n=0; n<nthrds; n++) tacc[n] = new double [3];

  userinfo();

  // Only turn on bar timing for extreme debugging levels
  if (VERBOSE>49) timing = true;

  gt = 0;
  gq = new LegeQuad(N);
}

UserEBarN::~UserEBarN()
{
  for (int n=0; n<nthrds; n++) delete [] tacc[n];
  delete [] tacc;
  delete ellip;
  delete gt;
  delete gq;
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
      cout << "using monopole with turn-on/off with fraction=" 
	   << monopole_frac << ", TmonoOn=" << TmonoOn
	   << ", TmonoOff=" << TmonoOff << ", DeltaMonoT=" << DeltaMonoT
	   << ", ";
    else
      cout << "using monopole, ";
    if (monopole_follow)
      cout << "self-consistent monopole centering, ";
    else
      cout << "monopole center fixed, ";
  }
  else
    cout << "without monopole, ";

  cout << "bar softness (alpha)=" << alpha << ", ";

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
  string val;

  if (get_value("ctrname", val))	ctr_name = val;
  if (get_value("angmname", val))	angm_name = val;
  if (get_value("length", val))		length = atof(val.c_str());
  if (get_value("bratio", val))		bratio = atof(val.c_str());
  if (get_value("cratio", val))		cratio = atof(val.c_str());
  if (get_value("amp", val))		amplitude = atof(val.c_str());
  if (get_value("angmomfac", val))	angmomfac = atof(val.c_str());
  if (get_value("barmass", val))	barmass = atof(val.c_str());
  if (get_value("Ton", val))		Ton = atof(val.c_str());
  if (get_value("Toff", val))		Toff = atof(val.c_str());
  if (get_value("TmonoOn", val))	TmonoOn = atof(val.c_str());
  if (get_value("TmonoOff", val))	TmonoOff = atof(val.c_str());
  if (get_value("DeltaT", val))		DeltaT = atof(val.c_str());
  if (get_value("DeltaMonoT", val))	DeltaMonoT = atof(val.c_str());
  if (get_value("DOmega", val))		DOmega = atof(val.c_str());
  if (get_value("tom0", val))     	tom0 = atof(val.c_str());
  if (get_value("dtom", val))     	dtom = atof(val.c_str());
  if (get_value("T0", val))		T0 = atof(val.c_str());
  if (get_value("Fcorot", val))		Fcorot = atof(val.c_str());
  if (get_value("omega", val))		omega0 = atof(val.c_str());
  if (get_value("fixed", val))		fixed = atoi(val.c_str()) ? true:false;
  if (get_value("self", val))		fixed = atoi(val.c_str()) ? false:true;
  if (get_value("alpha", val))		alpha = atof(val.c_str());
  if (get_value("monopole", val))	monopole = atoi(val.c_str()) ? true:false;  
  if (get_value("follow", val))		monopole_follow = atoi(val.c_str()) ? true:false;
  if (get_value("onoff", val))		monopole_onoff = atoi(val.c_str()) ? true:false;
  if (get_value("monofrac", val))	monopole_frac = atof(val.c_str());
  if (get_value("quadfrac", val))	quadrupole_frac = atof(val.c_str());
  if (get_value("filename", val))	filename = val;
  if (get_value("expon", val))          expon = atoi(val.c_str()) ? true : false;
  if (get_value("modelp", val))         modelp = atof(val.c_str());
  if (get_value("rmin", val))           rsmin = atof(val.c_str());
  if (get_value("rmax", val))           rsmax = atof(val.c_str());
  if (get_value("numt", val))           numt = atoi(val.c_str());
}


vector<double> a(3);

double find_fct(double u, vector<double> z)
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
  
  double ans;
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

double UserEBarN::Potential(vector<double> x)
{
  double mshell, ans=0.0, u, d, t, denom, m2;
  double tmin, tmax=0.5*M_PI;

  double ellip = 0.0;
  for (int k=0; k<3; k++) ellip += x[k]*x[k]/(a[k]*a[k]);
  ellip -= 1.0;

  if (ellip<0.0) {		// Inside
    tmin = 0.0;
  } else {			// Outside
    tmin = atan(solve(x, 1.0));
  }

  for (int i=1; i<=N; i++) {
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

    if (expon)
      mshell = 2.0*rho0*modelp/a[0]*(1.0 - exp(-a[0]*sqrt(m2)/modelp));
    else 
      mshell = rho0/(modelp+1.0)*(1.0 - pow(m2, modelp+1.0));
    
    ans += d*gq->weight(i) * mshell/sqrt(denom);
  }

  return ans*(tmax - tmin);
}


double UserEBarN::U22(double r) 
{
  const double numfac = 0.25*sqrt(15.0/(2.0*M_PI));
  vector<double> z(3);
  const int nphi = 20, ntheta = 40;
  const double dphi = 2.0*M_PI/nphi;

  double cosx, sinx, phi, ans=0.0;

  if (gt==0) gt = new LegeQuad(ntheta);

  for (int i=0; i<nphi; i++) {
    phi = dphi*i;
    for (int j=1; j<=ntheta; j++) {
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


bool UserEBarN::quadpot(double r, double& f, double& fr)
{
  if (r<rr.front() || r>=rr.back()) {
    f = fr = 0.0;
    return false;
  }

  double lr = log(r);
  int indx = (int)floor((log(r) - lrmin)/ldr);
  double a = (lrmin + ldr*(indx+1) - lr)/ldr;
  double b = (lr - lrmin - ldr*indx)/ldr;

  f = a*uu[indx] + b*uu[indx+1];
  fr = (uu[indx+1] - uu[indx])/ldr;

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
    
    ellip = new EllipForce(length, length*bratio, length*bratio*cratio,
			   barmass, 200, 200);

    list<Component*>::iterator cc;
    Component *c;

    if (omega0 < 0.0) {

      double R=length*Fcorot;
      double phi, theta=0.5*M_PI;
      double dens0, potl0, dens, potl, potr, pott, potp;
      double avg=0.0;
      
      for (int n=0; n<8; n++) {
	phi = 2.0*M_PI/8.0 * n;

	for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
	  c = *cc;
	
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

    mass = barmass * fabs(amplitude);

    a[0] = length;
    a[1] = bratio*a[0];
    a[2] = cratio*a[1];

    // *M_PI*a[0]*a[1]*a[2] sucked into the leading factor for the 
    // gravitational potential
    //
    if (expon)
      rho0 = a[0]*a[0]*mass/(4.0*modelp*modelp) / 
	( 1.0 - (1.0 + a[0]/modelp)*exp(-a[0]/modelp) ); 
    else 
      rho0 = (2.0*modelp + 3.0)*mass/4.0; 
    
    lrmin = log(rsmin*a[0]);
    lrmax = log(rsmax*a[0]);
    ldr = (lrmax - lrmin)/(numt-1);
    double r;

    for (int j=0; j<numt; j++) {
      r = exp(lrmin + ldr*j);
      rr.push_back(r);
      uu.push_back(U22(r));
    }

    if (myid==0) {

      if (VERBOSE>3) {
	string name = "UserEBarN." + runtag + ".debug";
	ofstream out(name.c_str());

	out << "# a1, a2, a3=" << a[0] << ", " << a[1] << ", " << a[2] << endl;
	out << "# Expon=" << expon << endl;
	out << "# Param=" << modelp << endl;
	out << "# Mass=" << barmass << endl;
	out << "# Ampl=" << amplitude << endl;
	out << "# I_3=" << 0.2*mass*(a[0]*a[0] + a[1]*a[1]) << endl;
	out << "#" << endl;
	for (int j=0; j<numt; j++)
	  out << setw(18) << rr[j] << setw(18) << uu[j] << endl;
      }

      cout << "====================================================\n";
      cout << "Computed quadrupole fit for an ellipsoid using\n";
      if (expon)
	cout << "an exponential density with scale length " << modelp << endl;
      else
	cout << "an power-law density with exponent " << modelp << endl;
      cout << "Rho0=" << rho0 << endl;
      cout << "Mass=" << barmass << endl;
      cout << "Ampl=" << amplitude << endl;
      cout << "I_3=" << 0.2*mass*(a[0]*a[0] + a[1]*a[1]) << endl;
      cout << "Omega(0)=" << omega0 << endl;
      cout << "====================================================\n";
      
      name = filename;
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

    Iz = 0.2*mass*(a[0]*a[0] + a[1]*a[1]) * angmomfac;
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
	string backupfile = name + ".bak";
	string command("cp ");
	command += name + " " + backupfile;
	system(command.c_str());

	// Open new output stream for writing
	ofstream out(name.c_str());
	if (!out) {
	  cout << "UserEBarN: error opening new log file <" 
	       << filename << "> for writing\n";
	  MPI_Abort(MPI_COMM_WORLD, 121);
	  exit(0);
	}
	
	// Open old file for reading
	ifstream in(backupfile.c_str());
	if (!in) {
	  cout << "UserEBarN: error opening original log file <" 
	       << backupfile << "> for reading\n";
	  MPI_Abort(MPI_COMM_WORLD, 122);
	  exit(0);
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

  } else {

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
  if (monopole && monopole_follow) {
    for (int k=0; k<3; k++) {
      bps[k] += vel[k] * (tnow - teval[mlevel]);
      vel[k] += acc[k] * (tnow - teval[mlevel]);
    }
    for (unsigned m=mlevel; m<=multistep; m++) teval[m] = tnow;
  }

  if (myid==0 && update) 
    {
      ofstream out(name.c_str(), ios::out | ios::app);
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
	 << setw(18) << 1.0e-6*timer_tot.getTime().getRealTime() << endl
	 << setw(20) << "Bar threads: "
	 << setw(18) << 1.0e-6*timer_thrd.getTime().getRealTime() << endl;
    timer_tot.reset();
    timer_thrd.reset();
  }
}


void * UserEBarN::determine_acceleration_and_potential_thread(void * arg) 
{
  int id = *((int*)arg), nbodies, nbeg, nend, indx;
  double fac=0.0, ffac=0.0, dfac=0.0, amp=0.0, pot=0.0, dpot=0.0;
  double xx, yy, zz, rr, pp=0.0, extpot, M0=0.0;
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
    amp = afac * amplitude/fabs(amplitude) * quad_onoff;


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
      if (quadpot(rr, pot, dpot)) {
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

	M0 = ellip->getMass(rr);
	
	if (monopole_onoff) M0 *= mono_onoff;
	
	for (int k=0; k<3; k++) {
				// Add monopole acceleration
	  acct[k] += -M0*pos[k]/(rr*rr*rr);

				// Force on bar (via Newton's 3rd law)
	  tacc[id][k] += -cC->Mass(indx) * acct[k];
	}

				// Monopole potential
	extpot += ellip->getPot(rr);
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
  ExternalForce *makerEBarN(string& line)
  {
    return new UserEBarN(line);
  }
}

class proxyebarn { 
public:
  proxyebarn()
  {
    factory["userebarn"] = makerEBarN;
  }
};

proxyebarn p;
