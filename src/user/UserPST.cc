#include <math.h>
#include <sstream>

#include "expand.h"
#include <localmpi.h>
#include <UserPST.H>
#include <Timer.h>
#include <ZBrent.H>

static Timer timer_tot, timer_thrd;
static bool timing = false;

/*
  Param 0: r_max
  Param 1: rhoB
  Param 2: M_bulge
 */
double fit_mass(double rb, vector<double>& p)
{
  double y = p[0]/rb;
  return 4.0*M_PI*p[1]*rb*rb*rb*(asinh(y) - y/sqrt(1.0+y*y)) - p[2];
}

double UserPST::get_fid_bulge_dens()
{
  //
  // Fiducial parameters from Piner, Stone, & Teuben
  //

  double arat = 2.5;
  double Qm = 4.5e10;
  double rL = 6.0;
  double rhoC = 2.4e10;
  double nu = 1.0;

  //
  // Scaling
  //
  const double pc    = 3.08568025e18; // One parsec in cm
  const double Msun  = 1.98892e33;    // Solar mass in grams
  const double Grav  = 6.673e-8;      // In cgs
  const double Lunit = 300.0;	      // Virial length is 300 kpc
  const double Munit = 1.0e12;	// Virial mass is 1e12 solar masses
				// Time units in seconds 
				// (1 time unit is 2.45 Gyr)
  const double Tunit = 1.0/sqrt(Grav*Munit*Msun/pow(Lunit*1.0e3*pc, 3.0));
				// 1 km/s in system units (0.00835)
  const double Vunit = Tunit*1.0e5/(Lunit*1.0e3*pc);

  /*
    First compute the constants from the input parameters
  
    1) arat = a/b (major over minor axis ratio)
    2) Quad moment, Q_m
    3) Corotation radius, r_L
    4) \rho_c (central density) = \rho_0 + \rho_b (disk + bulge density)

    Fixed parameters: a=5 kpc (semi-major axis), nu=1,2 (Ferrer's index)

  */
				// Input parameters to system units
  Qm   /= Lunit*Lunit*Munit;
  rhoC /= Munit/(Lunit*Lunit*Lunit);
  rL   /= Lunit;

				// Implied bar parameters
  double a = 5.0/Lunit;
  double b = a/arat;
  double Mbar = Qm*(5.0 + 2.0*nu)/(a*a*(1.0 - 1.0/(arat*arat)));

  double rho0 = Mbar/( 2.0*M_PI*a*b*b *
		       exp(lgamma(nu+1.0)+lgamma(1.5)-lgamma(nu+2.5)) );

  double rhoB = rhoC - rho0;

  //
  // Disk surface density coef: v_o^2/(2*pi*G*r_o) where v_o = 200 km/s
  // Disk surface density scaling: r_o = 14.1 kpc
  //
  double v_o = 200.0 * Vunit;
  double r_o = 14.1 / Lunit;

  //
  // Mass within 16kpc
  //
  double r_max = 16.0/Lunit;
  double Mencl = 1.54e11/Munit;
  // double Mdisk = v_o*v_o*r_o*r_o*(1.0 - 1.0/sqrt(1.0 + r_max*r_max/(r_o*r_o)));
  double Mdisk = v_o*v_o*r_max/sqrt(1.0 + r_max*r_max/(r_o*r_o));
  double Mleft = Mencl - Mdisk - Mbar;

  if (Mleft<0) {
    cerr << "Impossible parameters: Mencl (" << Mencl << ") < Mdisk + Mbar ("
	 << Mdisk + Mbar << ")" << endl;
    exit(-1);
  }
  
  Mscale = Mencl/(Mdisk + Mleft);

  double r_b;
  ZBrent< vector<double> > zb;
  vector<double> param;
/*
  Param 0: r_max
  Param 1: rhoB
  Param 2: M_bulge
 */
  param.push_back(r_max);
  param.push_back(rhoB);
  param.push_back(Mleft);

  zb.find(fit_mass, param, 1.0e-5, 1.0, 1.0e-8, r_b);

  return r_b;
}

UserPST::UserPST(string &line) : ExternalForce(line)
{
  id = "Piner-Stone-Teuben";

                                // Parameters
				// ----------
  rmin = 0.0001;		// Minimum model radius
  rmax = 2.0;			// Maximum model radius
  numr = 2000;		    // Number of radial points for bulge model
  blog = false;		    // Use log spacing for bulge model
  dlog = false;		    // Use log spacing for disk model

  arat = 2.5;		      // Maximum to minimum axis ratio for bar
  Qm   = 4.5e10;	      // Quadrupole strength
  rL   = 6.0;		      // Corotation radius
  rhoC = 2.4e10;	      // Central density
  nu   = 1.0;		      // Ferrer's index

  Lmax = 32;	      // Maximum harmonic index for disk pot expansion
  Nmax = 16;	      // Maximum radial index for expansion
  numR = 10000;	      // Number of radial table elements
  numt = 400;	      // Integration knots for theta
  numg = 100;	      // 


  Ton    = 0.5e8;		// Turn on start time (in years)
  DeltaT = 0.25e8;		// Turn on width (in years)

				// Output file name
  filename = outdir + "BarRot." + runtag;

  initialize();


  //
  // Scaling
  //
  // const double mp = 1.67e-24;	      // Proton mass
  // const double kb  = 1.381e-16;	      // Boltzmann constant
  const double pc    = 3.08568025e18; // One parsec in cm
  const double Msun  = 1.98892e33;    // Solar mass in grams
  const double Grav  = 6.673e-8;      // In cgs
  const double Lunit = 300.0;	      // Virial length is 300 kpc
  const double Munit = 1.0e+12;	// Virial mass is 1e12 solar masses
				// Time units in seconds 
				// (1 time unit is 2.45 Gyr)
  const double Tunit = 1.0/sqrt(Grav*Munit*Msun/pow(Lunit*1.0e3*pc, 3.0));
				// 1 km/s in system units (0.00835)
  const double Vunit = Tunit*1.0e5/(Lunit*1.0e3*pc);

  // const double r_max = 16.0/Lunit;	       // Mass within 16kpc
  // const double c_s = 5.0;		       // Sounds speed: 5 km/s
  // const double T  = (mp*pow(c_s*1e5, 2.0))/kb; // Temperature

  Ton    *= 365.25*24*3600/Tunit;
  DeltaT *= 365.25*24*3600/Tunit;

  /*
    First compute the constants from the input parameters
  
    1) arat = a/b (major over minor axis ratio)
    2) Quad moment, Q_m
    3) Corotation radius, r_L
    4) \rho_c (central density) = \rho_0 + \rho_b (disk + bulge density)

    Fixed parameters: a=5 kpc (semi-major axis), nu=1,2 (Ferrer's index)

  */
				// Input parameters to system units
  Qm   /= Lunit*Lunit*Munit;
  rhoC /= Munit/(Lunit*Lunit*Lunit);
  rL   /= Lunit;
  
				// Implied bar and bulge parameters
  a    = 5.0/Lunit;
  b    = a/arat;
  Mbar = Qm*(5.0 + 2.0*nu)/(a*a*(1.0 - 1.0/(arat*arat)));
  rho0 = Mbar/( 2.0*M_PI*a*b*b *
		       exp(lgamma(nu+1.0)+lgamma(1.5)-lgamma(nu+2.5)) );
  rhoB = rhoC - rho0;

  //
  // Disk surface density coef: v_o^2/(2*pi*G*r_o) where v_o = 200 km/s
  // Disk surface density scaling: r_o = 14.1 kpc
  //
  double v_o = 200.0 * Vunit;
  double r_o = 14.1 / Lunit;
  // double coef = v_o*v_o/(2.0*M_PI*r_o);

  //
  // Compute r_b by enforcing enclosed mass within 10 kpc
  //

  double r_b = get_fid_bulge_dens();

  //
  // Now compute spherical bulge model
  //

  double trmin = rmin;
  double trmax = rmax;
  if (rmin > 0.0 && blog) {
    trmin = log(rmin);
    trmax = log(rmax);
  } else {
    blog = false;
  }

  vector<double> rr(numr), dd(numr), mm(numr, 0), pp(numr, 0), p2(numr, 0);
  double dr = (trmax - trmin)/(numr - 1);
  for (int i=0; i<numr; i++) {
    rr[i] = trmin + dr*i;
    if (blog) rr[i] = exp(rr[i]);
    dd[i] = rhoB*pow(1.0 + (rr[i]/r_b)*(rr[i]/r_b), -1.5);
  }

  for (int i=1; i<numr; i++) {

    if (blog) {
      dr = log(rr[i]) - log(rr[i-1]);
      mm[i] = mm[i-1] + 
	2.0*M_PI*(dd[i-1]*rr[i-1]*rr[i-1]*rr[i-1] + 
		  dd[i  ]*rr[i  ]*rr[i  ]*rr[i  ] ) * dr;
      p2[i] = p2[i-1] + 
	2.0*M_PI*(dd[i-1]*rr[i-1]*rr[i-1] + 
		  dd[i  ]*rr[i  ]*rr[i  ] ) * dr;
    } else {
      dr = rr[i] - rr[i-1];
      mm[i] = mm[i-1] + 
	2.0*M_PI*(dd[i-1]*rr[i-1]*rr[i-1] + 
		  dd[i  ]*rr[i  ]*rr[i  ] ) * dr;
      p2[i] = p2[i-1] + 
	2.0*M_PI*(dd[i-1]*rr[i-1] + 
		  dd[i  ]*rr[i  ] ) * dr;
    }
  }

  for (int i=0; i<numr; i++) 
    pp[i] = -mm[i]/(rr[i]+1.0e-14) + p2[i] -  p2[numr-1];
  

  if (myid==0 && VERBOSE > 5) {

    ofstream out(string(outdir + "bulgemod.dat").c_str());
    for (int i=0; i<numr; i++) 
      out << setw(16) << rr[i]
	  << setw(16) << dd[i]
	  << setw(16) << mm[i]
	  << setw(16) << pp[i]
	  << endl;
  }

  bulge = new SphericalModelTable(numr, 
				  &rr[0]-1, &dd[0]-1, &mm[0]-1, &pp[0]-1,
				  0, 0.0, 0, "Bulge model");

  if (myid==0 && VERBOSE > 5) cout << "Bulge created" << endl;

  //
  // Bar model
  //
  EllipsoidForce::BarType bartype = EllipsoidForce::ferrers;
  int num  = 400;
  int numT = 200;
  
  bar = new EllipsoidForce(a, b, b, Mbar, rmin, 5.0*a, 
			   bartype, false, nu, num, numr);

  if (myid==0 && VERBOSE>5) cout << "Bar created" << endl;
  if (myid==0 && VERBOSE>5) cout << "Table creation . . . " << flush;
  bar->MakeTable(numT, numT, numT);
  if (myid==0 && VERBOSE>5) cout << "done" << endl;

  //
  // Disk model
  //
  vector<double> mparam;
  mparam.push_back(v_o);
  mparam.push_back(r_o);
  mparam.push_back(0.1/Lunit);

  int Lmax = 32, Nmax = 16;
  int numR = 10000;
  int numt = 400, numg = 100;

  disk.Initialize(rmin, rmax, dlog, Nmax, Lmax, numR, numt, numg, mparam);
    
  if (myid==0 && VERBOSE>5) cout << "Disk created" << endl;


  //
  // Compute pattern speed
  //

  double bfrc = -bulge->get_mass(rL)/(rL*rL), fr, fz;
  disk.force_eval(rL, 0.0, 0.0, fr, fz);

  omega = sqrt(-(bfrc+fr)*Mscale/rL);
  
  userinfo();

  //
  // Only turn on bar timing for extreme debugging levels
  //
  if (VERBOSE>49) timing = true;
}

UserPST::~UserPST()
{
  delete bulge;
  delete bar;
}

void UserPST::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine Piner-Stone-Teuben disk & bar model, " ;
  cout << "prescribed pattern speed with"
       << " omega=" <<  omega << ", t_on=" << Ton << ", delta_t=" << DeltaT
       << ", a=" << a << ", b=c=" << b
       << ", Qm=" << Qm << ", rL=" << rL << ", Mbar=" << Mbar
       << ", rho(total)=" << rhoC << ", rho(bar)=" << rho0 
       << ", rho(bulge)=" << rhoB << ", Ferrer's index=" << nu;

  if (blog)
    cout << ", using log spacing for bulge model";

  if (dlog) 
    cout << ", using log spacing for disk model";

  cout << endl;

  print_divider();
}

void UserPST::initialize()
{
  string val;

  if (get_value("rmin", val))	        rmin     = atof(val.c_str());
  if (get_value("rmax", val))	        rmax     = atof(val.c_str());
  if (get_value("numr", val))	        numr     = atoi(val.c_str());
  if (get_value("blog", val))	        blog     = atol(val);
  if (get_value("dlog", val))	        dlog     = atol(val);
  if (get_value("arat", val))		arat     = atof(val.c_str());
  if (get_value("Qm", val))		Qm       = atof(val.c_str());
  if (get_value("rL", val))		rL       = atof(val.c_str());
  if (get_value("rhoC", val))		rhoC     = atof(val.c_str());
  if (get_value("nu", val))		nu       = atof(val.c_str());
  if (get_value("Lmax", val))		Lmax     = atoi(val.c_str());
  if (get_value("Nmax", val))		Nmax     = atoi(val.c_str());
  if (get_value("numR", val))		numR     = atoi(val.c_str());
  if (get_value("numt", val))		numt     = atoi(val.c_str());
  if (get_value("numg", val))		numg     = atoi(val.c_str());
  if (get_value("Ton", val))		Ton      = atof(val.c_str());
  if (get_value("DeltaT", val))		DeltaT   = atof(val.c_str());
  if (get_value("filename", val))	filename = val;
}


void UserPST::determine_acceleration_and_potential(void)
{
  if (timing) timer_tot.start();
  
  if (timing) timer_thrd.start();
  exp_thread_fork(false);
  if (timing) timer_thrd.stop();

  if (timing) {
    timer_tot.stop();
    cout << setw(20) << "Bar total: "
	 << setw(18) << timer_tot.getTime()  << endl
	 << setw(20) << "Bar threads: "
	 << setw(18) << timer_thrd.getTime() << endl;
    timer_tot.reset();
    timer_thrd.reset();
  }

  print_timings("UserPST: acceleration timings");
}


void * UserPST::determine_acceleration_and_potential_thread(void * arg) 
{
  int id = *((int*)arg), nbodies, nbeg, nend, indx;
  double xx, yy, zz, rr, bfrc, fr, fz, extpot;
  vector<double> pos(3), pos1(3), acct(3), force(4), acc(3);

  double posang = omega*tnow;

  double cosp = cos(posang);
  double sinp = sin(posang);

  thread_timing_beg(id);

  double elip_frac = 0.5*(1.0 + erf( (tnow - Ton )/DeltaT ));
  double disk_frac = 1.0 + (Mscale - 1.0)*(1.0 - elip_frac);

  for (unsigned lev=mlevel; lev<=multistep; lev++) {

    nbodies = cC->levlist[lev].size();
    nbeg = nbodies*(id  )/nthrds;
    nend = nbodies*(id+1)/nthrds;

    for (int i=nbeg; i<nend; i++) {
      
      indx = cC->levlist[lev][i];
      
      for (int k=0; k<3; k++) pos[k] = cC->Pos(indx, k);

      xx = pos[0];
      yy = pos[1];
      zz = pos[2];
      rr = sqrt( xx*xx + yy*yy + zz*zz );

      //
      // Rotation for bar force
      //
      pos1[0] =  xx*cosp + yy*sinp;
      pos1[1] = -xx*sinp + yy*cosp;
      pos1[2] =  zz;

      bar->TableEval(pos1, force);

      acc[0] = force[1]*cosp - force[2]*sinp;
      acc[1] = force[1]*sinp + force[2]*cosp;
      acc[2] = force[3];
      
      //
      // Bulge force
      //
      bfrc = -bulge->get_mass(rr)/(rr*rr*rr + 1.0e-30);
      disk.force_eval(xx, yy, zz, fr, fz);
      

      for (int k=0; k<3; k++) 
	acct[k] = 
	  (fr/(rr+1.0e-10) + bfrc)*pos[k] * disk_frac + acc[k]*elip_frac;
    
      extpot = 
	(bulge->get_pot(rr) + disk.potential_eval(xx, yy, zz)) * disk_frac + 
	force[0]*elip_frac;
    
				// Add bar acceleration to particle
      cC->AddAcc(indx, acct);

				// Add external potential
      cC->AddPotExt(indx, extpot);

    }
  }

  thread_timing_end(id);

  return (NULL);
}


extern "C" {
  ExternalForce *makerPST(string& line)
  {
    return new UserPST(line);
  }
}

class proxypst { 
public:
  proxypst()
  {
    factory["userpst"] = makerPST;
  }
};

proxypst p;
