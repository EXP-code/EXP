                                // C++/STL headers
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>
#include <random>

                                // System libs
#include <getopt.h>

                                // MDW classes
#include <numerical.H>
#include <EllipsoidForce.H>
#include <CylindricalDisk.H>
#include <massmodel.H>
#include <interp.H>
#include <ZBrent.H>
#include <localmpi.H>
#include <libvars.H>

using namespace __EXP__;	// Reference to n-body globals

                                // For debugging
#ifdef DEBUG
#include <fpetrap.h>
#endif


// EXP library support
//
#include <libvars.H>
#include <localmpi.H>

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

double Mscale;

double get_fid_bulge_dens()
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

				// Bar semi-major axis
  double a    = 5.0/Lunit;
  double b    = a/arat;
  double Mbar = Qm*(5.0 + 2.0*nu)/(a*a*(1.0 - 1.0/(arat*arat)));

  double rho0 = Mbar/( 2.0*M_PI*a*b*b *
		       exp(lgamma(nu+1.0)+lgamma(1.5)-lgamma(nu+2.5)) );

  double rhoB = rhoC - rho0;

  //
  // Disk surface density coef: v_o^2/(2*pi*G*r_o) where v_o = 200 km/s
  // Disk surface density scaling: r_o = 14.1 kpc
  //
  double v_o  = 200.0 * Vunit;
  double r_o  = 14.1 / Lunit;

  //
  // Mass within 16kpc
  //
  double r_max = 16.0/Lunit;
  double Mencl = 1.54e11/Munit;
  double Mdisk = v_o*v_o*r_max/sqrt(1.0 + r_max*r_max/(r_o*r_o));
  double Mleft = Mencl - Mdisk - Mbar;

  if (Mleft<0) {
    cerr << "Impossible parameters: Mencl (" << Mencl << ") < Mdisk + Mbar ("
	 << Mdisk + Mbar << ")" << endl;
    exit(-1);
  }
  
  cout << "r_max = " << r_max << endl;
  cout << "Mdisk = " << Mdisk << endl;
  cout << "Rho_b = " << rhoB  << endl;
  cout << "M_b   = " << Mleft << endl;
  cout << "Mbar  = " << Mbar  << endl;
  cout << "Scale = " << (Mscale=Mencl/(Mdisk+Mleft)) << endl;

  double r_b;
  ZBrent< vector<double> > zb, zo;
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

  cout << "R_b   = " << r_b << endl;
  cout << "Fit   = " << fit_mass(r_b, param) << endl;

  return r_b;
}

void epi2(double r, SphericalModelTable& bulge, CylindricalDisk& disk,
	  EllipsoidForce& bar, double& kx2, double& ky2)
{
  //
  // kappa^2 = 3.0*Omega^2 + d^2U/dr^2
  //
  // Use simple 3 point scheme
  
  const double dfrac = 0.1;
  double fr, fz, bfrc, r1, omg0, pot_p, pot_m, pot0;
  vector<double> x(3, 0), force(4);
  double dr = r*dfrac;
  
  //
  // X-axis
  //
  x[0] = r;
  x[1] = 0.0;
  bar.TableEval(x, force);
  pot0 = bulge.get_pot(r) + disk.potential_eval(r, 0.0, 0.0) + force[0];
  
  bfrc = -bulge.get_mass(r)/(r*r);
  disk.force_eval(r, 0.0, 0.0, fr, fz);

  omg0 = -(force[1]+bfrc+fr)/r;

  r1 = r + dr;
  x[0] = r1;
  bar.TableEval(x, force);
  pot_p = bulge.get_pot(r1) + disk.potential_eval(r1, 0.0, 0.0) + force[0];

  r1 = r - dr;
  x[0] = r1;
  bar.TableEval(x, force);
  pot_m = bulge.get_pot(r1) + disk.potential_eval(r1, 0.0, 0.0) + force[0];

  kx2 = 3.0*omg0 + (pot_p + pot_m - 2.0*pot0)/(dr*dr);

  //
  // Y-axis
  //
  x[0] = 0.0;
  x[1] = r;
  bar.TableEval(x, force);
  pot0 = bulge.get_pot(r) + disk.potential_eval(0.0, r, 0.0) + force[0];
  
  bfrc = -bulge.get_mass(r)/(r*r);
  disk.force_eval(0.0, r, 0.0, fr, fz);

  omg0 = -(force[2]+bfrc+fr)/r;

  r1 = r + dr;
  x[1] = r1;
  bar.TableEval(x, force);
  pot_p = bulge.get_pot(r1) + disk.potential_eval(0.0, r1, 0.0) + force[0];

  r1 = r - dr;
  x[1] = r1;
  bar.TableEval(x, force);
  pot_m = bulge.get_pot(r1) + disk.potential_eval(0.0, r1, 0.0) + force[0];

  ky2 = 3.0*omg0 + (pot_p + pot_m - 2.0*pot0)/(dr*dr);

}

int 
main(int argc, char **argv)
{
  //
  // Inialize MPI stuff
  //
  // local_init_mpi(argc, argv);
  
  //
  // Now parse the rest of the arguments
  //

                                // Parameters
  double rmin = 0.0001;
  double rmax = 2.0;

  int    numr = 2000;
  int  Number = 0;
  int    SEED = 11;

  bool   blog = false;
  bool   dlog = false;

  double H    = 0.1;
  double arat = 2.5;
  double Qm   = 4.5e10;
  double rL   = 6.0;
  double rhoC = 2.4e10;
  double nu   = 1.0;
  double gasD = 10.0;
  double spd  = 5.0;
  double Rgas = 30.0;

  string file = "gas.bods";

  //========================= Parse command line ==============================

  int c;
  // int digit_optind = 0;

  while (1)
    {
      // int this_option_optind = optind ? optind : 1;
      int option_index = 0;
      static struct option long_options[] = {
	{"arat",      1, 0, 0},
	{"rmin",      1, 0, 0},
	{"rmax",      1, 0, 0},
	{"numr",      1, 0, 0},
	{"Number",    1, 0, 0},
	{"Rgas",      1, 0, 0},
	{"H",         1, 0, 0},
	{"Qm",        1, 0, 0},
	{"rL",        1, 0, 0},
	{"rhoC",      1, 0, 0},
	{"nu",        1, 0, 0},
	{"SEED",      1, 0, 0},
	{"file",      1, 0, 0},
	{"logB",      0, 0, 0},
	{"logD",      0, 0, 0},
	{0, 0, 0, 0}
      };

      c = getopt_long (argc, argv, 
		       "a:r:R:n:N:Q:L:C:S:G:o:lh",
		       long_options, &option_index);
      if (c == -1)
        break;

      string optname;

      switch (c)
        {
	case 0:			// Long options
	  optname = string(long_options[option_index].name);
	  if (!optname.compare("arat"))      arat   = atof(optarg);
	  if (!optname.compare("rmin"))      rmin   = atof(optarg);
	  if (!optname.compare("rmax"))      rmax   = atof(optarg);
	  if (!optname.compare("numr"))      numr   = atoi(optarg);
	  if (!optname.compare("H"))         H      = atof(optarg);
	  if (!optname.compare("Qm"))        Qm     = atof(optarg);
	  if (!optname.compare("rL"))        rL     = atof(optarg);
	  if (!optname.compare("rhoC"))      rhoC   = atof(optarg);
	  if (!optname.compare("nexp"))      nu     = atof(optarg);
	  if (!optname.compare("Number"))    Number = atoi(optarg);
	  if (!optname.compare("Rgas"))      Rgas   = atoi(optarg);
	  if (!optname.compare("SEED"))      SEED   = atoi(optarg);
	  if (!optname.compare("file"))      file   = optarg;
	  if (!optname.compare("logD"))      dlog   = true;
	  if (!optname.compare("logB"))      blog   = true;
	  break;

        case 'a':
          arat = atof(optarg);
          break;

        case 'r':
          rmin = atof(optarg);
          break;

        case 'R':
          rmax = atof(optarg);
          break;

        case 'n':
          numr = atoi(optarg);
          break;

        case 'N':
          nu = atof(optarg);
          break;

        case 'Q':
          Qm = atof(optarg);
          break;

        case 'L':
          rL = atof(optarg);
          break;

        case 'H':
          H = atof(optarg);
          break;

        case 'G':
          Rgas = atof(optarg);
          break;

        case 'S':
          SEED = atoi(optarg);
          break;

        case 'o':
          file = optarg;
          break;

        case 'l':
          blog = true;
          break;

        case 'C':
          rhoC = atof(optarg);
          break;

        case 'h':
        case '?':
        default:
          cout << "\nGenerates a spherical phase space with an embedded disk\n"
               << "  using Jeans' equations assuming velocity isotropy for\n"
               << "  the spherical component and in-plane isotropy for the\n"
               << "  disk component\n\n"
               << "\nUsage : mpirun -v -np <# of nodes> " << argv[0] 
               << " -- [options]\n\n"
               << "  -a, --arat  f    semi-major/semi-minor ratio"
	       << endl
               << "  -H, --H     f    scale height"
	       << endl
               << "  -r, --rmin  f    minimum model radius"
	       << endl
               << "  -R, --rmax  f    maximum model radius"
	       << endl
               << "  -G, --Rgas  f    radius of realized disk"
	       << endl
               << "  -n, --numr  i    number of mesh points for bulge model"
	       << endl
               << "  -N, --nexp  f    Ferrer's model index"
	       << endl
               << "  -Q, --Qm    f    quadrupole amplitude"
	       << endl
               << "  -L, --rL    f    corotation radius"
	       << endl
               << "  -C, --rhoC  f    central density"
	       << endl
               << "  -o, --file  s     file name for gas particle ICs"
	       << endl
               << "  -l, --logB       use logarithmic grid for bulge model"
	       << endl
               << "  --logD           use logarithmic grid for disk model"
	       << endl
               << "  --Number         Number of gas particles"
	       << endl
	       << "  -h               this help message"
	       << endl
               << endl;
          exit(0);
        }
    }

  //===========================================================================

#ifdef DEBUG                    // For gdb . . . 
  sleep(20);
  set_fpu_handler();            // Make gdb trap FPU exceptions
#endif

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
  H    /= Lunit;
  Qm   /= Lunit*Lunit*Munit;
  rhoC /= Munit/(Lunit*Lunit*Lunit);
  rL   /= Lunit;
  gasD /= Munit/(Lunit*Lunit*1.0e6);
  
  
  double r_max = 16.0/Lunit;	// Mass within 16kpc
  double mp = 1.67e-24;
  double k  = 1.381e-16;
  double T  = (mp*pow(spd*1e5, 2.0))/k;

  double r_disk = Rgas/Lunit;   // Radius of disk for particle realization

  cout << "Dens  = " << gasD << endl;
  cout << "Temp  = " << T    << endl;
  cout << "GasD  = " << gasD << endl;
  cout << "MGas  = " << gasD*M_PI*r_max*r_max*Munit << endl;
  cout << "10^8  = " << 1.0e8*365.25*24*3600/Tunit << endl;

				// Bar semi-major axis
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
  double coef = v_o*v_o/(2.0*M_PI*r_o);

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
  

  ofstream out("bulgemod.dat");
  for (int i=0; i<numr; i++) 
    out << setw(16) << rr[i]
	<< setw(16) << dd[i]
	<< setw(16) << mm[i]
	<< setw(16) << pp[i]
	<< endl;

  SphericalModelTable bulge(numr, 
			    &rr[0]-1, &dd[0]-1, &mm[0]-1, &pp[0]-1,
			    0, 0.0, 0, "Bulge model");

  cout << "Bulge created" << endl;

  //
  // Bar model
  //
  EllipsoidForce::BarType bartype = EllipsoidForce::ferrers;
  int num  = 400;
  int numT = 200;
  
  EllipsoidForce bar(a, b, b, Mbar, rmin, 5.0*a, 
		     bartype, false, nu, num, numr);

  cout << "Bar created" << endl;
  cout << "Table creation . . . " << flush;
  bar.MakeTable(numT, numT, numT);
  cout << "done" << endl;

  //
  // Disk model
  //
  vector<double> mparam;
  mparam.push_back(v_o);
  mparam.push_back(r_o);
  mparam.push_back(H);

  int Lmax = 32, Nmax = 16;
  int numR = 10000;
  int numt = 400, numg = 100;

  KuzminDisk disk;
  disk.Initialize(rmin, rmax, dlog, Nmax, Lmax, numR, numt, numg, mparam);
    
  cout << "Disk created" << endl;

  cout << "PST output begin . . . " << flush;

  ofstream out2("pst.dat");
  int NUMR = 10000;
  double Rmin = rmin*Lunit;
  double Rmax = 20.0;
  double R, dR = (Rmax - Rmin)/(NUMR - 1), r;
  double pot1, pot2, vel1, vel2, bfrc, fr, fz, dens, kap1, kap2;
  vector<double> x(3, 0.0), force(4);
  double Dunit = Munit/(Lunit*Lunit*Lunit)/1e9;

  for (int i=0; i<NUMR; i++) {

    R = Rmin + dR*i;
    r = R/Lunit;
				// Major
    x[0] = r;
    x[1] = 0.0;
    bar.TableEval(x, force);

    pot1 = 
      bulge.get_pot(r) + 
      force[0]   + 
      disk.potential_eval(r, 0.0, 0.0);
    
    bfrc = -bulge.get_mass(r)/(r*r);
    disk.force_eval(r, 0.0, 0.0, fr, fz);

    vel1 = sqrt(-(force[1]+bfrc+fr)*r);

    dens = 
      bulge.get_density(r) + 
      bar.getDens(x) + 
      disk.density_eval(r, 0.0, 0.0);

    x[0] = 0.0;
    x[1] = r;
    bar.TableEval(x, force);

    pot2 =
      bulge.get_pot(r) + 
      force[0]  + 
      disk.potential_eval(0.0, r, 0.0);

    vel2 = sqrt(-(force[2]+bfrc+fr)*r);

    x[0] = r;
    x[1] = 0.0;

    epi2(r, bulge, disk, bar, kap1, kap2);
    
    out2 << setw(16) << R                                              // 1
	 << setw(16) << pot1                           /(Vunit*Vunit)  // 2
	 << setw(16) << pot2                           /(Vunit*Vunit)  // 3
	 << setw(16) << vel1                           /Vunit          // 4
	 << setw(16) << vel2                           /Vunit          // 5
	 << setw(16) << dens                           *Dunit          // 6
	 << setw(16) << vel1/r                         /(Vunit*Lunit)  // 7
	 << setw(16) << vel2/r                         /(Vunit*Lunit)  // 8
	 << setw(16) << sqrt(kap1)                     /(Vunit*Lunit)  // 9
	 << setw(16) << sqrt(kap2)                     /(Vunit*Lunit)  // 10
	 << endl;
  }

  cout << "done" << endl;

  if (Number) {
    cout << "Particle output begin . . . " << flush;

    const int niatr = 0;
    const int ndatr = 3;
    ofstream nout(file.c_str());
    if (nout) {
      nout << setw(12) << Number << setw(12) << niatr << setw(12) << ndatr
	   << endl;

      random_gen.seed(SEED);
      std::uniform_real_distribution<> unit(0.0, 1.0);
      std::normal_distribution<> vel1(0.0, k*T/mp*Vunit*Vunit/1e10);

      double R, z, phi, bfrc, fr, fz, velc;
      double mass = M_PI*r_disk*r_disk*gasD/Number;

      for (int n=0; n<Number; n++) {
				// Pick a particle in the disk
	R = r_disk*sqrt(unit(random_gen));
	z = 1.0e-3*H/Lunit*atanh(2.0*(unit(random_gen) - 0.5));
	r = sqrt(R*R + z*z);
	phi = 2.0*M_PI*unit(random_gen);

				// Midplane force
	bfrc = -bulge.get_mass(R)/(R*R);
	disk.force_eval(R, 0.0, 0.0, fr, fz);
	velc = sqrt(-(bfrc+fr)*Mscale*R);

	nout << setw(16) << mass
	     << setw(16) << R*cos(phi)
	     << setw(16) << R*sin(phi)
	     << setw(16) << z
	     << setw(16) << -velc*sin(phi) + vel1(random_gen)
	     << setw(16) <<  velc*cos(phi) + vel1(random_gen)
	     << setw(16) << vel1(random_gen);
	for (int i=0; i<niatr; i++) nout << setw(10) << 0;
	for (int i=0; i<ndatr; i++) nout << setw(16) << 0.0;
	nout << endl;
      }
    } else {
      cerr << "Could not open <" << file << "> for gas particles"
	   << endl;
    }

    cout << "done" << endl;
  }

  // MPI_Finalize();

  return 0;
}

