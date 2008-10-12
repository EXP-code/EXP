/*
  Generates a Monte Carlo realization of a halo with an embedded
  disk using Jeans' equations.

  Assumptions:

  1) Spherical halo supplied in an input table

  2) Axisymmetric (but 3-dimensional) exponential disk with sech^2(Z/z)
     vertical profile

  3) Halo as spherical velocity ellipsoid

  4) Disk as axisymmetric velocity ellipsoid in the plane (that is,
     $\sigma_r = \sigma_\phi$ and $\sigma_z$ determined by solving
     Jeans' equations in cylindrical coordinates.

 Loosely based on Enrico Vesperini's initial.cc and diskANDhalo.cc
 (mostly rewritten)

 Added the basis expansion of the disk: 12/10/01. KHB

 Rewritten and debugged by MDW between 12/28/01-12/31/01.  

        Added command line parsing.  "gendisk -h" will list parameters.

        I removed everything but the "disk and halo" case.  Removed
        multiple inheritance.  Changed interface; expansions registered
        with DiskHalo class on construction

        Switched from biortho classes to expansion classes from the
        EXP code.

        Uses a vector of particle structures rather than a Matrix to
        store an pass phase space.

        Rewrote the particle component foliation code using a more
        general algorithm.

        Solution to Jeans' equations are now computed in parallel and
        tabulated.  

        Particles are stored on local nodes and written to disk by
        master.  

        Removed lots of other cruft.

 More debugging 03/05 by MDW

        Repaired EmpCylSL scaling

	Added additional debugging output

	Compared against expsl routines

	Removed orphaned parameters

 Updated to include gas disk using local Euler solution 04/08 by MDW

 Both constant scale height and isothermal gas disks 08/08 by MDW

*/
                                // System libs
#include <unistd.h>
#include <getopt.h>
#include <values.h>

                                // C++/STL headers
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

                                // MDW classes
#include <numerical.h>
#include <gaussQ.h>
#include <isothermal.h>
#include <hernquist.h>
#include <model3d.h>
#include <biorth.h>
#include <SphericalSL.h>
#include <interp.h>
#include <EmpOrth9thd.h>

#include <norminv.H>

#define M_SQRT1_3 (0.5773502691896257645091487)

                                // For debugging
#ifdef DEBUG
#include <fpetrap.h>
#endif

                                // Local headers
#include "SphericalSL.h"
#include "DiskHalo.h" 
#include "localmpi.h"
void local_init_mpi(int argc, char **argv);

                                // Parameters
int LMAX=4;
int NMAX=10;
int NUMR=2000;
double RMIN=0.005;
double RCYLMIN=0.001;
double RCYLMAX=20.0;
double SCSPH=1.0;
double RSPHSL=47.5;
double ASCALE=1.0;
double HSCALE=0.1;
double ZMAX=10.0;
double TIME=0.0;
double DMFAC=1.0;

double X0=0.0;
double Y0=0.0;
double Z0=0.0;
double U0=0.0;
double V0=0.0;
double W0=0.0;

int NDR=800;
int NDZ=200;
int NHR=800;
int NHT=200;
double SHFAC=16.0;

int NMAX2=8;
int LMAX2=36;
int MMAX=4;
int NUMX=128;
int NUMY=64;
int NORDER=16;

int DIVERGE=0;
double DIVERGE_RFAC=1.0;

int DF=0;
double R_DF=20.0;
double DR_DF=5.0;

double scale_height = 0.1;
double scale_length = 2.0;
double halo_mass = 1.0;
double disk_mass = 1.0;
double gas_mass = 1.0;
double ToomreQ = 1.2;
double Tmin = 500.0;

bool const_height = true;
bool images = false;

int SEED = 11;

string hbods = "halo.bods";
string dbods = "disk.bods";
string gbods = "gas.bods";
string centerfile = "center.dat";
  
// Hydrogen fraction
//
const double f_H = 0.76;


// Global variables

#include <Particle.H>

int nthrds = 1;
int this_step = 0;
unsigned multistep = 0;
unsigned maxlev = 100;
int mstep = 1;
int Mstep = 1;
vector<int> stepL(1, 0), stepN(1, 1);
char threading_on = 0;
pthread_mutex_t mem_lock;
pthread_mutex_t coef_lock;
double tpos = 0.0;
double tnow = 0.0;
  
int 
main(int argc, char **argv)
{
  //
  // Inialize MPI stuff
  //
  local_init_mpi(argc, argv);
  
  //
  // Now parse the rest of the arguments
  //

  static string halofile = "SLGridSph.model";
  string suffix;
  bool basis = false;
  bool zero = false;
  int nhalo = 1000;             // Halo particles
  int ndisk = 1000;             // Disk particles
  int ngas  = 1000;             // Gas particles

  //========================= Parse command line ==============================

  int c;
  // int digit_optind = 0;

  while (1)
    {
      // int this_option_optind = optind ? optind : 1;
      int option_index = 0;
      static struct option long_options[] = {
	{"halo_mass", 1, 0, 0},
	{"rcylmin", 1, 0, 0},
	{"rcylmax", 1, 0, 0},
	{"rmin", 1, 0, 0},
	{"scsph", 1, 0, 0},
	{"sccyl", 1, 0, 0},
	{"ascale", 1, 0, 0},
	{"hscale", 1, 0, 0},
	{"numr", 1, 0, 0},
	{"norder", 1, 0, 0},
	{"cfile", 1, 0, 0},
	{"seed", 1, 0, 0},
	{"tmin", 1, 0, 0},
	{"constant", 1, 0, 0},
	{0, 0, 0, 0}
      };

      c = getopt_long (argc, argv, 
		       "I:D:G:L:M:X:N:n:f:Q:A:Z:m:g:r:R:1:2:s:S:t:d:c:T:bBzih",
		       long_options, &option_index);
      if (c == -1)
        break;

      string optname;

      switch (c)
        {
	case 0:			// Long options
	  optname = string(long_options[option_index].name);
	  if (!optname.compare("halo_mass")) halo_mass = atof(optarg);
	  if (!optname.compare("rcylmin"))   RCYLMIN = atof(optarg);
	  if (!optname.compare("rcylmax"))   RCYLMAX = atof(optarg);
	  if (!optname.compare("rmin"))      RMIN = atof(optarg);
	  if (!optname.compare("scsph"))     SCSPH = atof(optarg);
	  if (!optname.compare("ascale"))    ASCALE = atof(optarg);
	  if (!optname.compare("hscale"))    HSCALE = atof(optarg);
	  if (!optname.compare("numr"))      NUMR = atoi(optarg);
	  if (!optname.compare("norder"))    NORDER = atoi(optarg);
	  if (!optname.compare("seed"))      SEED = atoi(optarg);
	  if (!optname.compare("tmin"))      Tmin = atof(optarg);
	  if (!optname.compare("cfile"))     centerfile = string(optarg);
	  if (!optname.compare("constant"))  const_height = atoi(optarg)?true:false;
	  break;

        case 'I':
          nhalo = atoi(optarg);
          break;

        case 'D':
          ndisk = atoi(optarg);
          break;

        case 'G':
          ngas = atoi(optarg);
          break;

        case 'L':
          LMAX = atoi(optarg);
          break;

        case 'X':
          LMAX2 = atoi(optarg);
          break;

        case 'M':
          MMAX = atoi(optarg);
          break;

        case 'N':
          NMAX = atoi(optarg);
          break;

        case 'n':
          NMAX2 = atoi(optarg);
          break;

        case 'Q':
          ToomreQ = atof(optarg);
          break;

        case 'A':
          scale_length = atof(optarg);
          break;

        case 'Z':
          scale_height = atof(optarg);
          break;

        case 'm':
          disk_mass = atof(optarg);
          break;

        case 'g':
          gas_mass = atof(optarg);
          break;

        case 's':
          SCSPH = atof(optarg);
          break;

	case 'r':
	  RSPHSL = atof(optarg);
	  break;

        case 'R':
          DIVERGE = 1;
          DIVERGE_RFAC = atof(optarg);
          break;

        case '1':
	  DF = 1;
          R_DF= atof(optarg);
          break;

        case '2':
	  DF = 1;
          DR_DF = atof(optarg);
          break;

        case 'c':
          const_height = atoi(optarg) ? true : false;
          break;

        case 'T':
          Tmin = atof(optarg);
          break;

        case 'b':
          basis = true;
          break;

        case 'z':
          zero = true;
          break;

        case 't':
          suffix = optarg;
          break;

        case 'd':
          DMFAC = atof(optarg);
          break;

        case 'i':
          images = true;
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
               << "  -H num     number of halo particles (1000)\n"
               << "  -D num     number of disk particles (1000)\n"
               << "  -G num     number of gas particles (1000)\n"
               << "  -L lmax    spherical harmonic order (4)\n"
               << "  -M mmax    cylindrical harmonic order (4)\n"
               << "  -X lmax    maximum l for cylindrical expansion (26)\n"
               << "  -N nmax    spherical radial order (10)\n"
               << "  -n nmax2   cylindrical radial order (8)\n"
               << "  -Q value   Toomre Q for radial dispersion (1.2)\n"
               << "  -l a       scale length (2.0)\n"
               << "  -Z h       scale height (0.1)\n"
               << "  -m mass    disk mass (1.0)\n"
               << "  -r rsphsl  edge for SL expansion (47.5)\n"
               << "  -R expon   power law divergence exponent (unset)\n"
               << "  -s scale   halo coordinate scale\n"
               << "  -1 rmin    minimum radius for change over to DF\n"
               << "  -2 rmax    maximum radius for change over to DF\n"
	       << "  -c bool    assume constant gas disk scale height if true, isothermal with T=10000 if false\n"
               << "  -b         print out basis images (false)\n"
               << "  -z         zero center of mass and velocity (false)\n"
               << "  -t tag     suffix string appended to output files\n"
               << endl;
          exit(0);
        }
    }

  //===========================================================================

#ifdef DEBUG                    // For gdb . . . 
  sleep(20);
  set_fpu_handler();            // Make gdb trap FPU exceptions
#endif

  int n_particlesH, n_particlesD, n_particlesG;

  if (suffix.size()>0) {
    hbods = hbods + "." + suffix;
    dbods = dbods + "." + suffix;
    gbods = gbods + "." + suffix;
  }

  MPI_Bcast(&nhalo,  1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ndisk,  1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ngas,   1, MPI_INT, 0, MPI_COMM_WORLD);

                                // Divvy up the particles
  n_particlesH = nhalo/numprocs;
  if (myid==0) n_particlesH = nhalo - n_particlesH*(numprocs-1);
  
  n_particlesD = ndisk/numprocs;
  if (myid==0) n_particlesD = ndisk - n_particlesD*(numprocs-1);

  n_particlesG = ngas/numprocs;
  if (myid==0) n_particlesG = ngas  - n_particlesG*(numprocs-1);


#ifdef DEBUG  
  cout << "Processor " << myid << ": n_particlesH=" << n_particlesH << "\n";
  cout << "Processor " << myid << ": n_particlesD=" << n_particlesD << "\n";
  cout << "Processor " << myid << ": n_particlesG=" << n_particlesG << "\n";
#endif

  if (n_particlesH + n_particlesD + n_particlesG <= 0) {
    if (myid==0) cout << "You have specified zero particles!\n";
    MPI_Abort(MPI_COMM_WORLD, 3);
    exit(0);
  }

                                // Vectors to contain phase space
                                // Particle structure is defined in
                                // Particle.h
  vector<Particle> dparticles, hparticles;
  

                                // Disk halo grid parameters
  DiskHalo::RDMIN = RCYLMIN*scale_length;
  DiskHalo::RHMIN = RMIN;
  DiskHalo::RHMAX = RSPHSL;
  DiskHalo::RDMAX = RCYLMAX*scale_length;
  DiskHalo::NDR = NDR;
  DiskHalo::NDZ = NDZ;
  DiskHalo::NHR = NHR;
  DiskHalo::NHT = NHT;
  DiskHalo::SHFACTOR = SHFAC;
  DiskHalo::DMFACTOR = DMFAC;
  DiskHalo::Q = ToomreQ;        // Toomre Q
  DiskHalo::R_DF = R_DF;
  DiskHalo::DR_DF = DR_DF;
  DiskHalo::SEED = SEED;

  AddDisk::use_mpi = true;
  AddDisk::Rmin = RMIN;

  //===========================Spherical expansion=============================

  // SLGridSph::diverge = DIVERGE;
  // SLGridSph::divergexp = DIVERGE_RFAC;
    
  SphericalSL::RMIN = RMIN;
  SphericalSL::RMAX = RSPHSL;
  SphericalSL::NUMR = NUMR;
                                // Create expansion only if needed . . .
  SphericalSL *expandh = NULL;
  if (n_particlesH) {
    expandh = new SphericalSL(LMAX, NMAX, SCSPH);
#ifdef DEBUG
    string dumpname("debug");
    expandh->dump_basis(dumpname);
#endif
  }

  //===========================Cylindrical expansion===========================


  EmpCylSL::RMIN = RCYLMIN;
  EmpCylSL::RMAX = RCYLMAX;
  EmpCylSL::NUMX = NUMX;
  EmpCylSL::NUMY = NUMY;
  EmpCylSL::CMAP = true;
  EmpCylSL::logarithmic = true;

  if (basis)
    EmpCylSL::DENS = true;
  else
    EmpCylSL::DENS = false;

                                // Create expansion only if needed . . .
  EmpCylSL* expandd = NULL;
  if (n_particlesD) 
    expandd = new EmpCylSL(NMAX2, LMAX2, MMAX, NORDER, ASCALE, HSCALE);
  cout << "Proccess " << myid << ": "
       << " rmin=" << EmpCylSL::RMIN
       << " rmax=" << EmpCylSL::RMAX
       << " a=" << ASCALE
       << " h=" << HSCALE
       << " nmax2=" << NMAX2
       << " lmax2=" << LMAX2
       << " mmax=" << MMAX
       << " nordz=" << NORDER
       << endl << flush;

  //====================Create the disk & halo model===========================


  DiskHalo diskhalo(expandh, expandd,
                    scale_height, scale_length, 
		    halo_mass, disk_mass, halofile,
                    DF, DIVERGE, DIVERGE_RFAC);

  ifstream center(centerfile.c_str());
  if (center) {

    bool ok = true;

    center >> X0;
    if (center.fail()) ok = false;

    center >> Y0;
    if (center.fail()) ok = false;

    center >> Z0;
    if (center.fail()) ok = false;

    if (ok) {
      diskhalo.set_pos_origin(X0, Y0, Z0);
      if (myid==0) cout << "Using position origin: " 
			<< X0 << ", " << Y0 << ", " << Z0 << endl;
    }

    center >> U0;
    if (center.fail()) ok = false;

    center >> V0;
    if (center.fail()) ok = false;

    center >> W0;
    if (center.fail()) ok = false;

    if (ok) {
      diskhalo.set_vel_origin(U0, V0, W0);
      if (myid==0) cout << "Using velocity origin: " 
			<< U0 << ", " << V0 << ", " << W0 << endl;
    }
  }

                                // Make zero center of mass and
                                // center of velocity
  diskhalo.zero_com_cov(zero);
  
  //===========================================================================

                                // Open output file (make sure it exists
                                // before realizing a large phase space)
  ofstream out_halo, out_disk;
  if (myid==0) {
    out_halo.open(hbods.c_str());
    if (!out_halo) {
      cout << "Could not open <" << hbods << "> for output\n";
      MPI_Abort(MPI_COMM_WORLD, 4);
      exit(0);
    }

    out_disk.open(dbods.c_str());
    if (!out_disk) {
      cout << "Could not open <" << dbods << "> for output\n";
      MPI_Abort(MPI_COMM_WORLD, 4);
      exit(0);
    }
  }

  //=================Make the phase space coordinates==========================

  if (n_particlesH) {
    if (myid==0) cout << "Generating halo coordinates . . . " << flush;
    diskhalo.set_halo_coordinates(hparticles, nhalo, n_particlesH);
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done\n";
  }

  if (n_particlesD) {
    if (myid==0) cout << "Generating disk coordinates . . . " << flush;
    diskhalo.set_disk_coordinates(dparticles, ndisk, n_particlesD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done\n";
  }

  if (n_particlesH) {
    if (myid==0) cout << "Beginning halo accumulation . . . " << flush;
    expandh->accumulate(hparticles);
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done\n";
  }
  
  if (n_particlesD) {
    if (myid==0) cout << "Beginning disk accumulation . . . " << flush;
    expandd->setup_eof();
    expandd->setup_accumulation();
    expandd->accumulate_eof(dparticles);
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done\n";
  
    if (myid==0) cout << "Making the EOF . . . " << flush;
    expandd->make_eof();
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done\n";
  
    if (myid==0) cout << "Making disk coefficients . . . " << flush;
    expandd->make_coefficients();
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done\n";

    if (myid==0) cout << "Reexpand . . . " << flush;
    expandd->accumulate(dparticles);
    expandd->make_coefficients();
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << "done\n";

    if (images && myid==0) {
      cout << "Images . . . " << flush;
      ostringstream dumpname;
      dumpname << "images.0";
      expandd->dump_images(dumpname.str(), 5.0*ASCALE, 5.0*HSCALE, 64, 64, true);
      cout << "done\n";
    }
  }
  

  //===========================Diagnostics=====================================

                                // For examining the coverage, etc.
                                // Images can be contoured in SM using
                                // the "ch" file type
  if (myid==0 && basis) {
    
    cout << "Dumping basis images . . . " << flush;
    
    if (n_particlesD) {
      int nout = 200;
      char dumpname[] = "basis.dump";
      expandd->dump_basis(dumpname, 0);
      string prefix = "gendisk2";
      expandd->dump_images(prefix, 5.0*scale_length, 5.0*scale_height,
			   nout, nout, false);
      expandd->dump_images_basis(prefix, 5.0*scale_length, 5.0*scale_height,
				 nout, nout, false, 0, MMAX, 0, NORDER-1);
    }


    if (n_particlesH) {
      string extn("test");
      expandh->dump_basis(extn);
    }
    
    if (n_particlesH) {
      
      const int nstr = 5;
      const char *names[nstr] = {".dens", ".potl", ".potr", ".pott", ".potp"};
      ofstream *out = new ofstream [nstr];
      
      int nout = 200;
      double rmax = 6.0*scale_length;
      double x, y, dr = 2.0*rmax/(nout-1);
      float f;
    
      for (int i=0; i<nstr; i++) {
        string name("halo");
        name += names[i];
        out[i].open(name.c_str());
        
        out[i].write((char *)&nout, sizeof(int));
        out[i].write((char *)&nout, sizeof(int));
        out[i].write((char *)&(f=-rmax), sizeof(float));
        out[i].write((char *)&(f= rmax), sizeof(float));
        out[i].write((char *)&(f=-rmax), sizeof(float));
        out[i].write((char *)&(f= rmax), sizeof(float));
      }
      
      double r, theta, phi;
      double dens, potl, potr, pott, potp;
    
      for (int j=0; j<nout; j++) {
        y = -rmax + dr*j;
      
        for (int i=0; i<nout; i++) {
          x = -rmax + dr*i;
        
          r = sqrt(x*x + y*y);
          theta = 0.5*M_PI;
          phi = atan2(y, x);
        
          expandh->determine_fields_at_point(r, theta, phi,
                                             &dens, &potl, 
                                             &potr, &pott, &potp);
        
          out[0].write((char *)&(f=dens), sizeof(float));
          out[1].write((char *)&(f=potl), sizeof(float));
          out[2].write((char *)&(f=potr), sizeof(float));
          out[3].write((char *)&(f=pott), sizeof(float));
          out[4].write((char *)&(f=potp), sizeof(float));
        }
        
      }
    
      for (int i=0; i<nstr; i++) out[i].close();
      delete [] out;
    }

    if (n_particlesD) {

      const int nstr = 5;
      const char *names[nstr] = {".dens", ".pot", ".fr", ".fz", ".fp"};
      ofstream *out = new ofstream [nstr];
    
      int nout = 200;
      double rmax = DiskHalo::RDMAX;
      double x, y, dr = 2.0*rmax/(nout-1);
      float f;
    
      for (int i=0; i<nstr; i++) {
        string name("disk");
        name += names[i];
        out[i].open(name.c_str());
        
        out[i].write((char *)&nout, sizeof(int));
        out[i].write((char *)&nout, sizeof(int));
        out[i].write((char *)&(f=-rmax), sizeof(float));
        out[i].write((char *)&(f= rmax), sizeof(float));
        out[i].write((char *)&(f=-rmax), sizeof(float));
        out[i].write((char *)&(f= rmax), sizeof(float));
      }
    
      double z = 0.0, d0, p0, d, p, fr, fz, fp;
    
      for (int j=0; j<nout; j++) {
        y = -rmax + dr*j;
      
        for (int i=0; i<nout; i++) {
          x = -rmax + dr*i;
        
	  if (x<0.0)
	    expandd->accumulated_eval(fabs(x), y, M_PI, p0, p, fr, fz, fp);
	  else
	    expandd->accumulated_eval(x, y,  0.0, p0, p, fr, fz, fp);


          d = expandd->accumulated_dens_eval(sqrt(x*x + y*y), z, atan2(y, x), d0);
        
          
          out[0].write((char *)&(f=d ), sizeof(float));
          out[1].write((char *)&(f=p ), sizeof(float));
          out[2].write((char *)&(f=fr), sizeof(float));
          out[3].write((char *)&(f=fz), sizeof(float));
          out[4].write((char *)&(f=fp), sizeof(float));
        }
        
      }
    
      for (int i=0; i<nstr; i++) out[i].close();
      delete [] out;
    }
    
    cout << "done\n";
  }
  
  MPI_Barrier(MPI_COMM_WORLD);

  //====================Make the phase space velocities========================

  if (myid==0) cout << "Generating halo velocities . . . " << flush;
  diskhalo.set_vel_halo(hparticles);
  if (myid==0) cout << "done\n";
  
  if (myid==0) cout << "Generating disk velocities . . . " << flush;
  diskhalo.set_vel_disk(dparticles);
  if (myid==0) cout << "done\n";
  

  //====================All done: write it out=================================

  if (myid==0) cout << "Writing phase space file . . . " << flush;

  diskhalo.write_file(out_halo, out_disk, hparticles, dparticles);
  if (myid==0) cout << "done\n";

  out_halo.close();
  out_disk.close();
                                // Diagnostic . . .
  diskhalo.virial_ratio(hparticles, dparticles);

  //====================Compute gas particles==================================

  if (myid==0 && n_particlesG) {
    cout << "Computing gas particles . . . " << endl;

				// UNITS
				// -------------------

				// cm
    const double pc = 3.08568025e18;
				// proton mass
    const double m_p = 1.67262158e-24;
				// g
    const double msun = 1.98892e33; //
				// cgs
    const double G = 6.67300e-08;
				// cgs
    const double boltz = 1.3806503e-16;

    double T = 10000;

    
    double Lunit = 3.0e5*pc;	// Virial radius
    double Munit = 1.0e12*msun;	// Virial mass
    double Tunit = sqrt(Lunit*Lunit*Lunit/(Munit*G));
    double Vunit = Lunit/Tunit;

    // Fac = kT*R_vir/(G*m_p*M_vir)
    // where M_vir = 1e12 Msun, R_vir=300 kpc
    //
    double fac = T/(G*m_p*Munit/(Lunit*boltz));

    // Thermal velocity in system units
    //
    double mm   = f_H*m_p + (1.0-f_H)*4.0*m_p;
    double vthermal = sqrt( (boltz*T)/mm ) / Vunit;
    double vmin2 = (boltz*Tmin/mm) / (Vunit*Vunit);

    // Compute using Jeans theorem
    //
    double rmin = RMIN;
    double rmax = 10.0*scale_length;
    double zmin = 0.001*scale_height;
    int nrint = 200;
    int nzint = 400;
    vector< vector<double> > zrho, zmas, vcir;
    double r, R, dR = (rmax - rmin)/(nrint-1);
    double z, dz = (log(rmax) - log(zmin))/(nzint-1);

    double p0, p, fr, fz, fp, dens, potl, potr, pott, potp;

    if (const_height) {

      for (int i=0; i<nrint; i++) {
	R = rmin + dR*i;

	vector<double> lrho(nzint), trho(nzint), tcir(nzint), tmas(nzint, 0);

	for (int j=0; j<nzint; j++) {
	  z = zmin*exp(dz*j);
	  r = sqrt(R*R + z*z);
	  
	  double pot=0.0, frt0=0.0, fzt0=0.0;
	  if (expandd) {
	    expandd->accumulated_eval(R, z, 0, p0, p, fr, fz, fp);
	    frt0 += -fr;
	    fzt0 += -fz;
	    pot += p0;
	  }
	  if (expandh) {
	    expandh->determine_fields_at_point(r, acos(z/(r+1.0e-8)), 0.0,
					       &dens, &potl, 
					       &potr, &pott, &potp);
	    
	    frt0 += potr;
	    fzt0 += (potr*z + pott*R*R/(r*r))/r;
	    pot += potl;
	  }
	  
	  trho[j] = fzt0*scale_height;
	  tcir[j] = sqrt(max<double>(R*frt0-R*trho[j]/scale_length, 0.0));
	}
	
	for (int j=0; j<nzint; j++) 
	  tmas[j] = 1.0 - exp(-zmin*exp(dz*j)/scale_height);
	
	zrho.push_back(trho);
	zmas.push_back(tmas);
	vcir.push_back(tcir);
      }

      //
      // Vertical table
      //
      string ztable("ztable.dat");
      cout << "Writing " << setw(15) << right << ztable
	   << " [gas] . . . " << flush;
      ofstream ztest(ztable.c_str());
      for (int i=0; i<nrint; i++) {
	for (int j=0; j<nzint; j++) {
	  ztest << setw(15) << rmin + dR*i
		<< setw(15) << zmin*exp(dz*j)
		<< setw(15) << zrho[i][j]
		<< setw(15) << zrho[i][j]*Vunit*Vunit*mm/boltz
		<< setw(15) << zmas[i][j]
		<< setw(15) << vcir[i][j]
		<< endl;
	}
	ztest << endl;
      }
      ztest.close();
      cout << "done" << endl;
      
    } else {

      for (int i=0; i<nrint; i++) {
	R = rmin + dR*i;



	vector<double> lrho(nzint), trho(nzint), tcir(nzint), tmas(nzint, 0);

	for (int j=0; j<nzint; j++) {
	  z = zmin*exp(dz*j);
	  r = sqrt(R*R + z*z);
	  
	  double frt0=0.0, fzt0=0.0;
	  if (expandd) {
	    expandd->accumulated_eval(R, z, 0, p0, p, fr, fz, fp);
	    frt0 += -fr;
	    fzt0 += -fz;
	  }
	  if (expandh) {
	    expandh->determine_fields_at_point(r, acos(z/(r+1.0e-8)), 0.0,
					       &dens, &potl, 
					       &potr, &pott, &potp);
	    frt0 += potr;
	    fzt0 += (potr*z + pott*R*R/(r*r))/r;
	  }
	  
	  trho[j] = -fzt0/(vthermal*vthermal);
	  tcir[j] = sqrt(max<double>(R*frt0-R*vthermal*vthermal/scale_length, 0.0));
	}
	
	double mass = 0.0;
	double zfac = 1.0 - exp(-dz);
				    
	lrho[0] = 0.0;
	for (int j=1; j<nzint; j++) 
	  lrho[j] = lrho[j-1] + 0.5*(trho[j-1] + trho[j]) * zmin*exp(dz*j)*zfac;
	
	for (int j=1; j<nzint; j++) 
	  tmas[j] = tmas[j-1] + 0.5*(exp(lrho[j-1]) + exp(lrho[j])) * zmin*exp(dz*j)*zfac;
	
	for (int j=0; j<nzint; j++) {
	  if (tmas[nzint-1]>0.0 && !isnan(tmas[nzint-1])) {
	    trho[j]  = exp(lrho[j])/tmas[nzint-1];
	    tmas[j] /= tmas[nzint-1];
	  } else {
	    trho[j] = 0.0;
	    if (j==0) tmas[j] = 0.0;
	    else      tmas[j] = 1.0;
	  }
	}
	zrho.push_back(trho);
	zmas.push_back(tmas);
	vcir.push_back(tcir);
      }


      //
      // Vertical table
      //
      cout << "Writing ztable.dat [gas] . . . " << flush;
      ofstream ztest("ztable.dat");
      for (int i=0; i<nrint; i++) {
	for (int j=0; j<nzint; j++) {
	  ztest << setw(15) << rmin + dR*i
		<< setw(15) << zmin*exp(dz*j)
		<< setw(15) << zrho[i][j]
		<< setw(15) << zmas[i][j]
		<< setw(15) << vcir[i][j]
		<< endl;
	}
	ztest << endl;
      }
      ztest.close();
      cout << "done" << endl;
      
    }

    // 
    // Prepare output stream
    //
    ofstream outps("gas.bods");
    if (!outps) {
      cerr << "Couldn't open <" << "gas.bods" << "> for output\n";
      exit (-1);
    }

    const int ITMAX=1000;
    const int NREPORT=1000;
    const int nparam = 3;
    
    //
    // Maximum enclosed disk mass given rmax
    //
    double rmx2 = 1.5*rmax;
    double mmx2 = 1.0 - (1.0 + rmx2/scale_length)*exp(-rmx2/scale_length);
    double mmax = 1.0 - (1.0 + rmax/scale_length)*exp(-rmax/scale_length);

    //
    // Random generators
    //
    ACG gen(10, 20);
    Uniform unit(0.0, 1.0, &gen);

    //
    // Trimmed Gaussian
    //
    double minK=0.0, maxK=1.0, sigma = 3.0;
    if (sigma>0) {
      minK = 0.5*(1.0+erf(-0.5*sigma));
      maxK = 0.5*(1.0+erf( 0.5*sigma));
    }
    Uniform unitN(minK, maxK, &gen);


    double gmass = gas_mass/ngas;
    double KE=0.0, VC=0.0;
    vector<double> mc2(nzint);

    fr = fz = potr = 0.0;

    outps << setw(8) << ngas
	  << setw(6) << 0 << setw(6) << nparam << endl;

    for (int n=0; n<ngas; n++) {

      double F, dF, M=mmax*unit(), Z=unit();
      double R = M*rmax, phi=2.0*M_PI*unit(), x, y, z, rr, vc;
      double ax, ay, az;

				// Narrow with bisection
      double rm = 0.0, rp = rmx2;
      double fm = -M, fp = mmx2 - M;
      for (int j=0; j<15; j++) {
	R = 0.5*(rm + rp);
	F = 1.0 - M - (1.0 + R/scale_length)*exp(-R/scale_length);
	if (fm*F<0.0) {
	  rp = R;
	  fp = F;
	} else {
	  rm = R;
	  fm = F;
	}
      }
				// Polish with Newton-Raphson
      for (int j=0; j<ITMAX; j++) {
	F = 1.0 - M - (1.0 + R/scale_length)*exp(-R/scale_length);
	dF = R/(scale_length*scale_length)*exp(-R/scale_length);
	R += -F/dF;
	if (fabs(F/dF)<1.0e-12) break;
      }
    
      int indr = static_cast<int>(floor(R/dR));
      if (indr<0) indr=0;
      if (indr>nrint-2) indr=nrint-2;
      double a = (dR*(indr+1) - R)/dR;
      double b = (R - indr*dR)/dR;

      vector<double> mz(nzint), vz(nzint);
      for (int j=0; j<nzint; j++) {
	mz[j] = a*zmas[indr][j] + b*zmas[indr+1][j];
	vz[j] = a*vcir[indr][j] + b*vcir[indr+1][j];
      }
      for (int j=0; j<nzint; j++) mz[j] /= mz[nzint-1];
      
      if (const_height) {
	for (int j=0; j<nzint; j++) 
	  mc2[j] = max<double>(a*zrho[indr][j] + b*zrho[indr+1][j], vmin2);
      }

      int indz = max<int>(0, min<int>(nzint-2, Vlocate(Z, mz)));

      a = (mz[indz+1] - Z)/(mz[indz+1] - mz[indz]);
      b = (Z - mz[indz  ])/(mz[indz+1] - mz[indz]);

      vc = fabs(a*vz[indr] + b*vz[indr+1]);

      z = zmin*exp(dz*(a*indz + b*(indz+1)));
      if (unit()<0.5) z *= -1.0;
      rr = sqrt(R*R + z*z);

      if (const_height) {
	vthermal = a*mc2[indz] + b*mc2[indz+1];
	vthermal = sqrt(max<double>(vmin2, vthermal));
      }

      double sinp = sin(phi), cosp = cos(phi);
      x = R*cosp;
      y = R*sinp;

      double u = -vc*sinp + vthermal*norminv(unitN());
      double v =  vc*cosp + vthermal*norminv(unitN());
      double w =  vthermal*norminv(unitN());
      
      outps << setw(18) << gmass
	    << setw(18) << R*cos(phi)
	    << setw(18) << R*sin(phi)
	    << setw(18) << z
	    << setw(18) << u
	    << setw(18) << v
	    << setw(18) << w;
      for (int k=0; k<nparam; k++) outps << setw(18) << 0.0;
      outps << endl;
    
      if (expandd)
	expandd->accumulated_eval(R, z, phi, p0, p, fr, fz, fp);

      if (expandh)
	expandh->determine_fields_at_point(rr, acos(z/(rr+1.0e-8)), 0.0,
					   &dens, &potl, 
					   &potr, &pott, &potp);
      KE += 0.5*gmass*(u*u + v*v + w*w);

      VC += gmass*(-rr*potr + R*fr + z*fz);

      if (!((n+1)%NREPORT)) cout << "\r." << n+1 << flush;
    }

    cout << endl << "Done!" << endl;

    cout << "****************************" << endl
	 << "  Gas disk" << endl
	 << "----------------------------" << endl
	 << "  KE       = " << KE << endl
	 << "  VC       = " << VC << endl;
    if (VC<0.0)
      cout << " -2T/W     = " << -2.0*KE/VC << endl;
    cout << "****************************" << endl;
  }

  //===========================================================================

  MPI_Barrier(MPI_COMM_WORLD);

  delete expandh;
  delete expandd;

  MPI_Finalize();

  return 0;
}

