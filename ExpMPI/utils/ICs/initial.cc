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
*/

                                // System libs
#include <math.h>
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <values.h>

                                // C++/STL headers
#include <iostream>
#include <iomanip>
#include <fstream>
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
double RDISK=10.0;
double ZMAX=10.0;
double TIME=0.0;

double X0=0.0;
double Y0=0.0;
double Z0=0.0;
double U0=0.0;
double V0=0.0;
double W0=0.0;

int NDR=200;
int NDZ=200;
double SHFACTOR=16.0;

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
double disk_mass = 1.0;
double ToomreQ = 1.2;

int SEED = 11;

string centerfile = "center.dat";
  
// Global variables

#include <Particle.H>

int nthrds = 1;
int this_step = 0;
int multistep = 0;
int mstep = 1;
int Mstep = 1;
vector<int> stepL(1, 0), stepN(1, 1);
char threading_on = 0;
pthread_mutex_t mem_lock;
double tpos = 0.0;
double tnow = 0.0;
  
int 
main(int argc, char **argv)
{
  static string halofile = "SLGridSph.model";
  bool basis = false;
  bool zero = false;
  int nhalo = 1000;             // Halo particles
  int ndisk = 1000;             // Disk particles

  //========================= Parse command line ==============================

  int c;
  // int digit_optind = 0;

  while (1)
    {
      // int this_option_optind = optind ? optind : 1;
      int option_index = 0;
      static struct option long_options[] = {
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
	{0, 0, 0, 0}
      };

      c = getopt_long (argc, argv, 
		       "H:D:L:M:X:N:n:f:Q:A:Z:m:r:R:1:2:s:S:bBzh",
		       long_options, &option_index);
      if (c == -1)
        break;

      string optname;

      switch (c)
        {
	case 0:			// Long options
	  optname = string(long_options[option_index].name);
	  if (!optname.compare("rcylmin")) RCYLMIN = atof(optarg);
	  if (!optname.compare("rcylmax")) RCYLMAX = atof(optarg);
	  if (!optname.compare("rmin"))    RMIN = atof(optarg);
	  if (!optname.compare("scsph"))   SCSPH = atof(optarg);
	  if (!optname.compare("ascale"))  ASCALE = atof(optarg);
	  if (!optname.compare("hscale"))  HSCALE = atof(optarg);
	  if (!optname.compare("numr"))    NUMR = atoi(optarg);
	  if (!optname.compare("norder"))  NORDER = atoi(optarg);
	  if (!optname.compare("seed"))    SEED = atoi(optarg);
	  if (!optname.compare("cfile"))   centerfile = string(optarg);
	  break;

        case 'H':
          nhalo = atoi(optarg);
          break;

        case 'D':
          ndisk = atoi(optarg);
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

        case 'b':
          basis = true;
          break;

        case 'z':
          zero = true;
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
               << "  -b         print out basis images (false)\n"
               << "  -z         zero center of mass and velocity (false)\n"
               << endl;
          exit(0);
        }
    }

  //===========================================================================

#ifdef DEBUG                    // For gdb . . . 
  sleep(20);
  set_fpu_handler();            // Make gdb trap FPU exceptions
#endif

  local_init_mpi(argc, argv);   // Inialize MPI stuff

  int n_particlesH, n_particlesD;

  MPI_Bcast(&nhalo,  1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ndisk,  1, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);

                                // Divvy up the particles
  n_particlesH = nhalo/numprocs;
  if (myid==0) n_particlesH = nhalo - n_particlesH*(numprocs-1);
  
  n_particlesD = ndisk/numprocs;
  if (myid==0) n_particlesD = ndisk - n_particlesD*(numprocs-1);


#ifdef DEBUG  
  cout << "Processor " << myid << ": n_particlesH=" << n_particlesH << "\n";
  cout << "Processor " << myid << ": n_particlesD=" << n_particlesD << "\n";
#endif

  if (n_particlesH + n_particlesD <= 0) {
    if (myid==0) cout << "You have specified zero particles!\n";
    MPI_Abort(MPI_COMM_WORLD, 3);
    exit(0);
  }

                                // Vectors to contain phase space
                                // Particle structure is defined in
                                // Particle.h
  vector<Particle> dparticles, hparticles;
  

                                // Disk halo grid parameters
  DiskHalo::RDMIN = RMIN;
  DiskHalo::RHMIN = RMIN;
  DiskHalo::RHMAX = RSPHSL;
  DiskHalo::RDMAX = RDISK*scale_length;
  DiskHalo::NDR = NDR;
  DiskHalo::NDZ = NDZ;
  DiskHalo::SHFACTOR = SHFACTOR;
  DiskHalo::Q = ToomreQ;        // Toomre Q
  DiskHalo::R_DF = R_DF;
  DiskHalo::DR_DF = DR_DF;
  DiskHalo::SEED = SEED;

  AddDisk::use_mpi = true;
  AddDisk::Rmin = RMIN;

  //===========================Spherical expansion=============================

  // SLGridSph::diverge = DIVERGE;
  // SLGridSph::divergexp = DIVERGE_RFAC;
  // SLGridSph::cmap = 1;
    
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
                    scale_height, scale_length, disk_mass, halofile,
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
    out_halo.open("out.halo");
    if (!out_halo) {
      cout << "Could not open <out.halo> for output\n";
      MPI_Abort(MPI_COMM_WORLD, 4);
      exit(0);
    }

    out_disk.open("out.disk");
    if (!out_disk) {
      cout << "Could not open <out.disk> for output\n";
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
    expandd->accumulate_eof(dparticles);
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
  }
  

  //===========================Diagnostics=====================================

                                // For examining the coverage, etc.
                                // Images can be contoured in SM using
                                // the "ch" file type
  if (myid==0 && basis) {
    
    cout << "Dumping basis images . . . " << flush;
    
    if (n_particlesD) {
      char dumpname[] = "basis.dump";
      expandd->dump_basis(dumpname, 0);
    }


    if (n_particlesH) {
      string extn("test");
      expandh->dump_basis(extn);
    }
    
    if (n_particlesH) {
      
      const int nstr = 5;
      char *names[nstr] = {".dens", ".potl", ".potr", ".pott", ".potp"};
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
      char *names[nstr] = {".dens", ".pot", ".fr", ".fz", ".fp"};
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

  //===========================================================================

  delete expandh;
  delete expandd;

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  return 0;
}

