// This is -*- C++ -*-

#ifndef _DiskHalo_H
#define _DiskHalo_H

#include <math.h>
#include <stdlib.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <stdexcept>

#include <ACG.h>
#include <Uniform.h>
#include <Normal.h>
#include <Random.h>

#include <exponential.h>
#include <massmodel.h>

#include <gaussQ.h>
#include <localmpi.h>
#include <Vector.h>

#include "AddDisk.h"
#include "SphericalSL.h"
#include "EmpCylSL.h"
#include "SParticle.H"

class DiskHaloException : public std::runtime_error
{
public:

  DiskHaloException(const string message, 
		    const string sourcefilename, const int sourcelinenumber) 
    : runtime_error(message) 
  {
    cerr << "Thrown from source file <"  << message << "> at line number <"
	 << sourcelinenumber << endl;
  }
};

class DiskHalo
{
 private:
  AddDisk *newmod;
  SphericalModelTable *halo, *halo2, *halo3;
  SphericalModelMulti *multi;

  ExponentialDisk *disk;
  double scalelength, scaleheight, dmass;
  double center_pos[3], center_vel[3];

  SphericalSLptr expandh;
  EmpCylSLptr    expandd;

  Matrix *disktableP, *disktableN, epitable, dv2table, asytable;
  double dP, dR, dZ;

  Matrix halotable;
  double dr, dc;

  Uniform *rndU;
  Normal *rndN;
  ACG *gen;

  bool DF;
  bool MULTI;
  bool com;
  bool cov;

  double disk_density(double R, double z);
  double disk_surface_density(double R);
  void write_record(ostream &out, SParticle &p);
  void write_record(ostream &out, Particle &p)
  {
    SParticle P(p);
    write_record(out, P);
  }
  void table_halo_disp();

  // For frequency computation
  //
				// Number of mass and number bins, log spaced
  static const unsigned nh = 200;
				// Minimum number on grid before computing
  static const unsigned mh = 5;
				// Boundaries and delta
  double hDmin, hDmax, dRh;
  vector<double>   nhD;		// Mass in bins
  vector<unsigned> nhN;		// Number in bins
  int nzero;			// Computed minimum mesh point

  // Use Spline integration rather than Trapezoidal rule
  //
  static const int use_spline = 0;

 public:
  enum DiskGenType { Jeans, Asymmetric, Epicyclic };

  static int NDP;		// Number of knots in disk table phi grid
				// Default: 16

  static int NDZ;		// Number of knots in disk table z grid
				// Default: 200

  static int NDR;		// Number of knots in disk table R grid
				// Default: 40

  static int NHR;		// Number of knots in halo table r grid
				// Default: 40

  static int NHT;		// Number of knots in halo table theta grid
				// Default: 40
  
  static double RHMIN;		// Smallest radius in halo grid

  static double RHMAX;		// Largest radius in halo grid

  static double RDMIN;		// Smallest radius in disk grid

  static double RDMAX;		// Largest radius in disk grid & model

  static double COMPRESSION;	// Extra effective mass due to disk
				// Default: 1.0

  static double SHFACTOR;	// Minimum vertical integration size
				// in units of scale height

  static double TOLE;		// Energy offset for particle selection
				// Default: 0.003

  static double Q;		// Toomre Q

  static double R_DF;		// Change over points for DF and Jeans
  static double DR_DF;

  static int LOGSCALE;		// Log DF in SphericalModelTable

  static unsigned NBUF;		// Number of particles in MPI buffer
				// default: 8192

  static bool LOGR;		// Radial grid for Eddington inversion

  static bool CHEBY;		// Use Cheybshev fit for epicylic derivatives

  static int  NCHEB;		// Order for Cheybshev fit

  static int SEED;		// Initial seed

  static double RA;		// Anisotropy radius (default: 1e20)

  static int NUMDF;		// Number of DF grid points (default: 1200)

  static int RNUM;		// Number of model grid points for added
				// component (default: 1000)

  static unsigned VFLAG;	// Verbose debugging flags
				//
				// Bit		Action
				// ---		------
				//  1		Informational
				//  2		Model diagnostics
				//  4		Table diagnostics

				// ID string for debugging output
  static string RUNTAG;		// Default: "debug"

  //! Dummy Constructor
  DiskHalo();

  //! Constructor: equal mass
  DiskHalo(SphericalSLptr expandh, EmpCylSLptr expandd,
	   double H, double A, double DMass,
	   string& filename, int DF=0, int DIVERGE=0, double DIVERGE_RFAC=1.0,
	   DiskGenType type=Asymmetric);

  //! Constructor: multi mass
  DiskHalo(SphericalSLptr haloexp, EmpCylSLptr diskexp,
	   double H, double A, double DMass, 
	   string& filename1, int DIVERGE,  double DIVERGE_RFAC,
	   string& filename2, int DIVERGE2, double DIVERGE_RFAC2,
	   DiskGenType type=Asymmetric);

  //! Copy constructor
  DiskHalo(const DiskHalo &);

  //! Destructor
  ~DiskHalo();

  void set_halo(vector<Particle>& phalo, int nhalo, int npart);
  void set_halo_coordinates(vector<Particle>& phalo, int nhalo, int npart);
  void set_disk_coordinates(vector<Particle>& pdisk, int ndisk, int npart);

  void set_halo_table_single(vector<Particle>& phalo);
  void set_halo_table_multi (vector<Particle>& phalo);

  void set_pos_origin(double& x, double& y, double& z) 
    { 
      center_pos[0] = x;
      center_pos[1] = y;
      center_pos[2] = z;
    }

  void set_vel_origin(double& x, double& y, double& z) 
    { 
      center_vel[0] = x;
      center_vel[1] = y;
      center_vel[2] = z;
    }

  double get_hpot(double xp, double yp, double zp);

  double get_dpot(double xp, double yp, double zp);

  double get_hdpot(double xp, double yp, double zp);

  double get_ddpot(double xp, double yp, double zp);

  double deri_pot(double xp, double yp, double zp, int n);

  double deri_pot_disk(double xp, double yp, double zp, int n);

  double deri_pot_halo(double xp, double yp, double zp, int n);

  void disk_eval(double R, double z, double phi,
		 double &p, double &fr, double &fz, double &fp);

  double epi(double xp, double yp, double zp);

  void table_disk(vector<Particle>& part);

  double vz_disp2(double xp,double yp,double zp);

  double vr_disp2(double xp, double yp,double zp);

  double a_drift(double xp, double yp,double zp);

  double vp_disp2(double xp, double yp, double zp);

  double v_circ(double xp, double yp, double zp);

  void set_vel_disk(vector<Particle>& part);

  void table_halo(vector<Particle>& part);

  double get_disp(double xp, double yp,double zp);

  void set_vel_halo(vector<Particle>& part);

  void write_file(ostream &fou, vector<Particle>& part);

  void virial_ratio(vector<Particle>& hpart, vector<Particle>& dpart);

  void virial_ratio(const char *, const char *);

  void zero_com(bool val) { com = val; }

  void zero_cov(bool val) { cov = val; }

  void profile(ostream &rotfile, vector<Particle>& dpart,
	       double rmin, double rmax, int numr);

private:
  DiskGenType type;

};

#endif

