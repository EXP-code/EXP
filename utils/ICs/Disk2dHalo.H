#ifndef _Disk2dHalo_H
#define _Disk2dHalo_H

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <string>
#include <random>
#include <cmath>

#include "exponential.H"
#include "massmodel.H"

#include "gaussQ.H"
#include "localmpi.H"

#include "AddDisk.H"
#include "SphericalSL.H"
#include "Disk2d.H"
#include "SParticle.H"
#include "libvars.H"

class Disk2dHaloException : public std::runtime_error
{
public:

  Disk2dHaloException(const string message, 
		      const string sourcefilename, const int sourcelinenumber) 
    : runtime_error(message) 
  {
    cerr << "Thrown from source file <"  << message << "> at line number <"
	 << sourcelinenumber << endl;
  }
};

class Disk2dHalo
{
 private:
  std::shared_ptr<AddDisk> newmod;
  std::shared_ptr<SphericalModelTable> halo, halo2, halo3;
  std::shared_ptr<SphericalModelMulti> multi;

  std::shared_ptr<ExponentialDisk> disk;
  double scalelength, dmass;
  double center_pos[3], center_vel[3];

  SphericalSLptr expandh;
  std::shared_ptr<Disk2d> expandd;

  Eigen::MatrixXd epitable, dv2table, asytable;
  double dP, dR, dZ, sigma0;

  Eigen::MatrixXd halotable;
  double dr, dc;

  std::mt19937 gen;
  std::uniform_real_distribution<> rndU;
  std::normal_distribution<> rndN;

  bool DF;
  bool MULTI;
  bool com;
  bool cov;

  const int bunchcnt = 16384;
  int bufcnt;
  std::ostringstream bufout;

  double disk_surface_density(double R);
  void write_record(ostream &out, SParticle &p);
  void write_record(ostream &out, Particle &p)
  {
    SParticle P(p);
    write_record(out, P);
  }
  void flush_buffer(ostream &out);
  void table_halo_disp();

  // For frequency computation
  //
				// Number of mass and number bins, log spaced
  static const unsigned nh = 50;
				// Minimum number on grid before computing
  static const unsigned mh = 5;
				// Boundaries and delta
  double hDmin, hDmax, dRh;
  vector<double>   nhD;		// Mass in bins
  vector<double>   nhM;		// Cumulative mass in bins
  vector<unsigned> nhN;		// Number in bins
  int nzero;			// Computed minimum mesh point
  int nzepi;			// Computed minimum epi value

  // Use Spline integration rather than Trapezoidal rule
  //
  static const int use_spline = 0;

  class CMap
  {
  public:
    virtual double x_to_r(double x) = 0;
    virtual double r_to_x(double r) = 0;
  };

  class CLin : public CMap
  {
  public:

    virtual double x_to_r(double x) { return x; }
    virtual double r_to_x(double r) { return r; }
  };

  class CLog : public CMap
  {
  public:

    virtual double x_to_r(double x) { return exp(x); }
    virtual double r_to_x(double r) { return log(r); }
  };

  class CRat : public CMap
  {
  public:
    double scale;

  public:

    CRat(double scale)  : scale(scale)   {}
    CRat(const CRat& p) : scale(p.scale) {}

    virtual double x_to_r(double x) { return scale*(x + 1.0)/(x - 1.0); }
    virtual double r_to_x(double r) { return (r + scale)/(r - scale);   }
  };

  std::shared_ptr<CMap> cmap;
  double Xmin, Xmax;
  std::string mtype;

  std::shared_ptr<CMap> make_cmap(std::string type, double scale)
  {
    std::for_each(type.begin(), type.end(), [](char & c){ c = ::tolower(c); });

    std::shared_ptr<CMap> cmap;

    if (type.find("linear") == 0)   cmap = std::make_shared<CLin>();
    else if (type.find("log") == 0) cmap = std::make_shared<CLog>();
    else if (type.find("rat") == 0) cmap = std::make_shared<CRat>(scale);
    else                            cmap = std::make_shared<CLog>();

    return cmap;
  }

  class GenPhi
  {
  private:

    double scale, baramp;
    unsigned numtot, numbar;
    int m;

    std::random_device dev;		   // Random device
    std::mt19937 rng;			   // Random generator
    std::uniform_real_distribution<> dist; // Random distribution

    static constexpr int numcnt = 32;
    std::vector<unsigned> counts;

  public:

    GenPhi(double scale, double amp, int m) :
      scale(scale), baramp(amp), numtot(0), numbar(0), m(m)
    {
      rng = std::mt19937(dev()); // Set generator
      dist = std::uniform_real_distribution<>(0.0, 1.0); // Set distribution
      counts.resize(numcnt, 0);
    }
    
    double Angle(double r)
    {
      numtot++;
      // Relative bar amplitude factor
      //
      double f = baramp*exp(-(r - scale)*(r - scale)/(2.0*scale*scale*0.09));

      // Do we have a bar particle? (Factor of 2 accounts for [-1, 1]
      // range of cos(m*phi) amplitude)
      //
      if (dist(rng) < f*2.0) {
	numbar++;
	// Compute quadrupole angle
	auto dp = [this](double phi)
	{ return 0.5*(1.0 + cos(phi*this->m)); };

	double phi = 2.0*M_PI*dist(rng);
	if (dp(phi) < dist(rng)) {
	  phi = 2.0*M_PI*dist(rng);
	}
	// Debug
	int cnt = floor(phi*numcnt/(2.0*M_PI));
	if (cnt>=0 and cnt<numcnt) counts[cnt]++;

	return phi;
      } else {
	// Random azimuthal angle
	return 2.0*M_PI*dist(rng);
      }
    }

    std::tuple<unsigned, unsigned, std::vector<unsigned>> getStats()
    {
      MPI_Allreduce(MPI_IN_PLACE, &numtot, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &numbar, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, counts.data(), numcnt, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
      return {numtot, numbar, counts};
    }

  };
  
public:

  enum DiskGenType { Jeans, Asymmetric, Epicyclic };
  static std::map<std::string, DiskGenType> getDiskGenType;

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

  static double Q;		// Toomre Q.  Default: 0.0 (use
				// exponential)

  static double XI;             // Ratio of radial to tangential
				// velocity dispersion scale
				// length. Default: 1.0

  static double SIG0;           // Velocity dispersion fraction at a
				// scale length. Default: 0.1

  static int MPERT;		// Harmonic perturbation order. Default: 2

  static double AMPL;		// Harmonic amplitude.  Default: 0.0

  static double R_DF;		// Change over points for DF and Jeans
  static double DR_DF;

  static unsigned NBUF;		// Number of particles in MPI buffer
				// default: 8192

  static bool LOGR;		// Radial grid for Eddington inversion

  static bool CHEBY;		// Use Cheybshev fit for epicylic derivatives

  static int  NCHEB;		// Order for Cheybshev fit

  static bool ALLOW;            // Turn off negative mass requeueing in
				// multimass model

  static int SEED;		// Initial seed


  static int ITMAX;		// Maximum number of DF iterations

  static bool use_mono;	        // Use monopole approximation for
				// computing total d(phi)/dr

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
				//  8           Velocity diagnostics

				// ID string for debugging output
  static string RUNTAG;		// Default: "debug"

  //! Dummy Constructor
  Disk2dHalo();

  //! Constructor: equal mass
  Disk2dHalo(SphericalSLptr expandh, std::shared_ptr<Disk2d> expandd,
	     double A, double DMass, std::string maptype,
	     string& filename, int DF=0, int DIVERGE=0, double DIVERGE_RFAC=1.0,
	     DiskGenType type=Asymmetric);

  //! Constructor: multi mass
  Disk2dHalo(SphericalSLptr haloexp, std::shared_ptr<Disk2d> diskexp,
	     double A, double DMass, std::string maptype,
	     string& filename1, int DIVERGE,  double DIVERGE_RFAC,
	     string& filename2, int DIVERGE2, double DIVERGE_RFAC2,
	     DiskGenType type=Asymmetric);

  //! Copy constructor
  Disk2dHalo(const Disk2dHalo &);

  //! Destructor
  ~Disk2dHalo();

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

  void profile(std::ostream &rotfile, vector<Particle>& dpart,
	       double rmin, double rmax, int numr);

  void disk_model(const std::string &modfile);

private:
  DiskGenType type;

};

#endif

