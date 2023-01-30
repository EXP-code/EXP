// This may look like C code, but it is really -*- C++ -*-

#ifndef _RespMat_h
#define _RespMat_h

#include <cmath>
#include <complex>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>

#include <numerical.H>
#include <massmodel.H>
#include <biorth.H>
#include <orbit.H>
#include <gaussQ.H>
#include <Timer.H>

enum IntegrationType {ratint, jacoint, rombint};
// static const char* SIType[] = {"Rational", "Jacobi", "Rombe"};

class OrbitTable;
class OrbResTable;
class ResTable;
class DiskHalo;

class RespMat : public std::enable_shared_from_this<RespMat>
{
  friend class OrbitTable;
  friend class OrbResTable;
  friend class ResTable;
  friend class DiskHalo;

protected:
				// Parameters
  Eigen::MatrixXcd mr;
  Eigen::MatrixXcd mr2;
  std::complex<double> omp;

  int l, m, lmax, nmax, nptsK, nptsE;
  AxiSymModPtr model;
  AxiSymBioPtr biorth;

  int OLD, VERBOSE;
  IntegrationType SITYPE;

  int pv;
  int dof;
  int isotropic;

				// Variables

  int alpha, beta, ikap, l1, l2;
  double E, kappa, E_r, Emodmin, Emodmax, Emax;
  std::vector<Eigen::VectorXcd> mab, hold;
  Eigen::MatrixXcd holdR;
  int num_E, num_K;
  Eigen::VectorXcd pp_x, pp_y;
  Eigen::VectorXd pp_x2;

				// Get previously recorded CHECK POINT set
  bool get_chkpnt();
				// Store CHECK POINT set (after every kappa)
  void put_chkpnt();

				// return maximum E for a given kappa
  void get_emax(void); 
  double E_max(double E);

  int matrix_computed;
  int matrix_computed_pv;

  void integrand_unroll(void);
  void integrand_rolled(void);
  void integrand_orig(void);
  void integrand_pv(void);
  double Ylm02(int ll, int mm);

  string ID;
  string CHK_NAME;
  double ESTEP;

  double time_in_inner_loop;
  double time_in_biorth_gen;

  int RATINT;
  double DELTA;
  int RECS;
  int PVALUE;

  int disp_computed;
  std::complex<double> disp;

  std::shared_ptr<OrbitTable>  Orb;
  std::shared_ptr<OrbResTable> ORes;

  int CPREAD, CHECKPT;

  void read_in(istream& in);
  void read_in(string& file);

  void bomb(string msg) {
    cerr << msg << endl;
    exit(-1);
  }

  static const int defaultRECS = 400;
  static const int defaultPVALUE = 1;

  Timer timer;

  std::shared_ptr<RespMat> getptr() {
    return shared_from_this();
  }


 public:

				// Constructor

  RespMat();

  RespMat(std::complex<double> omp0, int ll, int mm,
	  int llmax, int nnmax, int nnptsK, int nnptsE, 
	  AxiSymModPtr MODEL, AxiSymBioPtr BIORTH,
	  int old, int cpread, int checkpt, int verbose, 
	  IntegrationType type, 
	  int principal=1, int ratint=0, int isotropic=0);

				// Transformed response version
  RespMat(double omp0, int ll, int mm,
	  int llmax, int nnmax, int nnptsK, int nnptsE, 
	  AxiSymModPtr MODEL, AxiSymBioPtr BIORTH,
	  std::string& chkpt, int verbose);

  RespMat(std::istream& in, AxiSymModPtr MODEL, AxiSymBioPtr BIORTH) {
    model  = MODEL;
    biorth = BIORTH;
    read_in(in);
  }

  RespMat(std::string& from, AxiSymModPtr MODEL, AxiSymBioPtr BIORTH) {
    model = MODEL;
    biorth = BIORTH;
    read_in(from);
  }

  RespMat(const RespMat&);

  //  RespMat &operator=(RespMat &);

				// Parameters

  void set_params(double delta, int recs=200) { DELTA = delta; RECS = recs; }
  
  void set_ratint(int ratint) { RATINT = ratint; }

  void set_principal_value(int pvalue);


				// Do the work!
  void make_matrix(void);
  void make_matrix_pv(void);

				// Input, Output

  RespMat(string& file) {read_in(file);}
  RespMat(istream& in) {read_in(in);}

  void write_out(ostream& out);
  void write_out(string& file);

  void MPISynchronize(int id);

				// Access
  Eigen::MatrixXcd& mat() {
    if (!matrix_computed) {
      if (pv) make_matrix_pv();
      else make_matrix();
    }
    return mr;}

  Eigen::MatrixXcd& mat_asymptotic() {
    if (!matrix_computed_pv) make_matrix_pv();
    return mr2;}

  AxiSymBioPtr bio()  { return biorth; }
  AxiSymModPtr modl() { return model;  }

  // Utilities

  std::complex<double> disper(void);
  std::complex<double> omega(void) {return omp;}
  const int& get_lmax(void) {return lmax;}
  const string& id(void) {return ID;}
  const int& L(void) {return l;}
  const int& M(void) {return m;}

				// Response matrix calcs

// <b>none</b>    input vector is used directly to generate wake
// <b>noself</b>  input vector used as non-self-gravitating perturbation 
//			 to compute response
// <b>self</b>    self-gravitating respones computed

  enum gravity {none, noself, self};

// type of reponse:
  enum response {density, potential};

  Eigen::VectorXcd get_response (Eigen::VectorXcd& ext, gravity grav=self);

  Eigen::MatrixXd get_wake (Eigen::VectorXcd& ext,
		   double xmin, double xmax, double ymin, double ymax,
		   int numx, int numy,
		   gravity grav=self, response type=density,
		   double phi=0, double theta=0);

  void get_wake (Eigen::MatrixXd& wake, Eigen::VectorXcd& ext,
		   double xmin, double xmax, double ymin, double ymax,
		   int numx, int numy,
		   gravity grav=self, response type=density,
		   double phi=0, double theta=0);

  Eigen::MatrixXd get_bkgrnd (
		   double xmin, double xmax, double ymin, double ymax,
		   int numx, int numy, response type=density,
		   double phi=0, double theta=0);

  void print_wake (ostream& out, Eigen::VectorXcd& ext,
		   double xmin, double xmax, double ymin, double ymax,
		   int numx, int numy,
		   gravity grav=self, response type=density,
		   double phi=0, double theta=0);

  void print_wake_volume (ostream& out, Eigen::VectorXcd& ext,
			  double xmin, double xmax, 
			  double ymin, double ymax,
			  double zmin, double zmax,
			  int numx, int numy, int numz,
			  bool bkgrnd=true,
			  gravity grav=self, response type=density);

  void get_wake_volume (std::vector<Eigen::MatrixXd>& mat,
			Eigen::VectorXcd& ext,
			double xmin, double xmax, 
			double ymin, double ymax,
			double zmin, double zmax,
			int numx, int numy, int numz,
			bool bkgrnd=true,
			gravity grav=self, response type=density);


  void print_wake_volume_hips (ostream& out, Eigen::VectorXcd& ext,
			       double xmin, double xmax, 
			       double ymin, double ymax,
			       double zmin, double zmax,
			       int numx, int numy, int numz,
			       bool bkgrnd=true,
			       gravity grav=self, response type=density);


  void print_wake_volume_cubes (ostream& out, Eigen::VectorXcd& ext,
			       double xmin, double xmax, 
			       double ymin, double ymax,
			       double zmin, double zmax,
			       int numx, int numy, int numz,
			       bool bkgrnd=true,
			       gravity grav=self, response type=density);


};

				// Orbit list stuff

class OrbitTable
{
private:
  std::shared_ptr<RespMat> p;
  Timer timer;

public:
				// Constructor
  OrbitTable(std::shared_ptr<RespMat> p, int recs=200);
				// Destructor
  ~OrbitTable(void);
    
				// Data
  std::vector<SphericalOrbit>  orbits;
  std::vector<Eigen::VectorXd> PotTrans;
  Eigen::VectorXd dfqE, dfqL;
  Eigen::VectorXd norm;
  Eigen::VectorXcd EInt;
      
  void setup(void);
  void setup_p_array(void);
};


				// Orbit table and resonance location
class OrbResTable
{
private:
  std::shared_ptr<RespMat> p;
  Timer timer;

public:
  static double derivative_step;
				// Constructor
  

  OrbResTable(std::shared_ptr<RespMat> p, int recs=200);
				// Destructor
  ~OrbResTable(void);
    
				// Data
  std::vector<std::shared_ptr<SphericalOrbit>> orbits;
  SphericalOrbit torbit;
  std::vector<Eigen::VectorXd> PotTrans;
  std::vector<Eigen::VectorXd> PotTrans2;
  Eigen::VectorXd dfqE, dfqL;
  Eigen::VectorXd norm;
  Eigen::VectorXd EInt;
  Eigen::VectorXd ELoc;
  Eigen::VectorXd ERes;
  Eigen::VectorXd EJac;
  int number;
    
  void setup(void);
  void setup_E_array(void);
};


#endif

