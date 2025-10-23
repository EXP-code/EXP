#ifndef _SphericalSL_H
#define _SphericalSL_H

#include <vector>

#include "Particle.H"
#include "Particle.H"
#include "SLGridMP2.H"

class SphericalSL
{
private:
  int NMAX, LMAX, n, nthrds;

  double *MPIin, *MPIout;
  bool MPIset;

  std::shared_ptr<SLGridSph> ortho;

  Eigen::MatrixXd expcoef;
  std::vector<Eigen::MatrixXd> cc;
  std::vector<Eigen::MatrixXd> expcoef1;
  std::vector<std::vector<Eigen::MatrixXd>> cc1;

  Eigen::MatrixXd factorial;
  Eigen::MatrixXd normM, krnl;
  std::vector<Eigen::MatrixXd> potd, legs, dpot, dlegs, dend;
  std::vector<Eigen::VectorXd> cosm, sinm;

  bool coefs_made;

  bool initialized;
  void initialize(void);

  std::vector<int> use;
  int compute, used;


  void compute_coefficients_single(vector<Particle> &part);
  void compute_coefficients_thread(vector<Particle> &part);
  void compute_coefficients_thread_call(int id, vector<Particle>* p);

  void parallel_gather_coefficients
  (Eigen::MatrixXd& expcoef,
   Eigen::MatrixXd& expcoef1,
   std::vector<Eigen::MatrixXd>& cc,
   std::vector<Eigen::MatrixXd>& cc1,
   int lmax);

  void parallel_distribute_coefficients(Eigen::MatrixXd& expcoef);

  void pca_hall(int compute, Eigen::MatrixXd& expcoef,
		std::vector<Eigen::MatrixXd>& cc, Eigen::MatrixXd& normM);

				// Error function

  void bomb(const char *s) {
    cerr << "SphericalSL ERROR: " << s << '\n';
    exit(-1);
  }


public:
				// Global parameters

  static double RMIN;
  static double RMAX;
  static int NUMR;

  static int selector;		// For PCA methods
  static int tk_type;
  static double tksmooth;
  static double tkcum;



				// Constructors
  SphericalSL(void);
  SphericalSL(std::string file, int Nth, int lmax, int nmax, int CMAP, double SCALE=1.0);
  ~SphericalSL(void);

  void reset(std::string file, int Nth, int lmax, int nmax, int CMAP, double SCALE);

				// Parameter access
  int get_maxNR(void) {return NMAX;}
  int get_maxNL(void) {return LMAX;}
				// Pointer to orthgonal function instance
  std::shared_ptr<SLGridSph> SL(void) {return ortho;}

				// Main member functions

  void accumulate(vector<Particle>& particles);

  
  void determine_fields_at_point(double r, double theta, double phi,
				 double *tdens, double *tpotl, 
				 double *tpotr, double *tpott, 
				 double *tpotp, int id=0);
  

  void get_potl_dens(int l, int n, double r, int id=0)
  {
    ortho->get_dens (dend[id], r);
    ortho->get_pot  (potd[id], r);
    ortho->get_force(dpot[id], r);
  }

  void get_dens_coefs(int l, const Eigen::VectorXd& coef, double& p, int id=0);
  void get_pot_coefs (int l, const Eigen::VectorXd& coef, double& p, double& dp, int id=0);

  void dump_coefs(ofstream& out, bool binary);
  void dump_basis(string& dumpname);

  void set_compute()   { compute = 1;}
  void unset_compute() { compute = 0;}

};

typedef std::shared_ptr<SphericalSL> SphericalSLptr;

#endif
