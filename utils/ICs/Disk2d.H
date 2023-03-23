#ifndef _Disk2d_H
#define _Disk2d_H

#include <memory>
#include <tuple>
#include <map>

#include <yaml-cpp/yaml.h>

#include <BiorthCyl.H>
#include <Particle.H>

typedef std::shared_ptr<BiorthCyl> CylPtr;

/** Computes the potential, acceleration and density using an empirical
    function analysis for a flat Bessel basis and a target density

    \param acyltbl is the scale length of target density (default is
    0.6)

    \param rcylmin is the minimum value in the table for the 
    radial basis functions (default is 0.0)

    \param rcylmax is the maximum value in the table for the 
    radial basis functions in acyltbl units (default is 10.0)

    \param numr is the number radial coordinate knots in the table 
    (default is 1000)

    \param knots is the number quadrature points for integrals in
    basis construction

    \param logr scales the EmpCyl2d grid logarithmically

    \param model name for EmpCyl2d (default: expon)

    \param biorth set for EmpCyl2d (default: bess)

    \param Mmax is the azimuthal order

    \param nmaxfid is the radial order for the Bessel basis used to
    construct the EOF

    \param nmax is the radial order

    \param numx is the number of grid points in the scaled radial direction

    \param numy is the number of grid points in the scaled vertical direction

    \param NQDHT is the number of Bessel knots in the underlying
    Hankel transform used to get the vertical potential

    \param cmapR is the scaling method for radius (0=none, 1=rational)
    in the grid mapping

    \param cmapZ is the scaling method for vertical height (0=none, 1=sinh, 2=rational)
    in the grid mapping

    \param scale is the expansion factor to get the physical scale from
    the internal scale (default is 0.01)

    \param cachename is the file cache for the 2d EOF basis

    \param verbose set to true turns on some additional diagnostic
    info in BiorthCyl

    \param mmin changes the minimum harmonic order for evaluation to
    some value >=0

    \param mlim changes the maximum harmonic order for evaluation to
    some value <=mmax

    \param nmin changes the minimum radial order for evaluation to
    some value >=0

    \param nlim changes the maximum radial order for evaluation to
    some value <=nmax

    \param EVEN_M evaluates the expansion for even azimuthal harmonics
    only
*/
class Disk2d
{

private:

  //! Valid key set for parameter name checking
  static std::set<std::string> valid_keys;

  //! Config db
  YAML::Node conf;

  //! Basis instance
  CylPtr ortho;

  //! Initialize parameters
  void initialize(void);

  //@{
  //! Field evaluation members
  
  void get_dpotl(double r, double z,
		 Eigen::MatrixXd& p, Eigen::MatrixXd& dpr, Eigen::MatrixXd& dpz, int tid);

  void get_potl(double r, double z, Eigen::MatrixXd& p, int tid);

  virtual double getRtable() { return ortho->getRtable(); }

  void get_dens(double r, double z, Eigen::MatrixXd& p, int tid);

  void get_potl_dens(double r, double z,
		     Eigen::MatrixXd& p, Eigen::MatrixXd& d, int tid);

  //@}

  //@{
  //! Parameters
  double acyltbl, scale, rcylmin, rcylmax;
  int    numr;
  int    mmax;
  int    nmax;
  int    knots;
  string model;
  string biorth;
  bool   logr;
  //@}

  //@{
  //! Internal state variables
  int    nthrds;
  bool   accum_setup, coefs_made;
  //@}

protected:

  //! Use BiorthCyl to evaluate the basis at r, z
  void get_pot(Eigen::MatrixXd& Vc, Eigen::MatrixXd& Vs, double r, double z)
  { ortho->get_pot(Vc, Vs, r, z); }

  //! Matrices per thread for obtaining expansion coefficients
  std::vector<Eigen::MatrixXd> expcoef0;

  //! Potential coefficient computation
  std::vector<Eigen::MatrixXd> potd;

  //! Field evaluation
  Eigen::MatrixXd dend, potl, potr, potz;

  //! Final coefficiet set
  Eigen::MatrixXd expcoef;

  //! Mass tally
  std::vector<double> cylmass1;
  double cylmass;

  //! Setup accumulation arrays
  void setup_accum();

public:
				// Global parameters
  //! Constructor
  Disk2d(const YAML::Node& conf);

  //! Destructor
  virtual ~Disk2d();

  //! Set up for accumulations by zeroing counters, allocated storage, etc.
  void setup_accumulation();

  //! Finish making the coefficients after accumulation
  void make_coefficients();

  //! Add single particle to coefficients
  void accumulate(double r, double phi, double mass);

  //! Accumulate coefficients from particle distribution
  void accumulate(std::vector<Particle>& part);

  //! Evaluation of field. Returns: den0, den1, pot0, pot1, Fr, Fz, Fp
  std::tuple<double, double, double, double, double, double, double>
  accumulated_eval(double R, double z, double phi);
};



#endif


