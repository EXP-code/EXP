#ifndef _ShearSL_H
#define _ShearSL_H

#include <complex>
#include <SLGridMP2.H>
#include <biorth1d.H>
#include <PotAccel.H>

/**
   This routine computes the potential, acceleration and density
   using expansion periodic in X & Y and outgoing vacuum boundary
   condtions in Z for a shearing slab
*/
class ShearSL : public PotAccel
{
private:

  SLGridSlab *grid;

  std::vector<Eigen::VectorXcd> expccof;
  std::vector<std::vector<OneDTrig>> trig;

  int nminx, nminy;
  int nmaxx, nmaxy, nmaxz;
  double zmax, hslab, Omega, Kappa, R;

  int imx, imy, imz, jmax, nnmax;
  double  dfac;
  std::complex<double> kfac;

  std::vector<Eigen::VectorXd> zfrc, zpot;

  int NGRID;

  // Usual evaluation interface

  void determine_coefficients(void);
  void get_acceleration_and_potential(Component*);

  // Threading

  void * determine_coefficients_thread(void * arg);
  void * determine_acceleration_and_potential_thread(void * arg);

  // Biorth ID
  static const int ID=1;

protected:

  //! Parse parameters and initialize on first call
  void initialize(void);

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Id string
  string id;

  //! Constructor
  ShearSL(Component* c0, const YAML::Node& conf);

  //! Destructor
  virtual ~ShearSL();

  //! Print coefficients to output stream
  void dump_coefs(std::ostream& out);
};


#endif
