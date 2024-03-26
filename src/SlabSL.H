#ifndef _SlabSL_H
#define _SlabSL_H

#include <complex>

#include <Eigen/Eigen>
#include <unsupported/Eigen/CXX11/Tensor>

#include <Coefficients.H>
#include <SLGridMP2.H>
#include <biorth1d.H>
#include <PotAccel.H>

/*! This routine computes the potential, acceleration and density
  using expansion periodic in X & Y and outgoing vacuum boundary
  condtions in Z */
class SlabSL : public PotAccel
{

//! Header structure for Sturm-Liouville slab expansion
//! Used for deprecated native coefficient files
struct SlabSLCoefHeader {
  double time;
  double zmax;
  double h;
  int type;
  int nmaxx, nmaxy, nmaxz;
  int jmax;
};

private:

  std::shared_ptr<SLGridSlab> grid;

  //! Coefficients are a 3-tensor
  using coefType = Eigen::Tensor<std::complex<double>, 3>;

  //! Current coefficient tensor
  std::vector<coefType> expccof;

  int nminx, nminy;
  int nmaxx, nmaxy, nmaxz;
  double zmax, hslab;

  int imx, imy, imz, jmax, nnmax;
  double  dfac;
  std::complex<double> kfac;

  std::vector<Eigen::VectorXd> zfrc, zpot;

  SlabSLCoefHeader coefheader;

  int NGRID = 100;

  // Usual evaluation interface

  void determine_coefficients(void);
  void get_acceleration_and_potential(Component*);

  // Threading

  void * determine_coefficients_thread(void * arg);
  void * determine_acceleration_and_potential_thread(void * arg);

  //! Coefficient container instance for writing HDF5
  CoefClasses::SlabCoefs slabCoefs;

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
  SlabSL(Component* c0, const YAML::Node& conf);

  //! Destructor
  virtual ~SlabSL();

  //! Coefficient output
  void dump_coefs_h5(const std::string& file);

  //! Print coefficients to output stream
  void dump_coefs(ostream& out);
};


#endif
