#ifndef _ShearSL_H
#define _ShearSL_H

#include <complex>
#include <memory>

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

  std::shared_ptr<SLGridSlab> grid;

  std::vector<Eigen::VectorXcd> expccof;
  std::vector<std::vector<OneDTrig>> trig;

  int nminx, nminy;
  int nmaxx, nmaxy, nmaxz;
  double zmax, hslab, Omega, Kappa, R, Lx, Ly;

  int imx, imy, imz, jmax, nnmax;

  double dfacx, dfacy;
  double normx, normy;
  std::complex<double> kfacx, kfacy;
  
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

  //! Coefficient container instance for writing HDF5
  CoefClasses::SlabCoefs slabCoefs;

public:

  //! Id string
  string id;

  //! Constructor
  ShearSL(Component* c0, const YAML::Node& conf);

  //! Destructor
  virtual ~ShearSL();

  //! Dump current coefficients into named HDF5 file
  void dump_coefs_h5(const std::string& file);
};


#endif
