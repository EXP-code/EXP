#ifndef _Slab_H
#define _Slab_H

#include <complex>
#include <biorth1d.H>

#include <PotAccel.H>

/*! This routine computes the potential, acceleration and density 
  using periodic box expansion in X & Y and slab in Z */
class Slab : public PotAccel
{

//! Header structure for slab expansion
struct SlabCoefHeader {
  double time;
  double zmax;
  double h;
  int type;
  int nmaxx, nmaxy, nmaxz;
  int jmax;
};

private:

  std::vector<Eigen::VectorXcd> expccof;
  std::vector<std::vector<OneDTrig>> trig;

  int nmaxx, nmaxy, nmaxz;
  int nminx, nminy;
  double zmax;

  int imx, imy, imz, jmax, nnmax;
  double  dfac;
  std::complex<double> kfac;

  std::vector<Eigen::VectorXd> zfrc, zpot;

  SlabCoefHeader coefheader;

  void determine_coefficients(void);
  void get_acceleration_and_potential(Component*);

  // Threading

  void * determine_coefficients_thread(void * arg);
  void * determine_acceleration_and_potential_thread(void * arg);

  // Biorth ID
  static const int ID=0;

protected:

  //! Parse parameters and initialize on first call
  void initialize(void);

public:

  //! Id string
  string id;

  //! Constructor
  Slab(Component* c0, const YAML::Node& conf);

  //! Destructor
  virtual ~Slab();

  //! Print coefficients to output stream
  void dump_coefs(ostream& out);
};


#endif
