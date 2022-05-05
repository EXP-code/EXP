#ifndef _COEFS_H
#define _COEFS_H

#include <yaml-cpp/yaml.h>	// YAML support

#include "coef.H"

struct CylCoefs
{
  double time;
  int nmax, mmax;
  std::vector< std::vector<double> > cos_c, sin_c;

  //! Read from EXP coefficient file
  bool read(std::istream& in, bool verbose=false);

  //! Write an EXP coefficient file with the current data
  void write(std::ostream& out);

  //! Clear the data fields
  void clear()
  {
    for (auto & v : cos_c) std::fill(v.begin(), v.end(), 0.0);
    for (auto & v : sin_c) std::fill(v.begin(), v.end(), 0.0);
  }

  //! Set data from another pointer
  void clone(std::shared_ptr<CylCoefs> ptr)
  {
    time = ptr->time;
    nmax = ptr->nmax;
    mmax = ptr->mmax;

    cos_c.resize(mmax+1);
    sin_c.resize(mmax+1);

    for (int mm=0; mm<=mmax; mm++) {
      cos_c[mm].resize(nmax);
      if (mm) sin_c[mm].resize(nmax);
    }
  }

  //! Zero out coefficients
  void zero()
  {
    for (auto & val : cos_c) std::fill(val.begin(), val.end(), 0.0);
    for (auto & val : sin_c) std::fill(val.begin(), val.end(), 0.0);
  }

};

typedef std::shared_ptr<CylCoefs> CylCoefsPtr;

typedef std::tuple<unsigned, unsigned> LMkey;

struct SphCoefs
{
  SphCoefHeader header;
  std::vector< std::vector<double> > coefs;
  std::map<LMkey, std::vector<double>> cos_c, sin_c;
  bool normed;

  //! Write the EXP coefficient pack to vectors
  void coeftomap();

  //! Write the vectors to the EXp coefficient pack
  void maptocoef();

  //! Read the EXP coefficient file
  bool read(std::istream& in);

  //! Write an EXP coefficient file with the current data
  void write(std::ostream& out);

  //! Set data from another pointer
  void clone(std::shared_ptr<SphCoefs> ptr)
  {
    header = ptr->header;
    coefs  = ptr->coefs;
    normed = ptr->normed;

    for (auto val : ptr->cos_c)
      cos_c[val.first] = val.second;

    for (auto val : ptr->sin_c)
      sin_c[val.first] = val.second;
  }

  //! Zero out coefficients
  void zero()
  {
    for (auto & val : cos_c) std::fill(val.second.begin(), val.second.end(), 0.0);
    for (auto & val : sin_c) std::fill(val.second.begin(), val.second.end(), 0.0);
  }

};

typedef std::shared_ptr<SphCoefs> SphCoefsPtr;

#endif
