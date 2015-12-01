#ifndef _KARZAS_LATTER_H
#define _KARZAS_LATTER_H

class chdata;
class Ion;

#include "spline.H"

/** 
    Read free-bound gaunt factors for n=1-6 from Karzas and Latter,
    1961, ApJSS, 6, 167
*/
class KLGFdata
{
  //! Pair of unsigned shorts for indexing (Z, C) index pairs
  typedef std::pair<unsigned short, unsigned short> lQ;

private:

  chdata *ch;

  //! Number of levels
  int ngfb;

  //! Rank of energy vector
  int nume;

  //! Energy vector
  std::vector<double> pe;

  //! GF data
  std::map<std::pair<int, int>, std::vector<double> > gfb;

  //! A spline instance
  tk::spline spl;

public:

  //! Constructor
  KLGFdata() {}

  //! Initialize the data from the CHIANTI data file
  void initialize(chdata* ch);


  //! Evaluate the Gaunt factor
  double operator()(double scaledE, int n, int l)
  {
    std::pair<int, int> key(n, l);

    if (gfb.find(key) == gfb.end()) {
      return 0.0;
    }

    spl.set_points(pe, gfb[key]);
    return exp(spl(scaledE));
  }

};

#endif
