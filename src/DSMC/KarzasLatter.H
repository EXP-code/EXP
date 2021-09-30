#ifndef _KARZAS_LATTER_H
#define _KARZAS_LATTER_H

class atomicData;
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

  atomicData *ad;

  //! Number of levels
  int ngfb;

  //! Rank of energy vector
  int nume;

  //! Energy vector
  std::vector<double> pe;

  //! Spline pointer
  typedef std::shared_ptr<tk::spline> tksplPtr;

  //! GF data
  std::map<std::pair<int, int>, std::vector<double> > gfb;

  //@{
  //! Spline instances
  typedef std::map<std::pair<int, int>, tksplPtr> splMap;
  splMap spl;
  //@}
  
public:

  //! Constructor
  KLGFdata() {}

  //! Initialize the data from the CHIANTI data file
  void initialize(atomicData* ad);


  //! Evaluate the Gaunt factor
  double operator()(double scaledE, int n, int l)
  {
    std::pair<int, int> key(n, l);

    splMap::iterator s = spl.find(key);

    if (s == spl.end()) {
      return 0.0;
    }

    return exp((*(s->second))(scaledE));
  }

};

#endif
