#ifndef _VERNER_H
#define _VERNER_H

class chdata;
class Ion;

/** 
   Reads the Verner & Yakovlev (A&AS 109, 125, 1995) photoionization
   cross-section data
*/
class VernerData
{
  //! Pair of unsigned shorts for indexing (Z, C) index pairs
  typedef std::pair<unsigned short, unsigned short> lQ;

private:

  chdata *ch;

public:

  struct VernerRec 
  {
    int pql, l;
    double eth, e0, sig0, ya, p, yw;

    void sync(int myid=0);
  };
  
  typedef boost::shared_ptr<VernerRec> vrPtr;

  std::map<lQ, vrPtr> data;

  //! Constructor
  VernerData() {}

  //! Initialize the data from the CHIANTI data file
  void initialize(chdata* ch);


  //! Calculates the photoionization cross section.
  double cross(const lQ& Q, double EeV);
};

#endif
