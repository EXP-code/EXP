#ifndef _VERNER_H
#define _VERNER_H

class chdata;
class Ion;

extern "C" void phfit2_(int* nz, int* ne, int* is, float* e, float* s);

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


  //! Calculates the recombination cross section.
  double cross(const lQ& Q, double EeV);

  //! Calculates the photoionization cross section.
  double crossPhotoIon(const lQ& Q, double EeV);

  //! Call Verner-Korista-Yakolev fortran subroutine.  Returns cross
  //! section in nm^2
  double crossPhotoIon_VKY(const lQ& Q, double EeV)
  {
    int ZZ  = Q.first;
    int Nel = ZZ - Q.second + 1;
    float ee = EeV, cs, csum=0.0;

    for (int S=1; S<=7; S++) {
      phfit2_(&ZZ, &Nel, &S, &ee, &cs);
      if (S>0) csum += cs;
    }
    
    return csum * 1.0e-4;
  }
};

#endif
