#ifndef _VERNER_H
#define _VERNER_H

class atomicData;
class Ion;

extern "C" void phfit2_(int* nz, int* ne, int* is, float* e, float* s);

/** 
   Reads the Verner & Yakovlev (A&AS 109, 125, 1995) photoionization
   cross-section data
*/
class VernerData
{
  //! Pair of unsigned shorts for indexing (Z, C) index pairs
  using lQ = std::pair<unsigned short, unsigned short>;

private:

  atomicData *ad;

public:

  struct VernerRec 
  {
    int n, io, pql, l;
    double eth, e0, sig0, ya, p, yw;

    void sync(int myid=0);
  };
  
  using vrPtr = std::shared_ptr<VernerRec>;
  using vrVec = std::vector<vrPtr>;

  std::map<lQ, vrVec> data;

  //! Constructor
  VernerData() {}

  //! Initialize the data from the CHIANTI data file
  void initialize(atomicData* ad);

  //! Calculates the summed recombination cross section
  double cross(const lQ& Q, double EeV);

  //! Calculates the photoionization cross section.
  double crossPhotoIon(vrPtr V, double EeV);

  //! Call Verner-Ferlind-Korista-Yakolev fortran subroutine.  Returns
  //! cross section in nm^2
  double crossPhotoIon_VFKY(vrPtr V, double EeV)
  {
    int ZZ   = V->n;
    int Nel  = V->io;
    int S    = V->l + 1;
    float ee = EeV, cs;

    phfit2_(&ZZ, &Nel, &S, &ee, &cs);
    
    return cs * 1.0e-4;		// Return cross section in nm^2
  }
};

#endif
