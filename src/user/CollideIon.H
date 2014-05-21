#ifndef _CollideIon_H
#define _CollideIon_H

#include "Collide.H"
#include "interactSelect.H"
#include "Elastic.H"
#include "Ion.H"


//
// Define operator()+ for std::accumulate
//
namespace std {
  template <typename T1, typename T2>
  std::pair<T1, T2> operator+(const std::pair<T1, T2> &p1, const std::pair<T1, T2> &p2) 
  {
    return std::pair<T1, T2>(p1.first + p2.first, p1.second + p2.second);
  } 
}

class collDiag;

class CollisionTypeDiag 
{
  friend class collDiag;

public:

  typedef std::pair<unsigned int, double> udPr;
  typedef std::vector<udPr>               IDVP;
  typedef std::vector<double>             IDVD;
  typedef std::vector<unsigned int>       IDVI;

  IDVP ff, CE, CI, RR, dv;
  IDVD eV_av, eV_max, eV_min;
  IDVI eV_N, eV_10;
  
protected:

  udPr   ff_s, CE_s, CI_s, RR_s, dv_s;
  double eV_av_s, eV_max_s, eV_min_s;
  int    eV_N_s, eV_10_s;
  
public:

  CollisionTypeDiag() 
  {
    const udPr zero(0, 0);

    ff.resize(nthrds, zero);
    CE.resize(nthrds, zero);
    CI.resize(nthrds, zero);
    RR.resize(nthrds, zero);
    dv.resize(nthrds, zero);
    
    eV_av .resize(nthrds, 0);
    eV_N  .resize(nthrds, 0);
    eV_min.resize(nthrds, 999999);
    eV_max.resize(nthrds, 0);
    eV_10 .resize(nthrds, 0);
  }

  void reset() 
  {
    const udPr zero(0, 0);

    std::fill(ff.begin(), ff.end(), zero);
    std::fill(CE.begin(), CE.end(), zero);
    std::fill(CI.begin(), CI.end(), zero);
    std::fill(RR.begin(), RR.end(), zero);
    std::fill(dv.begin(), dv.end(), zero);
    
    std::fill(eV_av .begin(), eV_av .end(), 0);
    std::fill(eV_N  .begin(), eV_N  .end(), 0);
    std::fill(eV_min.begin(), eV_min.end(), 999999);
    std::fill(eV_max.begin(), eV_max.end(), 0);
    std::fill(eV_10 .begin(), eV_10 .end(), 0);
  }

  void sumUp()
  {
    const udPr zero(0, 0);
    
    ff_s     = std::accumulate(ff.begin(), ff.end(), zero);
    CE_s     = std::accumulate(CE.begin(), CE.end(), zero);
    CI_s     = std::accumulate(CI.begin(), CI.end(), zero);
    RR_s     = std::accumulate(RR.begin(), RR.end(), zero);
    dv_s     = std::accumulate(dv.begin(), dv.end(), zero);

    eV_av_s  = std::accumulate(eV_av .begin(), eV_av .end(), 0);
    eV_N_s   = std::accumulate(eV_N  .begin(), eV_N  .end(), 0);
    eV_max_s = std::accumulate(eV_max.begin(), eV_max.end(), 0);
    eV_min_s = std::accumulate(eV_min.begin(), eV_min.end(), 0);
    eV_10_s  = std::accumulate(eV_10 .begin(), eV_10 .end(), 0);
  }

};

typedef boost::shared_ptr<CollisionTypeDiag> collTDPtr;
typedef std::map<speciesKey, collTDPtr> sKeyCollTD;

class collDiag : public sKeyCollTD
{
private:
  std::string coll_file_debug, energy_file_debug;

  void initialize();

public:
  collDiag();
  void gather();
  void print();
  void reset();
};

class CollideIon : public Collide
{
private:
  //! Global variables for the calculations with the individual ions
  std::map <int, std::map<int, Ion> > IonList;

  //! Not sure what this does . . . looks like cruft.  Remove??
  std::vector<int> Ni;

  chdata ch;
  Elastic elastic;
  Geometric geometric;
  InteractSelect IS;
  unsigned int NUM;

  std::map<unsigned short, double> atomic_weights;
  std::map<unsigned short, std::pair<unsigned short, double> > interInfo;


  ACG *gen;
  Uniform *unit;
  
  //@{
  //! Inelastic work space
  typedef std::vector< std::pair<double, double> > CEvector;
  std::map< int, std::vector<double> > dCrossMap;
  std::map< int, std::vector<int>    > dInterMap;
  std::map< int, CEvector  > CE1, CE2;
  double kEi, kEe1, kEe2, Ein1, Ein2;
  //@}

protected:

  //! List of cross-sections (separate for each thread)
  std::vector<sKey2Dmap> csections;
  
  boost::shared_ptr<collDiag> collD;

protected:
  
  double deltaE;
  
  void initialize_cell(pHOT* tree, pCell* c, double rvmax, int id);

  void initialize_cell_dsmc(pHOT* tree, pCell* c, sKey2Umap& nsel, 
			    double rvmax, double tau, int id) {}

  void initialize_cell_epsm(pHOT* tree, pCell* c, sKey2Umap& nsel, 
			    double rvmax, double tau, int id) 
  {
    std::cout << "NO EPSM implemented!" << std::endl;
    exit(-1);
  }
  
  void finalize_cell(pHOT* tree, pCell* c, double kedsp, int id);

  
  sKey2Dmap& totalCrossSections(int id)
  {    
    return csections[id];
  }

  sKey2Dmap& totalScatteringCrossSections(double crm, pCell* c, int id);

  double crossSection(pHOT* tree, Particle* p1, Particle* p2, 
		      double crm, int id=0);
  int inelastic(pHOT *tree, Particle* p1, Particle* p2, double *crm, int id=0);

  bool hasHeatCool() { return true; }

  double getCoolingRate(int id) { return 0; }
  
  
  //@{
  //! Diagnostic output
  double totalSoFar, massSoFar;
  vector<double> lostSoFar;
  
  void list_sizes_proc(ostream*);
  //@}
  
public:
  static double Nmin;
  static double Nmax;
  static double Tmin;
  static double Tmax;
  static double TolV;
  static unsigned Nnum;
  static unsigned Tnum;
  static string cache;
  static bool frost_warning;

  //! Constructor
  CollideIon(ExternalForce *force, double hsDiam, double diamfac, int nth=1);

  //! Destructor
  ~CollideIon();
  
  virtual void set_timestep(int DTpos) { use_delt=DTpos; }
  
  virtual double Etotal();
  virtual double Mtotal();

  virtual void *timestep_thread(void *arg);
  
  virtual void Elost(double* collide, double* epsm);
  
  void Debug(double t);

  void printCollGather()  { collD->gather(); }
  void printCollSummary() { collD->print();  }
  void resetColls()       { collD->reset();  }

};

#endif
