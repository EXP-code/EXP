#ifndef _CollideIon_H
#define _CollideIon_H

#include <tuple>

#include "Collide.H"
#include "interactSelect.H"
#include "Elastic.H"
#include "Ion.H"
#include "atomic_constants.H"

//
// Define operator()+ for std::accumulate
//
namespace std {
  template <typename T1, typename T2>
  std::pair<T1, T2> operator+(const std::pair<T1, T2> &p1, 
			      const std::pair<T1, T2> &p2) 
  {
    return std::pair<T1, T2>(p1.first + p2.first, p1.second + p2.second);
  } 
}

class collDiag;

/**
   This helper class defines a data structure and manipulators for the
   diagnostic collDiag class.
*/
class CollisionTypeDiag 
{
  friend class collDiag;

public:

  //@{
  //! Container definition
  typedef std::tuple<unsigned, double, double>  UDDT;
  typedef std::vector<UDDT>                     IDVP;
  typedef std::vector<double>                   IDVD;
  typedef std::vector<unsigned long>            IDVI;
  //@}

  //@{
  //! Vectors for each type of diagnostic by thread id
  IDVP ff, CE, CI, RR, dv;
  IDVD eV_av, eV_max, eV_min;
  IDVI eV_N, eV_10;
  //@}
  
protected:

  //@{
  //! Variables to hold processes collective reductions
  UDDT          ff_s, CE_s, CI_s, RR_s, dv_s;
  double        eV_av_s, eV_max_s, eV_min_s;
  unsigned long eV_N_s, eV_10_s;
  //@}
  
  //! Synchronize a UDDT type (tuple)
  void syncUDDT(UDDT& value) 
  {
    unsigned int u;
    double z;

    MPI_Reduce(&(u=std::get<0>(value)), &std::get<0>(value), 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&(z=std::get<1>(value)), &std::get<1>(value), 1, MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&(z=std::get<2>(value)), &std::get<2>(value), 1, MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD);
  }

private:

  UDDT zero;

  UDDT accum(IDVP& v)
  {
    UDDT ret(zero);

    for (auto t : v) {
      std::get<0>(ret) += std::get<0>(t);
      std::get<1>(ret) += std::get<1>(t);
      std::get<2>(ret) += std::get<2>(t);
    }

    return ret;
  }

public:

  //! Construtor initializes the vectors
  CollisionTypeDiag() 
  {
    zero = UDDT(0, 0, 0);

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

  //! Reset the accumulation vectors
  void reset() 
  {
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

  //! Sum up over all threads for processing data by the collDiag
  //! instance
  void sumUp()
  {
    ff_s     =  accum(ff);
    CE_s     =  accum(CE);
    CI_s     =  accum(CI);
    RR_s     =  accum(RR);
    for (auto v : dv) {
      if (std::isnan(std::get<2>(v)))
	std::cout << "Error in sumUp (before), t=" << tnow 
		  << ", size=" << dv.size() << ", id=" << myid << std::endl;
    }
    dv_s     =  accum(dv);
    if (std::isnan(std::get<2>(dv_s)))
      std::cout << "Error in sumUp (before), t=" << tnow 
		<< ", id=" << myid << std::endl;

    eV_av_s  =  std::accumulate (eV_av .begin(), eV_av .end(), 0);
    eV_N_s   =  std::accumulate (eV_N  .begin(), eV_N  .end(), 0);
    eV_max_s = *std::max_element(eV_max.begin(), eV_max.end()   );
    eV_min_s = *std::min_element(eV_min.begin(), eV_min.end()   );
    eV_10_s  =  std::accumulate (eV_10 .begin(), eV_10 .end(), 0);
  }

  void sync()
  {
    syncUDDT(ff_s);
    syncUDDT(CE_s);
    syncUDDT(CI_s);
    syncUDDT(RR_s);
    syncUDDT(dv_s);

    double z;
    unsigned long u;
    MPI_Reduce(&(z=eV_av_s),  &eV_av_s,  1, MPI_DOUBLE,        MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&(u=eV_N_s),   &eV_N_s,   1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&(z=eV_min_s), &eV_min_s, 1, MPI_DOUBLE,        MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&(z=eV_max_s), &eV_max_s, 1, MPI_DOUBLE,        MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&(u=eV_10_s),  &eV_10_s,  1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  }
};

typedef boost::shared_ptr<CollisionTypeDiag> collTDPtr;
typedef std::map<speciesKey, collTDPtr> sKeyCollTD;

/**
 This helper class maintains database of diagnostics on the CollideIon
 collisions.
*/
class collDiag : public sKeyCollTD
{
private:
  CollideIon *p;
  std::string coll_file_debug, energy_file_debug;
  std::vector<double> Esum;
  double Esum_s;

  //! Make the file header the first time
  void initialize();

public:
  //! Constructor
  collDiag(CollideIon *caller);

  //! All processes send info to the root process
  void gather();

  //! Print the current data to the diagnostic file
  void print();

  //! Zero out the structures maintained by the CollisionTypeDiag maps
  void reset();

  //! Add cell energy
  void addCell(double E, int id) { Esum[id] += E; }

};

/**
   Implementation of a DSMC collision algorithm based on the CHIANTI
   atomic data base.  We have implemented two distinct collision algorithms:
   <ol>
   <li> A standard DSMC algorithm that assigns each element a distinct
   super particle with a dynamic ionization state.
   <li> A trace-species DSMC algorithm that assigns each element the
   same super-particle mass and keeps track of the multiple elements
   and their ionization states using a vector of mass weights.
   </ol>

   The parseSpecies member function reads a specification file whose
   first line must include the tag string "direct", "weight", or
   "trace" to select the algorithm.  The "direct" file contains lines
   specifying the atomic number, Z, of the species represented in the
   simulation, one value per line.  The "weight" file contains lines
   specifying the atomic number, Z, of the species represented in the
   simulation followed by the relative particle weight, two values per
   line.  The "trace" file contains line with triples of numbers
   containing Z, the charge state C (C=1 means neutral, C=2 is the
   first ionization state, etc.) and the index position of the mass
   fraction in the Particle dattrib vector.
 */
class CollideIon : public Collide
{
  friend class collDiag;

private:

  chdata ch;

  //! Neutral atom--electron cross-section for hydrogen and helium
  //! interpolated from experimentally calibrated theortical data (see
  //! Elastic.H)
  Elastic elastic;

  //! Geometric cross section for each atom based on theoretical
  //! comptuation of atomic radius in units of nm^2 for consistency
  //! with CHIANTI.
  Geometric geometric;

  InteractSelect IS;

  ACG     *gen;
  Uniform *unit;
  
  //@{
  //! Inelastic work space
  typedef std::vector< std::pair<double, double> >   CEvector;
  typedef std::pair<speciesKey, speciesKey>          dKey;
  typedef std::pair<dKey, double>                    dKeyD;
  typedef std::map<speciesKey, std::vector<double> > keyCross;
  typedef std::map<speciesKey, std::vector<int>    > keyInter;
  typedef std::map<speciesKey, CEvector>             keyColEx;
  typedef std::map<speciesKey, double>               keyWghts;
  typedef std::list<dKeyD>                           keyWexsList;

  std::vector< std::vector<double> >                 dCross;
  std::vector< std::vector<int>    >                 dInter;
  std::vector< keyCross >                            sCross; 
  std::vector< keyWghts >                            meanF;
  std::vector< keyWghts >                            sCrsTot1, sCrsTot2;
  std::vector< keyWexsList >                         excessW;
  std::vector< keyInter >                            sInter;
  std::vector< CEvector >                            CE1, CE2;
  std::vector< keyColEx >                            kCE;
  std::vector< double >                              meanE, neutF, meanR, meanM;
  std::vector< double >                              kEi, kEe1, kEe2;
  std::vector< double >                              Ein1, Ein2;
  std::vector< double >                              spCrm, spTau, spProb;
  std::vector< unsigned >                            spNsel;
  //@}

  //! This defines the various types of interactions used to tag
  //! non-zero cross sections and subsequent interaction processing
  enum InterTypes { 
    neut_neut   = 0,
    neut_elec   = 1,
    ion_elec    = 2,
    free_free   = 3,
    colexcite   = 4,
    recomb      = 5,
    ionize      = 6,
    neut_neut_1 = 100,
    neut_elec_1 = 101,
    ion_elec_1  = 102,
    free_free_1 = 103,
    colexcite_1 = 104,
    ionize_1    = 105,
    recomb_1    = 106,
    neut_neut_2 = 200,
    neut_elec_2 = 201,
    ion_elec_2  = 202,
    free_free_2 = 203,
    colexcite_2 = 204,
    ionize_2    = 205,
    recomb_2    = 206
  };    

  std::map<int, std::string> labels;

  //! Mnemonic for the type of algorithm selected by the user
  enum AlgorithmType { Direct, Weight, Trace };
  const char* AlgorithmLabels[3] = {"Direct", "Weight", "Trace"};

public:

  typedef std::set<unsigned short> ZLtype;
  typedef std::map<unsigned short, double> ZWtype;

protected:

  //! List of cross sections (separate for each thread)
  std::vector<sKey2Dmap> csections;

  boost::shared_ptr<collDiag> collD;

  //! Set of the elements (atomic number) for direct multi-species
  //! implementation
  ZLtype ZList;

  //! Set of the elements (atomic number) for weighted multi-species
  //! implementation
  ZWtype ZWList;

  //! Designated non-trace element
  unsigned short sFid;

  //! Map of species to attribute position for trace multi-species
  //! implementation
  std::map<speciesKey, int> SpList;

  //! Attribute location for energy excess in "weight" algorithm
  int use_cons;

  //! The algorithm type
  AlgorithmType aType;

protected:
  
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

  double crossSectionDirect(pHOT* tree, Particle* p1, Particle* p2, 
		      double crm, int id=0);
  
  int inelasticDirect(pHOT *tree, Particle* p1, Particle* p2, double *crm, int id=0);

  double crossSectionWeight(pHOT* tree, Particle* p1, Particle* p2, 
			    double crm, int id=0);
  
  int inelasticWeight(pHOT *tree, Particle* p1, Particle* p2, double *crm, int id=0);

  double crossSectionTrace(pHOT* tree, Particle* p1, Particle* p2, 
			   double crm, int id=0);

  int inelasticTrace(pHOT *tree, Particle* p1, Particle* p2, double *crm, int id=0);

  double crossSection(pHOT* tree, Particle* p1, Particle* p2, 
		      double crm, int id=0)
  {
    if (aType == Direct)
      return crossSectionDirect(tree, p1, p2, crm, id);

    if (aType == Weight)
      return crossSectionWeight(tree, p1, p2, crm, id);

    if (aType == Trace)
      return crossSectionTrace(tree, p1, p2, crm, id);

    return 0.0;
  }


  int inelastic(pHOT *tree, Particle* p1, Particle* p2, double *crm, int id=0)
  {
    if (aType == Direct)
      return inelasticDirect(tree, p1, p2, crm, id);

    if (aType == Weight)
      return inelasticDirect(tree, p1, p2, crm, id);

    if (aType == Trace)
      return inelasticTrace(tree, p1, p2, crm, id);

    return 1;
  }


  bool hasHeatCool() { return true; }

  double getCoolingRate(int id) { return 0; }
  
  void parseSpecies(const std::string& map);

  
  //@{
  //! Diagnostic output
  double totalSoFar, massSoFar;
  vector<double> lostSoFar;
  
  void list_sizes_proc(ostream*);
  //@}
  
  //@{
  //! Cross-section debugging
  std::string cross_debug;
  std::vector<double> cross1_dbg, cross2_dbg;
  double nextTime_dbg;
  int nCnt_dbg;
  static constexpr double delTime_dbg = 0.001;
  static constexpr int nCel_dbg = 10;
  void write_cross_debug();
  //@}


  virtual sKey2Umap generateSelectionDirect(pCell* c, sKeyDmap* Fn, double crm, double tau, int id,
					    double& meanLambda, double& meanCollP, 
					    double& totalNsel);

  virtual sKey2Umap generateSelectionWeight(pCell* c, sKeyDmap* Fn, double crm, double tau, int id,
					    double& meanLambda, double& meanCollP, 
					    double& totalNsel);

  virtual sKey2Umap generateSelectionTrace(pCell* c, sKeyDmap* Fn, double crm, double tau, int id,
					   double& meanLambda, double& meanCollP, 
					   double& totalNsel);


  //! Cache instance of current tree
  pHOT* curTree;

  //@{
  //! Trace species maps for diagnostics
  typedef std::map<speciesKey, int> spMap;
  typedef spMap::iterator spItr;
  void printSpeciesTrace();
  void gatherSpecies();

  typedef std::map<speciesKey, double> spDMap;
  typedef spDMap::iterator spDItr;
  double tempM, molwM;
  spDMap specM;
  //@}

  // For debugging
  void debugDeltaE(double delE, unsigned short Z, unsigned short C,
		   double dE, double prob, int interFlag);

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
  CollideIon(ExternalForce *force, Component *comp, 
	     double hsDiam, double diamfac, const std::string& spcMap, int nth=1);

  //! Destructor
  ~CollideIon();
  
  //! Sets timestep processing for the EXP multiple time-step algorithm
  virtual void set_timestep(int DTpos) { use_delt=DTpos; }
  
  virtual double Etotal();
  virtual double Mtotal();

  //! Times the per-thread processing time
  virtual void *timestep_thread(void *arg);
  
  virtual void Elost(double* collide, double* epsm);
  
  void Debug(double t);

  //! Print out species counts
  virtual void printSpecies(std::map<speciesKey, unsigned long>& spec, double T=0);

  //@{
  //! Galther and print collision statistics
  void printCollGather()  { collD->gather(); }
  void printCollSummary() { collD->print();  }
  void resetColls()       { collD->reset();  }
  //@}

  //! Select particle numbers for interaction
  virtual sKey2Umap generateSelection(pCell* c, sKeyDmap* Fn, 
				      double crm, double tau, int id,
				      double& meanLambda, double& meanCollP, 
				      double& totalNsel);

  //! Compute/retrieve mean molecular weight for this gas
  virtual double molWeight();

};

#endif
