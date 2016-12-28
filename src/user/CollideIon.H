#ifndef _CollideIon_H
#define _CollideIon_H

#include <bitset>
#include <tuple>
#include <array>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include <boost/circular_buffer.hpp>
#include <boost/assign.hpp>
#include <boost/bimap.hpp>

#include "Collide.H"
#include "interactSelect.H"
#include "Elastic.H"
#include "Ion.H"
#include "atomic_constants.H"
#include "InitContainer.H"
#include "AsciiHisto.H"
#include "EnumBitset.H"

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

// STL container pretty-print for std::array
//
template <class T, std::size_t N>
std::ostream& operator<<(std::ostream& o, const std::array<T, N>& arr)
{
  o << '(';
  std::copy(arr.cbegin(), arr.cend(), std::ostream_iterator<T>(o, " "));
  o << "\b\b)";
  return o;
}

// Swap values
//
template <class T>
void zswap(T & x, T & y)
{
  T t = x;
  x   = y;
  y   = t;
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
  typedef std::array<double, 3>        UDDA;
  typedef std::vector<UDDA>            IDVP;
  typedef std::vector<double>          IDVD;
  typedef std::vector<unsigned long>   IDVI;
  //@}

  //@{
  //! Vectors for each type of diagnostic by thread id
  IDVP nn, ne, np, ie, ff, CE, CI, RR, dv;
  IDVD eV_av, eV_max, eV_min;
  IDVI eV_N, eV_10;
  //@}
  
protected:

  //@{
  //! Variables to hold processes collective reductions
  UDDA          nn_s, ne_s, np_s, ie_s, ff_s, CE_s, CI_s, RR_s, dv_s;
  double        eV_av_s, eV_max_s, eV_min_s, eV_N_s, eV_10_s;
  //@}
  
  //! Synchronize a UDDT type (tuple)
  void syncUDDA(UDDA& value) 
  {
    UDDA tmp(value);
    MPI_Reduce(&tmp[0], &value[0], 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }

private:

  UDDA zero;

  UDDA accum(IDVP& v)
  {
    UDDA ret(zero);
    for (auto t : v) {
      for (size_t k=0; k<3; k++) ret[k] += t[k];
    }
    return ret;
  }

public:

  //! Construtor initializes the vectors
  CollisionTypeDiag() 
  {
    zero = {0, 0, 0};

    nn.resize(nthrds, zero);
    ne.resize(nthrds, zero);
    np.resize(nthrds, zero);
    ie.resize(nthrds, zero);
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
    std::fill(nn.begin(), nn.end(), zero);
    std::fill(ne.begin(), ne.end(), zero);
    std::fill(np.begin(), np.end(), zero);
    std::fill(ie.begin(), ie.end(), zero);
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
    nn_s     =  accum(nn);
    ne_s     =  accum(ne);
    np_s     =  accum(np);
    ie_s     =  accum(ie);
    ff_s     =  accum(ff);
    CE_s     =  accum(CE);
    CI_s     =  accum(CI);
    RR_s     =  accum(RR);
    for (auto v : dv) {
      if (std::isnan(v[2]))
	std::cout << "Error in sumUp (before), t=" << tnow 
		  << ", size=" << dv.size() << ", id=" << myid << std::endl;
    }
    dv_s     =  accum(dv);
    if (std::isnan(dv_s[2]))
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
    syncUDDA(nn_s);
    syncUDDA(ne_s);
    syncUDDA(np_s);
    syncUDDA(ie_s);
    syncUDDA(ff_s);
    syncUDDA(CE_s);
    syncUDDA(CI_s);
    syncUDDA(RR_s);
    syncUDDA(dv_s);

    double z;
    MPI_Reduce(&(z=eV_av_s),  &eV_av_s,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&(z=eV_N_s),   &eV_N_s,   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&(z=eV_min_s), &eV_min_s, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&(z=eV_max_s), &eV_max_s, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&(z=eV_10_s),  &eV_10_s,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
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
  std::vector<double> Esum, Elec, Edsp, Efrc, Emas, Epot, delI, delE, clrE, misE;
  double Esum_s, Efrc_s, Emas_s, Elec_s, Edsp_s, Etot_c, Epot_s, delI_s, delE_s, clrE_s, misE_s;

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

  //! Add cell energy
  void addCellEclr(double Eclr, double Emis, int id)
  { clrE[id] += Eclr; misE[id] += Emis; }

  //! Electron energy diagnostics (returns electron KE in momentum frame)
  double addCellElec(pCell* cell, int ue, int id);

  //! Add ionization potential energy
  void addCellPotl(pCell* cell, int id);

};

enum KE_Flags {
  KEpos,			//!< [positive KE]
  ExQ,				//!< [Exact q<1]
  AlgO,				//!< [Orthogonal exact]
  AlgW,				//!< [Weighted exact]
  StdE,				//!< [Standard exact]
  negE,				//!< [negative KE]
  zeroKE,			//!< [zero COM KE]
  momC,				//!< [Momentum cons]
  PreLoss,			//!< [Preloss check]
  Vfac,				//!< [Vfac set]
  MaxVal			//!< Mark the end of the enums
};

template<>
struct EnumTraits<KE_Flags>
{
  static const KE_Flags max = KE_Flags::MaxVal;
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

   <h2> Direct Algorithm </h2>

   Each super particle represents a particular element in a particular
   state with the same total number of atoms or molecules.  Therefore,
   there will be very few "metal" particles and this algorithm is not
   suitable for a simulation were trace elements are important either
   as a diagnostic or energetically (e.g.) as a cooling source.

   <h2> Weight Algorithm </h2>

   In this scheme, each particle carries a statistical weight.  This
   permits a simulation where trace elements and dominant elements may
   have the same number of super particles.  Of course, each super
   particle represents a different number of actual atoms and
   molecules.  Therefore, the collsions are modified with a splitting
   algorithm that breaks the dominant element into a non-interacting
   particle and an interacting particle with the same number of atoms
   or molecules as the trace particle.  The split particle is
   recombined after the interaction to conserve momentum.  The energy
   deficit (because momentum and energy cannot be simultaneously
   conserved in this algorithm) is added to the next non-trace
   interaction.

   <h2> Trace Algorithm </h2>

   In this algorithm, each particle represents all possible atoms,
   molecules and their various ionization and excitation states.  The
   basic strategy is as follows:
   <ol>
   <li> Each superparticle interacts according to the mean scattering
   cross section of all particles in its interaction cell
   <li> The "interesting" inelastic interactions are electron +
   neutral or electron + ion.  Therefore, when two particles interact,
   their internal excitation states are assumed to change in
   proportion to total probability of that interaction divided by the
   number of superparticle interactions that will occur in each cell. 
   </ol>

   <h2> Hybrid Algorithm </h2>

   <ol>
   <li> Each particle is divided into subspecies according to its
   ionization state.
   <li> The electron velocity is separately carried as 3 real
   attributes.  That is, the inertial transport is that of the ion,
   but the interaction energy is computed from the
   electron-ion/neutral interaction, where appropriate.
   <li> Since various electron-ion cross-sections are tail sensitive,
   the cell initialization computes the minimum, mean, and maximum
   ion-electron velocities and energies.  The per-cell cross sections
   are weighed by the mean fraction of each subspecies per cell (meanF
   variable) and the maximum <sigma*vel> for the min, mean, and max
   electron-ion interaction velocity.
   <li> Interaction selection counts are computed between the
   neutrals/ions and electrons for each subspecies from the mean
   super-particle number densities, the (possibly upweighted by
   electron to ion velocity ratio) cross sections, the mean velocities,
   and the time step.
   <li> The interactions are selected as follows:
   <ul>
   <li> The peak subspecies is included for full interaction.  The
   energy loss is weighted by the produce of subspecies ion/neutral
   fraction and electron fraction.  The particles are scattered in
   this interaction, unweighted by subspecies fraction, under the
   assumption that the peak cross section will tend to dominate.
   <li> The same number of non-peak interactions are chosen for
   non-peak subspecies.  The energy lost by inelastic interactions is
   computed with no scattering, downweighted by the ratio of predicted
   probabilities of interaction.
   </ul>

   </ol>

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
  typedef std::pair<double, double>                  CFreturn;
  typedef std::vector<CFreturn>                      CEvector;
  typedef std::pair<speciesKey, speciesKey>          dKey;
  typedef std::pair<dKey, double>                    dKeyD;
  typedef std::map<speciesKey, std::vector<double> > keyCross;
  typedef std::map<speciesKey, std::vector<int>    > keyInter;
  typedef std::map<speciesKey, CEvector>             keyColEx;
  typedef std::list<dKeyD>                           keyWexsList;
  typedef std::map<lQ, CEvector>                     CEmap;
  typedef std::map<lQ, CFreturn>                     CFmap;
  typedef Icont<std::map, int, std::pair<double, unsigned> >
                                                     TypeMap2;
  typedef Icont<std::map, unsigned short, double>    TypeMap;
  typedef Icont<std::map, speciesKey, double>        keyWghts;
  typedef std::pair< std::array<double, 3>,
		     std::array<double, 3> >         CVels;
  typedef std::tuple<double, CEvector, CFreturn>     XStup;
  typedef std::map<NTC::Interact::T, XStup>          InterXS;

  std::vector< std::vector<dKey>   >                 kInter;
  std::vector< std::vector<double> >                 dCross;
  std::vector< InterXS >                             hCross;
  std::vector< std::vector<double> >                 dCfrac;
  std::vector< std::vector<int>    >                 dInter;
  std::vector< keyCross >                            sCross; 
  std::vector< keyWghts >                            meanF, densE;
  std::vector< keyWexsList >                         excessW;
  std::vector< keyInter >                            sInter;
  std::vector< CEvector >                            CE1, CE2;
  std::vector< CFreturn >                            FF1, FF2;
  std::vector< CEmap >                               CEm;
  std::vector< CFmap >                               FFm;

  std::vector< keyColEx >                            kCE;
  std::vector< double >                              meanE, neutF, meanR, meanM;
  std::vector< double >                              numEf, colSc, elecDen, elecDn2;
  std::vector< double >                              kEi, kEe1, kEe2, kEee, clrE, misE;
  std::vector< double >                              Ein1, Ein2, Evel, testKE, nselRat;
  std::vector< double >                              spCrm, spTau, spProb, spEdel;
  std::vector< unsigned >                            spNsel, elecCnt, spNcol;
  std::vector< unsigned >                            elecOvr, elecAcc, elecTot, testCnt;
  std::vector< TypeMap >                             Escat, Etotl;
  std::vector< std::map<unsigned, TypeMap2> >        Italy;
  std::vector< CVels >                               cVels;
  //@}

  //! This defines the various types of interactions used to tag
  //! non-zero cross sections and subsequent interaction processing
  enum InterTypes { 
    any_type    = 0,
    neut_neut   = 1,
    neut_elec   = 2,
    neut_prot   = 3,
    ion_elec    = 4,
    ion_ion     = 5,
    free_free   = 6,
    colexcite   = 7,
    ionize      = 8,
    recomb      = 9,
    elec_elec   = 10,
    neut_neut_1 = 101,
    neut_elec_1 = 102,
    neut_prot_1 = 103,
    ion_elec_1  = 104,
    ion_ion_1   = 105,
    free_free_1 = 106,
    colexcite_1 = 107,
    ionize_1    = 108,
    recomb_1    = 109,
    neut_neut_2 = 201,
    neut_elec_2 = 202,
    neut_prot_2 = 203,
    ion_elec_2  = 204,
    ion_ion_2   = 205,
    free_free_2 = 206,
    colexcite_2 = 207,
    ionize_2    = 208,
    recomb_2    = 209
  };    

  //! Default electron-electron interaction type
  static NTC::Interact::T elecElec;

  //@{
  //! Specify the electron scattering algorithm
  enum ElectronScatter {none, always, classical, limited, fixed};
  typedef std::map<std::string, ElectronScatter> esMapType;
  static ElectronScatter esType;
  static esMapType esMap;
  //@}

  //! Mnemonic for the type of algorithm selected by the user
  enum AlgorithmType { Direct, Weight, Trace, Hybrid };
  const char* AlgorithmLabels[4] = {"Direct", "Weight", "Trace", "Hybrid"};

  //@{
  //! eVEL quantiles
  const std::vector<double> qnt = {0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95};
  std::vector<double> elecV, ionV;
  const size_t bufCap = 16384;
  std::vector< boost::circular_buffer<double> > velER, keER, keIR, momD, crsD;
  ahistoDPtr elecT, elecH, ionsT, ionH, lossH, keEH, keIH, momH, crsH;
  std::map<unsigned short, ahistoDPtr> elecZH, ionZH;
  unsigned Ovr_s, Acc_s, Tot_s, CntE;
  double RhoE, Rho2;
  //@}

  //@{
  //! For histogram electron interaction energy distribution in eV.
  //! Set to false for production (default: false)
  static bool distDiag;
  static bool elecDist;
  std::vector< boost::circular_buffer<double> > elecEV;
  std::vector< boost::circular_buffer<double> > elecEVmin;
  std::vector< boost::circular_buffer<double> > elecEVavg;
  std::vector< boost::circular_buffer<double> > elecEVmax;
  std::vector< boost::circular_buffer<double> > elecEVsub;
  ahistoDPtr elecEVH, elecEVHmin, elecEVHavg, elecEVHmax;
  ahistoDPtr elecEVHsub;
  //@}

  //@{
  //! Electron NTC diagnostic histograms
  typedef std::map<float, ahistoDPtr> NTCeeHist;
  typedef std::map<double, std::vector<double> > eeMap;

  //! Set to false for production (default: false)
  static bool ntcDist;

  //! Desired quantiles
  const std::vector<double> qv = {0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95};

  //! Store NTCdb data
  eeMap ee;
  
  //! Histogram of NTCdb data
  NTCeeHist eeHisto;
  //@}

  //@{
  //! Ionization-recombination check
  const bool IonRecombChk = false;
  std::vector< std::map<speciesKey, double> > ionCHK, recombCHK;
  //@}

  //@{
  //! Collision-limit counter
  typedef std::tuple<unsigned, double, double> clampDat;
  clampDat const clamp0 {0, 0.0, 0.0};
  std::vector<clampDat> clampdat;
  clampDat clampStat;
  //@}

  //@{
  //! Count electron-electron and ion-electron scattering in EPSM
  //! limit for each cell
  std::vector<unsigned> totlES, epsmES;
  unsigned totlES0, epsmES0;
  //@}

  //@{
  //! Count ion-electron scattering in EPSM limit for each cell
  std::vector<unsigned> totlIE, epsmIE;
  unsigned totlIE0, epsmIE0;
  //@}

  //@{
  //! Photon spectrum
  bool use_spectrum, logSpect, wvlSpect;
  double flrSpect, minSpect, maxSpect, delSpect;
  int numSpect;
  
  typedef std::map<int, std::vector<double> >   dNamedVec;
  typedef std::map<int, std::vector<unsigned> > iNamedVec;

  std::set<int>            spectTypes;
  std::vector<dNamedVec>   dSpect, eSpect;
  std::vector<iNamedVec>   nSpect;
  dNamedVec                tSpect, fSpect;
  iNamedVec                mSpect;

  void spectrumAdd(int id, int type, double energy, double weight);
  void spectrumSetup();
  void spectrumGather();
  void spectrumPrint();
  //@}

  //! Compute electron dispersion KE per cell
  std::pair<double, double> computeEdsp(pCell* cell);

public:

  typedef std::set<unsigned short> ZLtype;
  typedef std::map<unsigned short, double> ZWtype;

  //! Maximum number electron scatterings using the "limited"
  //! algorithm
  static unsigned esNum;

  //! Minimum q ratio for electron scatterings using the "limited"
  //! algorithm
  static double esThr;

  //! Coulombic log Lambda value
  static double logL;

  //! Artificially bost inelastic scattering cross section by factor
  static double ieBoost;

  //! Relative error tolerance for warnings (default: 1.0e-6)
  static double tolE;

  //! Multiplicative factor for choosing cooling time step (default: 0.05)
  static double TSCOOL;

  //! Floor KE/deltaE for choosing cooling time step (default: 0.001)
  static double TSFLOOR;

  //! Use EPSM algorithm for electron scattering (default: false)
  static bool ElectronEPSM;


protected:

  //! Process config file
  void processConfig();

  //! List of cross sections (separate for each thread)
  std::vector<sKey2Amap> csections, csectionsH;

  boost::shared_ptr<collDiag> collD;

  //! Set of the elements (atomic number) for direct multi-species
  //! implementation
  ZLtype ZList;

  //@{
  //! Set of the elements (atomic number) for weighted multi-species
  //! implementation
  ZWtype ZWList, ZMList;
  //@}

  //! Designated non-trace element
  unsigned short sFid;

  //! Map of species to attribute position for trace multi-species
  //! implementation
  std::map<speciesKey, int> SpList;

  //! Attribute location for energy excess in "weight" algorithm
  int use_cons;

  //! Electron velocity position in attribute vector
  int use_elec;

  //! Initial position in dattrib of ion-species fractions
  int hybrid_pos;

  //! The algorithm type
  AlgorithmType aType;

protected:
  
  //! Computes (min, mean, max) velocities in each cell
  virtual void cellMinMax(pHOT* const tree, pCell* const cell, int id);

  //@{
  //! Diagnostic for entire cell iteration loop
  std::vector<unsigned> cellEg, cellEb;
  virtual void post_cell_loop(int id);
  //@}

  virtual void 
  initialize_cell(pHOT* const tree, pCell* const c, double rvmax, int id);

  virtual void 
  initialize_cell_dsmc(pHOT* const tree, pCell* const c, sKey2Amap& nsel, 
		       double rvmax, double tau, int id) {}

  virtual void 
  initialize_cell_epsm(pHOT* const tree, pCell* const c, sKey2Amap& nsel, 
		       double rvmax, double tau, int id) 
  {
    std::cout << "NO EPSM implemented!" << std::endl;
    exit(-1);
  }
  
  virtual void 
  finalize_cell(pHOT* const tree, pCell* const c, 
		sKeyDmap* const Fn, double kedsp, double tau, int id);

  
  virtual sKey2Amap& totalCrossSections(int id)
  {    
    return csections[id];
  }

  //! Currently, only used to estimate necessary time step
  sKey2Amap& totalScatteringCrossSections(double crm, pCell* const c, int id);

  double crossSectionDirect(int id, pCell* const c, 
			    Particle* const p1, Particle* const p2, 
			    double crm);
  
  int inelasticDirect(int id, pCell* const c, 
		      Particle* const p1, Particle* const p2, 
		      double *crm);

  double crossSectionWeight(int id, pCell* const c, 
			    Particle* const p1, Particle* const p2, 
			    double crm);
  
  int inelasticWeight(int id, pCell* const c, 
		      Particle* const p1, Particle* const p2, 
		      double *crm);

  double crossSectionTrace(int id, pCell* const c, 
			   Particle* const p1, Particle* const p2, 
			   double crm);

  int inelasticTrace(int id, pCell* const c, 
		     Particle* const p1, Particle* const p2, 
		     double *crm);

  void secondaryScatter(Particle *p);

  double crossSectionHybrid(int id, pCell* const c, 
			    Particle* const p1, Particle* const p2, 
			    double crm, const NTC::Interact::T& ityp);

  int inelasticHybrid(int id, pCell* const c, 
		      Particle* const p1, Particle* const p2, 
		      double *crm, const NTC::Interact::T& itype, double weight);

  virtual double crossSection(int id, pCell* const c, 
			      Particle* const p1, Particle* const p2, 
			      double crm, const NTC::Interact::T& itype)
  {
    if (aType == Direct)
      return crossSectionDirect(id, c, p1, p2, crm);

    if (aType == Weight)
      return crossSectionWeight(id, c, p1, p2, crm);

    if (aType == Trace)
      return crossSectionTrace (id, c, p1, p2, crm);

    if (aType == Hybrid)
      return crossSectionHybrid(id, c, p1, p2, crm, itype);

    return 0.0;
  }


  virtual int 
  inelastic(int id, pCell* const c, Particle* const p1, Particle* const p2, 
	    double *crm, const NTC::Interact::T& itype, double prob)
  {
    if (aType == Direct)
      return inelasticDirect(id, c, p1, p2, crm);

    if (aType == Weight)
      return inelasticWeight(id, c, p1, p2, crm);

    if (aType == Trace)
      return inelasticTrace (id, c, p1, p2, crm);

    if (aType == Hybrid)
      return inelasticHybrid(id, c, p1, p2, crm, itype, prob);

    spNcol[id] += 1;

    return 1;
  }

  virtual void 
  velocityUpdate(Particle* const p1, Particle* const p2, double cr)
  {
    if (aType == Direct)
      Collide::velocityUpdate(p1, p2, cr);

    return;
  }

  virtual bool hasHeatCool() { return true; }

  virtual double getCoolingRate(int id) { return 0; }
  
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
  double nextTime_dbg, vavg_dbg;
  int nCnt_dbg;
  static constexpr double delTime_dbg = 0.001;
  static constexpr int nCel_dbg = 10;
  void write_cross_debug();
  //@}


  virtual sKey2Amap generateSelectionDirect(pCell* const c, sKeyDmap* const Fn,
					    double crm, double tau, int id,
					    double& meanLambda, double& meanCollP, 
					    double& totalNsel);

  virtual sKey2Amap generateSelectionWeight(pCell* const c, sKeyDmap* const Fn,
					    double crm, double tau, int id,
					    double& meanLambda, double& meanCollP, 
					    double& totalNsel);

  virtual sKey2Amap generateSelectionTrace(pCell* const c, sKeyDmap* const Fn, 
					   double crm, double tau, int id,
					   double& meanLambda, double& meanCollP, 
					   double& totalNsel);

  virtual sKey2Amap generateSelectionHybrid(pCell* const c, sKeyDmap* const Fn, 
					    double crm, double tau, int id,
					    double& meanLambda, double& meanCollP, 
					    double& totalNsel);

  virtual NTC::Interact generateSelectionHybridSub(int id, Particle* const p1, Particle* const p2, 
						   NTC::Interact::T& maxT, sKeyDmap* const Fn,
						   double *cr, double tau);

  //! Cache instance of current tree
  pHOT* curTree;

  //@{
  //! Trace species maps for diagnostics
  typedef std::map<speciesKey, int> spMap;
  typedef spMap::iterator spItr;
  void printSpeciesElectrons(std::map<speciesKey, unsigned long>& spec,
			     double T);
  void printSpeciesTrace();
  void printSpeciesColl();
  void gatherSpecies();

  typedef std::map<speciesKey, double> spDMap;
  typedef spDMap::iterator spDItr;
  double tempM, tempE, massE, elecE, molwM, consE, consG, totlE;
  spDMap specM;

  typedef std::tuple<double, double> DTup;
  typedef std::tuple<std::array<DTup, 3>, double> ZTup;

  typedef std::map<unsigned short, ZTup> ZMap;
  ZMap specE, specI;
  std::set<unsigned short> specZ;

  typedef std::array<unsigned long, 2> CollCounts;
  typedef std::map<dKey, CollCounts> CollCountMap;
  const CollCounts ccZ = {{0, 0}};
  std::vector<CollCountMap> collCount;
  //@}

  //! For debugging
  void debugDeltaE(double delE, unsigned short Z, unsigned short C,
		   double dE, double prob, int interFlag);

  //@{
  //! Expensive diagnostics, only use these if you must
  void electronGather();
  void electronPrint(std::ostream&out );
  //@}

  //! Consider commenting these tests out for production
  virtual void auxGather() 
  { 
    mfpCLGather(); 
    electronGather(); 
    spectrumGather();
  }

  //! Consider commenting these tests out for production
  virtual void auxPrint(std::ostream& out) 
  { 
    mfpCLPrint(out); 
    electronPrint(out);
    spectrumPrint();
  }

  //! Compute electron energy in a cell
  double electronEnergy(pCell* const cell, int dbg=-1);

  //! Debugging routine
  void eEdbg();

  const std::string tpaths[4] = {"./test.0", "./test.1", 
				 "./test.2", "./test.3"};
  boost::filesystem::path getcwd;
  std::string hostname;
  unsigned long itp;
  std::map<unsigned long, std::tuple<double, double> > data[2];
  double tot1[2], tot2[2];

  //! MFP algorithm type
  enum class MFP_t
  {
    Ncoll,			//!< [Collision number]
    Direct			//!< [Estimate from xsection]
  };

  //! Map MFP_t to string id for convenience in retrieving the
  //! association
  typedef boost::bimap<MFP_t, std::string> mfp_type;

  //! Assign enum values
  mfp_type MFP_s = boost::assign::list_of<mfp_type::relation>
					   (MFP_t::Ncoll,  "Ncoll")
					   (MFP_t::Direct, "Direct");
  //! Current MFP type
  MFP_t mfptype;
  
  //@{
  //! Hybrid algorithm support routines

  //! Maintain intermediate energy conservation quantities
  class KE_
  {
  private:
    double _i[2], _f[2];

    const std::vector<std::string> names = {
      "[positive KE]",
      "[Exact q<1]",
      "[Orthogonal exact]",
      "[Weighted exact]",
      "[Standard exact]",
      "[negative KE]",
      "[zero COM KE]",
      "[Momentum cons]",
      "[Preloss check]",
      "[Vfac set]",
    };

  public:
    double miss, delta, delE, gamma;
    double totE, kE, dKE, vfac, vrat, o1, o2;


    EnumBitset<KE_Flags> bs;

    KE_() 
    {
      for (auto & v : _i) v = 0.0;
      for (auto & v : _f) v = 0.0;
      miss = delta = delE = gamma = 0.0;
      totE = kE = dKE = vfac = vrat = o1 = o2 = 0.0;
    }

    KE_(double dE) : delE(dE)
    {
      for (auto & v : _i) v = 0.0;
      for (auto & v : _f) v = 0.0;
      miss = delta = gamma = 0.0;
      totE = kE = dKE = vfac = vrat = o1 = o2 = 0.0;
    }
    
    double& i(int i) { return _i[i-1]; }
    double& f(int i) { return _f[i-1]; }

    std::string decode()
    {
      std::string ret;
      for (size_t i=0; i<bs.size(); i++) {
	if (bs[i]) ret += names[i];
      }
      return ret;
    }
  };
  
  enum HybridColl {Ion1=1, Ion2=2, Neutral=4, Scatter=8};

  //! Interaction data structure for simplifying variable passing
  struct InteractData 
  {
    double m1, m2, W1, W2, q, vi2;
    unsigned short Z1, Z2;
    Particle *p1, *p2;

    InteractData(double m1, double m2, double W1, double W2, double q, 
		 unsigned short Z1, unsigned short Z2, Particle* p1, Particle* p2) :
      
      m1(m1), m2(m2), W1(W1), W2(W2), q(q), vi2(0.0), Z1(Z1), Z2(Z2), p1(p1), p2(p2) {}

  };

  //! Compute the interaction for the hybrid weighted method
  void scatterHybrid
  (const InteractData& d, KE_& KE,
   std::vector<double>& v1, std::vector<double>& v2);

  //! Add excess energy to conservation counter
  void deferredEnergyHybrid(const InteractData& d, const double E);

  //! Check energy convervation for hybrid method
  void checkEnergyHybrid
  (const InteractData& d, KE_& KE,
   const std::vector<double>& v1, const std::vector<double>& v2,
   unsigned iType, int id);

  //! Print cross-section weights for debugging
  void meanFdump(int id);

  //! Normalize the ionization state weights for the Hybrid method and
  //! test for correctness
  void normTest(Particle* const p, const std::string& lab);

  //@}

public:
  static double   Nmin;
  static double   Nmax;
  static double   Tmin;
  static double   Tmax;
  static unsigned Nnum;
  static unsigned Tnum;
  static unsigned maxSelA;
  static unsigned maxSelB;
  static unsigned NoDelC;
  static unsigned maxCoul;
  static string   cache;
  static bool     equiptn;
  static bool     scatter;
  static bool     ExactE;
  static bool     AlgOrth;
  static bool     AlgWght;
  static bool     DebugE;
  static bool     collLim;
  static bool     E_split;
  static string   config0;

  //! Constructor
  CollideIon(ExternalForce *force, Component *comp, 
	     double hsDiam, double crossfac, 
	     const std::string& spcMap, int nth=1);

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
  void printCollSummary() { collD->print (); }
  void resetColls()       { collD->reset (); }
  //@}

  //! Select particle numbers for interaction
  virtual sKey2Amap generateSelection(pCell* c, sKeyDmap* Fn, 
				      double crm, double tau, int id,
				      double& meanLambda, double& meanCollP, 
				      double& totalNsel);

  //! Generate subdominate interaction rates
  virtual NTC::Interact generateSelectionSub(int id, Particle* const p1, Particle* const p2, NTC::Interact::T& maxT,
					sKeyDmap* const Fn, double *cr, double tau);


  //! Check whether subdominant interactions are defined (currently
  //! only implemented for hybrid method)
  virtual bool selectSub() 
  { 
    // if (aType == Hybrid) return true;
    return false; 
  }

  //! Compute/retrieve mean molecular weight for this cell
  virtual double molWeight(sCell *c);

  //! Set electron scattering type
  static void setEStype(const std::string& type) 
  {
    esMapType::const_iterator it = esMap.find(type);
    if (it != esMap.end())
      esType = it->second;
    else {
      if (myid==0) {
	std::cout << "CollideIon::setEStype: could not identify type <"
		  << type << ">, defaulting to <always>" << std::endl;
				// Be verbose here for clarity
	std::cout << "--- Available types ---" << std::endl;
	for (auto s : esMap) std::cout << "  " << s.first << std::endl;
	std::cout << "-----------------------" << std::endl;
	esMapType::const_iterator it = esMap.find(type);
      }
      esType = always;
    }
  }

  //! Get electron scattering type
  static std::string getEStype()
  {
    for (auto i : esMap) 
      if (esType == i.second) return i.first;
    return "unknown";
  }

};

#endif
