#ifndef _Ion_H
#define _Ion_H

#include <string>
#include <vector>
#include <set>
#include <map>

#include <boost/shared_ptr.hpp>
#include <boost/multi_array.hpp>

#include <atomic_constants.H>
#include <TopBase.H>
#include <Cspline.H>
#include <Verner.H>
#include <KarzasLatter.H>

#include <gaussQ.h>

// Structures for caching the Cspline instances for each Ion instance
//
typedef Cspline<double, double>      CsplineD2;
typedef boost::shared_ptr<CsplineD2> CsplD2Ptr;
typedef std::map<int, CsplD2Ptr>     CacheSplList;


static std::string eleNameList[] = 
  {"h","he","li","be","b","c","n","o","f","ne","na", 
   "mg","al","si","p","s","cl","ar","k","ca","sc","ti", 
   "v","cr","mn","fe","co","ni","cu","zn",
   "ga","ge","as","se","br","kr"};

const int numEle = 36;

static std::string chargeList[] = 
  {"I", "II", "III", "IV", " V", "VI", "VII", "VIII", "IX", "X", 
   "XI", "XII", "XIII", "XIV", "XV"};

const int numCharge = 15;

const int numMaster = 288;

/** Base class for reading CHIANTI data */
class chianti_data
{
protected:
  
  void sync_string(std::string &s);
  void sync_vector(std::vector<double> &v);

public:

  virtual void synchronize() = 0;

};

/** Use a struct for each line of the data. Then, if it needs to be
    switched to multi index it can be easily done so. This also done so
    that there can be a 1D vector with multiple data structures in each
    one.  container to store all the data in 1 line that is relevant from
    the elvlc files. */

class elvlc_data : public chianti_data
{
public:

  int level;
  int conf;
  std::string designation;
  int spin;
  int l;
  std::string l_str;
  double J;
  int mult;
  double encm;
  double enry;
  double encmth;
  double enryth;

  void synchronize();
};

class wgfa_data  : public chianti_data
{
public:

  int lvl1;			// lower level
  int lvl2;			// upper level
  double wvl;			// wavelength in Angstroms
  double gf;			// weighted oscillator strength
  double avalue;		// Einstein A coefficient
  std::string pretty1;		// Description (lower level)
  std::string pretty2;		// Description (upper level)
  std::string ref;		// References

  void synchronize();
};

class pe_data : public chianti_data
{
  int ngfb;
  int nphot;
  std::vector< double > pe;

  void synchronize();
};


class fblvl_data : public chianti_data
{
public:

  int lvl;
  std::string conf;
  int pqn;
  int l;
  std::string l_str;
  int mult;
  double encm;
  double encmth;

  void synchronize();
};

class klgfb_data : public chianti_data
{
public:

  int n;
  int l;
  std::vector<double> factors;

  void synchronize();
};

class gffint_data : public chianti_data
{
public:

  int ngffint;
  std::vector<double> g2;
  std::vector<double> gffint;
  std::vector<double> s1;
  std::vector<double> s2;
  std::vector<double> s3;

  void synchronize();
};

class splups_data : public chianti_data
{
public:

  unsigned char Z;
  unsigned char C;
  int i;
  int j;
  int type;
  double gf;
  double delERyd;
  double Const;
  std::vector<double> spline;

  void synchronize();
};

class di_data : public chianti_data
{
public:

  double btf;
  std::vector<double> xspline;
  double ev;
  std::vector<double> yspline;

  void synchronize();
};


//! Pair of unsigned shorts for indexing (Z, C) and (l1, l2) index pairs
typedef std::pair<unsigned short, unsigned short> lQ;

class Ion;

/**
   Store database of atomic data
*/
class chdata 
{

protected:

  //! Gauss-Laguerre structure for integration of Maxwell-Boltzmann
  //! distribution
  boost::shared_ptr<LaguQuad> Lagu;

  //! Gauss-Legendre structure for integration of Maxwell-Boltzmann
  //! distribution with boundaries
  boost::shared_ptr<LegeQuad> Lege;

public:

  //! Constructor
  chdata();
  
  //! Create database of requested ions
  typedef std::set<unsigned short> ZLtype;
  void createIonList(const ZLtype& ZList);

  std::map<lQ, double> ipdata;
  std::set<std::string> masterNames;
  double abundanceAll[numEle];

  //! TopBase data
  boost::shared_ptr<TopBase> tb;

  //! Access to data for individual Ions
  std::map <lQ,  boost::shared_ptr<Ion> > IonList;

  //! Not sure what this does . . . looks like cruft.  Remove??
  std::map< int, double > Ni;

  //! Verner-Yakovlev cross-section data
  VernerData VernerXC;
  
  //! Gaunt factor data
  KLGFdata radGF;
  
  /** Read in the master list to store to be able to check if elements
      are in it */
  void readMaster();

  /** Get the ionization potential data

      If you want to get the ip of any Z, C, you call it as
      ipdata[lQ(Z, C-(int)die)] 
  */
  void readIp();

  /** read in the abundance file, in this situation, just for test using
      the cosmic.abund file. Can later put in a multidimensional array
      to allow for all the abundance files */
  void readAbundanceAll();

  //! List the ionization potential data
  void printIp();

  //! List names of all species to stdout
  void printMaster();

  //! Read the radiative cross section data
  void readVerner();

  //! Read the radiative recombination Gaunt factor data
  void readRadGF();

  //@{
  /** Compute ionization-recombination thermal equilibrium (for ICs
      and diagnostics */
  std::map<unsigned short, double> 
  fraction(unsigned short Z, double Temperature, int norder=20);

  std::map<unsigned short, double> 
  fraction(unsigned short Z, double Temperature, double Emin, double Emax, 
	   int norder=40);

  std::map<unsigned short, std::vector<double> >
  recombEquil(unsigned short Z, double Temperature, int norder=20);

  std::map<unsigned short, std::vector<double> >
  recombEquil(unsigned short Z, double Temperature, double Emin, double Emax, 
	      int norder=40, bool use_log=false);
  //@}

};

/**
   Structure containing all info about each ion type

   For each species you want in your gas you can declare it as:

        Ion speciesA("He_1", 2, 1);
        Ion speciesB("C_4", 6, 4);

   Then when you are doing the cross sections instead of calculating the
   hard sphere in collide.cc:

        cross1 = speciesA.getDICross();
        cross2 = speciesB.getDICross();

        totalCross = factor1*cross1 + factor2*cross2;

   same with the rates . . .
*/
class Ion 
{
private:

  chdata* ch;

  //@{
  //! Cache the spline instances for efficiency
  CacheSplList splUps;
  CacheSplList diSpln;
  //@}
  
protected:

  //
  // Physical constants
  //
				// Prefactor A / pi^2 (from Greene, 1959)
  static constexpr double A       = 5.72791733e-08/(M_PI*M_PI);

				// Classical radius of the electron in nm
  static constexpr double r0      = 2.81794044e-13 * 1.0e+07;

				// Bohr radius in nm
  static constexpr double a0      = 0.0529177211;
  
				// Rydberg energy in eV
  static constexpr double RydtoeV = 13.60569253;

				// eV in Rydberg energy
  static constexpr double eVtoRyd = 1.0/13.60569253;

				// 1 eV in ergs
  static constexpr double eVtoErg = 1.60217733000001e-12;

				// value of h-bar * c in eV*nm
  static constexpr double hbc     = 197.327;

				// Planck's constant in eV*s
  static constexpr double heV     = 4.135667222637445e-15;

				// Fine structure constant
  static constexpr double alpha0  = 7.2973525698e-03;


				// Thomson scattering cross-section in nm^2
  static constexpr double thom    = 6.652459263719353e-11;

public:

  //@{
  //! k-grid spacing in inverse nm (log base 10)
  static double kmin;		// default: -9.0
  static double kmax;		// default:  2.0
  static double kdel;		// default:  1.0
  //@}

  //@{
  //! nu-grid spacing relative to ground state level for photo ionization
  static double numin;		// default:  1.0
  static double numax;		// default:  3.0
  static double nudel;		// default:  0.2
  //@}

  //@{
  //! Structures for user-friendly specification of recombination types
  enum    RR_Type { mewe, topbase, kramers, spitzer, verner };
  typedef std::map<std::string, RR_Type> RR_Map;
  typedef std::map<RR_Type, std::string> RR_Lab;
  //@}

  //@{
  /** Structures for user-friendly specification of photoionization
      background spectra */
  enum    IB_Type { none, uvIGM };
  typedef std::map<std::string, IB_Type> IB_Map;
  typedef std::map<IB_Type, std::string> IB_Lab;
  //@}

  //! Use Verner et al. fortran subroutine (default: true)
  static bool use_VFKY;

  //! Recombine into ground state only for Verner cross sections
  static bool gs_only;

protected:

  //@{
  //! Static varibles for setting radiative recombination cross-section type
  static RR_Map  rr_map;
  static RR_Lab  rr_lab;
  static RR_Type rr_type;
  //@}

  //@{
  //! Static varibles for setting radiative ionization background type
  static IB_Map  ib_map;
  static IB_Lab  ib_lab;
  static IB_Type ib_type;
  //@}

  //@{
  //! For photoionizing background
  
  //! Grid created
  bool IBinit;

  //! Creation routine
  void IBcreate();

  //! Partial rate grid
  std::vector<double> IBrate;

  //! Total rate
  double IBtotl;

  //! Normalized cumulative rate grid for selecting recoil energy
  std::vector<double> IBcum;
  //@}

public:

  //@{
  //! Constructors
  Ion (std::string name, chdata* ch);
  Ion (unsigned short Z1, unsigned short C1, chdata* ch);
  //@}

  //@{
  //! Copy constructors
  Ion (const Ion &);
  Ion ();
  //@}
  
  //! element name, i.e. "He"
  std::string eleName; 

  //! nuclear charge, e.g. 2 = He
  unsigned short Z;

  //! ionization state, e.g. 1 = He I
  unsigned short C; 
  
  //! ionization potential in eV
  double ip;
  bool d;

  //! Store by level
  typedef std::map<int, elvlc_data> elvlcType; 
  elvlcType elvlc; 

  //! Store by level pair
  typedef std::map<lQ, wgfa_data> wgfaType;
  wgfaType wgfa;

  //! Store by level
  typedef std::map<int, fblvl_data> fblvlType; 
  fblvlType fblvl; 

  // fblvlType ifblvl;

  typedef std::vector< splups_data > splupsType;	
  splupsType splups;	

  typedef std::vector< di_data > diSplineType;
  diSplineType diSpline;

  class di_head
  {
  public:

    unsigned short Z;
    unsigned short C;
    int nspline;
    int nfac;
    int neav;

    void synchronize();

  } di_header;		
  
  //! include in other variables here as needed from the files
  //! could include variables for the elvlc and wgfa files
  // double fbcross, fbloss, ffcross, ffloss, bbcross, bbloss;

  //! Free-free cross section
  double ffcross;
  
  //! Element name
  std::string MasterName;

  //@{
  //! Grid values for the various cross sections
  std::vector< double > egrid;	// Currently initialized but unused (remove?)
  std::vector< double > kgrid;
  std::vector< double > kgr10;
  std::vector< double > nugrid;
  int kffsteps;
  int effsteps;
  int nusteps;
  //@}
  
  void convertName();
  void readelvlc();
  void readwgfa();
  void readfblvl();
  void readSplups();
  void readDi();

  //! Is a species is in the database?
  bool isInMasterList(const std::string name) 
  { 
    return ch->masterNames.find(name) != ch->masterNames.end();
  }

  /** Calculate the differential free-free cross section and return
      the cumulative cross section vector The formula used to
      calculate the cross section is 3BN(a) from Koch & Motz, 1959 

      ** Non-relativistic limit

      First value is cross section, second is photon energy
  */
  std::pair<double, double> freeFreeCross(double E, int id);

  //! Return type for collisional excitation cross-section data
  typedef std::vector< std::pair<double, double> > collType;

  //! Get the collisional excitation cross section
  collType collExciteCross(double E, int id);

  //! Ionization cross section
  double directIonCross(double E, int id);

  //@{

  /** Calculates the differential radiative recombination cross section
      as a function of incoming electron impact energy, and returns the
      vector cumulative cross section array. */

  std::vector<double> radRecombCross       (double E, int id);
  std::vector<double> radRecombCrossMewe   (double E, int id);
  std::vector<double> radRecombCrossSpitzer(double E, int id);
  std::vector<double> radRecombCrossKramers(double E, int id);
  std::vector<double> radRecombCrossTopBase(double E, int id);
  std::vector<double> radRecombCrossVerner (double E, int id);
  //@}

  /** Ground state photoionization cross section using Verner (for
      debugging */
  std::vector<double> photoIonizationCross (double E, int id);

  //! Photoionizing background flux (Enu in eV)
  double photoBackground(double Enu);

  //! Photoionization rate
  std::pair<double, double> photoIonizationRate();

  /** Radiative recombination cross section using Verner (for
      debugging */
  std::vector< std::tuple<int, double> > recombCrossV (double E, int id);

  //! Calculate the QRP value as in Fontes, Sampson, Zhang 1999
  double qrp(double u);
  
  //@{

  //!print functions mainly used in debugging

  void printInfo();
  void printelvlc();
  void printwgfa();
  void printfblvl();
  //@}
  
  //@{
  // Return constant data
  unsigned short getZ()             const { return Z; }
  unsigned short getC()             const { return C; }
  double         getIP(const lQ& Q) const { return ch->ipdata[Q]; }
  //@}

  //! Set cross-section type for all instances
  static int setRRtype(const std::string& rr) 
  {
    if (rr_map.find(rr) == rr_map.end()) {
      rr_type = Ion::mewe;	// Default type
      return 1;			// Failure
    } else { 			// Set request
      rr_type = rr_map[rr];
      return 0;			// Success
    }
  }

  //! Get cross-section type
  static const std::string& getRRtype()
  {
    return rr_lab[rr_type];
  }

  //! Set ionizing background type for all instances
  static int setIBtype(const std::string& ib) 
  {
    if (ib_map.find(ib) == ib_map.end()) {
      ib_type = Ion::none;	// Default type
      return 1;			// Failure
    } else { 			// Set request
      ib_type = ib_map[ib];
      return 0;			// Success
    }
  }

  //! Get cross-section type
  static const std::string& getIBtype()
  {
    return ib_lab[ib_type];
  }

};

typedef boost::shared_ptr<Ion> IonPtr;

#endif
