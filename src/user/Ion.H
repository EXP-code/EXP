#ifndef _Ion_H
#define _Ion_H

#include <string>
#include <vector>
#include <set>
#include <map>

#include <boost/multi_array.hpp>

#include <atomic_constants.H>

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

/** Use a struct for each line of the data. Then, if it needs to be
    switched to multi index it can be easily done so. This also done so
    that there can be a 1D vector with multiple data structures in each
    one.  container to store all the data in 1 line that is relevant from
    the elvlc files. */
struct elvlc_data 
{
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
};

struct wgfa_data 
{
  int lvl1;			// lower level
  int lvl2;			// upper level
  double wvl;			// wavelength in Angstroms
  double gf;			// weighted oscillator strength
  double avalue;		// Einstein A coefficient
  std::string pretty1;		// Description (lower level)
  std::string pretty2;		// Description (upper level)
  std::string ref;		// References
};

struct pe_data
{
  int ngfb;
  int nphot;
  std::vector< double > pe;
};


struct fblvl_data
{
  int lvl;
  std::string conf;
  int pqn;
  int l;
  std::string l_str;
  int mult;
  double encm;
  double encmth;
};

struct klgfb_data
{
  int n;
  int l;
  std::vector<double> factors;
};

struct gffint_data
{
  int ngffint;
  std::vector<double> g2;
  std::vector<double> gffint;
  std::vector<double> s1;
  std::vector<double> s2;
  std::vector<double> s3;
};

struct splups_data
{
  unsigned char Z;
  unsigned char C;
  int i;
  int j;
  int type;
  double gf;
  double delERyd;
  double Const;
  std::vector<double> spline;
};

struct di_data
{
  double btf;
  std::vector<double> xspline;
  double ev;
  std::vector<double> yspline;
};

class Ion;

/**
   Store database of atomic data
*/
class chdata 
{

public:

  //! Constructor
  chdata();
  
  double ipdata[30][30];
  std::set<std::string> masterNames;
  double abundanceAll[numEle];

  //! Access to data for individual Ions
  std::map < std::pair<int, int>,  Ion > IonList;

  //! Access to interaction info
  std::map< std::pair<unsigned short, unsigned short>, double > interInfo;

  //! Not sure what this does . . . looks like cruft.  Remove??
  std::map< int, double > Ni;

  /** Read in the master list to store to be able to check if elements
      are in it */
  void readMaster();

  /** Get the ionization potential data

      If you want to get the ip of any Z, C, you call it as
      ipdata[Z-1][C-1-(int)die] 
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
public:

  //@{
  //! Constructors
  Ion (std::string name, chdata ch);
  Ion (unsigned char Z1, unsigned char C1, chdata ch);
  //@}

  //@{
  //! Copy constructors
  Ion (const Ion &);
  Ion ();
  //@}
  
  
  //! element name, i.e. "He"
  std::string eleName; 

  //! nuclear charge, e.g. 2 = He
  unsigned char Z;

  //! ionization state, e.g. 1 = He I
  unsigned char C; 
  
  //! ionization potential in eV
  double ip;
  bool d;

  //double abundance; 
  //abundance of the element in the system (for test, using cosmic.abund)

  //! Store by level
  typedef std::map<int, elvlc_data> elvlcType; 
  elvlcType elvlc; 

  //! Store by level pair
  typedef std::pair<int, int> ZC;
  typedef std::map<ZC, wgfa_data> wgfaType;
  wgfaType wgfa;

  //! Store by level
  typedef std::map<int, fblvl_data> fblvlType; 
  fblvlType fblvl; 

  // fblvlType ifblvl;

  typedef std::vector< splups_data > splupsType;	
  splupsType splups;	

  typedef std::vector< di_data > diSplineType;
  diSplineType diSpline;

  struct di_head
  {
    unsigned char Z;
    unsigned char C;
    int nspline;
    int nfac;
    int neav;
  } di_header;		
  
  //! include in other variables here as needed from the files
  //! could include variables for the elvlc and wgfa files
  // double fbcross, fbloss, ffcross, ffloss, bbcross, bbloss;

  //! Free-free cross section
  double ffcross;
  
  //! Element name
  std::string MasterName;

  //@{
  //! Cross sections that are dynamic
  std::vector< std::pair<double, double> > CEcross1j;
  std::vector< std::pair<double, double> > CEcrossCum;
  std::vector< double > dEij_CE;
  std::vector< std::vector < double> > ffDiffCross;
  std::vector< std::vector < double> > ffCumCross;
  std::vector< double > radRecCross;
  std::vector< double > radRecCrossCum;
  double diCross;
  //@}
  
  //@{
  //! Grid values for the various cross sections
  std::vector< double > egrid;
  std::vector< double > kgrid;
  int kffsteps;
  int effsteps;
  //@}
  
  
  void convertName();
  void readelvlc();
  void readwgfa();
  void readfblvl();
  void readSplups();
  void readDi();
  //! Is a species is in the database?
  bool isInMasterList(chdata& ch, const std::string name) 
  { 
    return ch.masterNames.find(name) != ch.masterNames.end();
  }
  double getIp(chdata ch);
  double freeFreeCross(chdata ch, double E);

  /** Calculate the differential free-free cross section and return
      the cumulative cross section vector The formula used to
      calculate the cross section is 3BS(a) from Koch&Motz 1959 */
  void freeFreeDifferential(chdata ch);


  //! Return type for collisional excitation cross-section data
  typedef std::vector< std::pair<double, double > > collType;

  //! Get the collisional excitation cross section
  collType collExciteCross(chdata ch, double E);

  //! Ionization cross section
  double directIonCross(chdata ch, double E);

  //@{
  /** Calculates the differential radiative recombination cross section
      as a function of incoming electron impact energy, and returns the
      vector cumulative cross section array. */
  std::vector<double> radRecombCross       (chdata ch, double E);
  std::vector<double> radRecombCrossMewe   (chdata ch, double E);
  std::vector<double> radRecombCrossSpitzer(chdata ch, double E);
  std::vector<double> radRecombCrossKramers(chdata ch, double E);
  std::vector<double> radRecombCrossMilne  (chdata ch, double E);
  //@}

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
  // get functions
  unsigned char getZ() const { return Z; }
  unsigned char getC() const { return C; }
  //@}
};

#endif
