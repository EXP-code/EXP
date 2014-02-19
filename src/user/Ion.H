/**
GOAL OF OBJECT USAGE:

For each species you want in your gas you can declare it as:

	Ion speciesA("He_1", 2, 1);
	Ion speciesB("C_4", 6, 4);

Then when you are doing the cross sections instead of calculatin the
hard sphere in collide.cc:

	cross1 = speciesA.getDICross();
	cross2 = speciesB.getDICross();

	totalCross = factor1*cross1 + factor2*cross2;

same with the rates
*/


#ifndef _Ion_H
#define _Ion_H

#include <string>
#include <vector>
#include "boost/multi_array.hpp"

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

#define light 29979245800.
#define boltzmann 1.3806504e-16


/** Use a struct for each line of the data. Then, if it needs to be
    switched to multi index it can be easily done so. This also done so
    that there can be a 1D vector with multiple data structures in each
    one.  container to store all the data in 1 line that is relevant from
    the elvlc files. */
struct elvlc_data {
  /*elvlc_data( int lvl, int cnf, std::string des, int s, int l_, std::string l_s, double j,
    int mul, double ecm, double ery, double ecth, double erth): level(lvl), conf(cnf), 
    designation(des), spin(s), l(l_), l_str(l_s), J(j), mult(mul), encm(ecm),
    enry(ery), encmth(ecth), enryth(erth)
    {}*/
  
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

//! container to store all the data in 1 line that is relevant from the elvlc files
/*
  struct wgfa_data{
  //wgfa_data( int l1, int l2, int lambda, int g, int a): lvl1(l1), lvl2(l2), wavelength(lambda), gf(g), A(a)
  //{}

  int lvl1;
  int lvl2;
  double wavelength;
  double gf;
  double A;
  };
**/

struct pe_data
{
  int ngfb;
  int nphot;
  std::vector< double > pe;
};


struct fblvl_data
{
  /*fblvl_data(int lv, std::string cnf, int pn, int L, std::string l_st, int m, double ecm, double ecmt): lvl(lv), conf(cnf), pqn(pn), l(L), l_str(l_st), mult(m), encm(ecm), encmth(ecmt)
    {}*/	
  
  int lvl;
  std::string conf;
  int pqn;
  int l;
  std::string l_str;
  int mult;
  double encm;
  double encmth;
};

struct klgfb_data{
  int n;
  int l;
  std::vector<double> factors;
};

struct gffint_data{
  int ngffint;
  std::vector<double> g2;
  std::vector<double> gffint;
  std::vector<double> s1;
  std::vector<double> s2;
  std::vector<double> s3;
};

struct splups_data{
  int Z;
  int C;
  int i;
  int j;
  int type;
  double gf;
  double delERyd;
  double Const;
  std::vector<double> spline;
};

struct di_data{
  double btf;
  std::vector<double> xspline;
  double ev;
  std::vector<double> yspline;
};

/**
   What is this class for?
*/
class chdata 
{
public:

  //! Constructor
  chdata();
  
  double ipdata[30][30];
  std::vector< std::string > masterNames;
  double abundanceAll[numEle];

  //std::vector< klgfb_data > klgfb;
  //pe_data peData;
  //int nklgfb; //length of the file
  
  //verner cross section data
  //static const unsigned nVern = 465;
  //static const unsigned maxZ = 31; //maxZ = 30 + 1
  //static const unsigned maxNel = 31; //maximum number of electrons
  /*struct verner_data {
    unsigned pqn[maxZ][maxNel];
    unsigned l[maxZ][maxNel];
    double eth[maxZ][maxNel];
    double e0[maxZ][maxNel];
    double sig0[maxZ][maxNel];
    double y_a[maxZ][maxNel];
    double p[maxZ][maxNel];
    double y_w[maxZ][maxNel];
    } verner;
    
    gffint_data gffint; */
  //std::vector<double > ionEqAll;
  
  /*struct ionT_data{
    int ntemp;
    int nele;
    std::vector< double > ionT;
    };*/
  
  // itoh data for free-free relativistic
  // double itoh[30][11][11];
  
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

  //void readKlgfb();
  //void readVern();
  //void readIoneq();
  //void readItoh();
  //void readgffint();
  //void readVerner();
  
  //! Print the ipdata set
  void printIp();

  //!  print functions for the chdata types
  void printMaster();

  //! No idea . . . 
  void printKlgfb();
};

/**
   Structure containing all info about each ion type
 */
class Ion 
{
public:

  //@{
  //! Constructors
  Ion(std::string name, chdata ch);
  Ion(int Z1, int C1,   chdata ch);
  //@}

  //@{
  //! Copy constructors
  Ion(const Ion &);
  Ion();
  //@}
  
  
  std::string eleName; //element name, i.e. "He"
  int Z; //2 = He
  int C; //1 = He I
  
  //! ionization potential in eV
  double ip;
  bool d;

  //double abundance; 
  //abundance of the element in the system (for test, using cosmic.abund)

  //! store as vector since the file contains string and float information
  std::vector< elvlc_data > elvlc; 

  //std::vector< wgfa_data > wgfa; //store as vector since the file contains string and float information
  std::vector< fblvl_data > fblvl; //store as vector since the file contains string and float information
  //std::vector< fblvl_data > ifblvl;
  std::vector< splups_data > splups;	
  std::vector< di_data > diSpline;
  struct di_head{
    int Z;
    int C;
    int nspline;
    int nfac;
    int neav;
  } di_header;		
  
  int nelvlc; //lengths of the files
  int nsplups;
  //int nwgfa; 
  int nfblvl;
  int ndispline;
  //include in other variables here as needed from the files
  //could include variables for the elvlc and wgfa files
  //double fbcross, fbloss, ffcross, ffloss, bbcross, bbloss;
  double ffcross;
  
  std::string MasterName;
  //! Cross sections that are dynamic
  std::vector< std::pair<double, double> > CEcross1j;
  std::vector< std::pair<double, double> > CEcrossCum;
  std::vector< double > dEij_CE;
  std::vector< std::vector < double> > ffDiffCross;
  std::vector< std::vector < double> > ffCumCross;
  std::vector< double > radRecCross;
  std::vector< double > radRecCrossCum;
  double diCross;
  
  //! Grid values for the various cross sections
  std::vector< double > egrid;
  std::vector< double > kgrid;
  int kffsteps;
  int effsteps;
  
  
  void convertName();
  void readelvlc();
  void readwgfa();
  void readfblvl();
  void readSplups();
  void readDi();
  bool isInMasterList(chdata ch, std::string name);
  double getIp(chdata ch);
  double freeFreeCross(chdata ch, double E);

  /** Calculate the differential free-free cross section and return
      the cumulative cross section vector The formula used to
      calculate the cross section is 3BS(a) from Koch&Motz 1959 */
  void freeFreeDifferential(chdata ch);

  //! Get the collisional excitation cross section
  std::vector<std::pair<double, double> > 
  collExciteCross(chdata ch, double p, double Eth, double m);

  //! Ionization cross section
  double directIonCross(chdata ch, double E);

  /** Calculates the differential radiative recombination cross section
      as a function of incoming electron impact energy, and returns the
      vector cumulative cross section array. */
  std::vector<double> radRecombCross(chdata ch, double E);

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
  int getZ(){ return Z; }
  int getC(){ return C; }
  //@}
};

#endif


/****
     Ion::convertName() {
     char ele[2];
     char charge[2];
     
     sscanf(MasterName, "%s_%s", ele, charge);
     for(int i = 0; i < numEle; i++) {
     if(strcmp(ele, eleNameList[i]) == 0) 
     Z = i+1;
     }
     for(int i = 0; i < numCharge; i++) {
     if(strcmp(charge, chargeList[i]) == 0) 
     C = i+1;
     }
     }
     
     you give it C_V, and you get Z = 6, C = 5. need to know how Z and C are related to see if the atom is bare or not
     
****/
