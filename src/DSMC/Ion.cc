#include <boost/algorithm/string.hpp>
#include <string>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <tuple>
#include <map>

#include <localmpi.h>

#include "Ion.H"
#include "interactSelect.H"

// For setting cross-section type (c++-11 initialization style)
//
Ion::RR_Map Ion::rr_map = {
  { "Mewe",    Ion::mewe    },
  { "TopBase", Ion::topbase },
  { "Kramers", Ion::kramers },
  { "Spitzer", Ion::spitzer },
  { "Verner",  Ion::verner  }
};

// For printing cross-section type (c++-11 initialization style)
//
Ion::RR_Lab Ion::rr_lab = {
    { Ion::mewe,    "Mewe",   },
    { Ion::topbase, "TopBase" },
    { Ion::kramers, "Kramers" },
    { Ion::spitzer, "Spitzer" },
    { Ion::verner,  "Verner"  }
};

Ion::RR_Type Ion::rr_type = Ion::verner;

// For setting ionization background type (c++-11 initialization style)
//
Ion::IB_Map Ion::ib_map = {
  { "None",    Ion::none    },
  { "uvIGM",   Ion::uvIGM   }
};

// For printing ionization background type (c++-11 initialization style)
//
Ion::IB_Lab Ion::ib_lab = {
    { Ion::none,    "None",   },
    { Ion::uvIGM,   "uvIGM"   }
};

Ion::IB_Type Ion::ib_type = Ion::none;

bool Ion::use_VFKY  = true;
bool Ion::no_verner = true;

// Free-free grid
//
double Ion::kmin    = -10.0;
double Ion::kmax    =   8.0;
double Ion::kdel    =   0.01;

// Photo-ionization grid
//
double Ion::numin   = 1.0;
double Ion::numax   = 3.0;
double Ion::nudel   = 0.05;

// Ionization background
//
double Ion::HandM_coef  = 1.5e-22;
double Ion::HandM_expon = -0.5;

// Energy grids (in eV)
//
bool   Ion::useFreeFreeGrid  = true;
bool   Ion::useRadRecombGrid = false;
bool   Ion::useExciteGrid    = true;
bool   Ion::useIonizeGrid    = true;
bool   Ion::GridDebug        = false; // Set to true for debugging
int    Ion::GridReport       = 10000; // Used for debugging only
double Ion::EminGrid         = 0.001; // eV
double Ion::EmaxGrid         = 2000.; // eV
double Ion::DeltaEGrid       = 0.5;  // eV

// Chianti element list
//
std::map<unsigned short, std::string> chElems {{1, "h"}, {2, "he"}, {3, "li"}, {4, "be"}, {5, "b"}, {6, "c"}, {7, "n"}, {8, "o"}, {9, "f"}, {10, "ne"}, {11, "na"}, {12, "mg"}, {13, "al"}, {14, "si"}, {15, "p"}, {16, "s"}, {17, "cl"}, {18, "ar"}, {19, "k"}, {20, "ca"}, {21, "sc"}, {22, "ti"}, {23, "v"}, {24, "cr"}, {25, "mn"}, {26, "fe"}, {27, "co"},	{28, "ni"}, {29, "cu"}, {30, "zn"} };




//
// Convert the master element name to a (Z, C) pair
//
void Ion::convertName() 
{
  std::string ele;
  std::string charge;
  
  std::vector<std::string> v;
  std::string die = "d";	// Set to dielectronic
  size_t      isd;
  
  // split the name up into its element ab. and charge
  // In C, this would be: sscanf(MasterName, "%s_%s", ele, charge);
  //
  boost::split(v, MasterName, boost::is_any_of("_") );

  eleName = v[0];
  
  // get the Z value for the element by looking up through table
  for (int i = 0; i < numEle; i++) {
    if (boost::iequals(v[0], eleNameList[i])) {
      Z = i+1; break;
    }
  }
  
  // filter out the d after the charge if it is a dielectronic
  isd = v[1].find(die);
  if (isd!=string::npos) {
    std::cout << "FOUND" << std::endl;
    d = true;
				// expecting it to be at the end
    v[1].replace(v[1].find(die), die.length(), "\0");
  }
  else {
    d = false;
  }
  // get the charge value
  C = atoi(v[1].c_str());
}

//
// Convert a given Z,C pair into a master name string
//
std::string Ion::ZCtoName(unsigned char Z, unsigned char C) 
{
  std::stringstream ss;
  ss << eleNameList[Z-1] << "_" << static_cast<unsigned>(C);
  return ss.str();
}

Ion::cFile Ion::chiantiFile(std::string suffix)
{
  std::pair<unsigned char, std::string> ret {0, ""};

  if (myid==0) {

    static bool first = true;
    static char * val;

    if (first) {
      if ( (val = getenv("CHIANTI_DATA")) == 0x0) {
	std::cout << "Could not find CHIANTI_DATA environment variable"
		  << " . . . exiting" << std::endl;
	ret.first = 2;
      }
    }

    if (ret.first == 0) {

      std::string MasterNameT = ZCtoName(Z, C);

      std::string fileName(val);

      fileName.append("/");
      fileName.append(eleName); 
      fileName.append("/");
      fileName.append(MasterName); 
      fileName.append("/"); 
      fileName.append(MasterName);
      fileName.append("." + suffix);
      
      std::string inLine;
      std::ifstream file(fileName.c_str());
      
      if (file.is_open()) {
	ret.second = fileName;
      } else {
	ret.first = 1;
      }
    }
  }

  return ret;
}


/** 
    Functions to read in the elvlc and wgfa files if they are found
    since they contain both string and float information, all
    information will be stored as strings for now could eventually
    just reject string data and switch to floats if string data isn't
    wanted
*/
void Ion::readelvlc() 
{
  unsigned char nOK = 0;
  
  if (myid==0) {

    cFile ret = chiantiFile("elvlc");
    
    nOK = ret.first;

    if (nOK == 0) {

      elvlc_data e;
      std::string inLine;
      std::ifstream elvlcFile(ret.second);
      
      if (elvlcFile.is_open()) {

	while (elvlcFile.good()) {
	  
	  std::vector <std::string> v(8);
	  getline(elvlcFile, inLine);

	  /*
	    1 i7    LVL   Level index
	    2 a30   CONF  Configuration description
	    3 a5    LABEL Level label string
	    4 i5    2S+1  Spin multiplicity
	    5 a5    L     Orbital angular momentum
	    6 f5.1  J     Total angular momentum
	    7 f15.3 Eo    Observed energy (cm−1)
	    8 f15.3 Eb    ‘Best-guess’ theoretical energy (cm−1)
	  */


	  // Split according to fortran formats field widths
	  //
	  int c = 0;
	  v[0] = inLine.substr(c,  7); c += 7;
	  if (atoi(v[0].c_str()) == -1) break;
	  v[1] = inLine.substr(c, 30); c += 30;
	  v[2] = inLine.substr(c,  5); c += 5;
	  v[3] = inLine.substr(c,  5); c += 5;
	  v[4] = inLine.substr(c,  5); c += 5;
	  v[5] = inLine.substr(c,  5); c += 5;
	  v[6] = inLine.substr(c, 15); c += 15;
	  v[7] = inLine.substr(c, 15); c += 15;

	  e.level       = atoi(v[0].c_str());
	  e.designation = v[1];
	  e.spin        = atoi(v[2].c_str());
	  e.l_str       = v[3];
	  e.J           = atof(v[4].c_str());
	  e.encm        = atof(v[5].c_str());
	  e.enry        = atof(v[6].c_str());
	  e.encmth      = atof(v[7].c_str());
	  
	  elvlc[e.level] = e;
	}
	elvlcFile.close();
      }
      else {
	if (nOK==1)
	  std::cout << "Cannot find file: " << ret.second << std::endl;
	if (nOK==2)
	  std::cout << "Cannot find CHIANTI database"     << std::endl;
      }
    }
  }

  MPI_Bcast(&nOK, 1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  if (nOK)  {
    if (myid == 0) {
      std::cerr << "Ion::readelvlc: problem reading CHIANTI files "
		<< "for Z=" << Z << " C=" << C << std::endl;
    }
    MPI_Finalize();
    exit(44);
  } 

  unsigned number = elvlc.size();
  MPI_Bcast(&number, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

  if (myid==0) {
    elvlcType::iterator it = elvlc.begin();
    for (unsigned i=0; i<number; i++) {
      it->second.synchronize();
      it++;
    }
  } else {
    elvlc.clear();
    for (unsigned i=0; i<number; i++) {
      elvlc_data e;
      e.synchronize();
      elvlc[e.level] = e;
    }
  }
}

//
// Reads the CHIANTI .wgfa file into a structure.
//
// CHIANTI wgfa files contain radiative decay rates and, in some cases
// autoionization rates and/or two photon transitions.
//
void Ion::readwgfa() 
{
  unsigned char nOK = 0;

  if (myid==0) {

    cFile ret = chiantiFile("wgfa");
    
    nOK = ret.first;

    if (nOK == 0) {

      std::string inLine;
      wgfa_data w;

      std::ifstream wgfaFile(ret.second);
      
      if (wgfaFile.is_open()) {
	
	while (wgfaFile.good()) {
	  
	  std::vector <std::string> v(6);
	  getline(wgfaFile, inLine);

	  /*
	    1. index of the lower energy level, format=i5 (consistent
	       with the ordering in the .elvlc file)

	    2. index of the upper energy level, format=i5

	    3. wavelength in Angstroms, format=f15. If the wavelength
	       does not connect 2 observed energy levels, the
	       wavelength is given as a negative number. Generally
	       wavelengths are given to three decimal places, with
	       exceptions for very short and very long wavelength
	       transitions. Two-photon transitions are given a zero
	       wavelength. Also, the corrections to the inner-shell
	       levels above the ionization limit due to autoionization
	       are included as radiationless transitions to the ground
	       state.

	    4. gf value (weighted oscillator strength), format=e15

	    5. A-value, format=e15. Generally this will be the
	       radiative decay rate (or A-value) in units of s-1. It
	       can also be the autoionization rate or the two-photon
	       decay rate, both in units of s-1.

	    6. In some cases additional columns have extra information
	       on the transition.
	  */

	  // Split according to fortran formats field widths
	  //
	  int c = 0;
	  v[0] = inLine.substr(c,  5); c += 5;
	  if (atoi(v[0].c_str()) == -1) break;
	  v[1] = inLine.substr(c,  5); c += 5;
	  v[2] = inLine.substr(c, 15); c += 15;
	  v[3] = inLine.substr(c, 15); c += 15;
	  v[4] = inLine.substr(c, 15); c += 15;
	  v[5] = inLine.substr(c, 15); c += 15;

	  w.lvl1    = atoi(v[0].c_str());
	  w.lvl2    = atoi(v[1].c_str());
	  w.wvl     = atof(v[2].c_str());
	  w.gf      = atof(v[3].c_str());
	  w.avalue  = atof(v[4].c_str());
	  
	  wgfa[lQ(w.lvl1, w.lvl2)] = w;
	}
	wgfaFile.close();
      }
      else {
	if (nOK==1)
	  std::cout << "Cannot find file: " << ret.second << std::endl;
	else
	  std::cout << "Cannot find CHIANTI database"     << std::endl;
      }
    }
  }

  MPI_Bcast(&nOK, 1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  if (nOK)  {
    if (myid == 0) {
      std::cerr << "Ion::readwgfa: problem reading CHIANTI files "
		<< "for Z=" << Z << " C=" << C << std::endl;
    }
    MPI_Finalize();
    exit(43);
  } 
   
   unsigned number = wgfa.size();
   MPI_Bcast(&number, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

  if (myid==0) {
    wgfaType::iterator it = wgfa.begin();
    for (unsigned i=0; i<number; i++) {
      it->second.synchronize();
      it++;
    }
  } else {
    wgfa.clear();
    for (unsigned i=0; i<number; i++) {
      wgfa_data w;
      w.synchronize();
      wgfa[lQ(w.lvl1, w.lvl2)] = w;
    }
  }

}

//
// Read energy levels and free-bound data found in the CHIANTI
// database
//
void Ion::readfblvl() 
{
  unsigned char nOK = 0;

  if (myid==0) {

    cFile ret = chiantiFile("fblvl");
    
    nOK = ret.first;

    if (nOK == 0) {
  
      std::string inLine;
      std::ifstream fblvlFile(ret.second);
    
      fblvl_data f;
      if (fblvlFile.is_open()) {
	
	while (fblvlFile.good()) {

	  std::vector <std::string> v(8);
	  getline(fblvlFile, inLine);

	  // Split according to fortran formats field widths
	  //
	  int c = 0;
	  v[0] = inLine.substr(c,  5); c += 5;
	  if (atoi(v[0].c_str()) == -1) break;
	  v[1] = inLine.substr(c, 20); c += 20;
	  v[2] = inLine.substr(c,  5); c += 5;
	  v[3] = inLine.substr(c,  5); c += 5;
	  v[4] = inLine.substr(c,  3); c += 3;
	  v[5] = inLine.substr(c, 20); c += 5;
	  v[6] = inLine.substr(c, 20); c += 15;
	  v[7] = inLine.substr(c, 15); c += 15;

	  f.lvl      = atoi(v[0].c_str());
	  f.conf     = v[1];
	  f.pqn      = atoi(v[2].c_str());
	  f.l        = atoi(v[3].c_str());
	  f.l_str    = v[4];
	  f.mult     = atoi(v[5].c_str());
	  f.encm     = atof(v[6].c_str());
	  f.encmth   = atof(v[7].c_str());
	  
	  fblvl[f.lvl] = f;
	}
	fblvlFile.close();
      }
      else {
	if (nOK==1)
	  std::cout << "Cannot find file: " << ret.second << std::endl;
	else
	  std::cout << "Cannot find CHIANTI database"     << std::endl;
      }
    }
  }

  MPI_Bcast(&nOK, 1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  if (nOK)  {
    if (myid == 0) {
      std::cerr << "Ion::readfblvl: problem reading CHIANTI files "
		<< "for Z=" << Z << " C=" << C << std::endl;
    }
    MPI_Finalize();
    exit(44);
  } 

  unsigned number = fblvl.size();
  MPI_Bcast(&number, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

  if (myid==0) {
    fblvlType::iterator it = fblvl.begin();
    for (unsigned i=0; i<number; i++) {
      it->second.synchronize();
      it++;
    }
  } else {
    fblvl.clear();
    for (unsigned i=0; i<number; i++) {
      fblvl_data f;
      f.synchronize();
      fblvl[f.lvl] = f;
    }
  }

}

//
// Read in the spline fits to the collision strengths from the CHIANTI
// database
//
void Ion::readSplups() 
{
  unsigned char nOK = 0;

  if (myid==0) {

    cFile ret = chiantiFile("scups");
    
    nOK = ret.first;

    if (nOK == 0) {

      std::string inLine;
      splups_data s;

      std::ifstream sFile(ret.second);

      if (sFile.is_open()) {

	while(sFile.good()) {

	  std::vector <std::string> v(8);
	  getline(sFile, inLine);
	  istringstream iss(inLine);

	  /*
	    1 i7 L1 Lower level of transition (integer)
	    2 i7 L2 Upper level of transition (integer)
	    3 e12.3 DE Energy of transition, Rydberg (float)
	    4 e12.3 GF Oscillator strength (float)
	    5 e12.3 LIM High-temperature limit value (float)
	    6 i5 NT Number of scaled temperatures
	    7 i3 T TYPE BT92 Transition type (integer)
	    8 e12.3 C VAL BT92 scaling parameter (float)
	  */

	  // Split according to fortran formats field widths
	  //
	  int c = 0;
	  v[0] = inLine.substr(c,  7); c +=  7;
	  if(atoi(v[0].c_str()) == -1) break;
	  v[1] = inLine.substr(c,  7); c +=  7;
	  v[2] = inLine.substr(c, 12); c += 12;
	  v[3] = inLine.substr(c, 12); c += 12;
	  v[4] = inLine.substr(c, 12); c += 12;
	  v[5] = inLine.substr(c,  5); c +=  5;
	  v[6] = inLine.substr(c,  5); c +=  5;
	  v[7] = inLine.substr(c, 12); c += 12;
	  
	  s.Z       = Z;
	  s.C       = C;
	  s.i       = atoi(v[0].c_str());
	  s.j       = atoi(v[1].c_str());
	  s.delERyd = atof(v[2].c_str());
	  s.gf      = atof(v[3].c_str());
	  s.Lim     = atof(v[4].c_str());
	  s.nt      = atoi(v[5].c_str());
	  s.type    = atoi(v[6].c_str());
	  s.Const   = atof(v[7].c_str());

				// Discard scaled temp grid
	  getline(sFile, inLine);

				// Get spline coefficients
	  getline(sFile, inLine);

	  iss.clear();		// Restock input string stream
	  iss.str(inLine);
				// Erase and repopulate string vector
	  v.clear();
	  copy(istream_iterator<std::string>(iss), 
	       istream_iterator<std::string>(), 
	       back_inserter<vector<std::string> >(v));

	  for (int i=0; i<s.nt; i++) {
	    s.spline.push_back(atof(v[i].c_str()));
	  }
	
	  splups.push_back(s);
	  s.spline.clear();
	}
	sFile.close();
      }
      else {
	if (nOK==1)
	  std::cout << "Cannot find file: " << ret.second << std::endl;
	else
	  std::cout << "Cannot find CHIANTI database"     << std::endl;
      }
    }
  }

  MPI_Bcast(&nOK, 1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  if (nOK)  {
    if (myid == 0) {
      std::cerr << "Ion::readSplups: problem reading CHIANTI files "
		<< "for Z=" << Z << " C=" << C << std::endl;
    }
    MPI_Finalize();
    exit(45);
  } 
  
  unsigned number = splups.size();
  MPI_Bcast(&number, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

  if (myid==0) {
    splupsType::iterator it = splups.begin();
    for (unsigned i=0; i<number; i++) {
      it->synchronize();
      it++;
    }
  } else {
    splups.clear();
    for (unsigned i=0; i<number; i++) {
      splups_data s;
      s.synchronize();
      splups.push_back(s);
    }
  }
}

//
// Read in the direct ionization cross section splines from the
// CHIANTI database files
//
void Ion::readDi() 
{
  unsigned char nOK = 0;

  if (myid==0) {

    cFile ret = chiantiFile("diparams");
    
    nOK = ret.first;

    if (nOK == 0) {

      std::string inLine;
      di_data s;
      ifstream sFile(ret.second);
      int i = 0;
      int i_fac = 0;
      if (sFile.is_open()) {
	while (sFile.good()) {
	  std::vector <std::string> v;
	  getline(sFile, inLine);
	  istringstream iss(inLine);
	  copy(istream_iterator<std::string>(iss), istream_iterator<std::string>(), 
	       back_inserter<vector<std::string> >(v));
	  
	  if(atoi(v[0].c_str()) == -1) break;
	  
	  if (i == 0) {
	    di_header.Z       = atoi(v[0].c_str());
	    di_header.C       = atoi(v[1].c_str());
	    di_header.nspline = atoi(v[2].c_str());
	    di_header.nfac    = atoi(v[3].c_str());
	    di_header.neav    = atoi(v[4].c_str());
	  }
	  else if (i%2 == 1) {
	    if (i_fac < di_header.nfac) {
	      s.btf = atof(v[0].c_str());
	      for(int i = 0; i < di_header.nspline; i++) {
		s.xspline.push_back(atof(v[i+1].c_str()));
	      }
	    }
	  }
	  else {
	    if (i_fac < di_header.nfac) {
	      s.ev = atof(v[0].c_str());
	      for(int i = 0; i < di_header.nspline; i++) {
		s.yspline.push_back(atof(v[i+1].c_str()));
	      }
	      diSpline.push_back(s);
	      s.xspline.clear();
	      s.yspline.clear();
	      i_fac++;
	    }
	  }
	  i++;
	}
	sFile.close();
      }
      else {
	if (nOK==1)
	  std::cout << "Cannot find file: " << ret.second << std::endl;
	else
	  std::cout << "Cannot find CHIANTI database"     << std::endl;
      }
    }
  }

  MPI_Bcast(&nOK, 1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  if (nOK)  {
    if (myid == 0) {
      std::cerr << "Ion::readDi: problem reading CHIANTI files "
		<< "for Z=" << Z << " C=" << C << std::endl;

    }
    MPI_Finalize();
    exit(46);
  } 

  unsigned number = diSpline.size();
  MPI_Bcast(&number, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

  di_header.synchronize();

  if (myid==0) {
    diSplineType::iterator it = diSpline.begin();
    for (unsigned i=0; i<number; i++) {
      it->synchronize();
      it++;
    }
  } else {
    diSpline.clear();
    for (unsigned i=0; i<number; i++) {
      di_data s;
      s.synchronize();
      diSpline.push_back(s);
    }
  }
}


//
// Read in the direct ionization cross section splines from the
// CHIANTI database files
//

//
// Initialization function when the master name is given
//
Ion::Ion(std::string name, chdata* ch) : ch(ch)
{
  MasterName = name;

  std::vector<std::string> v;
  boost::split(v, MasterName, boost::is_any_of("_") );
  eleName = v[0];

  convertName();		// Sets Z and C . . . 
  ip = ch->ipdata[lQ(Z, C)];

  freeFreeGridComputed  = false;
  radRecombGridComputed = false;
  exciteGridComputed    = false;
  ionizeGridComputed    = false;

  freeFreeGridMiss      = 0;
  freeFreeGridTry       = 0;
  freeFreeMissMin       = std::numeric_limits<double>::max();
  freeFreeMissMax       = 0;

  radRecombGridMiss     = 0;
  radRecombGridTry      = 0;
  radRecombMissMin      = std::numeric_limits<double>::max();
  radRecombMissMax      = 0;

  directIonGridMiss     = 0;
  directIonGridTry      = 0;
  directIonMissMin      = std::numeric_limits<double>::max();
  directIonMissMax      = 0;

  if (isInMasterList(MasterName)) {
    readfblvl();
    readelvlc();
    readwgfa();
    readSplups();
    readDi();
  }
  // else std::cout << "NOT IN LIST" <<std::endl;
  
  // E = h * nu
  //
  // E = h/(2*pi) * c * 2*pi/lambda = h-bar * c * k
  //
  // lambda = 2*pi/k.  Example: 1216 Ang = 121.6 nm --> k = 2*pi/121.6 = 0.052
  //                  1 mu  =  10000 Ang = 1000. nm --> k = 2*pi/1000. = 0.006
  //                 10 mu  = 100000 Ang = 1e04  nm --> k = 2*pi/1e04  = 0.0006
  //         
  //

  // initialize the k-grid (in inverse nm) for ff and the energy grid (in eV)
  //
  for (double k = kmin; k < kmax; k += kdel) 
    kgrid.push_back(k);


  kffsteps = kgrid.size();

  kgr10.resize(kffsteps);
  for (int i=0; i<kffsteps; i++) {
    kgr10[i] = pow(10, kgrid[i]) * hbc;
  }

  // Initialize the nu-grid for photoionization (in units of the
  // ionization potential) 

  IBinit = false;
  if (ib_type != Ion::none) {
    for (double nu = numin; nu < numax; nu += nudel) 
      nugrid.push_back(nu);

    nusteps = nugrid.size();
  }
}

//
// Constructor when the Z, C pair is given
//
Ion::Ion(unsigned short Z, unsigned short C, chdata* ch) : ch(ch), Z(Z), C(C)
{
  d = false;
  MasterName = ZCtoName(Z, C);

  std::vector<std::string> v;
  boost::split(v, MasterName, boost::is_any_of("_") );
  eleName = v[0];

  freeFreeGridComputed  = false;
  radRecombGridComputed = false;
  exciteGridComputed    = false;
  ionizeGridComputed    = false;

  ip = 0.0;

  if (Z>=C) {

    ip = ch->ipdata[lQ(Z, C)];

    if (isInMasterList(MasterName)) {
      readfblvl();
      readSplups();
      readDi();
      readelvlc();
      readwgfa();
    } else {
      std::cerr << "MasterName [" << MasterName << "] not in master list" 
		<< std::endl;
      std::cerr << "Attempting to read fblvl and diparams files";

      std::string MasterNameT = ZCtoName(Z, C);

      char * val = getenv("CHIANTI_DATA");

      {
	std::string fileName(val);
    
	fileName.append("/");
	fileName.append(eleName); 
	fileName.append("/");
	fileName.append(MasterNameT); 
	fileName.append("/"); 
	fileName.append(MasterNameT);
	fileName.append(".fblvl");
	
	std::string inLine;
	ifstream file(fileName.c_str());
	
	if (file.is_open()) {
	  file.close();
	  std::cerr << "...fblvl ok";
	  readfblvl();
	} else {
	  std::cerr << "...no fblvl";
	}
      }

      {
	std::string fileName(val);
    
	fileName.append("/");
	fileName.append(eleName); 
	fileName.append("/");
	fileName.append(MasterNameT); 
	fileName.append("/"); 
	fileName.append(MasterNameT);
	fileName.append(".diparams");
	
	std::string inLine;
	ifstream file(fileName.c_str());
	
	if (file.is_open()) {
	  file.close();
	  std::cerr << "...diparams ok";
	  readDi();
	} else {
	  std::cerr << "...no diparams";
	}
      }
      cerr << std::endl;
    }
  }
  
  // Initialize the k-grid (in inverse nm) for ff and the energy grid
  // (in eV)
  double k = 0;
  for (k = kmin; k < kmax; k += kdel) {
    kgrid.push_back(k);
  }
  kffsteps = kgrid.size();

  kgr10.resize(kffsteps);
  for (int i=0; i<kffsteps; i++) {
    kgr10[i] = pow(10, kgrid[i]) * hbc;
  }

  // Initialize the nu-grid for photoionization (in units of the
  // ionization potential) 

  IBinit = false;
  if (ib_type != Ion::none) {
    for (double nu = numin; nu < numax; nu += nudel) 
      nugrid.push_back(nu);

    nusteps = nugrid.size();
  }
}

//
// Default constructor: NOT CURRENTLY USED
//
Ion::Ion() 
{
  Z = 1;
  C = 1;

  freeFreeGridComputed  = false;
  radRecombGridComputed = false;
  exciteGridComputed    = false;
  ionizeGridComputed    = false;
}

//
// Copy constructor
//
Ion::Ion(const Ion &I) 
{
  Z          = I.Z;
  C          = I.C;
  ip         = I.ip;

  kffsteps   = I.kffsteps;

  kgrid      = I.kgrid;
  kgr10      = I.kgr10;
  fblvl      = I.fblvl;

  diSpline   = I.diSpline;
  di_header  = I.di_header;

  MasterName = I.MasterName;
  eleName    = I.eleName;
  elvlc      = I.elvlc;
  wgfa       = I.wgfa;

  freeFreeGridComputed  = I.freeFreeGridComputed;
  radRecombGridComputed = I.radRecombGridComputed;

  exciteGridComputed    = I.exciteGridComputed;
  ionizeGridComputed    = I.ionizeGridComputed;

  NfreeFreeGrid         = I.NfreeFreeGrid;
  NradRecombGrid        = I.NradRecombGrid;

  freeFreeGrid          = I.freeFreeGrid;
  radRecombGrid         = I.radRecombGrid;

  freeFreeGridMiss      = I.freeFreeGridMiss;
  freeFreeGridTry       = I.freeFreeGridTry;
  freeFreeMissMin       = I.freeFreeMissMin;
  freeFreeMissMax       = I.freeFreeMissMax;

  radRecombGridMiss     = I.radRecombGridMiss;
  radRecombGridTry      = I.radRecombGridTry;
  radRecombMissMin      = I.radRecombMissMin;
  radRecombMissMax      = I.radRecombMissMax;

  directIonGridMiss     = I.directIonGridMiss;
  directIonGridTry      = I.directIonGridTry;
  directIonMissMin      = I.directIonMissMin;
  directIonMissMax      = I.directIonMissMax;

  collideEmin           = I.collideEmin;
  ionizeEmin            = I.ionizeEmin;
  collideEmax           = I.collideEmax;
  ionizeEmax            = I.ionizeEmax;
  delCollideE           = I.delCollideE;

  collideDataGrid       = I.collideDataGrid;
  ionizeDataGrid        = I.ionizeDataGrid;

}

/** 
    Calculate the collision excitation cross sections and return the
    cumulative cross section array.

    Returns a vector with a pair of (cross section, energy difference)
    since the file input, and thus array, are not in any specific order
*/
Ion::collType
Ion::collExciteCrossSingle(double E, int id)
{
				// This will contain the cumulative
				// cross section
  collType CEcum;
				// Zero-valued datum
  const std::pair<double,double> Null(0, 0);

				// If the data is missing, assume zero
				// cross section
  if (splups.size() == 0) {
    CEcum.push_back(Null);
    return CEcum;
  }

  // The collision strengths have already by divided by the
  // statistical weight of the ground level 2j+1.  Saving the ground
  // state multiplicity, mult0, allows the ratio to be used in
  // computing Omega

  double totalCross = 0.0, mult0 = 2*elvlc[1].J + 1;

  for (size_t i=0; i<splups.size(); i++) {

    double EijEv = splups[i].delERyd*RydtoeV;
    double Const = splups[i].Const;

    if (splups[i].i==1 and E >= EijEv) {

      assert(splups[i].i == 1);
      assert(splups[i].spline.size() != 0);
      
      // Following Burgess & Tully (BT), 1992, Section 3
      //
      double Ej = E - EijEv, x = 0, y = 0;
      int  type = splups[i].type;
      
				// BT eq. 5, eq. 13
      if (type==1 or type==4) {
	x = 1.0 - (log(Const)/(log((Ej/EijEv) + Const)));
      }
				// BT eq. 9, eq. 11
      if (type==2 or type==3) {
	x = (Ej/EijEv)/((Ej/EijEv) + Const);
      }
				// BT eq. 11

      // xmin is 0 and xmax is 1, so this if statement is to make sure
      // x is within the bounds of interpolation
      //
      if ( x < 0 or x > 1.0) {
	std::cout << "ERROR IN EXCITATION CROSS: Ej = " << Ej
		  << " Eij = " << EijEv << " x = " << x
		  << " [Node ID=" << myid << "]" << std::endl;
      }

      // An extra couple of sanity checks for the interpolation
      //
      if (type > 4) break;
      
      int nspl = splups[i].spline.size();

      // Need at least 3 points for a spline
      //
      if (nspl<3) break;

      CacheSplList::iterator it = splUps.find(i);
      CsplD2Ptr sp;

      if (it == splUps.end()) {

	if (x_array.find(nspl) == x_array.end()) {
	  std::vector<double> t(nspl);
	  double dx = 1.0/(nspl - 1);
	  t[0] = 0.0;
	  t[nspl-1] = 1.0;
	  for (int n=1; n<nspl-1; n++) t[n] = dx*n;
	  x_array[nspl] = t;
	}

	splUps[i] = sp = CsplD2Ptr(new CsplineD2(x_array[nspl], splups[i].spline));
      } else {
	sp = it->second;
      }

      y = (*sp)(x);

      // Calculate the collision strength from the interpolated value
      double CStrength = 0.0;
				// BT, eq. 6
      if (type == 1) {
	CStrength = y * log((Ej/EijEv) + M_E);
      }
				// BT, eq. 10
      if (type == 2) {
	CStrength = y;
      }
				// BT, eq. 12
      if (type == 3) {
	double fac = Ej/EijEv + 1.0;
	CStrength = y/(fac*fac);
      }
				// BT, eq. 14
      if (type == 4) {
	CStrength = y * log((Ej/EijEv) + C);
      }
      
      // From Dere et al. 1997 
      //
      elvlcType::iterator eit = elvlc.find(splups[i].j-1);
      if (eit != elvlc.end()) {
	double weight = 1.0;
	if (mult0>0.0) weight = (2*eit->second.J+1)/mult0;
	if (weight>0) {
	  double crs1 = (M_PI*a0*a0*(CStrength/weight))/(E*Ion::eVtoRyd);
	  if (std::isinf(crs1)) {
	    std::cout << "crs1 is Inf: weight=" << weight << ", E="
		      << E << std::endl;
	  } else {
	    totalCross += crs1;
	    std::pair<double, double> cumi(totalCross, EijEv);
	    CEcum.push_back(cumi); // Add to the cumulative tally
	  }
	} else {
	  std::cout << "Coll crs for level=" << splups[i].j-1
		    << " at (Z, C)=(" << Z << ", " << C << ")"
		    << " has zero weight" << std::endl;
	}
      }
    }
    
  }
  
  if (CEcum.size() == 0) CEcum.push_back(Null);
  
  return CEcum;
}

Ion::collType
Ion::collExciteCrossGrid(double E, int id)
{
				// This will contain the cumulative
				// cross section
  collType CEcum;

				// Zero-valued datum
  const std::pair<double,double> Null(0, 0);

				// If the data is missing, assume zero
				// cross section
  if (splups.size() == 0) {
    CEcum.push_back(Null);
    return CEcum;
  }

  if (not exciteGridComputed) collExciteMakeGrid(id);

  if (E > collideEmax or E < collideEmin) {
    CEcum.push_back(Null);
    return CEcum;
  }

  // Interpolate the values
  //
  size_t indx = std::floor( (E - collideEmin)/delCollideE );

  // Sanity check
  //
  indx = std::max<size_t>(0, std::min<size_t>(indx, NcollideGrid-2));

  double eA   = collideEmin + delCollideE*indx;
  double eB   = collideEmin + delCollideE*(indx+1);

  double A = (eB - E)/delCollideE;
  double B = (E - eA)/delCollideE;

  typedef std::pair<double, double> Elem;
  std::array<Elem, 2> v { collideDataGrid[indx+0].back(), collideDataGrid[indx+1].back()};
  
  std::vector<Elem> ret(1);
  ret[0] = {A*v[0].first  + B*v[1].first, A*v[0].second + B*v[1].second};

  // Test
  if (false) {
    collType tst = collExciteCrossSingle(E, id);
    double diffXS = fabs(tst.back().first  - ret[0].first );
    double diffCM = fabs(tst.back().second - ret[0].second);
    if (tst.back().first >0.0) diffXS /= tst.back().first;
    if (tst.back().second>0.0) diffCM /= tst.back().second;
    if (diffXS > 0.01 or diffCM > 0.01) {
      std::cout << "Diff [collide]: XS expected=" << tst.back().first  << ", got "
		<< ret[0].first  << std::endl
		<< "                CM expected=" << tst.back().second << ", got "
		<< ret[0].second << std::endl;
    }
  }

  return ret;
}


void Ion::collExciteMakeGrid(int id)
{
  exciteGridComputed = true;

  if (splups.size()) {

    // Get min/max energy
    //
    collideEmax = 0.0;
    collideEmin = std::numeric_limits<double>::max();
    for (size_t i=0; i<splups.size(); i++) {
      double Elev = splups[i].delERyd*RydtoeV;
      if (splups[i].i==1) {
	collideEmin = std::min<double>(collideEmin, Elev);
	collideEmax = std::max<double>(collideEmax, Elev);
      }
    }

    // Number of elements in energy grid
    //
    NcollideGrid = 1 + std::floor((collideEmax - collideEmin)/DeltaEGrid );
    // NcollideGrid = std::max<int>(NcollideGrid, 10000);
    NcollideGrid = std::max<int>(NcollideGrid, 100);
    NcollideGrid = std::min<int>(NcollideGrid, 2000);
    delCollideE  = (collideEmax - collideEmin)/(NcollideGrid-1);
    
    // Compute the grid
    //
    collideDataGrid.resize(NcollideGrid);
    for (int n=0; n<NcollideGrid; n++) 
      collideDataGrid[n] = collExciteCrossSingle(collideEmin + delCollideE*n, id);
  } else {
    collideEmax  = collideEmin = delCollideE = 0.0;
    NcollideGrid = 0;
  }
}

//
// Calculate the Qr-prime value as in Fontes, Sampson, Zhang 1999
//
double Ion::qrp(double u) 
{
  double A = 1.13;
  double D, C, d, c;
  if (Z >= 16) {
    c = -0.28394;
    d = 1.95270;
    C = 0.20594;
    D = 3.70590;
  }
  else {
    c = -0.80414;
    d = 2.32431;
    C = 0.14424;
    D = 3.82652;
  }

  if (Z > 20) {
    C += pow(((Z-20.0)/50.5), 1.11);
  }

  double z  = 1.0 - 1.0/u;
  double z2 = z*z;
  double z4 = z2*z2;
  double q  = (A*log(u) + D*z2 + C*u*z4 + (c/u + d/(u*u))*z) / u;

  return q;
}

/** 
    Calculate the direct ionization cross section from the spline,
    which is a function of the interaction energy of the electron

    See: Dere, K. P., 2007, A&A, 466, 771
    ADS ref:  http://adsabs.harvard.edu/abs/2007A%26A...466..771D
*/
double Ion::directIonCrossSingle(double E, int id) 
{
  // Classical Bohr cross section in nm^2
  //
  constexpr double bohr_cs = M_PI*a0*a0;

  // Test for hydrogen-like/helium-like ion
  //
  unsigned char I = Z - C + 1;

  if (C == (Z+1)) {
    return -1;
  }

  // Scaled energy
  //
  double u        = E/ip;

  if (u<1.0) return 0.0;

  // Ionization potential in Rydbergs
  //
  double ipRyd    = ip/RydtoeV;

  double F, qr, cross;
  
  // From Fontes, et al. Phys Rev A, 59, 1329 (eq. 2.11)
  //
  if (Z >= 20) {
    F = (140.0+pow((double(Z)/20.0),3.2))/141.;
  }
  else {
    F = 1.0;
  }

  qr = qrp(u)*F;

  // first two if statements are whether or not to use Fontes cross
  // sections
  //
  if (I == 1 && Z >= 6) {	// Hydrogenic
    cross = bohr_cs*qr/(ipRyd*ipRyd);
  }
  else if (I==2 && Z >= 10) {	// Helium like
    cross = 2.0*bohr_cs*qr/(ipRyd*ipRyd);
  }
  else {
    cross = 0;
    for (int i = 0; i < di_header.nfac; i++) {
      if (E >= diSpline[i].ev) {

	double u1  = E/diSpline[i].ev;
	double bte = 1.0 - log(diSpline[i].btf)/log(u1-1.0+diSpline[i].btf);

	CacheSplList::iterator it = diSpln.find(i);
	CsplD2Ptr sp;

	if (it == diSpln.end()) {
	  sp = CsplD2Ptr(new CsplineD2(diSpline[i].xspline, diSpline[i].yspline));
	  diSpln[i] = sp;
	} else {
	  sp = it->second;
	}

	double btcross = (*sp)(bte);
	double cross_i = (log(u1) + 1.0)*btcross/(u1*diSpline[i].ev*diSpline[i].ev);
				// convert to cross section in nm^2
	cross += cross_i;
      }
    }
  }

  return cross;
}

double Ion::directIonCrossGrid(double E, int id) 
{
  if (C == (Z+1)) {
    return -1;
  }

  if (E<ip) return 0.0;

  if (not ionizeGridComputed) directIonMakeGrid(id);

  double eMin = ionizeEmin, eMax = ionizeEmin + DeltaEGrid*(NionizeGrid-1);
  
  bool gridProc = false;
  if (GridDebug and myid==0 and id==0) {
    directIonGridTry++;
    gridProc = true;
  }
  
  if (E<eMin or E>eMax) {
    if (gridProc) {
      directIonGridMiss++;
      directIonMissMin = std::min<double>(directIonMissMin, E);
      directIonMissMax = std::max<double>(directIonMissMax, E);
    }
    return directIonCrossSingle(E, id);
  }

  if (gridProc and (directIonGridMiss % GridReport == 1)) {
    std::cout << "DirectIonGrid DEBUG: "
      << static_cast<double>(directIonGridMiss)/directIonGridTry
	      << " [" << directIonGridTry << "], (min, max)=("
	      << directIonMissMin << ", " << directIonMissMax << ")"
	      << std::endl;

    directIonMissMin = std::numeric_limits<double>::max();
    directIonMissMax = 0;
  }

  size_t indx = std::floor( (E - eMin)/DeltaEGrid );
  double eA   = eMin + DeltaEGrid*indx, eB = eMin + DeltaEGrid*(indx+1);

  double A = (eB - E)/DeltaEGrid;
  double B = (E - eA)/DeltaEGrid;

  double est = A*ionizeDataGrid[indx+0] + B*ionizeDataGrid[indx+1];

  // Test
  if (false) {
    double tst = directIonCrossSingle(E, id);
    double diffXS = fabs(tst - est);
    if (tst >0.0) diffXS /= tst;
    if (diffXS > 0.01) {
      std::cout << "Diff [ionize]: XS expected=" << tst << ", got "
		<< est << std::endl;
    }
  }

  return est;
}

void Ion::directIonMakeGrid(int id)
{
  ionizeGridComputed = true;

  // Get max energy
  ionizeEmin = ip;
  ionizeEmax = EmaxGrid;

  // Number of elements in energy grid
  NionizeGrid = 1 + std::floor( (ionizeEmax - ionizeEmin)/DeltaEGrid );
  NionizeGrid = std::max<int>(NionizeGrid, 256);

  // Make the grid
  ionizeDataGrid.resize(NionizeGrid);
  for (int n=0; n<NionizeGrid; n++) 
    ionizeDataGrid[n] = directIonCrossSingle(ionizeEmin + DeltaEGrid*n, id);
}


/** 
    Cross section is 3BN(a) from Koch & Motz 1959, with the low-energy
    Elwert factor (see Koch & Motz eq. II-6), nonrelativistic limit

    Using the parametrization by Greene (1959)
 */
std::pair<double, double> Ion::freeFreeCrossSingle(double Ei, int id) 
{
  // No free-free with a neutral
  //
  if (C==1) return std::pair<double, double>(0.0, 0.0);

  // Scaled inverse energy (initial)
  //
  double ni       = sqrt( RydtoeV*(C-1)*(C-1)/Ei );

  // Integration variables
  //
  double cum      = 0;
  double dk       = (kgrid[1] - kgrid[0])*log(10.0);

  std::vector<double> diff, cuml;

  for (int j = 0; j < kffsteps; j++) {
    //
    // Photon energy in eV
    //
    double k      = kgr10[j];

    //
    // Final kinetic energy
    //
    double Ef     = Ei - k;

    //
    // Can't emit a photon if not enough KE!
    //
    if (Ef <= 0.0) break;

    //
    // Scaled inverse energy (final)
    //
    double nf     = sqrt( RydtoeV*(C-1)*(C-1)/Ef );

    //
    // Elwert factor
    //
    double corr   = (1.0 - exp(-2.0*M_PI*ni))/(1.0 - exp(-2.0*M_PI*nf));

    //
    // Differential cross section contribution
    //
    double dsig   = A * ni*nf * log((nf + ni)/(nf - ni)) * corr * dk;

    cum = cum + dsig;

    diff.push_back(dsig/dk);
    cuml.push_back(cum);
  }


  double phi = 0.0, ffWaveCross = 0.0;

  // If cross section is offgrid, set values to zero
  //
  if (cum > 0.0) {

    // Location in cumulative cross section grid
    //
    double rn = cum * static_cast<double>(rand())/RAND_MAX;
  
    // Interpolate the cross section array
    //
    
    // Points to first element that is not < rn
    // but may be equal
    std::vector<double>::iterator lb = 
      std::lower_bound(cuml.begin(), cuml.end(), rn);
    
    // Assign upper end of range to the
    // found element
    //
    std::vector<double>::iterator ub = lb;
    //
    // If is the first element, increment
    // the upper boundary
    //
    if (lb == cuml.begin()) { if (cuml.size()>1) ub++; }
    //
    // Otherwise, decrement the lower boundary
    //
    else { lb--; }
    
    // Compute the associated indices
    //
    size_t ii = lb - cuml.begin();
    size_t jj = ub - cuml.begin();
    double  k = kgrid[ii];
	  
    // Linear interpolation
    //
    if (*ub > *lb) {
      double d = *ub - *lb;
      double a = (rn - *lb) / d;
      double b = (*ub - rn) / d;
      k  = a * kgrid[ii] + b * kgrid[jj];
    }
    
    // Assign the photon energy
    //
    ffWaveCross = pow(10, k) * hbc;

    // Use the integrated cross section from the differential grid
    //
    phi = cuml.back();
  }

  return std::pair<double, double>(phi, ffWaveCross);
}


void Ion::freeFreeMakeEvGrid(int id)
{
  if (C==1) return;

  freeFreeGridComputed = true;

  // Integration variables
  //
  double dk       = (kgrid[1] - kgrid[0])*log(10.0);

  // Energy grid
  NfreeFreeGrid = 1 + std::floor( (EmaxGrid - EminGrid)/DeltaEGrid );

  // Allocate storage for grid
  freeFreeGrid.resize(NfreeFreeGrid);

  for (size_t n = 0; n < NfreeFreeGrid; n++) {

    // Energy in eV
    //
    double Ei = EminGrid + DeltaEGrid*n;

    // Scaled inverse energy
    //
    double ni = sqrt( RydtoeV*(C-1)*(C-1)/Ei );

    freeFreeGrid[n].clear();

    double cum = 0;

    for (int j = 0; j < kffsteps; j++) {
      //
      // Photon energy in eV
      //
      double k  = kgr10[j];

      //
      // Final kinetic energy
      //
      double Ef = Ei - k;
      
      //
      // Can't emit a photon if not enough KE!
      //
      if (Ef <= 0.0) break;

      //
      // Scaled inverse energy (final)
      //
      double nf     = sqrt( RydtoeV*(C-1)*(C-1)/Ef );
      
      //
      // Elwert factor
      //
      double corr   = (1.0 - exp(-2.0*M_PI*ni))/(1.0 - exp(-2.0*M_PI*nf));
      
      //
      // Differential cross section contribution
      //
      double dsig   = A * ni*nf * log((nf + ni)/(nf - ni)) * corr * dk;
      
      cum = cum + dsig;
      
      freeFreeGrid[n].push_back(cum);
    }

  }

}


std::pair<double, double> Ion::freeFreeCrossEvGrid(double E, int id)
{
  // No free-free with a neutral
  //
  if (C==1) return std::pair<double, double>(0.0, 0.0);

  if (not freeFreeGridComputed) freeFreeMakeEvGrid(id);

  double eMin = EminGrid, eMax = EminGrid + DeltaEGrid*(NfreeFreeGrid-1);

  bool gridProc = false;
  if (GridDebug and myid==0 and id==0) {
    freeFreeGridTry++;
    gridProc = true;
  }

  if (E<eMin) E = eMin;		// Enforce minimum energy to prevent
				// off grid evaluatoin
  if (E<eMin or E>eMax) {
    if (gridProc) {
      freeFreeGridMiss++;
      freeFreeMissMin = std::min<double>(freeFreeMissMin, E);
      freeFreeMissMax = std::max<double>(freeFreeMissMax, E);
    }
    return freeFreeCrossSingle(E, id);
  }

  if (gridProc and (freeFreeGridMiss % GridReport == 1)) {
    std::cout << "FreeFreeGrid DEBUG: "
      << static_cast<double>(freeFreeGridMiss)/freeFreeGridTry
	      << " [" << freeFreeGridTry << "], (min, max)=("
	      << freeFreeMissMin << ", " << freeFreeMissMax << ")"
	      << std::endl;

    freeFreeMissMin = std::numeric_limits<double>::max();
    freeFreeMissMax = 0;
  }

  double phi = 0.0, ffWaveCross = 0.0;

  size_t indx = std::floor( (E - eMin)/DeltaEGrid );
  double eA   = eMin + DeltaEGrid*indx, eB = eMin + DeltaEGrid*(indx+1);

  double A = (eB - E)/DeltaEGrid;
  double B = (E - eA)/DeltaEGrid;

  std::array<double, 2> cum
  { freeFreeGrid[indx+0].back(), freeFreeGrid[indx+1].back()};

  // If cross section is offgrid, set values to zero
  //
  if (cum[0] > 0.0 and cum[1] > 0.0) {

    std::array<double,  2> k;

    // Location in cumulative cross section grid
    //
    double rn = static_cast<double>(rand())/RAND_MAX;

    for (int i=0; i<2; i++) {

      // Interpolate the cross section array
      //
    
      if (rn*cum[i] < freeFreeGrid[indx+i].front()) {
	
	k[i] = kgrid[0];
	

      } else {

	// Points to first element that is not < rn
	// but may be equal
	std::vector<double>::iterator lb = 
	  std::lower_bound(freeFreeGrid[indx+i].begin(), freeFreeGrid[indx+i].end(), rn*cum[i]);
    
	// Assign upper end of range to the
	// found element
	//
	std::vector<double>::iterator ub = lb;

	//
	// If is the first element, increment
	// the upper boundary
	//
	if (lb == freeFreeGrid[indx+i].begin()) {
	  if (freeFreeGrid[indx+i].size()>1) ub++;
	}
	//
	// Otherwise, decrement the lower boundary
	//
	else { lb--; }
      
	// Compute the associated indices
	//
	size_t ii = std::distance(freeFreeGrid[indx+i].begin(), lb);
	size_t jj = std::distance(freeFreeGrid[indx+i].begin(), ub);
	
	k[i] = kgrid[ii];
	
	// Linear interpolation
	//
	if (*ub > *lb) {
	  double d = *ub - *lb;
	  double a = (rn*cum[i] - *lb) / d;
	  double b = (*ub - rn*cum[i]) / d;
	  k[i] = a * kgrid[ii] + b * kgrid[jj];
	  if (fabs(a)>1.0 or fabs(b)>1.0) {
	    std::cout << "Oops" << std::endl;
	  }
	}
      }
    }

    double K = A*k[0] + B*k[1];


    // Assign the photon energy
    //
    ffWaveCross = pow(10, K) * hbc;

    // Use the integrated cross section from the differential grid
    //
    phi = A*cum[0] + B*cum[1];
  }

  return std::pair<double, double>(phi, ffWaveCross);
}

std::pair<double, double> Ion::freeFreeCrossEvGridTest(double E, double rn, int id)
{
  // No free-free with a neutral
  //
  if (C==1) return std::pair<double, double>(0.0, 0.0);

  if (not freeFreeGridComputed) freeFreeMakeEvGrid(id);

  double eMin = EminGrid, eMax = EminGrid + DeltaEGrid*(NfreeFreeGrid-1);

  bool gridProc = false;
  if (GridDebug and myid==0 and id==0) {
    freeFreeGridTry++;
    gridProc = true;
  }

  if (E<eMin) E = eMin;		// Enforce minimum energy to prevent
				// off grid evaluatoin
  if (E<eMin or E>eMax) {
    if (gridProc) {
      freeFreeGridMiss++;
      freeFreeMissMin = std::min<double>(freeFreeMissMin, E);
      freeFreeMissMax = std::max<double>(freeFreeMissMax, E);
    }
    return freeFreeCrossSingle(E, id);
  }

  if (gridProc and (freeFreeGridMiss % GridReport == 1)) {
    std::cout << "FreeFreeGrid DEBUG: "
      << static_cast<double>(freeFreeGridMiss)/freeFreeGridTry
	      << " [" << freeFreeGridTry << "], (min, max)=("
	      << freeFreeMissMin << ", " << freeFreeMissMax << ")"
	      << std::endl;

    freeFreeMissMin = std::numeric_limits<double>::max();
    freeFreeMissMax = 0;
  }

  double phi = 0.0, ffWaveCross = 0.0;

  size_t indx = std::floor( (E - eMin)/DeltaEGrid );
  double eA   = eMin + DeltaEGrid*indx, eB = eMin + DeltaEGrid*(indx+1);

  double A = (eB - E)/DeltaEGrid;
  double B = (E - eA)/DeltaEGrid;

  std::array<double, 2> cum
  { freeFreeGrid[indx+0].back(), freeFreeGrid[indx+1].back()};

  // If cross section is offgrid, set values to zero
  //
  if (cum[0] > 0.0 and cum[1] > 0.0) {

    std::array<double,  2> k;

    // Location in cumulative cross section grid
    //

    for (int i=0; i<2; i++) {

      // Interpolate the cross section array
      //
    
      if (rn*cum[i] < freeFreeGrid[indx+i].front()) {
	
	k[i] = kgrid[0];
	

      } else {

	// Points to first element that is not < rn
	// but may be equal
	std::vector<double>::iterator lb = 
	  std::lower_bound(freeFreeGrid[indx+i].begin(), freeFreeGrid[indx+i].end(), rn*cum[i]);
    
	// Assign upper end of range to the
	// found element
	//
	std::vector<double>::iterator ub = lb;

	//
	// If is the first element, increment
	// the upper boundary
	//
	if (lb == freeFreeGrid[indx+i].begin()) {
	  if (freeFreeGrid[indx+i].size()>1) ub++;
	}
	//
	// Otherwise, decrement the lower boundary
	//
	else { lb--; }
      
	// Compute the associated indices
	//
	size_t ii = std::distance(freeFreeGrid[indx+i].begin(), lb);
	size_t jj = std::distance(freeFreeGrid[indx+i].begin(), ub);
	
	k[i] = kgrid[ii];
	
	// Linear interpolation
	//
	if (*ub > *lb) {
	  double d = *ub - *lb;
	  double a = (rn*cum[i] - *lb) / d;
	  double b = (*ub - rn*cum[i]) / d;
	  k[i] = a * kgrid[ii] + b * kgrid[jj];
	  if (fabs(a)>1.0 or fabs(b)>1.0) {
	    std::cout << "Oops" << std::endl;
	  }
	}
      }
    }

    double K = A*k[0] + B*k[1];


    // Assign the photon energy
    //
    ffWaveCross = pow(10, K) * hbc;

    // Use the integrated cross section from the differential grid
    //
    phi = A*cum[0] + B*cum[1];
  }

  return std::pair<double, double>(phi, ffWaveCross);
}

void Ion::radRecombMakeEvGrid(int id)
{
  radRecombGridComputed = true;

  // Energy grid
  NradRecombGrid = 1 + std::floor( (EmaxGrid - EminGrid)/DeltaEGrid );

  // Allocate storage for grid
  radRecombGrid.resize(NradRecombGrid);

  // Compute grid
  for (size_t n = 0; n < NradRecombGrid; n++) {
    double Ei = EminGrid + DeltaEGrid*n;
    radRecombGrid[n] = radRecombCrossSingle(Ei, id).back();
  }
}

std::vector<double> Ion::radRecombCrossEvGrid(double E, int id)
{
  // No recombination with a neutral (this test probably not needed here)
  //
  if (C==1) return std::vector<double>(1, 0.0);

  if (not radRecombGridComputed) radRecombMakeEvGrid(id);

  double eMin = EminGrid, eMax = EminGrid + DeltaEGrid*(NradRecombGrid-1);

  bool gridProc = false;
  if (GridDebug and myid==0 and id==0) {
    radRecombGridTry++;
    gridProc = true;
  }
  
  if (E<eMin) E = eMin;		// Enforce minimum energy to prevent
				// off grid evaluation
  if (E<eMin or E>eMax) {
    if (gridProc) {
      radRecombGridMiss++;
      radRecombMissMin = std::min<double>(radRecombMissMin, E);
      radRecombMissMax = std::max<double>(radRecombMissMax, E);
    }
    return radRecombCrossSingle(E, id);
  }
  
  if (gridProc and (radRecombGridMiss % GridReport == 1)) {
    std::cout << "RadRecombGrid DEBUG: "
      << static_cast<double>(radRecombGridMiss)/radRecombGridTry
	      << " [" << radRecombGridTry << "], (min, max)=("
	      << radRecombMissMin << ", " << radRecombMissMax << ")"
	      << std::endl;

    radRecombMissMin = std::numeric_limits<double>::max();
    radRecombMissMax = 0;
  }

  size_t indx = std::floor( (E - eMin)/DeltaEGrid );
  double eA   = eMin + DeltaEGrid*indx, eB = eMin + DeltaEGrid*(indx+1);

  double A = (eB - E)/DeltaEGrid;
  double B = (E - eA)/DeltaEGrid;

  return std::vector<double>(1, A*radRecombGrid[indx] + B*radRecombGrid[indx+1]);
}

std::vector<double> Ion::radRecombCrossSingle(double E, int id)
{
  // For testing . . .
  //
  //  +--- True for verbose debug reporting for all cross-section types
  //  |
  //  v
  if (false) {
    std::vector<double> v1 = radRecombCrossMewe   (E, id);
    std::vector<double> v2 = radRecombCrossTopBase(E, id);
    std::vector<double> v3 = radRecombCrossKramers(E, id);
    std::vector<double> v4 = radRecombCrossSpitzer(E, id);
    std::vector<double> v5 = radRecombCrossVerner (E, id);

    std::cout << "  [Z, C] = [" << Z << ", " << C << "]"     << std::endl;
    std::cout << "  E (eV) = " << std::setw(16) << E         << std::endl;
    std::cout << "    Mewe = " << std::setw(16) << v1.back() << std::endl;
    std::cout << " TopBase = " << std::setw(16) << v2.back() << std::endl;
    std::cout << " Kramers = " << std::setw(16) << v3.back() << std::endl;
    std::cout << " Spitzer = " << std::setw(16) << v4.back() << std::endl;
    std::cout << "  Verner = " << std::setw(16) << v5.back() << std::endl;
    std::cout << std::string(60, '-')                        << std::endl;
    
    if      (rr_type == mewe)    return v1;
    else if (rr_type == topbase) return v2;
    else if (rr_type == kramers) return v3;
    else if (rr_type == spitzer) return v4;
    else if (rr_type == verner)  return v5;
    else                         return v5;

  } else {

    if      (rr_type == mewe)    return radRecombCrossMewe   (E, id);
    else if (rr_type == topbase) return radRecombCrossTopBase(E, id);
    else if (rr_type == kramers) return radRecombCrossKramers(E, id);
    else if (rr_type == spitzer) return radRecombCrossSpitzer(E, id);
    else if (rr_type == verner)  return radRecombCrossVerner (E, id);
    else                         return radRecombCrossVerner (E, id);

  }
}


/**
   Compute total recombination cross section using Kramers b-f cross
   section and the Milne relation to get the f-b cross section

   The recombination cross section is related to the absorption cross
   section using the Milne relation:

   \sigma_R = \frac{g_A}{2g^+_A} \frac{(h\nu)^2}{Em_ec^2} \sigma_P

   where g_A is the degeneracy of the target state and g^+_A is the
   degeneracy of the ion (which we assume to be in the ground state)

   Semi-classical cross section (no Gaunt factor)
*/
std::vector<double> Ion::radRecombCrossKramers(double E, int id) 
{
  const double nfac   = 8.0*M_PI/pow(3.0, 1.5);

  const double incmEv = light * planck / eV;

  // Return vector
  //
  std::vector<double> radRecCum;
  double cross = 0.0;

  // This is the target neutral
  //
  Ion* N = ch->IonList[lQ(Z, C-1)].get();

  // Ionization threshold (eV)
  //
  double Eph = N->ip;

  // Compute the effective charge
  //
  double zz   = Z;
  double ii   = C - 1;
  double Zeff = 0.0;
  if (    zz >= ii && ii >= 0.5*zz) Zeff = 0.5*(zz + ii);
  if (0.5*zz >= ii && ii >= 1.0   ) Zeff = sqrt(zz * ii);
  
  double aeff  = a0/Zeff;
  
  for (auto j : N->fblvl) {

    fblvl_data* f = &j.second;

    // Line energy (eV)
    //
    double Elv = Eph;

    if (f->encm == 0 and f->encmth!=0) Elv -= f->encmth * incmEv;
    else if (f->encm != 0)             Elv -= f->encmth * incmEv;

    // Photon energy (eV)
    //
    double Enu = E + Elv;

    // Kramers cross section
    //
    double Erat   = Elv/Enu;
    double sigmaP = nfac*alpha0*alpha0*alpha0*2.0/f->lvl*aeff*aeff*Erat*Eph/E;

    cross += sigmaP;

    if (cross == 0) {
      std::cout << "NULL IN RAD RECOMB: Chi=" << ip
		<< ", E="      << E 
		<< ", n="      << f->lvl
		<< ", Elv="    << Elv
		<< ", Enu="    << Enu
		<< ", Erat="   << Erat
		<< ", sigmaP=" << sigmaP
		<< std::endl;
    }
      
    if (std::isnan(cross)) {
      std::cout << "NAN IN RAD RECOMB: Chi=" << ip
		<< ", E="      << E 
		<< ", n="      << f->lvl
		<< ", Elv="    << Elv
		<< ", Enu="    << Enu
		<< ", Erat="   << Erat
		<< ", sigmaP=" << sigmaP
		<< std::endl;
    }
  }

  radRecCum.push_back(cross);

  return radRecCum;
}


/** Calculates the differential radiative recombination cross section
    as a function of incoming electron impact energy, and returns the
    vector cumulative cross section array. 

    Details of implementation from CHIANTI. See "The Free-Bound
    Continuum", P.R. Young, Ver. 1.1, 8-Sep-2009.

    Uses Milne relation.
*/
std::vector<double> Ion::radRecombCrossMewe(double E, int id) 
{
  double incmEv = 1.239842e-4; // 1 inverse cm = 1.239.. eV

  // constant infront of the photo-cross using the Mewe method
  //
  double D = 1.075812e-23;

  // Electron rest mass in keV
  //
  double mec2 = 510.998896;

  // Key of parent ion
  //
  lQ Q(Z, C-1);

  // Get pointers to Ion data
  //
  double IP = ch->ipdata[Q];
  Ion* N    = ch->IonList[Q].get();

  // Convert kinetic energy to keV
  //
  E *= 1.0e-3;

  // Return values
  //
  std::vector<double> radRecCum;
  double cross = 0.0;
  
  if (E>0) {
    
    double mult0 = (C<=Z ? fblvl[1].mult : 1);

    for (auto j : N->fblvl) {

      fblvl_data* f = &j.second;
      double I = IP;

      if (f->encm == 0 and f->encmth!=0) I -= f->encmth*incmEv;
      else if (f->encm != 0)             I -= f->encm  *incmEv;

      if (I<=0.0) {
	std::cout << "ERROR in energy level for radRecombCrossMewe!" 
		  << "  ip=" << IP
		  << ", En=" << f->encmth
		  << ", Em=" << f->encm
		  <<std::endl;
      }

      I *= 1.0e-3;		// convert the energy to keV
      
      double mult = f->mult;	// Level multiplicity
      double n    = f->lvl ;	// Principle quantum number

      if (I >= 0 and n>0) {
				// Total radiated photon energy (line
				// + KE) in keV
	double hnu    = E + I;
				// Assume that g_{bf} = 1
				// 
	double sigmaP = D * I*I*pow(hnu, -3.0)/n; // in m^2
	
				// Apply Milne relation
				//
	double Erat   = (hnu*hnu)/(2.0*mec2*E);
	double crossi = mult/mult0 * Erat * sigmaP;

	cross += crossi;

	if (cross == 0) {
	  std::cout << "NULL in radRecombCrossMewe:" 
		    << "  Chi="   << ip
		    << ", I="     << I
		    << ", h*nu="  << hnu
		    << ", Erat="  << Erat
		    << ", mult="  << mult
		    << std::endl;
	}
	if (std::isnan(cross)) {
	  std::cout << "NaN in radRecombCrossMewe:" 
		    << "  Chi="   << ip
		    << ", I="     << I
		    << ", h*nu="  << hnu
		    << ", Erat="  << Erat
		    << ", mult="  << mult
		    << std::endl;
	}
      }
    }
  }

  // Convert to nm^2 -------+
  //                        |
  //                        v
  radRecCum.push_back(cross*1.0e18);

  return radRecCum;
}

std::vector<double> Ion::radRecombCrossSpitzer(double E, int id) 
{
				// 1 inverse cm = 1.239.. eV
  const double incmEv = 1.239842e-4;

				// Cross-section prefactor in nm^2
  const double coef   = 2.105310889751809e-08;

				// This is the target neutral
  Ion* N = ch->IonList[lQ(Z, C-1)].get();

				// Ionization energy in eV
  double ionE         = ch->ipdata[lQ(Z, 1)];

  std::vector<double> radRecCum;
  double cross = 0.0;
  if (E > 0) {
    for (auto j : N->fblvl) {
      fblvl_data* f = &j.second;

      double Ej = 0.0;
      if (f->lvl==1) 
	Ej = ionE;
      else if (f->lvl>1 && f->encmth > 0) 
	Ej = ionE - f->encmth * incmEv;
      else if (f->lvl>1 && f->encm > 0) 
	Ej = ionE - f->encm * incmEv;
      else continue;
      //
      double mult = static_cast<double>(f->mult);
      double n    = static_cast<double>(f->lvl );
      double Ephot  = E + Ej;
      double Erat   = Ej / Ephot;
      double crossn = coef * (Ej / Ephot) * (0.5*Ephot/E) * (1.0/(mult*n));

      cross += crossn;

      if (cross == 0) {
	std::cout << "NULL IN RAD RECOMB: " << ip << "\t" << Ej << "\t" 
		  << Ephot << "\t" << Erat << "\t" << mult << "\t" 
		  << n <<std::endl;
      }
      if (std::isnan(cross)) {
	std::cout << cross << "\t" << Ej << "\t" << Ephot << "\t" 
		  << n << "\t" << Erat << std::endl;
      }
    }
  }
  radRecCum.push_back(cross);

  return radRecCum;
}

//
// Use the TOPbase photoionization cross sections to compute the
// recombination cross sections using the Milne relation
//
std::vector<double> Ion::radRecombCrossTopBase(double E, int id) 
{
  // Initialize TopBase data (once) if needed
  //
  if (ch->tb.get() == 0) {
    if (myid==0) {
      std::cerr << "Ion: creating new TopBase instance"
		<< std::endl;
    }
    ch->tb = boost::shared_ptr<TopBase>(new TopBase);
    // For debugging (set to 'false' for production)
    //  |
    //  v
    if (false && myid==0) ch->tb->printInfo();
  }

  // Call for the cross section
  //
  TopBase::iKey k(Z, C);
  std::vector<double> ret(1, ch->tb->sigmaFB(k, E));

  return ret;
}

// 
// Summed recombination cross sections using the Verner relation
//
std::vector<double> Ion::radRecombCrossVerner(double E, int id) 
{
  // Call for the cross section
  //
  lQ Q(Z, C);
  std::vector<double> ret(1, ch->VernerXC.cross(Q, E));

  return ret;
}

//
// Printe various internal databases for debugging
//
void Ion::printInfo() 
{
  std::cout << "Master list name: " << MasterName <<std::endl;
  std::cout << "\t" << "Element: " << eleName <<std::endl << "\tZ = " << Z << "\n" << "\tC = " << C <<std::endl;
  std::cout << "\td = " << d <<std::endl;
  // std::cout << "\tAdundance = " << abundance <<std::endl;
  std::cout << "\tip = " << ip <<std::endl;
}

void Ion::printelvlc() {
  std::cout << "elvlc file for element " << MasterName <<std::endl;
  for(size_t i = 0; i < elvlc.size(); i++) {
    std::cout << elvlc[i].level       << "\t" 
	      << elvlc[i].designation << "\t"
	      << elvlc[i].spin        << "\t" 
	      << elvlc[i].l           << "\t" 
	      << elvlc[i].l_str       << "\t"
	      << elvlc[i].J           << "\t" 
	      << elvlc[i].encm        << "\t"
	      << elvlc[i].enry        << "\t" 
	      << elvlc[i].encmth      << std::endl;
  }
}

void Ion::printfblvl() 
{
  std::cout << "fblvl file for element " << MasterName <<std::endl;
  for (size_t i = 0; i < fblvl.size(); i++) {
    std::cout << fblvl[i].lvl         << "\t" 
	      << fblvl[i].conf        << "\t" 
	      << fblvl[i].pqn         << "\t"
	      << fblvl[i].l           << "\t" 
	      << fblvl[i].l_str       << "\t" 
	      << fblvl[i].mult        << "\t"
	      << fblvl[i].encm        << "\t" 
	      << fblvl[i].encmth      << std::endl;
  }	
}


//------------------------------------------------------------
// chdata functions
//------------------------------------------------------------

//
// Read in the master list
// 
void chdata::readMaster() 
{
  char * val;
  if ( (val = getenv("CHIANTI_DATA")) == 0x0) {
    if (myid==0)
      std::cout << "chdata::readMaster: "
		<< "could not find CHIANTI_DATA environment variable"
		<< " . . . exiting" << std::endl;
    MPI_Finalize();
    exit(51);
  }

  std::string fileName(val);
  fileName.append("/masterlist/masterlist.ions");
  std::string line;
  ifstream masterFile(fileName.c_str());

  if (masterFile.is_open()) {
    while(masterFile.good()) {
      getline(masterFile, line);
      std::vector<std::string> v;
      // std::cout << line <<std::endl;
      boost::split(v, line, boost::is_any_of(" "));
      masterNames.insert(v[0]);			
    }
    masterFile.close();
  }
  else {
    if (myid==0) std::cout << "chdata:readMaster: master list file <"
			   << fileName << "> not found" << std::endl;
    MPI_Finalize();
    exit(46);
  }
  
}

//
// Get the ipdata set so that if you want to get the ip of any Z, C,
// you call it as ipdata[lQ(Z, C-(int)die)]
//
void chdata::readIp() 
{
  char * val;
  if ( (val = getenv("CHIANTI_DATA")) == 0x0) {
    if (myid==0)
      std::cout << "chdata::readIp: "
		<< "could not find CHIANTI_DATA environment variable"
		<< " . . . exiting" << std::endl;
    MPI_Finalize();
    exit(52);
  }

  std::string fileName(val);
  fileName.append("/ip/chianti.ip");
  ifstream ipFile(fileName.c_str());
  
  int count = 0;
  unsigned char Z, C;
  double ip;
  double convert = 1.239841875e-4;
  std::string lineIn;

  if (ipFile.is_open() ) {

    while (ipFile.good() and count < 365) {
      getline(ipFile, lineIn);
      
      std::vector <std::string> v;
      // split the string away from the white space to get the values
      istringstream iss(lineIn);
      copy(istream_iterator<std::string>(iss), istream_iterator<std::string>(), 
	   back_inserter<vector<std::string> >(v));
      // assign the values from the string
      Z  = atoi(v[0].c_str());
      C  = atoi(v[1].c_str());
      ip = atof(v[2].c_str())*convert; // Convert to eV
      // set up the array
      ipdata[lQ(Z, C)] = ip;
      count++;
    }
    ipFile.close();
  }
  else {
    if (myid==0) std::cout << "chdata:readIp: file <" 
			   << fileName << "> not found" << std::endl;
    MPI_Finalize();
    exit(47);
  }
}

//
// Read in the abundance file, in this situation, just for test using
// the cosmic.abund file. Can later put in a multidimensional array to
// allow for all the abundance files
//
void chdata::readAbundanceAll() 
{
  char * val;
  if ( (val = getenv("CHIANTI_DATA")) == 0x0) {
    if (myid==0)
      std::cout << "chdata::readAbundanceAll: "
		<< "could not find CHIANTI_DATA environment variable"
		<< " . . . exiting" << std::endl;
    MPI_Finalize();
    exit(53);
  }

  std::string fileName(val);
  fileName.append("/abundance/cosmic_1973_allen.abund");
  
  ifstream abFile(fileName.c_str());
  std::string lineIn;
  if (abFile.is_open() ) {
    while(abFile.good()) {
      getline(abFile, lineIn);
      std::vector <std::string> v;
      istringstream iss(lineIn);
      copy(istream_iterator<std::string>(iss), istream_iterator<std::string>(),
	   back_inserter<vector<std::string> >(v));
      if (atoi(v[0].c_str()) == -1) break;
      unsigned char Z = atoi(v[0].c_str());
      abundanceAll[Z-1] = atof(v[1].c_str());
    }
    abFile.close();
  }
  else {
    if (myid==0) std::cout << "chdata::readAbundanceAll: "
			   << "abundance file <" 
			   << fileName << "> not found! " << std::endl;
    MPI_Finalize();
    exit(48);
  }
  
}

//
// Read in the Verner-Yakovlev data for radiative cross section
// determination using the short table provided by CHIANTI
//
void chdata::readVerner() 
{
  VernerXC.initialize(this);
}

//
// Read in the Karzas-Latter gaunt factor data table provided by
// CHIANTI
//
void chdata::readRadGF() 
{
  radGF.initialize(this);
}

//
// list names of all species to stdout
//
void chdata::printMaster() 
{
  if (myid==0) {
    std::cout << "Elements in the master list: " << std::endl;
    for (std::set<std::string>::iterator 
	   i=masterNames.begin(); i!=masterNames.end(); i++) {
      std::cout << "\t" << *i << std::endl;
    }
  }
}

void chdata::printIp() 
{
  std::cout << std::string(60, '-') << std::endl
	    << std::setw( 3) << "Z" << std::setw( 3) << "C"
	    << std::setw(16) << "Energy (eV)" << std::endl;

  for (auto i : ipdata) {
    if (i.second != 0) {
      std::cout 
	<< std::setw( 3) << i.first.first
	<< std::setw( 3) << i.first.second
	<< std::setw(16) << i.second
	<< std::endl;
    }
  }

  std::cout << std::string(60, '-') << std::endl;
}

//
// chdata constructor
//
chdata::chdata() 
{
  for (int i = 0; i < numEle; i++) abundanceAll[i] = 0;
  
  // std::cout << "Reading ip file\n";
  readIp();

  // std::cout << "Reading master file\n";
  readMaster();

  // std::cout << "Reading abundance file\n";
  readAbundanceAll();
  
  // std::cout << "Reading radiative cross section file\n";
  readVerner();

  // std::cout << "Reading radiative recombination Gaunt factor file\n";
  readRadGF();

  // Done
}

void chdata::createIonList(const std::set<unsigned short>& ZList)
{
  // Fill the Chianti data base
  //
  for (auto i : ZList) {
    for (int j=1; j<i+2; j++) {
      lQ Q(i, j);
      IonList[Q] = IonPtr(new Ion(i, j, this));
      // IonList[Q].freeFreeUltrarel();
    }
    Ni[i] = 1.0;		// Not sure what this does . . . 
  }

#if HAVE_LIBCUDA==1
  cuda_initialize();
#endif
}


void chianti_data::sync_string(std::string &s)
{
  int ssz = s.size();
  MPI_Bcast(&ssz, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (myid != 0) s.resize(ssz);
  MPI_Bcast(&s[0], ssz, MPI_CHAR, 0, MPI_COMM_WORLD);
}

void chianti_data::sync_vector(std::vector<double> &v)
{
  int ssz = v.size();
  MPI_Bcast(&ssz, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (myid!=0) v.resize(ssz);
  MPI_Bcast(&v[0], ssz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}


void elvlc_data::synchronize()
{
  MPI_Bcast(&level,  1, MPI_INT,      0, MPI_COMM_WORLD);

  sync_string(designation);

  MPI_Bcast(&spin,   1, MPI_INT,      0, MPI_COMM_WORLD);
  MPI_Bcast(&l,      1, MPI_INT,      0, MPI_COMM_WORLD);

  sync_string(l_str);

  MPI_Bcast(&J,      1, MPI_DOUBLE,   0, MPI_COMM_WORLD);
  MPI_Bcast(&encm,   1, MPI_DOUBLE,   0, MPI_COMM_WORLD);
  MPI_Bcast(&enry,   1, MPI_DOUBLE,   0, MPI_COMM_WORLD);
  MPI_Bcast(&encmth, 1, MPI_DOUBLE,   0, MPI_COMM_WORLD);
};

void wgfa_data::synchronize()
{
  MPI_Bcast(&lvl1,   1, MPI_INT,      0, MPI_COMM_WORLD);
  MPI_Bcast(&lvl2,   1, MPI_INT,      0, MPI_COMM_WORLD);
  MPI_Bcast(&wvl,    1, MPI_DOUBLE,   0, MPI_COMM_WORLD);
  MPI_Bcast(&gf,     1, MPI_DOUBLE,   0, MPI_COMM_WORLD);
  MPI_Bcast(&avalue, 1, MPI_DOUBLE,   0, MPI_COMM_WORLD);
};

void pe_data::synchronize()
{
  MPI_Bcast(&ngfb,   1, MPI_INT,      0, MPI_COMM_WORLD);
  MPI_Bcast(&nphot,  1, MPI_INT,      0, MPI_COMM_WORLD);

  sync_vector(pe);
};


void fblvl_data::synchronize()
{
  MPI_Bcast(&lvl,    1, MPI_INT,      0, MPI_COMM_WORLD);

  sync_string(conf);

  MPI_Bcast(&pqn,    1, MPI_INT,      0, MPI_COMM_WORLD);
  MPI_Bcast(&l,      1, MPI_INT,      0, MPI_COMM_WORLD);

  sync_string(l_str);

  MPI_Bcast(&mult,   1, MPI_INT,      0, MPI_COMM_WORLD);
  MPI_Bcast(&encm,   1, MPI_DOUBLE,   0, MPI_COMM_WORLD);
  MPI_Bcast(&encmth, 1, MPI_DOUBLE,   0, MPI_COMM_WORLD);
};

void klgfb_data::synchronize()
{
  MPI_Bcast(&n,      1, MPI_INT,      0, MPI_COMM_WORLD);
  MPI_Bcast(&l,      1, MPI_INT,      0, MPI_COMM_WORLD);

  sync_vector(factors);
};

void gffint_data::synchronize()
{
  MPI_Bcast(&ngffint, 1, MPI_INT,      0, MPI_COMM_WORLD);

  sync_vector(g2);
  sync_vector(gffint);
  sync_vector(s1);
  sync_vector(s2);
  sync_vector(s3);
};

void splups_data::synchronize()
{
  MPI_Bcast(&Z,       1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(&C,       1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(&i,       1, MPI_INT,           0, MPI_COMM_WORLD);
  MPI_Bcast(&j,       1, MPI_INT,           0, MPI_COMM_WORLD);
  MPI_Bcast(&nt,      1, MPI_INT,           0, MPI_COMM_WORLD);
  MPI_Bcast(&type,    1, MPI_INT,           0, MPI_COMM_WORLD);
  MPI_Bcast(&gf,      1, MPI_DOUBLE,        0, MPI_COMM_WORLD);
  MPI_Bcast(&Lim,     1, MPI_DOUBLE,        0, MPI_COMM_WORLD);
  MPI_Bcast(&delERyd, 1, MPI_DOUBLE,        0, MPI_COMM_WORLD);
  MPI_Bcast(&Const,   1, MPI_DOUBLE,        0, MPI_COMM_WORLD);

  sync_vector(spline);
};

void di_data::synchronize()
{
  MPI_Bcast(&btf,     1, MPI_DOUBLE,        0, MPI_COMM_WORLD);
  MPI_Bcast(&ev,      1, MPI_DOUBLE,        0, MPI_COMM_WORLD);

  sync_vector(xspline);
  sync_vector(yspline);
};

void Ion::di_head::synchronize()
{
  MPI_Bcast(&Z,       1, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&C,       1, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nspline, 1, MPI_INT,            0, MPI_COMM_WORLD);
  MPI_Bcast(&nfac,    1, MPI_INT,            0, MPI_COMM_WORLD);
  MPI_Bcast(&neav,    1, MPI_INT,            0, MPI_COMM_WORLD);
}


/** 
   Reads the Verner & Yakovlev (A&AS 109, 125, 1995) photoionization
   cross-section data
*/
void VernerData::VernerRec::sync(int nid) 
{
  MPI_Bcast(&n,    1, MPI_INT,    nid, MPI_COMM_WORLD);
  MPI_Bcast(&io,   1, MPI_INT,    nid, MPI_COMM_WORLD);
  MPI_Bcast(&pql,  1, MPI_INT,    nid, MPI_COMM_WORLD);
  MPI_Bcast(&l,    1, MPI_INT,    nid, MPI_COMM_WORLD);
  MPI_Bcast(&eth,  1, MPI_DOUBLE, nid, MPI_COMM_WORLD);
  MPI_Bcast(&e0,   1, MPI_DOUBLE, nid, MPI_COMM_WORLD);
  MPI_Bcast(&sig0, 1, MPI_DOUBLE, nid, MPI_COMM_WORLD);
  MPI_Bcast(&ya,   1, MPI_DOUBLE, nid, MPI_COMM_WORLD);
  MPI_Bcast(&p,    1, MPI_DOUBLE, nid, MPI_COMM_WORLD);
  MPI_Bcast(&yw,   1, MPI_DOUBLE, nid, MPI_COMM_WORLD);
}

void VernerData::initialize(chdata* ch)
{
  this->ch = ch;
  int nOK = 0;
  
  // Only used by root process
  //@{
  unsigned nVern = 465, nVernExt = 1677;
  bool extended = false;
  //@}
  
  std::vector<vrPtr> indat;


  if (myid==0) {
    
    char * val;

    // Look for extended Verner file
    //
    if ( (val = getenv("VERNER_DATA")) != 0x0) {
      extended = true;
    } else {
      std::cout << "Could not find VERNER_DATA environment variable"
		<< "for the extended Verner-Yakovlev table . . . "
		<< "will look for the CHIANTI version" << std::endl;
    }

    // Use CHIANTI version of extended version is not found
    //
    if (not extended and (val = getenv("CHIANTI_DATA")) == 0x0) {
      std::cout << "Could not find CHIANTI_DATA environment variable"
		<< " . . . exiting" << std::endl;
      nOK = 1;
    } else {
      nOK = 0;
    }

    if (nOK == 0) {

      std::string inLine, fileName(val);
      if (not extended) fileName.append("/continuum/verner_short.txt");

      std::ifstream vdFile(fileName);
      
      if (vdFile.is_open()) {
	
	while (vdFile.good()) {
	  
	  std::vector<std::string> v;
	  std::getline(vdFile, inLine);
	  
	  if (extended) {
	    // Handle header
	    if (inLine.find("ELE") != std::string::npos)
	      std::getline(vdFile, inLine);
	    // Handle eof
	    if (inLine.size()>2) {
	      std::fill(&inLine[10], &inLine[29], ' ');
	      inLine = inLine.substr(2, std::string::npos);
	    }
	  }

	  std::istringstream iss(inLine);
	  std::copy(std::istream_iterator<std::string>(iss), 
		    std::istream_iterator<std::string>(), 
		    std::back_inserter<std::vector<std::string> >(v));
	  
	  if (v.size() < 10) break;
	  
	  vrPtr dat(new VernerRec);
	  dat->n    = atoi(v[0].c_str());
	  dat->io   = atoi(v[1].c_str());
	  dat->pql  = atoi(v[2].c_str());
	  dat->l    = atoi(v[3].c_str());
	  dat->eth  = atof(v[4].c_str());
	  dat->e0   = atof(v[5].c_str());
	  dat->sig0 = atof(v[6].c_str());
	  dat->ya   = atof(v[7].c_str());
	  dat->p    = atof(v[8].c_str());
	  dat->yw   = atof(v[9].c_str());
	  
	  indat.push_back(dat);
	}

	vdFile.close();
      } else {
	std::cout << "Could not open the Verner data file <"
		  << fileName << ">" << std::endl;
	nOK = 1;
      }
    }
  }

  
  MPI_Bcast(&nOK, 1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  if (nOK)  {
    if (myid == 0) {
      std::cerr << "VernerData: problem reading VERNER data file" << std::endl;
    }
    MPI_Finalize();
    exit(59);
  } 

  if (myid==0) {
    
    if (extended) {
      if (indat.size() != nVernExt) {
	std::cout << "Root node: Verner short data size=" << data.size() 
		  << ", expected: " << nVernExt << std::endl;
      }
    } else {
      if (indat.size() != nVern) {
	std::cout << "Root node: Verner extended data size=" << data.size() 
		  << ", expected: " << nVern << std::endl;
      }
    }
    
    for (auto v : indat) {

      v->sync();

      int Z = v->n;
      int C = v->n - v->io + 1;
	  
      // Stage: 1 is neutral, etc.

      lQ key(Z, C);

      data[key].push_back(v);
    }
    
    unsigned short done = 0;
    MPI_Bcast(&done,  1, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);
    
  } else {
    lQ key;
    
    while (1) {
      vrPtr dat(new VernerRec);
      dat->sync();

      int Z = dat->n;
      int C = dat->n - dat->io + 1;
	  
      lQ key(Z, C);

      data[key].push_back(dat);
    }
    
    if (data.size() != nVern)
      std::cout << "Node " << myid << ": "
		<< " Verner data size=" << data.size() 
		<< ", expected: " << nVern << std::endl;
  }
}


/**
   Compute the radiative recombination cross section from the
   photoionization cross section.
   
   The data are from Verner and Yakolev (1995).
   
   The cross section is evaluated for the next lower ionization
   stage.
*/

double VernerData::cross(const lQ& Q, double EeV)
{
				// Rydberg in eV
  constexpr double ryd = 27.2113845/2.0;

				// 1 inverse cm = 1.239.. eV
  constexpr double incmEv = 1.0/8.06554465e+03;

				// Electron rest mass in eV
  constexpr double mec2 = 510.998896 * 1.0e3;
  
  // The ion after recombination
  //
  lQ rQ(Q.first, Q.second-1);

  // No data for this ion
  //
  if (data.find(rQ) == data.end()) return 0.0;

  double mult0  = 2.0;
  if (Q.second > Q.first) mult0 = 1.0;
  
  double vCross = 0.0;
  
  if (Ion::no_verner and Q.first==1 and Q.second==2) {

    for (int n=1; n<11; n++) {
    
      double Eiz = ryd/(n*n);	// These energies are now in eV
      double Eph = EeV + Eiz;

      for (int l=0; l<n; l++) {
				// Milne relation
	double Milne = 0.5*Eph*Eph/(mec2*EeV)*2.0*(2.0*l+1.0)/mult0;

	double crossPh = 7.90706903681e-04 * pow(Eiz/ryd, 2.0)*
	  pow(ryd/Eph, 3.0)/n;

				// Gaunt factor
	double scaledE = log(Eph/Eiz);
	double gf = ch->radGF(scaledE, n, l);
				// Cross section x Gaunt factor
	double cross = crossPh * gf * Milne;
	
	vCross += cross;
      }
    }

  } else {

    for (auto u : data[rQ]) {
      double Eiz = u->eth;
      double Eph = EeV + Eiz;
				// Milne relation
      double Milne = 0.5*Eph*Eph/(mec2*EeV)*2.0*(2.0*u->l+1.0)/mult0;

      double crossPh = crossPhotoIon(u, Eph);
      double cross   = crossPh * Milne;

      vCross += cross;
    }
  }

  return vCross;
}

// 
// Ground-state photoionization cross section
//
std::vector<double> Ion::photoIonizationCross(double E, int id)
{
  double val = 0.0;

  // Call for the cross section
  //
  lQ Q(Z, C);
  auto ion = ch->VernerXC.data.find(Q);
  if (ion != ch->VernerXC.data.end()) {
    val = ch->VernerXC.crossPhotoIon(*(ion->second.begin()), E);
  }
  
  std::vector<double> ret(1, val);
  return ret;
}

void Ion::IBcreate()
{
  if (IBinit) return;

  // For conversion of nm^2 to cm^2
  //           |
  //           v
  double ff =  1.0e-14 * ip * nudel / (heV*eVtoErg) * 4.0 * M_PI;

  IBrate.resize(nusteps);
  IBtotl = 0.0;
  for (int n=0; n<nusteps; n++) {
    double E  = ip * nugrid[n];
    double x  = ff * photoBackground(E) * photoIonizationCross(E, 0)[0] / E;
    IBrate[n] = x;
    IBtotl   += x;
  }

  IBcum.resize(nusteps);
  IBcum[0] = IBrate[0]/IBtotl;
  for (int n=1; n<nusteps; n++) IBcum[n] = IBrate[n]/IBtotl + IBcum[n-1];

  IBinit = true;
}

//
// Photoionizing background (photon energy in eV)
//
// Quasar exponent and value following Haardt and Madau (2012)
//
double Ion::photoBackground(double Enu)
{
  if (ib_type==Ion::none)
    return 0.0;
  else
    return HandM_coef*pow(Enu/RydtoeV, HandM_expon);
}

//
// Photoionizaton with background
//
std::pair<double, double> Ion::photoIonizationRate()
{
  if (ib_type == Ion::none or C > Z) {
    return std::pair<double, double>(0.0, 0.0);
  } else {
    if (not IBinit) IBcreate();

    // Location in cumulative cross section grid
    //
    double rn = static_cast<double>(rand())/RAND_MAX;
  
    // Interpolate the cross section array
    //
    
    // Points to first element that is not < rn
    // but may be equal
    std::vector<double>::iterator lb = 
      std::lower_bound(IBcum.begin(), IBcum.end(), rn);
    
    // Assign upper end of range to the
    // found element
    //
    std::vector<double>::iterator ub = lb;
    //
    // If is the first element, increment
    // the upper boundary
    //
    if (lb == IBcum.begin()) { if (IBcum.size()>1) ub++; }
    //
    // Otherwise, decrement the lower boundary
    //
    else { lb--; }
    
    // Compute the associated indices
    //
    size_t ii = lb - IBcum.begin();
    size_t jj = ub - IBcum.begin();
    double nu = nugrid[ii];
	  
    // Linear interpolation
    //
    if (*ub > *lb) {
      double d = *ub - *lb;
      double a = (rn - *lb) / d;
      double b = (*ub - rn) / d;
      nu  = a * nugrid[ii] + b * nugrid[jj];
    }
    
    return std::pair<double, double>(IBtotl, (nu-1.0)*ip);
  }
}

double VernerData::crossPhotoIon(vrPtr vdata, double Eph)
{
  if (Ion::use_VFKY) return crossPhotoIon_VFKY(vdata, Eph);

  double y    = Eph/vdata->e0;
  double y1   = y - 1.0;
    
  // Verner and Yakolev, equation 1
  //
  double fy  = vdata->sig0*(y1*y1 + vdata->yw*vdata->yw) * 
    pow(y, -5.5 - vdata->l + 0.5*vdata->p) * 
    pow(1.0 + sqrt(y/vdata->ya), -vdata->p);
  
  // Convert from Mbarnes to nm^2
  //
  return fy * 1.0e-4;
}



// Read CHIANTI files file containing the free-bound gaunt factors for
// n=1-6 from Karzas and Latter, 1961, ApJSS, 6, 167, the photon
// energy and the free-bound gaunt factors
//
void KLGFdata::initialize(chdata* ch)
{
  this->ch = ch;

  int nOK = 0;

  if (myid==0) {
    
    char * val;
    if ( (val = getenv("CHIANTI_DATA")) == 0x0) {
      std::cout << "Could not find CHIANTI_DATA environment variable"
		<< " . . . exiting" << std::endl;
      nOK = 1;
    }
    
    if (nOK == 0) {
      
      std::string fileName(val);
      
      fileName.append("/continuum/klgfb.dat");
      
      std::string inLine;
      std::ifstream klgfFile(fileName.c_str());
      
      if (klgfFile.is_open()) {
	
	std::getline(klgfFile, inLine);
	{
	  std::istringstream sin(inLine);
	  sin >> ngfb;
	  sin >> nume;
	}

	pe.resize(nume);

	std::getline(klgfFile, inLine);
	{
	  std::istringstream sin(inLine);
	  for (auto & v : pe) sin >> v;
	}

	
	while (klgfFile.good()) {

	  std::getline(klgfFile, inLine);
	  std::istringstream sin(inLine);

	  int n, l;
	  sin >> n;
	  sin >> l;
	  
	  std::pair<int, int> key(n, l);

	  while (sin.good()) {
	    double V;
	    sin >> V;
	    if (sin.good() or sin.eof())
	      gfb[key].push_back(V);
	    else
	      break;
	  }
	}
      }
      klgfFile.close();
    }
  }
  
  MPI_Bcast(&nOK, 1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  if (nOK)  {
    if (myid == 0) {
      std::cerr << "KLGFdata: problem reading CHIANTI files" << std::endl;
    }
    MPI_Finalize();
    exit(60);
  } 

  if (myid==0) {
    
    int sz = gfb.size();

    MPI_Bcast(&ngfb, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nume, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sz,   1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Bcast(&pe[0], nume, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (auto m : gfb) {
      int n   = m.first.first;
      int l   = m.first.second;
      int gsz = m.second.size();

      MPI_Bcast(&n,   1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&l,   1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&gsz, 1, MPI_INT, 0, MPI_COMM_WORLD);

      MPI_Bcast(&m.second[0], gsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    
  } else {
  
    int sz, n, l, gsz;

    MPI_Bcast(&ngfb, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nume, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sz,   1, MPI_INT, 0, MPI_COMM_WORLD);
      
    pe.resize(nume);
    MPI_Bcast(&pe[0], nume, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int i=0; i<sz; i++) {

      MPI_Bcast(&n,   1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&l,   1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&gsz, 1, MPI_INT, 0, MPI_COMM_WORLD);
      
      std::pair<int, int> key(n, l);
      gfb[key].resize(gsz);

      MPI_Bcast(&gfb[key][0], gsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
  }

  // Ok, now make and store the splines
  //
  for (auto v : gfb) {
    tksplPtr s(new tk::spline);	 // New spline instance
    s->set_points(pe, v.second); // Initialize
    spl[v.first] = s;		 // Store
  }

  // We can safely clear the gfb structure now . . .
  //
  gfb.clear();
}

std::map<unsigned short, double> 
chdata::fraction(unsigned short Z, double T, int norder)
{
  if (Lagu.get() == 0) 
    Lagu = boost::shared_ptr<LaguQuad>(new LaguQuad(norder, 0.0));
  else if (Lagu->n != norder) 
    Lagu = boost::shared_ptr<LaguQuad>(new LaguQuad(norder, 0.0));

  std::map<unsigned short, double> ret;

  std::map<lQ, double> ionize, recomb;

  for (unsigned short C=1; C<=Z+0; C++) ionize[lQ(Z, C)] = 0.0;
  for (unsigned short C=2; C<=Z+1; C++) recomb[lQ(Z, C)] = 0.0;

  const double beta = 0.5*me/(boltz*T);
  const double K    = 2.0/sqrt(M_PI*beta);
  const double kTeV = boltzEv * T;

  for (auto & v : ionize) v.second = 0.0;
  for (auto & v : recomb) v.second = 0.0;

  for (int i=1; i<=norder; i++) {

    double y = Lagu->knot(i);
    double w = Lagu->weight(i);
    
    // Ionization
    //
    for (auto & v : ionize) {

      double ab   = ipdata[v.first]/kTeV; // Scaled ionization energy
      
      v.second += w * exp(-ab) * K * (y + ab) *
	IonList[v.first]->directIonCross((y + ab)*kTeV, 0);
    }

    // Recombination
    //
    for (auto & v : recomb) {
      v.second += w * K * y *
	IonList[v.first]->radRecombCross(y*kTeV, 0).back();
    }
  }
  
  std::vector<double> nn(Z+1, 1);
  double norm = 1.0;
  for (unsigned short C=1; C<Z+1; C++) {
    if (recomb[lQ(Z,C+1)]>0.0)
      nn[C] = nn[C-1] * ionize[lQ(Z,C)]/recomb[lQ(Z,C+1)];
    else
      nn[C] = 0.0;
    norm += nn[C];
  }
  for (unsigned short C=1; C<Z+2; C++) ret[C] = nn[C-1]/norm;
    
  return ret;
}


std::map<unsigned short, double> 
chdata::fraction(unsigned short Z, double T, 
		 double Emin, double Emax, int norder)
{
  if (Lege.get() == 0) 
    Lege = boost::shared_ptr<LegeQuad>(new LegeQuad(norder));
  else if (Lege->n != norder) 
    Lege = boost::shared_ptr<LegeQuad>(new LegeQuad(norder));

  std::map<unsigned short, double> ret;

  std::map<lQ, double> ionize, recomb;

  for (unsigned short C=1; C<=Z+0; C++) ionize[lQ(Z, C)] = 0.0;
  for (unsigned short C=2; C<=Z+1; C++) recomb[lQ(Z, C)] = 0.0;

  for (auto & v : ionize) v.second = 0.0;
  for (auto & v : recomb) v.second = 0.0;

  const double beta = 0.5*me/(boltz*T);
  const double K    = 2.0/sqrt(M_PI*beta);
  const double kTeV = boltzEv * T;

  // Ionization
  //
  for (auto & v : ionize) {
				// ionization energy in eV
    double emin = std::max<double>(Emin, ipdata[v.first]);

    if (emin < Emax) {

      double ymax = Emax/kTeV;	// Emin and Emax are in eV
      double ymin = emin/kTeV;
      double dy   = ymax - ymin;

      for (int i=1; i<=norder; i++) {
	double y = ymin + dy*Lege->knot(i);
	v.second += Lege->weight(i)*dy * K * y*exp(-y) *
	  IonList[v.first]->directIonCross(y*kTeV, 0);
      }
    }
  }

  // Recombination
  //
  if (Emin < Emax) {
    
    double ymax = Emax/kTeV;
    double ymin = Emin/kTeV;
    double dy   = ymax - ymin;

    for (auto & v : recomb) {
      for (int i=1; i<=norder; i++) {
	double y = ymin + dy*Lege->knot(i);
	v.second += Lege->weight(i)*dy * K * y*exp(-y) *
	  IonList[v.first]->radRecombCross(y*kTeV, 0).back();
      }
    }
  }
  
  std::vector<double> nn(Z+1, 1);
  double norm = 1.0;
  for (unsigned short C=1; C<Z+1; C++) {
    if (recomb[lQ(Z,C+1)]>0.0)
      nn[C] = nn[C-1] * ionize[lQ(Z,C)]/recomb[lQ(Z,C+1)];
    else
      nn[C] = 0.0;
    norm += nn[C];
  }
  for (unsigned short C=1; C<Z+2; C++) ret[C] = nn[C-1]/norm;
    
  return ret;
}

std::map<unsigned short, std::vector<double> > 
chdata::recombEquil(unsigned short Z, double T, int norder)
{
  if (Lagu.get() == 0) 
    Lagu = boost::shared_ptr<LaguQuad>(new LaguQuad(norder, 0.0));
  else if (Lagu->n != norder) 
    Lagu = boost::shared_ptr<LaguQuad>(new LaguQuad(norder, 0.0));

  std::map<unsigned short, std::vector<double> > ret;

  std::map<lQ, double> ionize, recomb;

  for (unsigned short C=1; C<=Z+0; C++) ionize[lQ(Z, C)] = 0.0;
  for (unsigned short C=2; C<=Z+1; C++) recomb[lQ(Z, C)] = 0.0;

  const double beta = 0.5*me/(boltz*T);
  const double K    = 2.0/sqrt(M_PI*beta);
  const double kTeV = boltzEv*T;

  for (auto & v : ionize) v.second = 0.0;
  for (auto & v : recomb) v.second = 0.0;

  for (int i=1; i<=norder; i++) {

    double y = Lagu->knot(i);
    double w = Lagu->weight(i);
    
    // Ionization
    //
    for (auto & v : ionize) {
				// scaled ionization energy
      double ab   = ipdata[v.first]/kTeV;
      
      v.second += w * exp(-ab) * K * (y + ab) *
	IonList[v.first]->directIonCross((y + ab)*kTeV, 0);
    }

    // Recombination
    //
    for (auto & v : recomb) {
      v.second += w * K * y *
	IonList[v.first]->radRecombCross(y*kTeV, 0).back();
    }
  }
  
  std::vector<double> nn(Z+1, 1);
  double norm = 1.0;
  for (unsigned short C=1; C<Z+1; C++) {
    if (recomb[lQ(Z,C+1)]>0.0)
      nn[C] = nn[C-1] * ionize[lQ(Z,C)]/recomb[lQ(Z,C+1)];
    else
      nn[C] = 0.0;
    norm += nn[C];
  }

  ret[1] = {nn[0]/norm, ionize[lQ(Z, 1)], 0.0};
  for (unsigned short C=2; C<Z+1; C++) 
    ret[C] = {nn[C-1]/norm, ionize[lQ(Z, C)], recomb[lQ(Z, C)]};
  ret[Z+1] = {nn[Z]/norm, 0.0, recomb[lQ(Z, Z+1)]};
    
  return ret;
}


std::map<unsigned short, std::vector<double> >
chdata::recombEquil(unsigned short Z, double T, 
		    double Emin, double Emax, int norder, bool use_log)
{
  if (Lege.get() == 0) 
    Lege = boost::shared_ptr<LegeQuad>(new LegeQuad(norder));
  else if (Lege->n != norder) 
    Lege = boost::shared_ptr<LegeQuad>(new LegeQuad(norder));

  std::map<unsigned short, std::vector<double> > ret;

  std::map<lQ, double> ionize, recomb;

  for (unsigned short C=1; C<=Z+0; C++) ionize[lQ(Z, C)] = 0.0;
  for (unsigned short C=2; C<=Z+1; C++) recomb[lQ(Z, C)] = 0.0;

  for (auto & v : ionize) v.second = 0.0;
  for (auto & v : recomb) v.second = 0.0;

  const double beta = 0.5*me/(boltz*T);
  const double K    = 2.0/sqrt(M_PI*beta);
  const double kTeV = boltzEv*T;

  if (use_log and Emin <= 0.0) use_log = false;

  // Ionization
  //
  for (auto & v : ionize) {
				// ionization energy in eV
    double emin = std::max<double>(Emin, ipdata[v.first]);

    if (emin < Emax) {

      double ymax = Emax/kTeV;	// Scaled min energy
      double ymin = emin/kTeV;	// Scaled max energy

      if (use_log) {
	ymin = log(ymin);
	ymax = log(ymax);
      }

      double dy = ymax - ymin;

      for (int i=1; i<=norder; i++) {
	if (use_log) {
	  double y = exp(ymin + dy*Lege->knot(i));
	  v.second += Lege->weight(i) * dy * K * y*y*exp(-y) *
	    IonList[v.first]->directIonCross(y*kTeV, 0);
	} else {
	  double y = ymin + dy*Lege->knot(i);
	  v.second += Lege->weight(i) * dy * K * y*exp(-y) *
	    IonList[v.first]->directIonCross(y*kTeV, 0);
	}
      }
    }
  }

  // Recombination
  //
  if (Emin < Emax) {

    double ymax = Emax/kTeV;
    double ymin = Emin/kTeV;
    double dy   = ymax - ymin;

    for (auto & v : recomb) {
      for (int i=1; i<=norder; i++) {
	double y = ymin + dy*Lege->knot(i);
	v.second += Lege->weight(i) * dy * K * y*exp(-y) *
	  IonList[v.first]->radRecombCross(y*kTeV, 0).back();
      }
    }
  }
  
  std::vector<double> nn(Z+1, 1);
  double norm = 1.0;
  for (unsigned short C=1; C<Z+1; C++) {
    if (recomb[lQ(Z,C+1)]>0.0)
      nn[C] = nn[C-1] * ionize[lQ(Z,C)]/recomb[lQ(Z,C+1)];
    else
      nn[C] = 0.0;
    norm += nn[C];
  }
    
  ret[1] = {nn[0]/norm, ionize[lQ(Z, 1)], 0.0};
  for (unsigned short C=2; C<Z+1; C++) 
    ret[C] = {nn[C-1]/norm, ionize[lQ(Z, C)], recomb[lQ(Z, C)]};
  ret[Z+1] = {nn[Z]/norm, 0.0, recomb[lQ(Z, Z+1)]};
  return ret;
}


/** CUDA strategy

(1) All cross sections with lookup tables to be copied to cuda texture objects

(2) Copy necessary constants to cuda symbol data

(3) Cross sections are then evaluated on the device so we need cuda
__device__ functions for each desired cross section

  Four grids
  ----------
  free free    2-dimensional
  rad recomb   1-dimensional
  col excite   1-dimensional (2 values)
  col ionize   1-dimensional
*/

#if HAVE_LIBCUDA==1

void chdata::cuda_initialize()
{
  // Remove any previous initialization data
  //
  destroy_cuda();
  cuIonElem.clear();

  for (auto v : IonList) {
    cuIonElement E;
				// Initialize cross-section data structures
				// 
    E.Z = v.first.first;
    E.C = v.first.second;
    E.I = E.NColl  = E.NIonz = 0;
    E.ceEmin = E.ceEmax = E.ceDelE = 0;
    E.ciEmin = E.ciEmax = E.ciDelE = 0;
    E.ff_0 = E.ff_d = E.rc_d = E.ce_d = E.ci_d = 0;
    E.IPval = ipdata[lQ(E.Z, E.C)];

    cuIonElem.push_back(E);
    
				// Make grids
				// 
    if (not v.second->freeFreeGridComputed)
      v.second->freeFreeMakeEvGrid(0);

    if (v.first.second>1 and not v.second->radRecombGridComputed)
      v.second->radRecombMakeEvGrid(0);

    if (v.first.second<=v.first.first) {
      if (not v.second->exciteGridComputed) v.second->collExciteMakeGrid(0);
      if (not v.second->ionizeGridComputed) v.second->directIonMakeGrid(0);
    }
  }

  cuda_initialize_grid_constants();
  cuda_initialize_textures();
}

void chdata::destroy_cuda()
{
  int i = 0;
  for (auto & u : cuIonElem) {
    for (auto & v : {u.ff_0, u.ff_d, u.rc_d, u.ce_d, u.ci_d} ) {
      std::ostringstream sout;
      sout << "trying to free TextureObject [" << i++
	   << "] with [Z, C]=[" << u.Z << ", " << u.C << "]";
      cuda_safe_call(cudaDestroyTextureObject(v),
		     __FILE__, __LINE__, sout.str());
    }
  }
}

#endif
