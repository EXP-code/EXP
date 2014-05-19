#include <boost/algorithm/string.hpp>
#include <string>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <map>

#include <localmpi.h>

#include "Ion.H"
#include "interactSelect.H"
#include "Cspline.H"

//! Convert the master element name to a (Z, C) pair
void Ion::convertName() 
{
  std::string ele;
  std::string charge;
  
  std::vector<std::string> v;
  std::string die = "d";	// Set to dielectronic
  size_t isd;
  
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

//! Convert a given Z,C pair into a master name string
std::string ZCtoName(unsigned char Z, unsigned char C) 
{
  std::string fileName;
  std::string ele;
  std::string c;
  
  // Compile the two parts of the name and "stitch" together
  ele = eleNameList[Z-1];
  c = static_cast<ostringstream*>( &(ostringstream() << C) )->str();
  std::stringstream ss;
  ss << ele << "_" << c;
  fileName = ss.str();
  
  return fileName;
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
  char * val;
  if ( (val = getenv("CHIANTI_DATA")) == 0x0) {
    if (myid==0)
      std::cout << "Could not find CHIANTI_DATA environment variable"
		<< " . . . exiting" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 47);
  }

  std::string fileName(val);

  fileName.append("/");
  fileName.append(eleName); 
  fileName.append("/");
  fileName.append(MasterName); 
  fileName.append("/"); 
  fileName.append(MasterName);
  fileName.append(".elvlc");
  
  std::string inLine;
  elvlc_data e;
  ifstream elvlcFile(fileName.c_str());

  if (elvlcFile.is_open()) {

    while (elvlcFile.good()) {

      std::vector <std::string> v;
      getline(elvlcFile, inLine);
      std::istringstream iss(inLine);
      copy(istream_iterator<std::string>(iss), istream_iterator<std::string>(), 
	   back_inserter<vector<std::string> >(v));

      if (atoi(v[0].c_str()) == -1) break;

      e.level = atoi(v[0].c_str());
      e.conf = atoi(v[1].c_str());
      e.designation = v[2];
      e.spin = atoi(v[3].c_str());
      e.l = atoi(v[4].c_str());
      e.l_str = v[5];
      e.J = atof(v[6].c_str());
      e.mult = atoi(v[7].c_str());
      e.encm = atof(v[8].c_str());
      e.enry = atof(v[9].c_str());
      e.encmth = atof(v[10].c_str());
      e.enryth = atof(v[11].c_str());
      
      elvlc.push_back(e);
    }
    elvlcFile.close();
    nelvlc = elvlc.size();
  }
  else {
    if (myid==0) std::cout << "Cannot find file: " << fileName << std::endl;
    nelvlc = 0;
    MPI_Abort(MPI_COMM_WORLD, 42);
  }
}

//! Read in the fblvl file found in the CHIANTI database
void Ion::readfblvl() 
{
  std::string MasterNameT = ZCtoName(Z, C-1);

  char * val;
  if ( (val = getenv("CHIANTI_DATA")) == 0x0) {
    if (myid==0)
      std::cout << "Could not find CHIANTI_DATA environment variable"
		<< " . . . exiting" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 48);
  }

  std::string fileName(val);

  fileName.append("/");
  fileName.append(eleName); 
  fileName.append("/");
  fileName.append(MasterNameT); 
  fileName.append("/"); 
  fileName.append(MasterNameT);
  fileName.append(".fblvl");
  
  std::string inLine;
  ifstream fblvlFile(fileName.c_str());
  
  fblvl_data f;
  if (fblvlFile.is_open()) {

    while (fblvlFile.good()) {

      std::vector <std::string> v;
      getline(fblvlFile, inLine);
      istringstream iss(inLine);
      copy(istream_iterator<std::string>(iss), istream_iterator<std::string>(), 
	   back_inserter<vector<std::string> >(v));
      if(atoi(v[0].c_str()) == -1) break;
      f.lvl = atoi(v[0].c_str());
      f.conf = v[1];
      f.pqn = atoi(v[2].c_str());
      f.l = atoi(v[3].c_str());
      f.l_str = v[4];
      f.mult = atoi(v[5].c_str());
      f.encm = atof(v[6].c_str());
      f.encmth = atof(v[7].c_str());
      
      fblvl.push_back(f);
    }
    fblvlFile.close();
    nfblvl = fblvl.size();
  }
  else {
    if (myid==0) std::cout << "Cannot find file: " << fileName << std::endl;
    nfblvl = 0;
    MPI_Abort(MPI_COMM_WORLD, 43);
  }
}

//! Read in the spline file from the CHIANTI database
void Ion::readSplups() 
{
  char * val;
  if ( (val = getenv("CHIANTI_DATA")) == 0x0) {
    if (myid==0)
      std::cout << "Could not find CHIANTI_DATA environment variable"
		<< " . . . exiting" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 49);
  }

  std::string fileName(val);

  fileName.append("/");
  fileName.append(eleName); 
  fileName.append("/");
  fileName.append(MasterName); 
  fileName.append("/"); 
  fileName.append(MasterName);
  fileName.append(".splups");
  
  std::string inLine;
  splups_data s;
  ifstream sFile(fileName.c_str());
  if (sFile.is_open()) {

    while(sFile.good()) {

      std::vector <std::string> v;
      getline(sFile, inLine);
      istringstream iss(inLine);
      copy(istream_iterator<std::string>(iss), istream_iterator<std::string>(), 
	   back_inserter<vector<std::string> >(v));

      if(atoi(v[0].c_str()) == -1) break;

      s.Z = atoi(v[0].c_str());
      s.C = atoi(v[1].c_str());
      s.i = atoi(v[2].c_str());
      s.j = atoi(v[3].c_str());
      s.type = atoi(v[4].c_str());
      s.gf = atof(v[5].c_str());
      s.delERyd = atof(v[6].c_str());
      s.Const = atof(v[7].c_str());
				// Spline coefficients
      for(unsigned i = 8; i < v.size(); i++) {
	s.spline.push_back(atof(v[i].c_str()));
      }
      
      splups.push_back(s);
      s.spline.erase(s.spline.begin(), s.spline.end());		
    }
    sFile.close();
    nsplups = splups.size();
  }
  else {
    if (myid==0) std::cout << "Cannot find file: " 
			   << fileName << std::endl;
    nsplups = 0;
    MPI_Abort(MPI_COMM_WORLD, 44);
  }
}

//! Read in the direct ionization cross section splines from the file
void Ion::readDi() 
{
  char * val;
  if ( (val = getenv("CHIANTI_DATA")) == 0x0) {
    if (myid==0)
      std::cout << "Could not find CHIANTI_DATA environment variable"
		<< " . . . exiting" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 50);
  }

  std::string fileName(val);

  fileName.append("/");
  fileName.append(eleName); 
  fileName.append("/");
  fileName.append(MasterName); 
  fileName.append("/"); 
  fileName.append(MasterName);
  fileName.append(".diparams");
  
  std::string inLine;
  di_data s;
  ifstream sFile(fileName.c_str());
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
	  s.xspline.erase(s.xspline.begin(), s.xspline.end());		
	  s.yspline.erase(s.yspline.begin(), s.yspline.end());		
	  i_fac++;
	}
      }
      i++;
    }
    sFile.close();
    ndispline = diSpline.size();
  }
  else {
    if (myid==0) std::cout << "Cannot find file: " 
			   << fileName << std::endl;
    ndispline = 0;
    MPI_Abort(MPI_COMM_WORLD, 45);
  }
}

//! Initialization function when the master name is given
Ion::Ion(std::string name, chdata ch) 
{
  MasterName = name;

  std::vector<std::string> v;
  boost::split(v, MasterName, boost::is_any_of("_") );
  eleName = v[0];

  convertName();		// Sets Z and C . . . 
  ip = ch.ipdata[Z-1][C-1];

  std::string MasterNameT = ZCtoName(Z, C-1);
  if (isInMasterList(ch, MasterNameT)) {
    readfblvl();
  }
  if (isInMasterList(ch, MasterName)) {
    // std::cout << "IN MASTER LIST" <<std::endl;
    readelvlc();
    readSplups();
    readDi();
  }
  // else std::cout << "NOT IN LIST" <<std::endl;
  
  // initialize the k-grid (in inverse nm) for ff and the energy grid (in eV)
  double e = 0.;
  double k = 0.;
  for (e = 0.00000001; e < 250 ; e += 0.25) {
    egrid.push_back(e);
  }
  for (k = -9.; k < -3.; k += 0.1) {
    // std::cout << k 
    kgrid.push_back(k);
  }
  kffsteps = kgrid.size();
  effsteps = egrid.size();
}

//! Constructor when the Z, C pair is given
Ion::Ion(unsigned char Z1, unsigned char C1, chdata ch) 
{
  d = false;
  Z = Z1;
  C = C1;
  MasterName = ZCtoName(Z1, C1);

  std::vector<std::string> v;
  boost::split(v, MasterName, boost::is_any_of("_") );
  eleName = v[0];
  ip = ch.ipdata[Z-1][C-1];

  std::string MasterNameT = ZCtoName(Z, C-1);
  if (isInMasterList(ch, MasterNameT)) {
    readfblvl();
  }

  if (isInMasterList(ch, MasterName)) {
    readSplups();
    readDi();
    readelvlc();
  }
  
  // Initialize the k-grid (in inverse nm) for ff and the energy grid
  // (in eV)
  double e = 0;
  double k = 0;
  for (e = 0.00000001; e < 250 ; e += 0.25) {
    egrid.push_back(e);
  }
  for (k = -9.; k < -3.; k += 0.1) {
    kgrid.push_back(k);
  }
  kffsteps = kgrid.size();
  effsteps = egrid.size();
}

//! Default constructor: NOT USED
Ion::Ion() 
{
  Z = 1;
  C = 1;
}

//! Copy constructor
Ion::Ion(const Ion &I) 
{
  Z          = I.Z;
  C          = I.C;
  ip         = I.ip;

  kffsteps   = I.kffsteps;
  effsteps   = I.effsteps;

  egrid      = I.egrid;
  kgrid      = I.kgrid;
  fblvl      = I.fblvl;

  diSpline   = I.diSpline;
  di_header  = I.di_header;
  nsplups    = I.nsplups;
  nfblvl     = I.nfblvl;
  ndispline  = I.ndispline;

  MasterName = I.MasterName;
  eleName    = I.eleName;
  elvlc      = I.elvlc;
  nelvlc     = I.nelvlc;
}

/** 
    Calculate the collision excitation cross sections and return the
    cumulative cross section array.

    Returns a vector with a pair of (cross section, energy difference)
    since the file input, and thus array, are not in any specific order
*/
std::vector< std::pair<double, double > > 
Ion::collExciteCross(chdata ch, double E, double Eth) 
{
  const double x_array5[5] = {0, 0.25, 0.5, 0.75, 1.0};
  const double x_array9[9] = {0, 0.125, 0.25 , 0.375, 0.5 , 
			      0.625 , 0.75, 0.875, 1.0};

  std::vector<double> x5(x_array5, x_array5+5);
  std::vector<double> x9(x_array9, x_array9+9);
  
  double eVtoRyd = 1.0/13.60569253;
  double RydtoeV = 1.0/eVtoRyd;

  double a0 = 0.0529177211; // Bohr radius in nm
  
  double totalCross = 0;
  std::vector<std::pair<double, double > > CEcum;
  std::pair<double,double> Null(0, 0);
  if (splups.size() == 0) {
    CEcum.push_back(Null);
    CEcrossCum = CEcum;
    return CEcum;
  }

  for (size_t i=0; i<splups.size(); i++) {

    double EijEv = splups[i].delERyd*RydtoeV;
    double Const = splups[i].Const;
    double cross = 0.;

    if (splups[i].i == 1  and E >= EijEv) {

      assert(splups[i].i == 1);
      assert(splups[i].spline.size() != 0);
      
      // double dE = E - EijEv;
      // Filter out the types
      
      int type = splups[i].type;
      double x=0, y=0;
      
      if (type==1) {
	x = 1.0 - (log(Const)/(log((Eth/EijEv) + Const)));
      }
      if (type == 2) {
	x = (Eth/EijEv)/((Eth/EijEv) + Const);
      }
      if (type == 3) {
	x = (Eth/EijEv)/((Eth/EijEv) + Const);	
      }
      if (type == 4) {
	x = 1.0 - (log(Const)/(log((Eth/EijEv) + Const)));
      }
      // xmin is 0 and xmax is 1, so this if statement is to make sure
      // x is within the bounds of interpolation
      if ( x <= 0 or x >= 1.0) {
	std::cout << "ERROR IN EXCITATION CROSS: Eth = " << Eth 
	     << " Eij = " << EijEv << " x = " << x <<std::endl;
	exit(-1);
      }

      // An extra couple of sanity checks for the interpolation
      assert(x >= 0 and x <= 1);
      assert(splups[i].spline.size() == 5 or splups[i].spline.size() == 9);
      if(type > 4) break;
      
      // Extra sanity check to make sure x is monotonic to make sure
      // the arrays point to the right values
      for(int j = 1; j < 9; j++) {
	if (j < 5) assert(x_array5[j] > x_array5[j-1]);
	assert(x_array9[j] > x_array9[j-1]);
      }
      
      
      if (splups[i].spline.size() == 5) {
	Cspline<double, double> sp(x5, splups[i].spline);
	y = sp(x);
      }
      
      if (splups[i].spline.size() == 9) {
	Cspline<double, double> sp(x9, splups[i].spline);
	y = sp(x);
      }
      
      // Calculate the collision strength from the interpolated value
      double CStrength = 0.0;
      if (type == 1) {
	CStrength = y*log((Eth/EijEv) + M_E);
      }
      if(type == 2) {
	CStrength = y;
      }
      if(type == 3) {
	CStrength = y/(((Eth/EijEv) + 1)*((Eth/EijEv) + 1));
      }
      if(type == 4) {
	CStrength = y*log((Eth/EijEv) + C);
      }
      
      // from Dere et al. 1997 
      int weight = elvlc[splups[i].j-1].mult;
      double Eryd = E*eVtoRyd;
      cross += (M_PI*a0*a0*(CStrength/weight))/(Eryd);
    }

    if (splups[i].i == 1) {
      totalCross += cross;
      std::pair<double, double> cumi(totalCross, EijEv);
      CEcum.push_back(cumi);
    }
  }
  if (CEcum.size() == 0) { 
    std::cout << "\nERROR IN CE CROSS!" << "\n\tSplups size: " << splups.size() 
	      << "\n\tEth = " << Eth << "\n\tZ = " << Z << "\n\tC = " 
	      << C << "fblvl size: " << fblvl.size() <<std::endl; 
    exit(-1);
  }
  CEcrossCum = CEcum;
  return CEcum;
}

//! Calculate the QRP value as in Fontes, Sampson, Zhang 1999
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
  double q;
  q = (A*log(u) + D*(1.0-(1.0/u))*(1.0-(1.0/u)) + C*u*(1.0-(1.0/u))*(1.0-(1.0/u))*(1.0-(1.0/u))*(1.0-(1.0/u)) + ((c/u)+((d/u)*(d/u))*(1.0-(1.0/u))))/u;

  return q;
  
}

//! Calculate the direct ionization cross section from the spline,
//! which is a function of the interaction energy of the electron
double Ion::directIonCross(chdata ch, double E) 
{
  double u        = E/ip;
  unsigned char I = Z - C + 1; //test if its hydrogen-like/helium-like
  double ryd      = 27.2113845/2.0;
  double ipRyd    = ip/ryd;
  double a0       = 0.0529177211; // Bohr radius in nm
  double bohr_cs  = M_PI*a0*a0;

  double F, qr, cross;
  
  if (C == (Z+1)) {
    diCross = 0;
    return -1;
  }

  if (Z >= 20) {
    F = (140.0+pow((double(Z)/20.0),3.2))/141.;
  }
  else {
    F = 1.0;
  }

  qr = qrp(u)*F;
  if (Z >=6 or Z >= 10) std::cout << "QR = " << qr <<std::endl;

  // first two if statements are whether or not to use Fontes cross
  // sections
  if (I == 1 and Z >= 6) {
    cross = bohr_cs*(qr/ipRyd)*(qr/ipRyd);
  }
  else if (I==2 and Z >= 10) {
    cross = 2.0*bohr_cs*(qr/ipRyd)*(qr/ipRyd);
  }
  else {
    cross = 0;
    for (int i = 0; i < di_header.nfac; i++) {
      if (E >= diSpline[i].ev) {

	double u1  = E/diSpline[i].ev;
	double bte = 1.0 - log(diSpline[i].btf)/log(u1-1.0+diSpline[i].btf);

	Cspline<double, double> sp(diSpline[i].xspline, diSpline[i].yspline);
	double btcross = sp(bte);
	double a = 1.0 - diSpline[i].btf + exp(log(diSpline[i].btf)/(1.0 - bte));
	double cross_i = (log(a) + 1.0)*btcross/(a*diSpline[i].ev*diSpline[i].ev);
	cross += cross_i;
	// std::cout << "cross_i = " << cross_i << std::endl;
      }
    }
  }
  diCross = cross; // recast the cross section in nm^2
  return diCross;
}

double Ion::freeFreeCross(chdata ch, double E) 
{
  double hbc      = 197.327;	     // value of h-bar * c in eV nm
  double r0       = 2.81794033e-6;   // classic electron radius in nm
  double factor   = (Z*Z*r0*r0)/(137.);
  double hb       = 1.054572e-27;    // h-bar in erg s
  double me       = 9.10938e-28;
  double inmtoicm = 1e7;	  // nm^(-1) per cm^(-1)
  double eV2erg   = 1.602177e-12; // ergs per eV
  double c        = 2.998e10;	  // cm/s
  
  double p0 = sqrt(2*me*E*eV2erg);
  double v0 = p0/me;
  double b0 = v0/c;
  
  double momi = b0/sqrt(1.-b0*b0);
  
  double cum = 0;
  double dk  = 0;
  for (int j = 0; j < kffsteps; j++) {

    if (j != static_cast<int>(kgrid.size())-1) 
      dk = fabs(pow(10, kgrid[j+1]) - pow(10, kgrid[j]));

    double k  = pow(10, kgrid[j]);
    double Ek = k*hbc;
    double pk = k*inmtoicm*hb;
    double x  = (E - Ek)/E;
    double p  = p0-pk;

    if (x >= 0 and x <= 1) {
      double vf = p/me;
      double bf = vf/c;
      double corr = 1.0;
      corr = (b0*(1.0 - exp((-2*M_PI*double(Z))/(137.*b0))))/(bf*(1.0 - exp((-2*M_PI*double(Z))/(137.*bf))));
      double momf = bf/sqrt(1.-bf*bf);
      double x = momi + momf;
      double y = momi - momf;
      double dsig = corr*(factor*(dk/k))*(16./3.)*(1.0/(momi*momi))*log(x/y);
      
      cum = cum + dsig;
    }
    
  }
  double y_tmp = cum;
  return y_tmp;
  
}

/** Calculate the differential free-free cross section and return the
    cumulative cross section vector The formula used to calculate the
    cross section is 3BS(a) from Koch & Motz 1959 */
void Ion::freeFreeDifferential(chdata ch) 
{
  // Value of h-bar * c in eV nm
  double hbc = 197.327; 
  
  // Classic electron radius in nm
  double r0 = 2.81794033e-6;
  
  // Go through a grid of all k up to the energy E0
  double factor = (4.*Z*Z*r0*r0)/(137.);
  
  for (int i = 0; i < effsteps; i++) {
    std::vector<double> temp;
    std::vector<double> cum_Temp;
    double cum = 0;
    double E0 = egrid[i];
    double dk = 0;
    for(int j = 0; j < kffsteps; j++) {

      if (j != static_cast<int>(kgrid.size())-1) 
	dk = fabs(pow(10, kgrid[j+1]) - pow(10, kgrid[j]));

      double k = pow(10, kgrid[j]);
      double Ek = k*hbc;
      double x = (E0 - Ek)/E0;
      double dsig = (factor*(dk/k))*((1+x*x-(2.0/3.0)*x)*log(183.0/(double)Z) + (1.0/9.0)*x);

      cum = cum + dsig;
      temp.push_back(dsig);
      cum_Temp.push_back(cum);
    }
    ffDiffCross.push_back(temp);
    ffCumCross.push_back(cum_Temp);
  }
}


std::vector<double> Ion::radRecombCross(chdata ch, double E)
{
  // For testing . . .
  if (0) {
    std::vector<double> v1 = radRecombCrossMewe   (ch, E);
    std::vector<double> v2 = radRecombCrossSpitzer(ch, E);

    if (v1.back() < v2.back()) {
      std::cout << "   Mewe = " << v1.back() << std::endl;
      std::cout << "Spitzer = " << v2.back() << std::endl;
    }
    
    return v1;
  } else {
    return radRecombCrossMewe(ch, E);
  }
}


/** Calculates the differential radiative recombination cross section
    as a function of incoming electron impact energy, and returns the
    vector cumulative cross section array. */
std::vector<double> Ion::radRecombCrossMewe(chdata ch, double E) 
{
  // double hbc = 197.327; //value of hbar * c in eV nm
  // double hbckev = 0.197327; //value of hbc in keV nm
  double incmEv = 1.239842e-4; //1 inverse cm = 1.239.. eV

  // constant infront of the photo-cross using the Mewe method
  double D = 1.075812e-23;
  double mec2 = 510998.9; //mass of electron*c^2
  // double cumCross = 0;
  // int count = 0;
  std::vector<double> radRecCum;
  // double dk = 0;
  
  // double IP = ch.ipdata[Z-1][C-2];
  
  double cross = 0.;
  if (E!=0) {
    for(int j = 1; j < nfblvl; j++) 
      {
	double I; double eTemp;
	if(fblvl[j].encm == 0 and fblvl[j].encmth!=0) {
	  eTemp = fblvl[j].encmth;
	}
	else if(fblvl[j].encm != 0) {
	  eTemp = fblvl[j].encm;
	}
	else {
	  eTemp = 0;
	  std::cout << "ERROR WITH ETEMP!" <<std::endl;
	}
	double mult = double(fblvl[j].mult);
	double n = double(fblvl[j].lvl);
	eTemp = eTemp*incmEv; //convert the energy to eV
	eTemp = eTemp/1000.0; //convert to keV
	// I = (IP/1000.0 - eTemp);
	I = eTemp;
	if (I >= 0) {
	  double ePhot = E/1000. + I;
	  double hnu = E + I*1000.;
	  double Erat = (hnu*hnu)/(2.*mec2*E);
	  double crossi = (Erat*mult*D*I*I*(1.0/ePhot)*(1.0/ePhot)*(1.0/ePhot)*(1.0/n));
	  cross += crossi;
	  //if (C > 1 and cross == 0) {
	  // std::cout << "IP: " << IP << "\t" << eTemp << "\t" << I << "\t" << ePhot << "\t" <<Erat << "\t" << mult << "\t" << n << "\t" << crossi*1e18 <<std::endl;
	  //}
	  if (cross == 0) {
	    std::cout << "NULL IN RAD RECOMB: " << ip << "\t" << eTemp << "\t" << I << "\t" << ePhot << "\t" <<Erat << "\t" << mult << "\t" << n <<std::endl;
	  }
	  if (isnan(cross)) {
	    std::cout << cross << "\t" << I << "\t" << ePhot << "\t" << (double)n << "\t" << Erat <<std::endl;
	  }
	}
      }
  }
  radRecCum.push_back(cross*1.e18);
  radRecCrossCum = radRecCum;
  return radRecCum;
  
}

std::vector<double> Ion::radRecombCrossSpitzer(chdata ch, double E) 
{
  double incmEv = 1.239842e-4;	// 1 inverse cm = 1.239.. eV
  double eVincm = 8065.54446;	// 1 eV = 8065.54446 cm^{-1}
  double Ryd    = 13.6056923;	// Rydberg in eV
				// Ionization energy in cm^{-1}
  double ionE   = ch.ipdata[Z-1][0];

				// Cross-section prefactor in nm^2
  double coef   = 2.105310889751809e-08;
  
  std::vector<double> radRecCum;
  double cross = 0.0;
  if (E > 0) {
    for (int j = 0; j < nfblvl; j++) {
      double Ej = 0.0;
      if (j==0) 
	Ej = ionE;
      else if (j>0 && fblvl[j].encmth > 0) 
	Ej = ionE - fblvl[j].encmth * incmEv;
      else if (j>0 && fblvl[j].encm > 0) 
	Ej = ionE - fblvl[j].encm * incmEv;
      else continue;
      //
      double mult = static_cast<double>(fblvl[j].mult);
      double n    = static_cast<double>(fblvl[j].lvl );
      double Ephot  = E + Ej;
      double Erat   = Ej / Ephot;
      double crossn = coef * (Ej / Ephot) * (0.5*Ephot/E) * (mult/n);
      cross += crossn;

      if (cross == 0) {
	std::cout << "NULL IN RAD RECOMB: " << ip << "\t" << Ej << "\t" 
		  << Ephot << "\t" << Erat << "\t" << mult << "\t" << n <<std::endl;
      }
      if (isnan(cross)) {
	std::cout << cross << "\t" << Ej << "\t" << Ephot << "\t" 
		  << n << "\t" << Erat << std::endl;
      }
    }
  }
  radRecCum.push_back(cross);
  radRecCrossCum = radRecCum;
  return radRecCum;
}

// Ion print functions
void Ion::printInfo() {
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
	      << elvlc[i].conf        << "\t" 
	      << elvlc[i].designation << "\t"
	      << elvlc[i].spin        << "\t" 
	      << elvlc[i].l           << "\t" 
	      << elvlc[i].l_str       << "\t"
	      << elvlc[i].J           << "\t" 
	      << elvlc[i].mult        << "\t" 
	      << elvlc[i].encm        << "\t"
	      << elvlc[i].enry        << "\t" 
	      << elvlc[i].encmth      << "\t" 
	      << elvlc[i].enryth      << std::endl;
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


//! Read in the master list to store to be able to check if elements
//! are in it
void chdata::readMaster() 
{
  char * val;
  if ( (val = getenv("CHIANTI_DATA")) == 0x0) {
    if (myid==0)
      std::cout << "Could not find CHIANTI_DATA environment variable"
		<< " . . . exiting" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 51);
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
    if (myid==0) std::cout << "MASTER LIST FILE: "
			   << fileName << " NOT FOUND" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 46);
  }
  
}

/** Get the ipdata set so that if you want to get the ip of any Z, C,
    you call it as ipdata[Z-1][C-1-(int)die]
*/
void chdata::readIp() 
{
  char * val;
  if ( (val = getenv("CHIANTI_DATA")) == 0x0) {
    if (myid==0)
      std::cout << "Could not find CHIANTI_DATA environment variable"
		<< " . . . exiting" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 52);
  }
  std::string fileName(val);
  fileName.append("/ip/chianti.ip");
  int count = 0;
  unsigned char Z, C;
  double ip;
  double convert = 1.239841875e-4;
  std::string lineIn;
  ifstream ipFile(fileName.c_str());
  
  if (ipFile.is_open() ) {
    while (ipFile.good() and count < 365) {
      getline(ipFile, lineIn);
      
      std::vector <std::string> v;
      // split the string away from the white space to get the values
      istringstream iss(lineIn);
      copy(istream_iterator<std::string>(iss), istream_iterator<std::string>(), 
	   back_inserter<vector<std::string> >(v));
      // assign the values from the string
      Z = atoi(v[0].c_str());
      C = atoi(v[1].c_str());
      ip = atof(v[2].c_str())*convert;
      // set up the array
      ipdata[Z-1][C-1] = ip;
      count++;
    }
    ipFile.close();
  }
  else {
    if (myid==0) std::cout << "IP FILE: " 
			   << fileName << " NOT FOUND" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 47);
  }
}

/** read in the abundance file, in this situation, just for test using
    the cosmic.abund file. Can later put in a multidimensional array
    to allow for all the abundance files */
void chdata::readAbundanceAll() 
{
  char * val;
  if ( (val = getenv("CHIANTI_DATA")) == 0x0) {
    if (myid==0)
      std::cout << "Could not find CHIANTI_DATA environment variable"
		<< " . . . exiting" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 53);
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
    if (myid==0) std::cout << "Abundance file: " 
			   << fileName << " not found! " << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 48);
  }
  
}

// list names of all species to stdout
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
  for(int i = 0; i < 30; i++) {
    for(int j = 0; j < 30; j++) {
      if (ipdata[i][j] != 0) {
	std::cout << ipdata[i][j] << "\t";
      }
    }
    std::cout << std::endl;
  }
}

// chdata constructor
chdata::chdata() 
{
  //nVern = 465;
  //maxZ = 31; // maxZ = 30 + 1
  //maxNel = 31; 
  
  for(int i = 0; i < numEle; i++) abundanceAll[i] = 0;
  
  // Start with zeroed array
  for (int i=0; i<30; i++) {
    for (int j=0; j<30; j++) {
      ipdata[i][j] = 0.0;
    }
  }
  
  // std::cout << "Reading ip file\n";
  readIp();

  // std::cout << "Reading master file\n";
  readMaster();

  // std::cout << "Reading abundance file\n";
  readAbundanceAll();
  
  // Done
}

