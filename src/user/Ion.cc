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

      e.level       = atoi(v[0].c_str());
      e.conf        = atoi(v[1].c_str());
      e.designation = v[2];
      e.spin        = atoi(v[3].c_str());
      e.l           = atoi(v[4].c_str());
      e.l_str       = v[5];
      e.J           = atof(v[6].c_str());
      e.mult        = atoi(v[7].c_str());
      e.encm        = atof(v[8].c_str());
      e.enry        = atof(v[9].c_str());
      e.encmth      = atof(v[10].c_str());
      e.enryth      = atof(v[11].c_str());
      
      elvlc[e.level] = e;
    }
    elvlcFile.close();
  }
  else {
    if (myid==0) std::cout << "Cannot find file: " << fileName << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 42);
  }
}

void Ion::readwgfa() 
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
  fileName.append(".wgfa");
  
  std::string inLine;
  wgfa_data w;
  ifstream wgfaFile(fileName.c_str());

  if (wgfaFile.is_open()) {

    while (wgfaFile.good()) {

      std::vector <std::string> v;
      getline(wgfaFile, inLine);
      std::istringstream iss(inLine);
      copy(istream_iterator<std::string>(iss), istream_iterator<std::string>(), 
	   back_inserter<vector<std::string> >(v));

      if (atoi(v[0].c_str()) == -1) break;

      w.lvl1    = atoi(v[0].c_str());
      w.lvl2    = atoi(v[1].c_str());
      w.wvl     = atof(v[2].c_str());
      w.gf      = atof(v[3].c_str());
      w.avalue  = atof(v[4].c_str());
      w.pretty1 = v[5].c_str();
      w.pretty2 = v[6].c_str();
      w.ref     = v[7].c_str();
      
      wgfa[lQ(w.lvl1, w.lvl2)] = w;
    }
    wgfaFile.close();
  }
  else {
    if (myid==0) std::cout << "Cannot find file: " << fileName << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 42);
  }
}

//! Read in the fblvl file found in the CHIANTI database
void Ion::readfblvl() 
{
  std::string MasterNameT = ZCtoName(Z, C);

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

      if (atoi(v[0].c_str()) == -1) break;

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
    if (myid==0) std::cout << "Cannot find file: " << fileName << std::endl;
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

      s.Z       = atoi(v[0].c_str());
      s.C       = atoi(v[1].c_str());
      s.i       = atoi(v[2].c_str());
      s.j       = atoi(v[3].c_str());
      s.type    = atoi(v[4].c_str());
      s.gf      = atof(v[5].c_str());
      s.delERyd = atof(v[6].c_str());
      s.Const   = atof(v[7].c_str());

				// Spline coefficients
      for(unsigned i = 8; i < v.size(); i++) {
	s.spline.push_back(atof(v[i].c_str()));
      }
      
      splups.push_back(s);
      s.spline.erase(s.spline.begin(), s.spline.end());		
    }
    sFile.close();
  }
  else {
    if (myid==0) std::cout << "Cannot find file: " 
			   << fileName << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 44);
  }
}

/**
   Read in the direct ionization cross section splines from the file
*/
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
  }
  else {
    if (myid==0) std::cout << "Cannot find file: " 
			   << fileName << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 45);
  }
}

//! Initialization function when the master name is given
Ion::Ion(std::string name, chdata* ch) : ch(ch)
{
  MasterName = name;

  std::vector<std::string> v;
  boost::split(v, MasterName, boost::is_any_of("_") );
  eleName = v[0];

  convertName();		// Sets Z and C . . . 
  ip = ch->ipdata[lQ(Z, C)];

  if (isInMasterList(MasterName)) {
    readfblvl();
    readelvlc();
    readwgfa();
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
Ion::Ion(unsigned short Z, unsigned short C, chdata* ch) : ch(ch), Z(Z), C(C)
{
  d = false;
  MasterName = ZCtoName(Z, C);

  std::vector<std::string> v;
  boost::split(v, MasterName, boost::is_any_of("_") );
  eleName = v[0];

  ip = 0.0;

  if (Z>=C) {

    ip = ch->ipdata[lQ(Z, C)];

    if (isInMasterList(MasterName)) {
      readfblvl();
      readSplups();
      readDi();
      readelvlc();
      readwgfa();
    }
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

  MasterName = I.MasterName;
  eleName    = I.eleName;
  elvlc      = I.elvlc;
  wgfa       = I.wgfa;
}

/** 
    Calculate the collision excitation cross sections and return the
    cumulative cross section array.

    Returns a vector with a pair of (cross section, energy difference)
    since the file input, and thus array, are not in any specific order
*/
Ion::collType
Ion::collExciteCross(double E)
{
  const double x_array5[5] = {0, 0.25, 0.5, 0.75, 1.0};
  const double x_array9[9] = {0, 0.125, 0.25 , 0.375, 0.5 , 
			      0.625 , 0.75, 0.875, 1.0};

  std::vector<double> x5(x_array5, x_array5+5);
  std::vector<double> x9(x_array9, x_array9+9);
  
  const double eVtoRyd = 1.0/13.60569253;
  const double RydtoeV = 13.60569253;

  const double a0 = 0.0529177211; // Bohr radius in nm
  
  collType CEcum;
  const std::pair<double,double> Null(0, 0);
  if (splups.size() == 0) {
    CEcum.push_back(Null);
    CEcrossCum = CEcum;
    return CEcum;
  }

  double totalCross = 0.0;

  for (size_t i=0; i<splups.size(); i++) {

    double EijEv = splups[i].delERyd*RydtoeV;
    double Const = splups[i].Const;

    if (splups[i].i==1 and E >= EijEv) {

      assert(splups[i].i == 1);
      assert(splups[i].spline.size() != 0);
      
      // Following Burgess & Tully (BT), 1992, Section 3
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
      if ( x <= 0 or x >= 1.0) {
	std::cout << "ERROR IN EXCITATION CROSS: Ej = " << Ej
	     << " Eij = " << EijEv << " x = " << x <<std::endl;
	exit(-1);
      }

      // An extra couple of sanity checks for the interpolation
      //
      assert(x >= 0 and x <= 1);
      assert(splups[i].spline.size() == 5 or splups[i].spline.size() == 9);
      if (type > 4) break;
      
      // Extra sanity check to make sure x is monotonic to make sure
      // the arrays point to the right values
      //
      for (int j = 1; j < 9; j++) {
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
				// BT, eq. 6
      if (type == 1) {
	CStrength = y * log((Ej/EijEv) + M_E);
      }
				// BT, eq. 10
      if(type == 2) {
	CStrength = y;
      }
				// BT, eq. 12
      if(type == 3) {
	double fac = Ej/EijEv + 1.0;
	CStrength = y/(fac*fac);
      }
				// BT, eq. 14
      if(type == 4) {
	CStrength = y * log((Ej/EijEv) + C);
      }
      
      // From Dere et al. 1997 
      int weight = elvlc[splups[i].j-1].mult;
      totalCross += (M_PI*a0*a0*(CStrength/weight))/(E*eVtoRyd);
      std::pair<double, double> cumi(totalCross, EijEv);
      CEcum.push_back(cumi);
    }

  }

  if (CEcum.size() == 0) CEcum.push_back(Null);
  
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
double Ion::directIonCross(double E) 
{
  double u        = E/ip;
				// Test for hydrogen-like/helium-like ion
  unsigned char I = Z - C + 1;
  double ryd      = 27.2113845/2.0;
  double ipRyd    = ip/ryd;
  double a0       = 0.0529177211; // Bohr radius in nm
  double bohr_cs  = M_PI*a0*a0;

  double F, qr, cross;
  
  if (C == (Z+1)) {
    diCross = 0;
    return -1;
  }

  // From Fontes, et al. Phys Rev A, 59, 1329 (eq. 2.11)
  if (Z >= 20) {
    F = (140.0+pow((double(Z)/20.0),3.2))/141.;
  }
  else {
    F = 1.0;
  }

  qr = qrp(u)*F;

  // first two if statements are whether or not to use Fontes cross
  // sections
  if (I == 1 && Z >= 6) {
    cross = bohr_cs*qr/(ipRyd*ipRyd);
  }
  else if (I==2 && Z >= 10) {
    cross = 2.0*bohr_cs*qr/(ipRyd*ipRyd);
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
      }
    }
  }
  diCross = cross; // recast the cross section in nm^2
  return diCross;
}

/** 
    Cross section is 3BN(a) from Koch & Motz 1959, with the low-energy
    Elwert factor (see Koch & Motz eq. II-6)
*/
double Ion::freeFreeCross(double E) 
{
  double hbc      = 197.327;	     // value of h-bar * c in eV nm
  double r0       = 2.81794033e-6;   // classic electron radius in nm
  double factor   = (Z*Z*r0*r0)/(137.);
  double hb       = 1.054572e-27;    // h-bar in erg s
  double me       = 9.10938e-28;     // electron mass in g
  double inmtoicm = 1e7;	     // nm^(-1) per cm^(-1)
  double eV2erg   = 1.602177e-12;    // ergs per eV
  double c        = 2.998e10;	     // cm/s
  
  double p0       = sqrt(2*me*E*eV2erg);
  double v0       = p0/me;
  double b0       = v0/c;
  
  double momi     = b0/sqrt(1.0 - b0*b0);
  
  double cum      = 0;
  double dk       = 0;

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
void Ion::freeFreeDifferential() 
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


std::vector<double> Ion::radRecombCross(double E)
{
  // For testing . . .
  if (0) {
    std::vector<double> v1 = radRecombCrossMewe   (E);
    std::vector<double> v2 = radRecombCrossTopBase(E);
    std::vector<double> v3 = radRecombCrossKramers(E);
    std::vector<double> v4 = radRecombCrossSpitzer(E);

    std::cout << "    Mewe = " << std::setw(16) << v1.back() << std::endl;
    std::cout << " TopBase = " << std::setw(16) << v2.back() << std::endl;
    std::cout << " Kramers = " << std::setw(16) << v3.back() << std::endl;
    std::cout << " Spitzer = " << std::setw(16) << v4.back() << std::endl;
    
    return v1;
  } else {
    return radRecombCrossTopBase(E);
    // return radRecombCrossMewe   (E);
    // return radRecombCrossKramers(E);
    // return radRecombCrossSpitzer(E);
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
*/
std::vector<double> Ion::radRecombCrossKramers(double E) 
{
  const double incmEv = light * planck / eV;

  // Bohr radius in nm
  //
  const double a0 = 0.0529177211;

  // Return vector
  //
  std::vector<double> radRecCum;
  double cross = 0.0;

  // This is the target neutral
  //
  Ion* N = &ch->IonList[lQ(Z, C-1)];

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
  
  for (fblvlType::iterator j=N->fblvl.begin(); j!= N->fblvl.end(); j++) {

    fblvl_data* f = &j->second;

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
    double sigmaP = 0.25*f->lvl*aeff*aeff*Erat*Erat*Erat;

    // Ion statistical weight
    //
    double mult0 = (C<=Z ? fblvl[1].mult : 1);
      
    // Recombined statistical weight
    //
    double mult1 = f->mult;
      
    // Milne relation
    //
    double sigmaR = 0.5*mult1/mult0 * Enu*Enu / (E*mec2*1.0e6) * sigmaP;
	  
    cross += sigmaR;

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
      
    if (isnan(cross)) {
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
  radRecCrossCum = radRecCum;
  return radRecCum;
}


/** Calculates the differential radiative recombination cross section
    as a function of incoming electron impact energy, and returns the
    vector cumulative cross section array. 

    Details of implementation from CHIANTI. See "The Free-Bound
    Continuum", P.R. Young, Ver. 1.1, 8-Sep-2009.

    Uses Milne relation.
*/
std::vector<double> Ion::radRecombCrossMewe(double E) 
{
  double incmEv = 1.239842e-4; //1 inverse cm = 1.239.. eV

  // constant infront of the photo-cross using the Mewe method
  double D = 1.075812e-23;

  // Electron rest mass in keV
  double mec2 = 510.998896;

  // Key of parent ion
  lQ Q(Z, C-1);

  // Get pointers to data
  double IP = ch->ipdata[Q];
  Ion* N    = &ch->IonList[Q];

  // Convert kinetic energy to keV
  E *= 1.0e-3;

  // Return values
  std::vector<double> radRecCum;
  double cross = 0.0;
  
  if (E!=0) {
    
    double mult0 = (C<=Z ? fblvl[1].mult : 1);

    for (fblvlType::iterator j=N->fblvl.begin(); j!=N->fblvl.end(); j++) {

      fblvl_data* f = &j->second;
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

      if (I >= 0) {
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
	if (isnan(cross)) {
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

  // Convert to nm ---------+
  //                        |
  //                        v
  radRecCum.push_back(cross*1.0e18);

  radRecCrossCum = radRecCum;

  return radRecCum;
}

std::vector<double> Ion::radRecombCrossSpitzer(double E) 
{
				// 1 inverse cm = 1.239.. eV
  const double incmEv = 1.239842e-4;

				// Cross-section prefactor in nm^2
  const double coef   = 2.105310889751809e-08;

				// Ionization energy in eV
  double ionE         = ch->ipdata[lQ(Z, 1)];

  std::vector<double> radRecCum;
  double cross = 0.0;
  if (E > 0) {
    for (fblvlType::iterator j=fblvl.begin(); j!=fblvl.end(); j++) {
      fblvl_data* f = &j->second;

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
      double crossn = coef * (Ej / Ephot) * (0.5*Ephot/E) * (mult/n);

      cross += crossn;

      if (cross == 0) {
	std::cout << "NULL IN RAD RECOMB: " << ip << "\t" << Ej << "\t" 
		  << Ephot << "\t" << Erat << "\t" << mult << "\t" 
		  << n <<std::endl;
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

std::vector<double> Ion::radRecombCrossTopBase(double E) 
{
  // Initialize TopBase data (once) if needed
  //
  if (ch->tb.get() == 0) 
    ch->tb = boost::shared_ptr<TopBase>(new TopBase);

  // Call for the cross section
  //
  TopBase::iKey k(Z, C);
  std::vector<double> ret(ch->tb->sigmaFB(k, E), 1);
  return ret;
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
    you call it as ipdata[lQ(Z, C-(int)die)]
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
  std::cout << std::string(60, '-') << std::endl
	    << std::setw( 3) << "Z" << std::setw( 3) << "C"
	    << std::setw(16) << "Energy (eV)" << std::endl;

  for (std::map<lQ, double>::iterator i=ipdata.begin(); i!=ipdata.end(); i++) {
    if (i->second != 0) {
      std::cout 
	<< std::setw( 3) << i->first.first
	<< std::setw( 3) << i->first.second
	<< std::setw(16) << i->second
	<< std::endl;
    }
  }

  std::cout << std::string(60, '-') << std::endl;
}

// chdata constructor
chdata::chdata() 
{
  // nVern = 465;
  // maxZ = 31; // maxZ = 30 + 1
  // maxNel = 31; 
  
  for(int i = 0; i < numEle; i++) abundanceAll[i] = 0;
  
  // std::cout << "Reading ip file\n";
  readIp();

  // std::cout << "Reading master file\n";
  readMaster();

  // std::cout << "Reading abundance file\n";
  readAbundanceAll();
  
  // Done
}

void chdata::createIonList(const std::set<unsigned short>& ZList)
{
  // Fill the Chianti data base
  //
  for (ZLtype::const_iterator i=ZList.begin(); i!=ZList.end(); i++) {
    for (int j=1; j<*i+2; j++) {
      lQ Q(*i, j);
      IonList[Q] = Ion(*i, j, this);
      IonList[Q].freeFreeDifferential();
    }
    Ni[*i] = 1.0;		// Not sure what this does . . . 
  }
}
