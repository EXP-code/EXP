/*
  Test radiative cross section 
*/

#include <filesystem>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <random>
#include <tuple>
#include <regex>
#include <map>
#include <set>

#include <errno.h>
#include <sys/stat.h>

#include "Ion.H"
#include <cxxopts.H>
#include <libvars.H>

#include <mpi.h>

bool file_exists(const std::string& fileName)
{
  std::ifstream infile(fileName);
  return infile.good();
}


// Write ChiantiPy script, if necessary
//
void writeScript()
{
  const char *py =
#include "recomb_py.h"
    ;

  const std::string file("recomb.py");
    
  if (not file_exists(file)) {
    std::ofstream out(file);
    out << py;
    if (chmod(file.c_str(), S_IWUSR | S_IRUSR | S_IXUSR | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH)) {
      perror("Error in chmod:");
    }
  }
}

class Nigel
{
protected:

  //! Read and create T vs rate tables for the desired Z value
  void initialize(unsigned short Z, bool resonance);

  //! A simple table structure
  struct XY
  {
    std::vector<double> X;
    std::vector<double> Y;
  };

  std::map<unsigned short, XY> dataRR;
  std::map<unsigned short, std::vector<XY>> dataDR;

  const std::string datapath = "./";
  const std::string coupling = "[h]*ic";

public:

  //! Constructor
  Nigel(unsigned short Z, bool resonance=false)
  {
    initialize(Z, resonance);
  }

  //! Interpolate
  double operator()(unsigned short C, double T);

};


double Nigel::operator()(unsigned short C, double T)
{
  double ans = 0.0;

  auto u = dataRR.find(C);

  if (u != dataRR.end()) {
    auto & v = u->second;
    if (T>=v.X.front() and T<=v.X.back()) {
      auto lb = std::lower_bound(v.X.begin(), v.X.end(), T);
      auto ub = lb;
      if (lb == v.X.begin()) ub++;
      else                   lb--;

      size_t li = std::distance(v.X.begin(), lb);
      size_t ui = std::distance(v.X.begin(), ub);

      double A = (log(*ub) - log(T))/(log(*ub) - log(*lb));
      double B = (log(T) - log(*lb))/(log(*ub) - log(*lb));
      ans += exp(A*log(v.Y[li]) + B*log(v.Y[ui]));
    }
  }

  if (dataDR.find(C) != dataDR.end()) {

    for (auto & v : dataDR[C]) {

      if (T>=v.X.front() and T<=v.X.back()) {
	auto lb = std::lower_bound(v.X.begin(), v.X.end(), T);
	auto ub = lb;
	if (lb == v.X.begin()) ub++;
	else                   lb--;

	size_t li = std::distance(v.X.begin(), lb);
	size_t ui = std::distance(v.X.begin(), ub);
	
	double A = (log(*ub) - log(T))/(log(*ub) - log(*lb));
	double B = (log(T) - log(*lb))/(log(*ub) - log(*lb));
	ans += exp(A*log(v.Y[li]) + B*log(v.Y[ui]));
      }
    }
  }

  return ans;
}

void Nigel::initialize(unsigned short Z, bool resonance)
{
  // Periodic table
  //
  PeriodicTable table;

  /** From Nigel

      The naming convention is:
      user,year-code,sequence,element,charge,coupling scheme,core
      excitation (n->n').

      The year code is when the approach was introduced (NOT when
      calculated).  So, "20" is a new year code for the binned cross
      sections.  And "00" is the old year code for j-resolved rate
      coeffs.
  */

  // Set default datapath
  //
  std::filesystem::path dpath{datapath};

  // Look for atomic data path in environment
  //
  if (const char* env_p = std::getenv("ATOMIC_DATA_PATH")) {
    dpath = std::filesystem::path( env_p );
  }
  
  // Filter strings
  //
  const std::string user("nrb"), year("05"); // RR rates


  // File regular expression for RR
  //
  std::ostringstream rstr;

  // Path--+               Sequence-----+     +---Separator
  //       |                            |     |
  //       |               Separator--+ |     |
  //       |                          | |     |
  //       v                          v v     v
  rstr << ".*[/]" << user << year << "#[a-z]+[_]"
       << "([a-z]+)([0-9]+)" << coupling << "[0-9]*[.]dat";
  //        ^       ^           ^            ^        ^
  //        |       |           |            |        |
  //        |       |           |            +--------|----- Core excitation
  //        |       |           |                     |      n->n' for
  // Ion name       Charge      Coupling     Suffix---+      dielectronic
  // (Group 1)      (Group 2)

  std::regex species(rstr.str());
  std::smatch matches;

  // Will contain matched file(s).  Should only be one.
  //
  std::map<lQ, std::string> names;

  try {
    for (auto const& dir_entry : std::filesystem::recursive_directory_iterator{dpath / "RR"}) {
      std::string name = dir_entry.path().string();
      if (std::regex_search(name, matches, species)) {
	auto it = table[matches[1].str()];
	if (it == 0) {
	  std::ostringstream sout;
	  sout << "Badnell::initialize: [0] could not locate element <"
	       << matches[1].str() << ">";
	  throw std::runtime_error(sout.str());
	}

	// The element and and charge
	//
	lQ ZC = {it->number(), stoi(matches[2].str()) + 1};

	// Add the file name
	//
	names[ZC] = name;
      }
    }
  }
  catch(std::filesystem::filesystem_error const& ex) {
    std::cout
      << "Directory scan error for RR..."
      << "what():  " << ex.what() << '\n'
      << "path1(): " << ex.path1() << '\n'
      << "path2(): " << ex.path2() << '\n'
      << "code().value():    " << ex.code().value() << '\n'
      << "code().message():  " << ex.code().message() << '\n'
      << "code().category(): " << ex.code().category().name() << '\n';
    
    exit(32);
  }
    
  // File regular expressions for DR
  //
  std::ostringstream rdr1, rdr2;
  std::string yr1("00"), yr2("12");

  // Path--+              Sequence-----+     +---Separator
  //       |                           |     |
  //       |              Separator--+ |     |
  //       |                         | |     |
  //       v                         v v     v
  rdr1 << ".*[/]" << user << yr1 << "#[a-z]+[_]"
       << "([a-z]+)([0-9]+)" << coupling << "([0-9]*)";
  rdr2 << ".*[/]" << user << yr2 << "#[a-z]+[_]"
       << "([a-z]+)([0-9]+)" << coupling << "([0-9]*)";
  //        ^       ^           ^             ^
  //        |       |           |             |
  //        |       |           |             +--- Core excitation
  //        |       |           |                  n->n' for
  // Ion name       Charge      Coupling           dielectronic
  // (Group 1)      (Group 2)                      (Group 3)

  // Add suffix
  //
  if (resonance) {
    rdr1 << ".dat";
    rdr2 << ".dat";
  } else {
    rdr1 << "_t3.dat";
    rdr2 << "_t3.dat";
  }

  std::regex spc1(rdr1.str());
  std::regex spc2(rdr2.str());

  std::smatch match1, match2;

  std::map<lQ, std::vector<std::pair<std::string, std::string>>> names2;

  try {
    for (auto const& dir_entry : std::filesystem::recursive_directory_iterator{dpath / "DR"}) {
      std::string name = dir_entry.path().string();
      if (std::regex_search(name, match1, spc1)) {
	auto it = table[match1[1].str()];
	if (it == 0) {
	  std::ostringstream sout;
	  sout << "Badnell::initialize: [1] could not locate element <"
	       << match1[1].str() << ">";
	  throw std::runtime_error(sout.str());
	}

	// The element and and charge
	//
	lQ ZC = {it->number(), stoi(match1[2].str()) + 1};

	// Add this excitaton to the db
	//
	names2[ZC].push_back({match1[3].str(), name});
      }
      if (std::regex_search(name, match2, spc2)) {
	auto it = table[match2[1].str()];
	if (it == 0) {
	  std::ostringstream sout;
	  sout << "Badnell::initialize: [2] could not locate element <"
	       << match2[1].str() << ">";
	  throw std::runtime_error(sout.str());
	}

	// The element and and charge
	//
	lQ ZC = {it->number(), stoi(match2[2].str()) + 1};

	// Add this excitaton to the db
	//
	names2[ZC].push_back({match1[3].str(), name});
      }
    }
  }
  catch(std::filesystem::filesystem_error const& ex) {
    std::cout
      << "Directory scan error for DR..."
      << "what():  " << ex.what() << '\n'
      << "path1(): " << ex.path1() << '\n'
      << "path2(): " << ex.path2() << '\n'
      << "code().value():    " << ex.code().value() << '\n'
      << "code().message():  " << ex.code().message() << '\n'
      << "code().category(): " << ex.code().category().name() << '\n';
    
    exit(34);
  }

  // Loop through each desired ion and process the files
  //
  for (unsigned short C=2; C<Z+2; C++) {

    lQ ZC{Z, C};

    std::string line;
    double T, R;
    
    // Begin with RR files
    //
    auto it = names.find(ZC);

    if (it != names.end()) {

      std::ifstream in{dpath / it->second};
    
      std::getline(in, line);

      bool looking = true;
      while (in) {
	if (looking) {
	  // Keep reading lines until we reach the totals
	  //
	  if (line.find("T(K)") != std::string::npos) {
	    looking = false;
	    std::getline(in, line); // Read header
	  }
	} else {
	  // We are done!
	  //
	  if (line[0] == 'C') break;
	  
	  // Read the line
	  //
	  std::istringstream ins(line);
	  
	  // Values (Energy in Rydbergs and Energy*Sigma in Mbarn*Rydberg)
	  //
	  ins >> T >> R;
	  dataRR[C].X.push_back(T);
	  dataRR[C].Y.push_back(R);
	}
	
	// Read next line
	std::getline(in, line);
      }
      // END: RR file read
    }

    // Now, on to DR files
    //
    auto it2 = names2.find(ZC);

    if (it2 != names2.end()) {
      
      for (auto v : it2->second) {

	XY xy;

	auto cr = v.first;
	// TEST
	if (cr.size()>2) continue;

	auto in = std::ifstream{dpath / "DR" / v.second};
	
	std::getline(in, line);
	
	bool looking = true;
	while (in) {
	  if (looking) {
	    // Keep reading lines until we reach the totals
	    //
	    if (line.find("T(K)") != std::string::npos) {
	      looking = false;
	      std::getline(in, line); // Read header
	    }
	  } else {
	    // We are done!
	    //
	    if (line[0] == 'C') break;
	    
	    // Read the line
	    //
	    std::istringstream ins(line);
	    
	    // Default values
	    //
	    T = R = 0.0;
	    
	    // Values (Energy in Rydberg and energy averaged cross
	    // section in Mbarn)
	    ins >> T >> R;
	    xy.X.push_back(T); // Energy is now in eV
	    xy.Y.push_back(R);
	  }
	
	  // Read next line
	  std::getline(in, line);
	}
	// END: DR file read
	
	dataDR[C].push_back(xy);
      }
      // END: DR stanza

      //       +--- Verbose reporting on DR files found
      //       |    (not for production mode)
      //       v
    } else if (false) {
      std::ostringstream sout;
      sout << "Badnell::initialize: could not locate DR file for element "
	   << table[ZC.first]->name() << ", ion (Z, C) = ("
	   << ZC.first << ", " << ZC.second << ")";
      std::cout << sout.str() << std::endl;
    }
    // END: file loop

  }

  // DONE
}


int main (int ac, char **av)
{
  writeScript();

  //===================
  // MPI preliminaries 
  //===================

  MPI_Init(&ac, &av);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  std::string cmd_line;
  for (int i=0; i<ac; i++) {
    cmd_line += av[i];
    cmd_line += " ";
  }

  unsigned short Z;
  double Tmin, Tmax, Emin, Emax;
  std::string RRtype;
  int numT, norder, numE;
  bool use_grid;

  cxxopts::Options options(av[0], "Test radiative cross section ratios");

  options.add_options()
    ("h,help", "produce help message")
    ("l,long", "long output: print each rate and ratio")
    ("resonance", "Use Nigel's resonance averaged Maxwellian files")
    ("Z,elem", "atomic number",
     cxxopts::value<unsigned short>(Z)->default_value("2"))
    ("t,Tmin", "minimum temperature",
     cxxopts::value<double>(Tmin)->default_value("1000.0"))
    ("T,Tmax", "maximum temperature",
     cxxopts::value<double>(Tmax)->default_value("10000000.0"))
    ("e,Emin", "minimum energy in eV",
     cxxopts::value<double>(Emin)->default_value("0.001"))
    ("E,Emax", "maximum energy in eV",
     cxxopts::value<double>(Emax)->default_value("1000.0"))
    ("N,NumT", "number of temperature points",
     cxxopts::value<int>(numT)->default_value("200"))
    ("M,NumE", "number of incremental energy points",
     cxxopts::value<int>(numE)->default_value("100"))
    ("n,norder", "Laguerre order",
     cxxopts::value<int>(norder)->default_value("20"))
    ("R,RRtype", "cross-section type",
     cxxopts::value<std::string>(RRtype)->default_value("Badnell"))
    ("g,grid", "use radiative recombination grid",
     cxxopts::value<bool>(use_grid)->default_value("true"))
    ("reweight", "use energy bin reweighting for dielectronic cross sections")
    ;

  cxxopts::ParseResult vm;

  try {
    vm = options.parse(ac, av);
  } catch (cxxopts::OptionException& e) {
    if (myid==0) std::cout << "Option error: " << e.what() << std::endl;
    MPI_Finalize();
    exit(-1);
  }

  if (vm.count("help")) {
    std::cout << options.help() << std::endl;
    std::cout << "Example: Helium ionization fractions at T=40,000 K" << std::endl;
    std::cout << "\t" << av[0]
	      << " -Z 2 -T 40000" << std::endl;
    MPI_Finalize();
    return 1;
  }

  if (vm.count("reweight")) BadnellData::reweight = true;

  std::string prefix("testNigel");
  std::string cmdFile = prefix + ".cmd_line";
  std::ofstream cmd(cmdFile.c_str());
  if (!cmd) {
    std::cerr << "testCrossSection: error opening <" << cmdFile
	      << "> for writing" << std::endl;
  } else {
    cmd << cmd_line << std::endl;
  }

  cmd.close();

  // Get recombationation rates from ChiantiPy
  //
  std::string inFile("testRatio.in");
  std::ostringstream sout;
  sout << "$PYTHON ./recomb.py"
       << " -Z " << Z
       << " -t " << Tmin
       << " -T " << Tmax
       << " -n " << numT
       << " -o " << inFile << '\0';

  std::cout << "Command: " << &sout.str()[0] << std::endl;
  int sret = system(&sout.str()[0]);

  std::ifstream in(inFile);

  // For temperature data
  //
  std::vector<double> temp(numT);

  // For recombination rate data
  //
  std::vector<std::vector<double>> cdata(Z);
  for (int C=2; C<=Z+1; C++) cdata[C-2].resize(numT);

  // Get data from generated file
  //
  for (int nt=0; nt<numT; nt++) {
    std::string line;
    if (std::getline(in, line)) {
      std::istringstream ins(line);
      ins >> temp[nt];
      for (int C=2; C<=Z+1; C++) ins >> cdata[C-2][nt];
    }
  }

  // Open output file
  //
  std::ofstream rcb(prefix + ".data", ios::out);
  if (rcb.bad()) {
    std::cerr << "testRatio: error opening <" << prefix + ".data"
	      << "> for writing" << std::endl;
    exit(-1);
  }

  // Initialize CHIANTI
  //
  std::set<unsigned short> ZList = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 15, 16, 19, 20};

  if (ZList.find(Z) == ZList.end()) {
    if (myid==0) {
      std::cout << "Z=" << Z 
		<< " is not in element list.  Current list contains:";
      for (auto z : ZList) std::cout << " " << z;
      std::cout << std::endl;
    }
    MPI_Finalize();
    exit(1);
  }

  PeriodicTable pt;

  bool resonance = false;
  if (vm.count("resonance")) resonance = true;

  Nigel nigel(Z, resonance);

  atomicData ad;

  ad.createIonList(ZList);

  if (myid) {
    MPI_Finalize();
    return 0;
  }

  Ion::setRRtype(RRtype);

  bool logE = false;
  if (Emin>0.0) {
    Emin = log(Emin);
    Emax = log(Emax);
    logE = true;
  }

  typedef std::map<unsigned short, std::vector<double> > rateMap;

  numE = std::max<int>(numE, 10);
  
  Ion::useRadRecombGrid = use_grid;

  double dE = (Emax - Emin)/numE;

  for (int nt=0; nt<numT; nt++) {

    double T = temp[nt];

    rateMap val0;
    std::vector<rateMap> val1(numE);
    
    for (int ne=0; ne<numE; ne++) {
      
      double Eb = Emin + dE*ne;
      double Ef = Emin + dE*(ne+1);
      if (logE) {
	Eb = exp(Eb);
	Ef = exp(Ef);
      }
      
      if (val0.size()) {
	std::map<unsigned short, std::vector<double> >
	  valT = ad.recombEquil(Z, T, Eb, Ef, norder);
	for (auto v : val0) {
	  unsigned short C = v.first;
	  size_t        sz = v.second.size();
	  for (size_t j=0; j<sz; j++) val0[C][j] += valT[C][j];
	}
      } else {
	val0 = ad.recombEquil(Z, T, Ef, Eb, norder);
      }
      val1[ne] = val0;
    }
      
    rcb << std::setw(16) << T;
    for (auto v : val1.back()) {
      unsigned short C = v.first;
      if (C>1) {
	if (cdata[C-2][nt]>0.0) {
	  rcb << std::setw(16) << v.second[2]*1.0e-14
	      << std::setw(16) << cdata[C-2][nt]
	      << std::setw(16) << nigel(C, T);
	} 
      }
    }
    rcb << std::endl;
  }

  MPI_Finalize();

  return 0;
}
