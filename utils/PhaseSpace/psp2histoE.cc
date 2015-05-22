/*
  Separate a psp structure and make a 1-d histogram

  MDWeinberg 08/26/11
*/

using namespace std;

#include <cstdlib>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <list>
#include <map>

#include <Species.H>

#include <StringTok.H>
#include <header.H>
#include <PSP.H>

#include <boost/program_options.hpp>

namespace po = boost::program_options;


				// Globals for exputil library
				// Unused here
int myid = 0;
char threading_on = 0;
pthread_mutex_t mem_lock;
string outdir, runtag;

int
main(int ac, char **av)
{
  char *prog = av[0];
  double time, Emax, Lunit, Tunit;
  bool verbose = false;
  std:: string cname;
  int numb, comp, sindx, eindx;

  // Parse command line

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",		"produce help message")
    ("verbose,v",       "verbose output")
    ("time,t",		po::value<double>(&time)->default_value(1.0e20),
     "find closest time slice to requested value")
    ("Lunit,L",		po::value<double>(&Lunit)->default_value(3.086e18),
     "System length in physical units (cgs)")
    ("Tunit,T",		po::value<double>(&Tunit)->default_value(3.15569e12),
     "System time in physical units (cgs)")
    ("Emax,E",		po::value<double>(&Emax)->default_value(200.0),
     "Maximum energy in eV")
    ("bins,b",	        po::value<int>(&numb)->default_value(40),
     "number of bins")
    ("species,s",	po::value<int>(&sindx)->default_value(0),
     "position of species index")
    ("electrons,e",	po::value<int>(&eindx)->default_value(6),
     "position of electron index")
    ("name,c",	        po::value<std::string>(&cname)->default_value("gas"),
     "component name")
    ("files,f",         po::value< std::vector<std::string> >(), 
     "input files")
    ;


  po::variables_map vm;

  try {
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);    
  } catch (po::error& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    std::cout << "Example: " << std::endl;
    std::cout << "\t" << av[0]
	      << " --temp=25000 --number=250000 --output=out.bod" << std::endl;
    return 1;
  }

  if (vm.count("verbose")) {
    verbose = true;
  }

  // Units
  //
  const double m_e = 9.10938215e-28; // electron mass in g
  const double eV  = 1.60217653e-12; // erg per eV
  double Vunit     = Lunit/Tunit;
  double factor    = m_e/eV * Vunit*Vunit;

  // Get file arguments
  //
  std::vector<std::string> files = vm["files"].as< std::vector<std::string> >();

  bool first = true;
  
  for (auto file : files ) {

    ifstream *in = new ifstream(file.c_str());
    if (!*in) {
      cerr << "Error opening file <" << file << "> for input\n";
      exit(-1);
    }

    if (verbose) cerr << "Using filename: " << file << endl;


				// Parse the PSP file
				// ------------------
    PSPDump psp(in);

    in->close();

				// Now write a summary
				// -------------------
    if (verbose) {
      
      psp.PrintSummary(in, cerr);
    
      cerr << "\nBest fit dump to <" << time << "> has time <" 
	   << psp.SetTime(time) << ">\n";
    } else 
      psp.SetTime(time);

				// Dump ascii for each component
				// -----------------------------
    delete in;
    in = new ifstream(file);
    
  
    vector<double> pos(3), vel(3);

    PSPstanza *stanza;
    SParticle* part;


				// Will contain array for each species
				//
    std::map< speciesKey, std::vector<float> > shist;
    double dkE   = Emax/numb;
    double total = 0.0;

    for (stanza=psp.GetStanza(); stanza!=0; stanza=psp.NextStanza()) {
    
      if (stanza->name != cname) continue;


				// Position to beginning of particles
				// 
      in->seekg(stanza->pspos);

      for (part=psp.GetParticle(in); part!=0; part=psp.NextParticle(in)) {

	double kE = 0.0;
	for (size_t i=0; i<3; i++) {
	  double v = part->datr(eindx+i);
	  kE += v*v;
	}
	kE *= factor;

	KeyConvert kc(part->iatr(sindx));
	speciesKey k = kc.getKey();

	if (shist.find(k) == shist.end()) shist[k].resize(numb+1, 0.0);
      
	size_t l = std::min<size_t>(kE/dkE, numb);

	double wgt = part->mass() * (k.second - 1);
	shist[k][l] += wgt;
	total += wgt;
      }
    }
    
    //
    // Output
    //
    const size_t fw = 14;
    const size_t sw =  9;
    double Time = psp.CurrentTime();

    if (first) {
      std::cout << setw(fw) << "Time"
		<< setw(fw) << "Energy";
      
      for (auto v : shist) {
	speciesKey k = v.first;
	ostringstream str;
	str << "(" << k.first << ", " << k.second << ")";
	cout << setw(fw) << str.str();
      }
      cout << std::endl;

      std::cout << setw(fw) << std::string(sw, '-')
		<< setw(fw) << std::string(sw, '-');

      for (auto v : shist) {
	cout << setw(fw) << std::string(sw, '-');
      }
      cout << std::endl;

      first = false;
    }

    for (int i=0; i<numb; i++) {
      double energy = dkE*(0.5+i);
      cout << setw(fw) << Time 
	   << setw(fw) << energy;
      for (auto v : shist) {
	double z = v.second[i];
	cout << setw(fw) << z/total;
      }
      cout << endl;
    }
    cout << setw(fw) << Time 
	 << setw(fw) << "Overflow";
    for (auto v : shist) {
      double z = v.second[numb];
      cout << setw(fw) << z/total;
    }
    cout << endl << endl;
  }

  return 0;
}
