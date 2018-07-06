/*
  Separate a psp structure and make a 1-d histogram.  Trace species version.

  MDWeinberg 01/20/16
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

bool readSpeciesFile(std::string file,
		     std::map<speciesKey, int>& SpList,
		     std::list<unsigned short>& ZList,
		     int& elec)
{
  std::ifstream in(file);
  bool status = true;
  int cons;

  if (in) {
    const int nline = 2048;
    char line[nline];
      
    in.getline(line, nline);
    std::string type(line);

    if (type.compare("trace")==0) {

	in.getline(line, nline);

	if (in.good()) {
	  std::istringstream sz(line);
	  sz >> cons;
	  if (sz.good()) {
	    sz >> elec;
	  }
	} else {
	  status = false; // Can't read electrons or use_cons value, fatal
	}

	if (status) {

	  speciesKey key;
	  int pos;
	  while (1) {
	    in.getline(line, nline);
	    if (in.good()) {
	      std::istringstream sz(line);
	      sz >> key.first;
	      sz >> key.second;
	      sz >> pos;
	      // Add to the species list
	      if (!sz.bad()) {
		SpList[key] = pos;
		ZList.push_back(key.first);
	      }
	    } else {
	      break;
	    }
	  }
	}
    } else {
      status = false;
    }
    
  } else {
    status = false;
  }

  return status;
}

int
main(int ac, char **av)
{
  char *prog = av[0];
  double time, Emin, Emax, Lunit, Tunit;
  bool verbose = false, logE = false;
  std:: string cname, spfile;
  int numb, comp;

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
    ("Emin,e",		po::value<double>(&Emin)->default_value(1.0),
     "Minimum energy in eV")
    ("Emax,E",		po::value<double>(&Emax)->default_value(200.0),
     "Maximum energy in eV")
    ("bins,b",	        po::value<int>(&numb)->default_value(40),
     "number of bins")
    ("species,s",	po::value<std::string>(&spfile)->default_value("species.spec"),
     "species definition file")
    ("name,c",	        po::value<std::string>(&cname)->default_value("gas"),
     "component name")
    ("logE,L",
     "use log scaling for energy range")
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

  if (vm.count("logE")) {
    logE = true;
  }

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    std::cout << "Example: " << std::endl;
    std::cout << "\t" << av[0]
	      << " -E 300 -n 100 -f OUT.run.00001" << std::endl;
    return 1;
  }

  if (vm.count("verbose")) {
    verbose = true;
  }

  // Units
  //
  const double amu = 1.66053892e-24; // atomic mass unit in g
  const double eV  = 1.60217653e-12; // erg per eV
  double Vunit     = Lunit/Tunit;
  double KEfac     = 0.5 * amu/eV * Vunit*Vunit;

  const std::vector<double> atomic_mass = {0.000549,  // 0  electron
					   1.00797,   // 1  H
					   4.00260,   // 2  He
					   6.941,     // 3  Li
					   9.01218,   // 4  Be
					   10.81,     // 5  B
					   12.011,    // 6  C
					   14.0067,   // 7  N
					   15.9994,   // 8  O
					   18.998403, // 9  F
					   20.179,    // 10 Ne
					   22.98977,  // 11 Na
					   24.305,    // 12 Mg
					   26.98154,  // 13 Al
					   28.0855,   // 14 Si
					   30.97376,  // 15 P
					   32.06,     // 16 S
					   35.453,    // 17 Cl
					   39.948,    // 18 Ar
					   39.0983,   // 19 K
					   40.08,     // 20 Ca
					   44.9559,   // 21 Sc
					   47.90,     // 22 Ti
					   50.9415,   // 23 V
					   51.996,    // 24 Cr
					   54.9380,   // 25 Mn
					   55.847,    // 26 Fe
					   58.9332,   // 27 Co
					   58.70,     // 28 Ni
					   63.546,    // 29 Cu
					   65.38 };   // 30 Zn

  std::map<speciesKey, int> SpList;
  std::list<unsigned short> ZList;
  int elec;

  if (not readSpeciesFile(spfile, SpList, ZList, elec)) {
    std::cout << "Error reading species file <" << spfile << ">" << std::endl;
    return -1;
  }
  

  // Compute grid parameters and set up structures
  //
  if (logE) {
    if (Emin==0.0 or Emax==0.0) {
      std::cerr << "Energy must be greater than zero for log scaling"
		<< std::endl;
      exit(-2);
    }
    Emin = log10(Emin);
    Emax = log10(Emax);
  }

  double dE = (Emax - Emin)/numb;
  int nEbin = floor((Emax - Emin)/dE+1.0e-8*(Emax - Emin));
  Emax = Emin + dE*nEbin;

  std::vector<double> E(nEbin), Eion(nEbin, 0.0), Eelc(nEbin, 0.0);
  for (int i=0; i<nEbin; i++) E[i] = Emin + dE*(0.5+i);

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

				// Reopen file for data input
				// --------------------------
    delete in;
    in = new ifstream(file);
    
  
    vector<double> pos(3), vel(3);

    PSPstanza *stanza;
    SParticle* part;


				// Will contain array for each species
				//
    
    typedef std::pair< std::vector<float>, std::vector<float> > pVec;
    const std::vector<float> zero(numb+1, 0.0);
    const pVec pZero(zero, zero);

    std::map<speciesKey, pVec> shist;
    pVec thist(pZero);

    std::map<speciesKey, double> sfrac;

    double dkE = Emax/numb;
    double etotal = 0.0, itotal = 0.0, mtotal = 0.0;
    int eIout = 0, eEout = 0, eIgrid = 0, eEgrid = 0, total = 0;

    for (stanza=psp.GetStanza(); stanza!=0; stanza=psp.NextStanza()) {
    
      if (stanza->name != cname) continue;

				// Position to beginning of particles
				// 
      in->seekg(stanza->pspos);

      for (part=psp.GetParticle(in); part!=0; part=psp.NextParticle(in)) {

	double kEe = 0.0, kEi = 0.0;
	for (size_t i=0; i<3; i++) {
	  double ve = part->datr(elec+i);
	  kEe += ve*ve;
	  double vi = part->vel(i);
	  kEi += vi*vi;
	}

	// Molecular weight
	//

	double efrac = 0.0, mu = 0.0;
	
	for (auto v : SpList) {
	  unsigned Z = v.first.first;
	  unsigned C = v.first.second;
	  double   W = part->datr(v.second);
	  mu += W/atomic_mass[Z];
	  efrac += W/atomic_mass[Z]*(C-1);
	}
	efrac /= mu;
	mu = 1.0/mu;

	kEe *= KEfac * atomic_mass[0];
	kEi *= KEfac * mu;

	if (logE) {
	  kEe = log10(kEe);
	  kEi = log10(kEi);
	}

	if (kEe >= Emin and kEe < Emax) {
	  int Eindx = floor( (kEe - Emin)/dE );
	  if (Eindx >= 0 and Eindx < nEbin) {
	    Eelc[Eindx] += efrac * part->mass()/mu;
	    eEgrid++;
	  }
	  else eEout++;
	}

	if (kEi >= Emin and kEi < Emax) {
	  int Eindx = floor( (kEi - Emin)/dE );
	  if (Eindx >= 0 and Eindx < nEbin) {
	    Eion[Eindx] += part->mass()/mu;
	    eIgrid++;
	  }
	  else eIout++;
	}

	total++;
	
      } // END: Particle loop

      std::cout << "File <" << file << ">: "
		<< eIgrid << "/" << total << " ions and "
		<< eEgrid << "/" << total << " electrons with "
		<< eEout << " electron oab, "
		<< eIout << " ion oab" << std::endl;

    } // END: stanza loop


    //
    // Normalize
    //
    double normI=0.0, normE=0.0;
    for (auto  v : Eion) normI += v;
    for (auto  v : Eelc) normE += v;
    for (auto& v : Eion) v /= normI;
    for (auto& v : Eelc) v /= normE;

    //
    // Output
    //
    const size_t fw = 14;
    const size_t sw =  9;
    double Time = psp.CurrentTime();

    if (first) {
      std::cout << setw(fw) << "Time"
		<< setw(fw) << "Energy"
		<< setw(fw) << "Ions"
		<< setw(fw) << "Electrons"
		<< std::endl;
      

      std::cout << setw(fw) << std::string(sw, '-')
		<< setw(fw) << std::string(sw, '-')
		<< setw(fw) << std::string(sw, '-')
		<< setw(fw) << std::string(sw, '-')
		<< std::endl;

      first = false;
    }

    for (int i=0; i<nEbin; i++)
      cout << setw(fw) << Time
	   << setw(fw) << (logE ? exp(E[i]) : E[i])
	   << setw(fw) << Eion[i]
	   << setw(fw) << Eelc[i]
	   << std::endl;

    cout << std::string(12+18*2, '-') << endl << endl;


    } // END: file loop


  return 0;
}