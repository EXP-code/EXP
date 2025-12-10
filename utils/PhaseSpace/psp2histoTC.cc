/*
  Separate a psp structure and make a 1-d histogram.  Trace species version.

  MDWeinberg 01/20/16
*/

#include <cstdlib>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <vector>
#include <string>
#include <list>
#include <map>

using namespace std;

#include <yaml-cpp/yaml.h>	// YAML support
#include "cxxopts.H"		// Command-line parsing
#include "libvars.H"		// EXP library globals

#include "Species.H"

#include "StringTok.H"
#include "header.H"
#include "PSP.H"

bool readSpeciesFileOld(std::string file,
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

bool readSpeciesFile(std::string file,
		     std::map<speciesKey, int>& SpList,
		     std::list<unsigned short>& ZList,
		     int& elec)
{
  YAML::Node conf;

  // Try to make a YAML config, fall back to old style species spec on
  // failure
  //
  try {
    conf = YAML::LoadFile(file);
  }
  catch (YAML::Exception & error) {
    std::cout << "Error parsing component config.  Trying old-style PSP"
	      << std::endl
	      << error.what() << std::endl;

    return readSpeciesFileOld(file, SpList, ZList, elec);
  }


  // Check for species map
  //
  if (conf["species_map"]) {

    YAML::Node cconf = conf["species_map"];
    std::string type = cconf["type"].as<std::string>();

    if (type.compare("trace")==0) {
      
      if (cconf["elec"]) {
	elec = cconf["elec"].as<int>();
      } else {
	std::cerr << "Error reading elec field for trace species" << std::endl;
	return false;
      }

      if (cconf["elements"]) {
	YAML::Node elems = cconf["elements"];

	for (YAML::const_iterator it=elems.begin(); it!=elems.end(); ++it) {
	  auto dd = it->as<std::vector<int>>();
	  speciesKey key(dd[0], dd[1]);
	  SpList[key] = dd[2];
	  ZList.push_back(key.first);
	}
      } else {
	std::cerr << "No <elements> key found . . . quitting" << std::endl;
	return false;
      }

      return true;

    } else {
      std::cerr << "Method/map type is <" << type << ">, expected <trace>" << std::endl;
      return false;
    }

  } else {
    std::cerr << "Could not locate species_map in YAML file <" << file << ">" << std::endl;
    return false;
  }
}

int
main(int ac, char **av)
{
  char *prog = av[0];
  double Emin, Emax, Lunit, Tunit, Temp;
  bool verbose = false, logE = false, flat = false, meanmass = false;
  std:: string cname, spfile, runtag;
  int numb, comp, ibeg, iend;

  // Parse command line
  //
  cxxopts::Options options(prog, "Separate a psp structure and make a 1-d histogram.  Trace species version.");

  options.add_options()
   ("h,help", "produce help message")
   ("v,verbose", "verbose output")
   ("OUT", "assume that PSP files are in original format")
   ("SPL", "assume that PSP files are in split format")
   ("L,Lunit", "System length in physical units (cgs)",
     cxxopts::value<double>(Lunit)->default_value("3.086e18"))
   ("T,Tunit", "System time in physical units (cgs)",
     cxxopts::value<double>(Tunit)->default_value("3.15569e10"))
   ("e,Emin", "Minimum energy in eV",
     cxxopts::value<double>(Emin)->default_value("1.0"))
   ("E,Emax", "Maximum energy in eV",
     cxxopts::value<double>(Emax)->default_value("200.0"))
   ("K,Temp", "Temperature in kelvins",
     cxxopts::value<double>(Temp)->default_value("3.0e4"))
   ("b,bins", "number of bins",
     cxxopts::value<int>(numb)->default_value("40"))
   ("s,species", "species definition file",
     cxxopts::value<std::string>(spfile)->default_value("species.yml"))
   ("c,name", "component name",
     cxxopts::value<std::string>(cname)->default_value("gas"))
   ("r,runtag", "runtag for using a range of PSP files",
     cxxopts::value<std::string>(runtag)->default_value("run"))
   ("beg", "initial value for PSP file sequence",
     cxxopts::value<int>(ibeg)->default_value("100"))
   ("end", "final value for PSP file sequence",
     cxxopts::value<int>(iend)->default_value("200"))
   ("meanMass", "use mean-mass algorithm for electron energy computation")
   ("logE", "use log scaling for energy range")
   ("flat", "use E^{3/2} scaling for energy range")
   ("f,files", "input files",
     cxxopts::value< std::vector<std::string> >())
    ;

  cxxopts::ParseResult vm;

  try {
    vm = options.parse(ac, av);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("logE")) {
    logE = true;
    flat = false;
  }

  if (vm.count("flat")) {
    logE = false;
    flat = true;
  }

  if (vm.count("meanMass")) {
    meanmass = true;
  }

  if (vm.count("help")) {
    std::cout << options.help() << std::endl;
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
  const double amu   = 1.66053892e-24; // atomic mass unit in g
  const double eV    = 1.60217653e-12; // erg per eV
  const double boltz = 1.3806504e-16;  // Boltzmann constant
  double Vunit       = Lunit/Tunit;    // system to cgs conversion
  double kT          = boltz*Temp/eV;  // kT in eV units
				// KE conversion from system to eV 
  double KEfac       = 0.5 * amu/eV * Vunit*Vunit;

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
    Emin = log(Emin);
    Emax = log(Emax);
  }

  if (flat) {
    Emin = pow(Emin, 1.5);
    Emax = pow(Emax, 1.5);
  }

  double dE = (Emax - Emin)/numb;
  int nEbin = floor((Emax - Emin)/dE+1.0e-8*(Emax - Emin));
  Emax = Emin + dE*nEbin;

  std::vector<double> E(nEbin), Eion(nEbin, 0.0), Eelc(nEbin, 0.0);
  std::vector<unsigned> Nion(nEbin, 0), Nelc(nEbin, 0);
  for (int i=0; i<nEbin; i++) E[i] = Emin + dE*(0.5+i);

  // Get file arguments
  //
  std::vector<std::string> files;
  if (vm.count("runtag")) {
    for (int i=ibeg; i<=iend; i++) {
      std::ostringstream str;
      str << "OUT." << runtag << "."
	  << std::setw(5) << std::setfill('0') << std::right << i;
      files.push_back(str.str());
    }
  } else {
    files = vm["files"].as< std::vector<std::string> >();
  }

  for (auto file : files ) {

    auto in = std::make_shared<std::ifstream>(file.c_str());
    if (!*in) {
      cerr << "Error opening file <" << file << "> for input\n";
      exit(-1);
    }

    if (verbose) cerr << "Using filename: " << file << endl;


				// Parse the PSP file
				// ------------------
    PSPptr psp;
    if (vm.count("SPL")) psp = std::make_shared<PSPspl>(file);
    else                 psp = std::make_shared<PSPout>(file);

				// Now write a summary
				// -------------------
    if (verbose) {
      
      psp->PrintSummary(cerr);
    
      std::cerr << std::endl << "Phase space has time <" 
		<< psp->CurrentTime() << ">" << std::endl;
    }

				// Reopen file for data input
				// --------------------------
    in = std::make_shared<std::ifstream>(file);
    
  
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

    for (stanza=psp->GetStanza(); stanza!=0; stanza=psp->NextStanza()) {
    
      if (stanza->name != cname) continue;

      // Position to beginning of particles
      // 
      in->seekg(stanza->pspos);

      for (part=psp->GetParticle(); part!=0; part=psp->NextParticle()) {

	double kEe = 0.0, kEi = 0.0;
	for (size_t i=0; i<3; i++) {
	  double ve = part->datr(elec+i);
	  kEe += ve*ve;
	  double vi = part->vel(i);
	  kEi += vi*vi;
	}

	// Molecular weight computation
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

	if (meanmass) kEe *= efrac;

	if (logE) {
	  kEe = log(kEe);
	  kEi = log(kEi);
	}

	if (flat) {
	  kEe = pow(kEe, 1.5);
	  kEi = pow(kEi, 1.5);
	}

	if (kEe >= Emin and kEe < Emax) {
	  int Eindx = floor( (kEe - Emin)/dE );
	  if (Eindx >= 0 and Eindx < nEbin) {
	    Eelc[Eindx] += efrac * part->mass()/mu;
	    Nelc[Eindx]++;
	    eEgrid++;
	  }
	  else eEout++;
	}

	if (kEi >= Emin and kEi < Emax) {
	  int Eindx = floor( (kEi - Emin)/dE );
	  if (Eindx >= 0 and Eindx < nEbin) {
	    Eion[Eindx] += part->mass()/mu;
	    Nion[Eindx]++;
	    eIgrid++;
	  }
	  else eIout++;
	}

	total++;
	
      } // END: Particle loop

      std::cerr << "File <" << file << ">: "
		<< eIgrid << "/" << total << " ions and "
		<< eEgrid << "/" << total << " electrons with "
		<< eEout << " electron oab, "
		<< eIout << " ion oab"
		<< std::endl;

    }
    // END: stanza loop


  }
  // END: file loop


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
  
  std::cout << "# "
	    << setw(fw) << "Energy"
	    << setw(fw) << "Ions"
	    << setw(sw) << "N(ion)"
	    << setw(fw) << "Electrons"
	    << setw(sw) << "N(elc)"
	    << setw(fw) << "Exact"
	    << std::endl;
    
    
  std::cout << "# "
	    << setw(fw) << std::string(fw-1, '-')
	    << setw(fw) << std::string(fw-1, '-')
	    << setw(sw) << std::string(sw-1, '-')
	    << setw(fw) << std::string(fw-1, '-')
	    << setw(sw) << std::string(sw-1, '-')
	    << setw(fw) << std::string(fw-1, '-')
	    << std::endl;
  
  double norm = M_2_SQRTPI*pow(kT, -1.5);
  for (int i=0; i<nEbin; i++) {
    double Energy = E[i];
    double exact  = norm * sqrt(Energy) * exp(-Energy/kT);

    if (logE) {
      Energy = exp(Energy);
      exact  = norm * pow(Energy, 1.5) * exp(-Energy/kT);
    }
    
    if (flat) {
      Energy = pow(Energy, 2.0/3.0);
      exact  = norm * 2.0/3.0 * exp(-Energy/kT);
    }
    
    cout << "  "
	 << setw(fw) << Energy
	 << setw(fw) << Eion[i]
	 << setw(sw) << Nion[i]
	 << setw(fw) << Eelc[i]
	 << setw(sw) << Nelc[i]
	 << setw(fw) << exact * dE
	 << std::endl;
  }
  
  cout << '#' << std::string(fw*5+sw*2-1, '-') << endl << endl;


  return 0;
}
