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
#include <cxxopts.H>
#include <header.H>
#include <PSP.H>

int
main(int ac, char **av)
{
  char *prog = av[0];
  double Emax, Lunit, Tunit;
  bool verbose = false;
  std:: string cname;
  int numb, comp, sindx, eindx;

  // Parse command line
  //
  cxxopts::Options options(prog, "Separate a psp structure and make a 1-d histogram.\n");

  options.add_options()
   ("h,help", "produce help message")
   ("v,verbose", "verbose output")
   ("OUT", "assume that PSP files are in original format")
   ("SPL", "assume that PSP files are in split format")
   ("L,Lunit", "System length in physical units (cgs)",
     cxxopts::value<double>(Lunit)->default_value("3.086e18"))
   ("T,Tunit", "System time in physical units (cgs)",
     cxxopts::value<double>(Tunit)->default_value("3.15569e12"))
   ("E,Emax", "Maximum energy in eV",
     cxxopts::value<double>(Emax)->default_value("200.0"))
   ("b,bins", "number of bins",
     cxxopts::value<int>(numb)->default_value("40"))
   ("s,species", "position of species index",
     cxxopts::value<int>(sindx)->default_value("0"))
   ("e,electrons", "position of electron index",
     cxxopts::value<int>(eindx)->default_value("7"))
   ("c,name", "component name",
     cxxopts::value<std::string>(cname)->default_value("gas"))
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

  // Get file arguments
  //
  std::vector<std::string> files = vm["files"].as< std::vector<std::string> >();

  bool first = true;
  
  for (auto file : files ) {

    std::ifstream in(file);
    if (!in) {
      std::cerr << "Error opening file <" << file << "> for input\n";
      exit(-1);
    }
    in.close();

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
    
      cerr << "\nBest fit dump to <" << time << "> has time <" 
	   << psp->CurrentTime() << ">\n";
    }

    vector<double> pos(3), vel(3);

    PSPstanza *stanza;
    SParticle* part;


				// Will contain array for each species
				//
    typedef std::pair< std::vector<float>, std::vector<float> > pVec;
    const std::vector<float> zero(numb+1, 0.0);
    const pVec pZero(zero, zero);

    std::map< speciesKey, pVec> shist;
    pVec thist(pZero);

    double dkE   = Emax/numb;
    double etotal = 0.0, itotal = 0.0;

    for (stanza=psp->GetStanza(); stanza!=0; stanza=psp->NextStanza()) {
    
      if (stanza->name != cname) continue;

      for (part=psp->GetParticle(); part!=0; part=psp->NextParticle()) {

	double kEe = 0.0, kEi = 0.0;
	for (size_t i=0; i<3; i++) {
	  double ve = part->datr(eindx+i);
	  kEe += ve*ve;
	  double vi = part->vel(i);
	  kEi += vi*vi;
	}
	KeyConvert kc(part->iatr(sindx));

	kEe *= KEfac * atomic_mass[0];
	kEi *= KEfac * atomic_mass[kc.Z()];

	speciesKey k = kc.getKey();
	if (shist.find(k) == shist.end()) shist[k] = pZero;
      
	size_t le   = std::min<size_t>(kEe/dkE, numb);
	size_t li   = std::min<size_t>(kEi/dkE, numb);
	double wgte = part->mass() * (kc.C() - 1);
	double wgti = part->mass() ;

	shist[k].first[le]  += wgte;
	thist   .first[le]  += wgte;
	etotal              += wgte;

	shist[k].second[li] += wgti;
	thist   .second[li] += wgti;
	itotal              += wgti;
      }
    }
    
    //
    // Output
    //
    const size_t fw = 14;
    const size_t sw =  9;
    double Time = psp->CurrentTime();

    if (first) {
      std::cout << setw(fw) << "Time"
		<< setw(fw) << "Energy";
      
      for (auto v : shist) {
	speciesKey k = v.first;
	ostringstream stre, stri;
	stre << "(" << k.first << "," << k.second << ")_e";
	stri << "(" << k.first << "," << k.second << ")_i";
	cout << setw(fw) << stre.str() << setw(fw) << stri.str();
      }
      cout << setw(fw) << "Total_e" << setw(fw) << "Total_i" << std::endl;

      std::cout << setw(fw) << std::string(sw, '-')
		<< setw(fw) << std::string(sw, '-');

      for (auto v : shist) {
	cout << setw(fw) << std::string(sw, '-')
	     << setw(fw) << std::string(sw, '-');
      }
      cout << setw(fw) << std::string(sw, '-')
	   << setw(fw) << std::string(sw, '-') << std::endl;

      first = false;
    }

    for (int i=0; i<numb; i++) {
      double energy = dkE*(0.5+i);
      cout << setw(fw) << Time 
	   << setw(fw) << energy;
      for (auto v : shist) {
	double ze = v.second.first [i];
	double zi = v.second.second[i];
	cout << setw(fw) << ze/etotal
	     << setw(fw) << zi/itotal;
      }
      cout << setw(fw) << thist.first [i]/etotal 
	   << setw(fw) << thist.second[i]/itotal 
	   << endl;
    }
    cout << setw(fw) << Time 
	 << setw(fw) << "Overflow";
    for (auto v : shist) {
      double ze = v.second.first [numb];
      double zi = v.second.second[numb];
      cout << setw(fw) << ze/etotal
	   << setw(fw) << zi/itotal;
    }
    cout << setw(fw) << thist.first [numb]/etotal 
	 << setw(fw) << thist.second[numb]/itotal 
	 << endl << endl;
  }

  return 0;
}
