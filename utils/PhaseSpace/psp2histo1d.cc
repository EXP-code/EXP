/*
  Separate a psp structure and make a 1-d histogram

  MDWeinberg 08/26/11
*/

#include <cstdlib>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <list>
#include <map>

using namespace std;

#include "Species.H"

#include "StringTok.H"
#include "cxxopts.H"
#include "libvars.H"
#include "header.H"
#include "PSP.H"


int
main(int ac, char **av)
{
  char *prog = av[0];
  double pmin, pmax;
  bool mweight = true;
  bool nweight = false;
  bool areal   = false;
  bool verbose = false;
  std:: string cname;
  int axis, numb, comp, sindx, eindx;

  // Parse command line
  //
  cxxopts::Options options(prog, "Separate a psp structure and make a 1-d histogram");
    
  options.add_options()
    ("h,help", "produce help message")
    ("m,mweight", "mass-weighted values")
    ("n,nweight", "number-weighted values")
    ("A,areal", "areal average")
    ("v,verbose", "verbose output")
    ("OUT", "assume that PSP files are in original format")
    ("SPL", "assume that PSP files are in split format")
    ("a,axis", "histogram along desired axis: x=1, y=2, z=3",
     cxxopts::value<int>(axis)->default_value("3"))
    ("p,pmin", "minimum position along axis",
     cxxopts::value<double>(pmin)->default_value("-100.0"))
    ("P,pmax", "maximum position along axis",
     cxxopts::value<double>(pmax)->default_value("100.0"))
    ("b,bins", "number of bins",
     cxxopts::value<int>(numb)->default_value("40"))
    ("i,comp", "index for extended value",
     cxxopts::value<int>(comp)->default_value("9"))
    ("s,species", "position of species index",
     cxxopts::value<int>(sindx)->default_value("-1"))
    ("e,electrons", "position of electron index",
     cxxopts::value<int>(eindx)->default_value("-1"))
    ("c,name", "component name",
     cxxopts::value<std::string>(cname)->default_value("comp"))
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
	      << " --temp=25000 --number=250000 --output=out.bod" << std::endl;
    return 1;
  }

  if (vm.count("mweight")) {
    mweight = true;
    nweight = false;
  }

  if (vm.count("nweight")) {
    mweight = false;
    nweight = true;
  }

  if (vm.count("areal")) {
    areal = true;
  }

  if (vm.count("verbose")) {
    verbose = true;
  }


				// Axis sanity check 
				// ------------------
  if (axis<1) axis = 1;
  if (axis>3) axis = 3;

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
    
      cerr << "\nPhase space has time <" << psp->CurrentTime() << ">\n";
    }

				// Dump ascii for each component
				// -----------------------------
    
  
    double rtmp, mass, fac, dp=(pmax - pmin)/numb;
    vector<double> pos(3), vel(3);
    int itmp, icnt, iv;

				// Make the array
				// --------------

    vector<float> value(numb, 0), bmass(numb, 0);

    PSPstanza *stanza;
    SParticle* part;

    std::map< speciesKey, std::vector<float> > shist;

    for (stanza=psp->GetStanza(); stanza!=0; stanza=psp->NextStanza()) {
    
      if (stanza->name != cname) continue;

      for (part=psp->GetParticle(); part!=0; part=psp->NextParticle()) {

	if (part->pos(axis-1)<pmin || part->pos(axis-1)>=pmax) continue;

	iv = static_cast<int>( floor( (part->pos(axis-1) - pmin)/dp ) );
	
	double mass = part->mass();

	if (mweight) {
	  bmass[iv] += mass;
	  fac = mass;
	} else {
	  bmass[iv] += 1.0;
	  fac = 1.0;
	}

	if (comp == 0)
	  value[iv] += fac*mass;
	else if (comp <= 3)
	  value[iv] += fac*part->pos(comp-1);
	else if (comp <= 6)
	  value[iv] += fac*part->vel(comp-4);
	else if (comp == 7)
	  value[iv] += fac*part->phi();
	else if (part->niatr() && comp <= 7 + part->niatr())
	  value[iv] += fac*part->iatr(comp-8);
	else if (part->ndatr())
	  value[iv] += fac*part->datr(comp-8-part->niatr());
	
	if (sindx >= 0) {
	  KeyConvert k(part->iatr(sindx));
	  if (shist.find(k.getKey()) == shist.end()) 
	    shist[k.getKey()].resize(numb, 0);
	  shist[k.getKey()][iv] += fac;
	}

      }
    
    }
    
    //
    // Output
    //
    const size_t fw = 12;
    const size_t sw =  9;
    double Time = psp->CurrentTime();
    float p, f, m=0.0;

    if (first) {
      std::cout << setw(fw) << "Time"
		<< setw(fw) << "Position"
		<< setw(fw) << "Value"
		<< setw(fw) << "Mass";
      
      for (auto v : shist) {
	speciesKey k = v.first;
	ostringstream str;
	str << "(" << k.first << ", " << k.second << ")";
	cout << setw(fw) << str.str();
      }
      cout << std::endl;

      std::cout << setw(fw) << std::string(sw, '-')
		<< setw(fw) << std::string(sw, '-')
		<< setw(fw) << std::string(sw, '-')
		<< setw(fw) << std::string(sw, '-');

      for (auto v : shist) {
	cout << setw(fw) << std::string(sw, '-');
      }
      cout << std::endl;

      first = false;
    }

    for (int i=0; i<numb; i++) {
      p  = pmin + dp*(0.5+i);
      f  = 0.0;
      m += bmass[i];
      if (areal)  {
	f = value[i]/dp;
      } else {
	if (bmass[i] > 0.0) f = value[i]/bmass[i];
      }
      cout << setw(fw) << Time 
	   << setw(fw) << p
	   << setw(fw) << f
	   << setw(fw) << m;
      if (sindx>=0) {
	for (auto v : shist) {
	  double z = v.second[i];
	  if (areal) 
	    z /= dp;
	  else if (bmass[i] > 0.0) 
	    z /= bmass[i];
	  cout << setw(fw) << z;
	}
      }
      cout << endl;
    }
    cout << endl;
  }

  return 0;
}
