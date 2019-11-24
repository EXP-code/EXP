/*
  Separate a psp structure and make a 1-d histogram

  MDWeinberg 11/24/19
*/

using namespace std;

#include <cstdlib>

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
#include <PSP2.H>

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
  bool verbose = false;
  bool areal   = false;
  bool use_sph = false;
  bool use_cyl = false;
  std:: string cname;
  int axis, numb, comp, sindx, eindx;
  double pmin, pmax;

  // Parse command line

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",		"produce help message")
    ("radial,r",        "use spherical radius")
    ("cylindrical,R",   "cylindrical radius")
    ("areal,A",         "areal average")
    ("verbose,v",       "verbose output")
    ("pmin,p",	        po::value<double>(&pmin)->default_value(-100.0),
     "minimum position along axis")
    ("pmax,P",	        po::value<double>(&pmax)->default_value(100.0),
     "maximum position along axis")
    ("bins,b",	        po::value<int>(&numb)->default_value(40),
     "number of bins")
    ("comp,i",		po::value<int>(&comp)->default_value(9),
     "index for extended value")
    ("name,c",	        po::value<std::string>(&cname)->default_value("comp"),
     "component name")
    ("axis,a",		po::value<int>(&axis)->default_value(3),
     "histogram along desired axis: x=1, y=2, z=3")
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
	      << " --output=out.bod" << std::endl;
    return 1;
  }

  if (vm.count("radial")) {
    use_sph = true;
    use_cyl = false;
  }

  if (vm.count("cylindrical")) {
    use_sph = false;
    use_cyl = true;
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
      std::cerr << "Error opening file <" << file << "> for input" << std::endl;
      exit(-1);
    }
    in.close();

    if (verbose) cerr << "Using filename: " << file << endl;

				// Parse the PSP file
				// ------------------
    PSPptr psp;
    if (file.find("SPL") != std::string::npos)
      psp = std::make_shared<PSPspl>(file);
    else
      psp = std::make_shared<PSPout>(file);


				// Now write a summary
				// -------------------
    if (verbose) {
      
      psp->PrintSummary(cerr);
    
      std::cerr << "\nPSP file <" << file << "> has time <" 
	   << psp->CurrentTime() << ">\n";
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

	double val = 0.0;
	if (use_sph) {
	  for (int k=0; k<3; k++) {
	    val += part->pos(k) * part->pos(k);
	  }
	  val = sqrt(val);
	}
	else if (use_cyl) {
	  for (int k=0; k<2; k++) {
	    val += part->pos(k) * part->pos(k);
	  }
	  val = sqrt(val);
	}
	else {
	  val = part->pos(axis-1);
	}

	if (val<pmin || val>=pmax) continue;
	iv = static_cast<int>( floor( (part->pos(axis-1) - pmin)/dp ) );
	
	double mass = part->mass();
	
	bmass[iv] += 1.0;
	fac = 1.0;

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
