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

#include <header.H>
#include <PSP.H>

#include <boost/program_options.hpp>
#include <boost/random/mersenne_twister.hpp>

namespace po = boost::program_options;


				// Globals for exputil library
				// Unused here
int myid = 0;
char threading_on = 0;
pthread_mutex_t mem_lock;
string outdir, runtag;
#include <boost/random/mersenne_twister.hpp>

int
main(int ac, char **av)
{
  char *prog = av[0];
  bool verbose = false;
  bool areal   = false;
  bool vnorm   = false;
  bool snorm   = false;
  bool use_sph = false;
  bool use_cyl = false;
  std::string cname, new_dir;
  int axis, numb, comp;
  double pmin, pmax;

  // Parse command line

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",		"produce help message")
    ("radial,r",        "use spherical radius")
    ("cylindrical,R",   "cylindrical radius")
    ("areal,A",         "areal average")
    ("vnorm,V",         "compute density for radial bins")
    ("snorm,S",         "compute surface density for cylindrical bins")
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
    ("dir,d",           po::value<std::string>(&new_dir)->default_value("./"),
     "rewrite directory location for SPL files")
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
    snorm   = false;
    if (vm.count("vnorm")) vnorm = true;
  }

  if (vm.count("cylindrical")) {
    use_sph = false;
    use_cyl = true;
    if (vm.count("snorm")) snorm = true;
    vnorm   = false;
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

    if (verbose) cerr << "Using filename: " << file << endl;

				// Parse the PSP file
				// ------------------
    PSPptr psp;
    if (file.find("SPL") != std::string::npos)
      psp = std::make_shared<PSPspl>(file, new_dir);
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
  
    double rtmp, mass, dp=(pmax - pmin)/numb;
    vector<double> pos(3), vel(3);
    int itmp, icnt, iv;

				// Make the array
				// --------------

    vector<float> value(numb, 0), bmass(numb, 0);

    PSPstanza *stanza;
    SParticle* part;

    for (stanza=psp->GetStanza(); stanza!=0; stanza=psp->NextStanza()) {
    
      if (stanza->name != cname) continue;

      for (part=psp->GetParticle(); part!=0; part=psp->NextParticle()) {

	double val = 0.0;
	if (use_sph) {
	  for (int k=0; k<3; k++) {
	    val += part->pos(k) * part->pos(k);
	  }
	  val = sqrt(val);
	  iv = static_cast<int>( floor( (val - pmin)/dp ) );
	}
	else if (use_cyl) {
	  for (int k=0; k<2; k++) {
	    val += part->pos(k) * part->pos(k);
	  }
	  val = sqrt(val);
	  iv = static_cast<int>( floor( (val - pmin)/dp ) );
	}
	else {
	  val = part->pos(axis-1);
	  iv = static_cast<int>( floor( (part->pos(axis-1) - pmin)/dp ) );
	}

	if (iv < 0 || iv >= numb) continue;

	bmass[iv] += 1.0;

	if (comp == 0)
	  value[iv] += part->mass();
	else if (comp <= 3)
	  value[iv] += part->pos(comp-1);
	else if (comp <= 6)
	  value[iv] += part->vel(comp-4);
	else if (comp == 7)
	  value[iv] += part->phi();
	else if (part->niatr() && comp <= 7 + part->niatr())
	  value[iv] += part->iatr(comp-8);
	else if (part->ndatr())
	  value[iv] += part->datr(comp-8-part->niatr());
      }
    
    }
    
    //
    // Output
    //
    const size_t fw = 12;
    const size_t sw =  9;
    float p, f, m=0.0;

    if (first) {
      std::cout << setw(fw) << "Position"
		<< setw(fw) << "Value"
		<< setw(fw) << "Mass"
		<< std::endl;

      std::cout << setw(fw) << std::string(sw, '-')
		<< setw(fw) << std::string(sw, '-')
		<< setw(fw) << std::string(sw, '-')
		<< std::endl;

      first = false;
    }

    for (int i=0; i<numb; i++) {
      p  = pmin + dp*(0.5+i);
      f  = 0.0;
      m += bmass[i];
      if (vnorm) {
	double rmax3 = pow(pmin+dp*(1.0+i), 3.0);
	double rmin3 = pow(pmin+dp*(0.0+i), 3.0);
	f = value[i]/(4.0*M_PI/3.0*(rmax3 - rmin3));
      } else if (snorm) {
	double rmax2 = pow(pmin+dp*(1.0+i), 2.0);
	double rmin2 = pow(pmin+dp*(0.0+i), 2.0);
	f = value[i]/(M_PI*(rmax2 - rmin2));
      } else if (areal) {
	f = value[i]/dp;
      } else {
	if (bmass[i] > 0.0) f = value[i]/bmass[i];
      }
      cout << setw(fw) << p
	   << setw(fw) << f
	   << setw(fw) << m
	   << endl;
    }
    cout << endl;
  }

  return 0;
}
