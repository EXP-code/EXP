/*
  Separate a psp structure and make a 1-d histogram

  MDWeinberg 08/26/11
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
  double time, pmin, pmax;
  bool mweight = true;
  bool nweight = false;
  bool areal   = false;
  bool verbose = false;
  std:: string cname;
  int axis, numb, comp, sindx;

  // Parse command line

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",		"produce help message")
    ("mweight,m",       "mass-weighted values")
    ("nweight,n",       "number-weighted values")
    ("areal,A",         "areal average")
    ("verbose,v",       "verbose output")
    ("time,t",		 po::value<double>(&time)->default_value(1.0e20),
     "find closest time slice to requested value")
    ("axis,a",		po::value<int>(&axis)->default_value(3),
     "histogram along desired axis: x=1, y=2, z=3")
    ("pmin,p",	        po::value<double>(&pmin)->default_value(-100.0),
     "minimum position along axis")
    ("pmax,P",	        po::value<double>(&pmax)->default_value(100.0),
     "maximum position along axis")
    ("bins,b",	        po::value<int>(&numb)->default_value(40),
     "number of bins")
    ("comp,i",		po::value<int>(&comp)->default_value(9),
     "index for extended value")
    ("species,s",	po::value<int>(&sindx)->default_value(-1),
     "position of species index")
    ("name,c",	        po::value<std::string>(&cname)->default_value("comp"),
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

  ifstream *in;

  std::vector<std::string> files = vm["files"].as< std::vector<std::string> >();

  if (vm.count("files")) {
    ifstream *in2 = new ifstream(files[0].c_str());
    if (!*in2) {
      cerr << "Error opening file <" << files[0] << "> for input\n";
      exit(-1);
    }

    if (verbose) cerr << "Using filename: " << files[0] << endl;

				// Assign file stream to input stream
    in = in2;

  } else {
    std::cout << desc << std::endl;
  }

				// Axis sanity check 
				// ------------------
  if (axis<1) axis = 1;
  if (axis>3) axis = 3;

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
  in = new ifstream(files[0]);

  
  double rtmp, mass, fac, dp=(pmax - pmin)/numb;
  vector<double> pos(3), vel(3);
  int itmp, icnt, iv;

				// Make the array
				// --------------

  vector<float> value(numb, 0), bmass(numb, 0);

  PSPstanza *stanza;
  SParticle* part;

  std::map< speciesKey, std::vector<float> > shist;

  for (stanza=psp.GetStanza(); stanza!=0; stanza=psp.NextStanza()) {
    
    if (stanza->name != cname) continue;


				// Position to beginning of particles
    in->seekg(stanza->pspos);

    for (part=psp.GetParticle(in); part!=0; part=psp.NextParticle(in)) {

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
  double Time = psp.CurrentTime();
  float p, f, m=0.0;

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

  return 0;
}
