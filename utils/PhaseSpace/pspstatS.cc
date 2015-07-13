/*
  Compute simple statistics from psp dump

  MDWeinberg 06/10/02
*/

using namespace std;

#include <unistd.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <tuple>
#include <list>

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

double atomic_masses[3] = {0.000548579909, 1.00794, 4.002602};

int
main(int ac, char **av)
{
  char *prog = av[0];
  double time;
  int sindx, eindx;
  std::string cname;
  bool verbose = false;

  // Parse command line

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",		"produce help message")
    ("verbose,v",       "verbose output")
    ("time,t",		po::value<double>(&time)->default_value(1.0e20),
     "find closest time slice to requested value")
    ("species,s",	po::value<int>(&sindx)->default_value(-1),
     "position of species index")
    ("electrons,e",	po::value<int>(&eindx)->default_value(7),
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
	      << " -c gas -s 0 -f OUT.run.00001" << std::endl;
    return 1;
  }

  if (vm.count("verbose")) {
    verbose = true;
  }


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

				// Will contain array for each gas species
				// ---------------------------------------
    typedef std::tuple<double, unsigned> shtup;
    typedef std::map<speciesKey, shtup> shist;
    const shtup tupzero(0.0, 0);

    PSPstanza *stanza;
    SParticle* part;
    double rtmp;

    for (stanza=psp.GetStanza(); stanza!=0; stanza=psp.NextStanza()) {
      
      if (stanza->name != cname) continue;


				// Setup stats for each component
				// -----------------------------

      double com1[3] = {0.0, 0.0, 0.0};
      double cov1[3] = {0.0, 0.0, 0.0};
      double ang1[3] = {0.0, 0.0, 0.0};
      double KE1     = 0.0;
      double PE1     = 0.0;
      double EE1     = 0.0;
      double mass1   = 0.0;

      shist  hist1;		// For gas only

				// Phase space stuff
				// -----------------
      double ms;
      double pos[3];
      double vel[3];
      double mom[3];
      double pot;

				// Print the header

      cout << "Comp name: " << stanza->name << endl
	   << "     Bodies:\t\t"
	   << setw(15) << stanza->comp.nbod 
	   << setw(10) << stanza->comp.niatr 
	   << setw(10) << stanza->comp.ndatr 
	   << endl;


				// Position to beginning of particles
      in->seekg(stanza->pspos);

      for (part=psp.GetParticle(in); part!=0; part=psp.NextParticle(in)) {

	mom[0] = part->pos(1)*part->vel(2) - part->pos(2)*part->vel(1);
	mom[1] = part->pos(2)*part->vel(0) - part->pos(0)*part->vel(2);
	mom[2] = part->pos(0)*part->vel(1) - part->pos(1)*part->vel(0);

				// Accumulate statistics
	double ms = part->mass();
	mass1 += ms;
	for (int i=0; i<3; i++) com1[i] += ms*part->pos(i);
	for (int i=0; i<3; i++) cov1[i] += ms*part->vel(i);
	for (int i=0; i<3; i++) ang1[i] += ms*mom[i];
	rtmp = 0.0;
	for (int i=0; i<3; i++) rtmp += part->vel(i)*part->vel(i);
	KE1 += 0.5*ms*rtmp;
	PE1 += 0.5*ms*part->phi();

	if (sindx >= 0) {
	  KeyConvert kc(part->iatr(sindx));
	  speciesKey k = kc.getKey();
	  if (hist1.find(k) == hist1.end()) hist1[k] = tupzero;
	  std::get<0>(hist1[k]) += ms;
	  std::get<1>(hist1[k]) ++;
	  double ke = 0.0;
	  for (int i=0; i<3; i++) {
	    double t = part->datr(eindx+i);
	    ke += t*t;
	  }
	  EE1 += ke * 0.5*ms*atomic_masses[0]/atomic_masses[k.first];
	}
      }
    
      cout  << "     COM:\t\t";
      for (int i=0; i<3; i++) cout << setw(15) << com1[i]/mass1;
      cout << endl;
      cout  << "     COV:\t\t";
      for (int i=0; i<3; i++) cout << setw(15) << cov1[i]/mass1;
      cout << endl;
      cout  << "     Ang mom:\t\t";
      for (int i=0; i<3; i++) cout << setw(15) << ang1[i];
      cout << endl;
      cout  << "     Stats:\t\tKE=" << KE1 << " PE=" << PE1 << " -2T/W=" << -2.0*KE1/PE1 << " E_e=" << EE1
	    << " Mass=" << mass1 << endl;

      if (sindx >= 0) {
	cout  << "     Species" << endl;
	for (auto v : hist1) {
	  speciesKey k = v.first;
	  cout  << "            <"
		<< setw(2) << k.first << "," << setw(2) << k.second << ">"
		<< "  :  "
		<< std::setw(16) << std::get<0>(v.second)
		<< std::setw(10) << std::get<1>(v.second)
		<< endl;
	}
      }
    }
  }
  
  return 0;
}
  
