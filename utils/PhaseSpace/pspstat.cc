/*
  Compute simple statistics from psp dump

  MDWeinberg 06/10/02, 11/24/19
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
#include <memory>
#include <list>

#include <StringTok.H>
#include <header.H>
#include <PSP2.H>

				// Globals for exputil library
				// Unused here
int myid = 0;
char threading_on = 0;
pthread_mutex_t mem_lock;
string outdir, runtag;

//-------------
// Help message
//-------------

void Usage(char* prog) {
  cerr << prog << ": [-t time -v -h] filename\n\n";
  cerr << "    -o name         prefix name for each component (default: comp)\n";
  cerr << "    -d dir          replacement SPL file directory\n";
  cerr << "    -h              print this help message\n";
  cerr << "    -v              verbose output\n\n";
  exit(0);
}


int
main(int argc, char **argv)
{
  char *prog = argv[0];
  double time=1e20;
  bool verbose = false;
  string cname("comp"), new_dir("");

  // Parse command line

  while (1) {

    int c = getopt(argc, argv, "o:vh");

    if (c == -1) break;

    switch (c) {

    case 'v':
      verbose = true;
      break;

    case 'o':
      cname.erase();
      cname = string(optarg);
      break;

    case 'd':
      new_dir.erase();
      new_dir = string(optarg);
      break;

    case '?':
    case 'h':
    default:
      Usage(prog);
    }

  }

  std::ifstream in;
  std::string filename;

  if (optind < argc) {

    filename = std::string(argv[optind]);

    std::ifstream in2(filename);
    if (!in2) {
      cerr << "Error opening file <" << filename << "> for input\n";
      exit(-1);
    }

    if (verbose) cerr << "Using filename: " << filename << endl;
  }

				// Parse the PSP file
				// ------------------
  PSPptr psp;
  if (filename.find("SPL") != std::string::npos)
    psp = std::make_shared<PSPspl>(filename, new_dir);
  else
    psp = std::make_shared<PSPout>(filename);
  

				// Now write a summary
				// -------------------
  if (verbose) {
    psp->PrintSummary(cerr);
    
    std::cout << std::endl << "PSP file <" << filename << "> has time <" 
	      << psp->CurrentTime() << ">" << std::endl;
  }

				// Setup stats for all components
				// ------------------------------

  double com[3] = {0.0, 0.0, 0.0};
  double cov[3] = {0.0, 0.0, 0.0};
  double ang[3] = {0.0, 0.0, 0.0};
  double KE     = 0.0;
  double PE     = 0.0;
  double mass   = 0.0;
  int   totbod  = 0;
  
  PSPstanza *stanza;
  SParticle* part;
  double rtmp;

  for (stanza=psp->GetStanza(); stanza!=0; stanza=psp->NextStanza()) {


				// Setup stats for each component
				// -----------------------------

    double com1[3] = {0.0, 0.0, 0.0};
    double cov1[3] = {0.0, 0.0, 0.0};
    double ang1[3] = {0.0, 0.0, 0.0};
    double KE1     = 0.0;
    double PE1     = 0.0;
    double mass1   = 0.0;

				// Phase space stuff
				// -----------------
    double ms;
    double pos[3];
    double vel[3];
    double mom[3];
    double pot;

				// Open an output file
				// -------------------

    ostringstream oname;
    oname << cname << "." << stanza->name << '\0';
    ofstream out(oname.str().c_str());
    out.setf(ios::scientific);
    out.precision(10);

    if (!out) {
      cerr << "Couldn't open output name <" << oname.str() << ">\n";
      exit(-1);
    }

				// Print the header

    cout << "Comp name: " << stanza->name << endl
	 << "     Bodies:\t\t"
	 << setw(15) << stanza->comp.nbod 
	 << setw(10) << stanza->comp.niatr 
	 << setw(10) << stanza->comp.ndatr 
	 << endl;

    totbod += stanza->comp.nbod;

    for (part=psp->GetParticle(); part!=0; part=psp->NextParticle()) {

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
    cout  << "     Stats:\t\tKE=" << KE1 << " PE=" << PE1 << " -2T/W=" << -2.0*KE1/PE1
	  << " Mass=" << mass1 << endl;

    mass += mass1;
    for (int i=0; i<3; i++) com[i] += com1[i];
    for (int i=0; i<3; i++) cov[i] += cov1[i];
    for (int i=0; i<3; i++) ang[i] += ang1[i];
    KE += KE1;
    PE += PE1;
    


  }
  
  cout << endl << "Total:" << endl
       << "     Bodies:\t\t"
       << setw(15) << totbod << endl; 
  cout  << "     COM:\t\t";
  for (int i=0; i<3; i++) cout << setw(15) << com[i]/mass;
  cout << endl;
  cout  << "     COV:\t\t";
  for (int i=0; i<3; i++) cout << setw(15) << cov[i]/mass;
  cout << endl;
  cout  << "     Ang mom:\t\t";
  for (int i=0; i<3; i++) cout << setw(15) << ang[i];
  cout << endl;
  cout  << "     Stats:\t\tKE=" << KE << " PE=" << PE << " -2T/W=" << -2.0*KE/PE
	<< " Mass=" << mass << endl;
  
  return 0;
}
  
