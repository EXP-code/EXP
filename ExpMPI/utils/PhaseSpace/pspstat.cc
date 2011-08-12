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
#include <list>

#include <StringTok.H>
#include <header.H>
#include <PSP.H>

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
  cerr << "    -t time         use dump closest to <time>\n";
  cerr << "    -o name         prefix name for each component (default: comp)\n";
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
  string cname("comp");

  // Parse command line

  while (1) {

    int c = getopt(argc, argv, "t:o:vh");

    if (c == -1) break;

    switch (c) {

    case 't':
      time = atof(optarg);
      break;

    case 'v':
      verbose = true;
      break;

    case 'o':
      cname.erase();
      cname = string(optarg);
      break;

    case '?':
    case 'h':
    default:
      Usage(prog);
    }

  }

  ifstream *in;

  if (optind < argc) {

    ifstream *in2 = new ifstream(argv[optind]);
    if (!*in2) {
      cerr << "Error opening file <" << argv[optind] << "> for input\n";
      exit(-1);
    }

    if (verbose) cerr << "Using filename: " << argv[optind] << endl;

				// Assign file stream to input stream
    in = in2;

  }


				// Parse the PSP file
				// ------------------
  PSPDump psp(in);

  in->close();
  delete in;

				// Now write a summary
				// -------------------
  if (verbose) {

    psp.PrintSummary(in, cerr);
    
    cerr << "\nBest fit dump to <" << time << "> has time <" 
	 << psp.SetTime(time) << ">\n";
  } else 
    psp.SetTime(time);



				// Setup stats for all components
				// ------------------------------

  double com[3] = {0.0, 0.0, 0.0};
  double cov[3] = {0.0, 0.0, 0.0};
  double ang[3] = {0.0, 0.0, 0.0};
  double KE     = 0.0;
  double PE     = 0.0;
  double mass   = 0.0;
  int   totbod  = 0;
  
  list<PSPstanza>::iterator its;
  double rtmp;
  int itmp;
  
  in = new ifstream(argv[optind]);

  for (its = psp.CurrentDump()->stanzas.begin(); 
       its != psp.CurrentDump()->stanzas.end(); its++) {


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
    oname << cname << "." << its->name << '\0';
    ofstream out(oname.str().c_str());
    out.setf(ios::scientific);
    out.precision(10);

    if (!out) {
      cerr << "Couldn't open output name <" << oname.str() << ">\n";
      exit(-1);
    }

				// Print the header

    cout << "Comp name: " << its->name << endl
	 << "     Bodies:\t\t"
	 << setw(15) << its->nbod 
	 << setw(10) << its->niatr 
	 << setw(10) << its->ndatr 
	 << endl;

    totbod += its->nbod;

				// Position to beginning of particles
    in->seekg(its->pspos);

    for (int j=0; j<its->nbod; j++) {
      in->read((char *)&ms, sizeof(double));
      for (int i=0; i<3; i++) in->read((char *)&pos[i], sizeof(double));
      for (int i=0; i<3; i++) in->read((char *)&vel[i], sizeof(double));
      in->read((char *)&pot, sizeof(double));
				// Dump these . . . 
      for (int i=0; i<its->niatr; i++) in->read((char *)&itmp, sizeof(int));
      for (int i=0; i<its->ndatr; i++) in->read((char *)&rtmp, sizeof(double));

      mom[0] = pos[1]*vel[2] - pos[2]*vel[1];
      mom[1] = pos[2]*vel[0] - pos[0]*vel[2];
      mom[2] = pos[0]*vel[1] - pos[1]*vel[0];

				// Accumulate statistics
      mass1 += ms;
      for (int i=0; i<3; i++) com1[i] += ms*pos[i];
      for (int i=0; i<3; i++) cov1[i] += ms*vel[i];
      for (int i=0; i<3; i++) ang1[i] += ms*mom[i];
      rtmp = 0.0;
      for (int i=0; i<3; i++) rtmp += vel[i]*vel[i];
      KE1 += 0.5*ms*rtmp;
      PE1 += 0.5*ms*pot;
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
  
