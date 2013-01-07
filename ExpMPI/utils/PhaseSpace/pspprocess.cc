/*
  Process gas statistics

  MDWeinberg 01/24/10
*/

using namespace std;

#include <values.h>
#include <cstring>
#include <cstdlib>
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
#include <FileUtils.H>

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
  cerr << "    -p name         prefix name for each component (default: 'OUT.run1')\n";
  cerr << "    -c name         component name (default: 'gas disk')\n";
  cerr << "    -R float        maximum inner cylindrical radius\n";
  cerr << "    -Z float        maximum inner cylindrical height\n";
  cerr << "    -I int          initial index for subset\n";
  cerr << "    -F int          final index for subset\n";
  cerr << "    -o name         output file name (default: 'gas stat')\n";
  cerr << "    -h              print this help message\n";
  cerr << "    -v              verbose output\n\n";
  exit(0);
}


int
main(int argc, char **argv)
{
  char *prog = argv[0];
  double rmax = 0.05;
  double zmax = 0.005;
  bool verbose = false;
  unsigned long initial=0, final=MAXLONG;
  string pname("OUT.run1");
  string cname("gas disk");
  string oname("gas.stat");

  // Parse command line

  while (1) {

    int c = getopt(argc, argv, "p:c:I:F:R:Z:o:vh");

    if (c == -1) break;

    switch (c) {

    case 'I':
      initial = atoi(optarg);
      break;

    case 'F':
      final = atoi(optarg);
      break;

    case 'R':
      rmax = atof(optarg);
      break;

    case 'Z':
      zmax = atof(optarg);
      break;

    case 'p':
      pname.erase();
      pname = string(optarg);
      break;

    case 'c':
      cname.erase();
      cname = string(optarg);
      break;

    case 'o':
      oname.erase();
      oname = string(optarg);
      break;

    case 'v':
      verbose = true;
      break;

    case '?':
    case 'h':
    default:
      Usage(prog);
    }

  }

				// Open an output file
				// -------------------

  ofstream out(oname.c_str());
  out.setf(ios::scientific);
  out.precision(10);
  
  if (!out) {
    cerr << "Couldn't open output name <" << oname << ">" << endl;
    exit(-1);
  }

  const int nlabels = 19;
  const int nwid = 20;
  const char *labels[] = {"Time", 
			  "Total mass", "Inner mass", "Outer mass",
			  "L_tot(x)", "L_tot(y)", "L_tot(z)", 
			  "L_in(x)", "L_in(y)", "L_in(z)", 
			  "V_in(x)", "V_in(y)", "V_in(z)", 
			  "L_out(x)", "L_out(y)", "L_out(z)", 
			  "V_out(x)", "V_out(y)", "V_out(z)"};
  
  for (int i=0; i<nlabels; i++) {
    if (i==0) out << setfill('-') << left << setw(nwid) << "#";
    else      out << setfill('-') << left << setw(nwid) << "+";
  }
  out << setfill(' ') << endl;
  for (int i=0; i<nlabels; i++) {
    if (i==0) out << "# ";
    else      out << "+ ";
    out << setw(nwid-2) << left << labels[i];
  }
  out << endl;
  for (int i=0; i<nlabels; i++) {
    if (i==0) out << "# "; 
    else      out << "+ ";
    ostringstream sout;
    sout << "[" << i+1 << "]";
    out << setw(nwid-2) << left << sout.str();
  }
  out << endl;
  for (int i=0; i<nlabels; i++) {
    if (i==0) out << setfill('-') << left << setw(nwid) << "#";
    else      out << setfill('-') << left << setw(nwid) << "+";
  }
  out << setfill(' ') << endl;

  ifstream in;
  char fname[128];
  unsigned ipsp=0;

  while(1) {
    
    snprintf(fname, 128, "%s.%05d", pname.c_str(), ipsp++);

    if (!FileExists(fname)) {
      cerr << "Error opening file <" << fname << "> for input\n";
      break;
    }
    
    if (verbose) cerr << "Using filename: " << fname << endl;

				// Parse the PSP file
				// ------------------
    in.close();
    in.open(fname);
    PSPDump psp(&in);

				// Now write a summary
				// -------------------
    if (verbose) {

      psp.PrintSummary(&in, cerr);
    
    }
    in.close();
				// Dump ascii for each component
				// -----------------------------
    in.open(fname);
    
    list<PSPstanza>::iterator its;
    double rtmp;
    int itmp;
    
    vector<double> pos(3), vel(3), angmom(3);
    vector<double> LTOT(3, 0.0), LIN(3, 0.0), LOUT(3, 0.0), VIN(3, 0.0), VOUT(3, 0.0);
    double mass, M=0.0, MIN=0.0, MOUT=0.0, r;
    
    for (its = psp.CurrentDump()->stanzas.begin();
	 its != psp.CurrentDump()->stanzas.end(); its++) {
      
      if (cname != its->name) continue;


				// Position to beginning of particles
      in.seekg(its->pspos);
      for (int j=0; j<its->comp.nbod; j++) {

	in.read((char *)&mass, sizeof(double));

	for (int i=0; i<3; i++) {
	  in.read((char *)&pos[i], sizeof(double));
	}
	r = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);

	for (int i=0; i<3; i++)
	  in.read((char *)&vel[i], sizeof(double));

	in.read((char *)&rtmp, sizeof(double));

	for (int i=0; i<its->comp.niatr; i++) {
	  in.read((char *)&itmp, sizeof(int));
	}

	for (int i=0; i<its->comp.ndatr; i++) {
	  in.read((char *)&rtmp, sizeof(double));
	}      
	
	// Compute angular momentum for this particle
	angmom[0] = mass*(pos[1]*vel[2] - pos[2]*vel[1]);
	angmom[1] = mass*(pos[2]*vel[0] - pos[0]*vel[2]);
	angmom[2] = mass*(pos[0]*vel[1] - pos[1]*vel[0]);
	
	// Add to total
	M += mass;
	for (int i=0; i<3; i++) LTOT[i] += angmom[i];

	// Subset
	if (j>=initial && j<final) {
				// Inner disk
	  if (r < rmax && fabs(pos[2]) < zmax) {
	    MIN += mass;
	    for (int i=0; i<3; i++) VIN[i] += mass*vel[i];
	    for (int i=0; i<3; i++) LIN[i] += angmom[i];
	  } else {		// Outer disk
	    MOUT += mass;
	    for (int i=0; i<3; i++) VOUT[i] += mass*vel[i];
	    for (int i=0; i<3; i++) LOUT[i] += angmom[i];
	  }
	}
      }
    }
      
    // Print the results for this dump
    out << setw(nwid) << psp.CurrentTime() 
	<< setw(nwid) << M << setw(nwid) << MIN << setw(nwid) << MOUT;

    if (M>0.0) {
      for (int i=0; i<3; i++) out << setw(nwid) << LTOT[i]/M;
    } else {
      for (int i=0; i<3; i++) out << setw(nwid) << 0.0;
    }
    
    if (MIN>0.0) {
      for (int i=0; i<3; i++) out << setw(nwid) << LIN[i]/MIN;
      for (int i=0; i<3; i++) out << setw(nwid) << VIN[i]/MIN;
    } else {
      for (int i=0; i<3; i++) out << setw(nwid) << 0.0;
      for (int i=0; i<3; i++) out << setw(nwid) << 0.0;
    }

    if (MOUT>0.0) {
      for (int i=0; i<3; i++) out << setw(nwid) << LOUT[i]/MOUT;
      for (int i=0; i<3; i++) out << setw(nwid) << VOUT[i]/MOUT;
    } else {
      for (int i=0; i<3; i++) out << setw(nwid) << 0.0;
      for (int i=0; i<3; i++) out << setw(nwid) << 0.0;
    }

    out << endl;
  }
  
  return 0;
}
  
