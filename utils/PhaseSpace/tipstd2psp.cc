/*
  Translate from standard tipsy components into PSP

  MDWeinberg 10/17/07
*/

using namespace std;

#include <cstdlib>
#include <cstring>

extern "C" {
#include "tipsydefs.h"
}

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>

#include <libvars.H>		// For library globals
#include <header.H>

#include <getopt.h>		// C-style option parsing

//-------------
// help message
//-------------

void Usage(char* prog) {
  cerr << prog << " < tipsy.file > psp.file\n\n";
  exit(0);
}


int
main(int argc, char **argv)
{
  char *prog = argv[0];
  double time = 0.0;
  bool verbose = false;

  // Parse command line

  while (1) {

    int c = getopt(argc, argv, "h");

    if (c == -1) break;

    switch (c) {

    case '?':
    case 'h':
    default:
      Usage(prog);
    }

  }

  int N = xdr_init();
  if (N==0) {
    cerr << "No components found!" << endl;
    exit(-1);
  }

  vector<ifstream*> in;
  vector<ComponentHeader> headers;

  int ntmp, ntot=0;
  ComponentHeader header_psp;

  if (header.nsph) {
    xdr_gas();

    header_psp.nbod = header.nsph;
    header_psp.niatr = 0;
    header_psp.ndatr = 4;

    ntot += header.nsph;

    string cname("gas");
    string cparam("");
    string fname("");
    string fparam("");

    ostringstream outs;
    outs << cname << " : " << fname << " : " << cparam << " : " << fparam << '\0';
    strncpy(header_psp.info.get(), outs.str().c_str(), header_psp.ninfochar);

    headers.push_back(header_psp);
  }

  if (header.ndark) {
    xdr_dark();

    header_psp.nbod = header.ndark;
    header_psp.niatr = 0;
    header_psp.ndatr = 1;

    ntot += header.ndark;

    string cname("dark halo");
    string cparam("");
    string fname("");
    string fparam("");

    ostringstream outs;
    outs << cname << " : " << fname << " : " << cparam << " : " << fparam << '\0';
    strncpy(header_psp.info.get(), outs.str().c_str(), header_psp.ninfochar);

    headers.push_back(header_psp);
  }

  if (header.nstar) {
    xdr_star();

    header_psp.nbod = header.nstar;
    header_psp.niatr = 0;
    header_psp.ndatr = 3;

    ntot += header.nstar;

    string cname("star");
    string cparam("");
    string fname("");
    string fparam("");

    ostringstream outs;
    outs << cname << " : " << fname << " : " << cparam << " : " << fparam << '\0';
    strncpy(header_psp.info.get(), outs.str().c_str(), header_psp.ninfochar);

    headers.push_back(header_psp);
  }

  MasterHeader master;
  master.time = header.time;
  master.ntot = ntot;
  master.ncomp = N;

  // Write master header
  cout.write((char *)&master, sizeof(MasterHeader));
  
  double mass, pos[3], vel[3], phi;

  int icomp = 0;

  if (header.nsph) {
    vector<int>     ivec(max<int>(1, headers[icomp].niatr));
    vector<double>  dvec(max<int>(1, headers[icomp].ndatr));

    headers[icomp].write(&cout);

    for (int k=0; k<headers[icomp].nbod; k++) {

      mass = gas_particles[k].mass;
      for (int j=0; j<3; j++) {
	pos[j] = gas_particles[k].pos[j];
	vel[j] = gas_particles[k].vel[j];
      }

      dvec[0] = gas_particles[k].rho;
      dvec[1] = gas_particles[k].temp;
      dvec[2] = gas_particles[k].hsmooth;
      dvec[3] = gas_particles[k].metals;
      phi     = gas_particles[k].phi;

      // Write phase space
      cout.write((char *)&mass, sizeof(double));
      for (int j=0; j<3; j++) cout.write((char *)&(pos[j]), sizeof(double));
      for (int j=0; j<3; j++) cout.write((char *)&(vel[j]), sizeof(double));
      cout.write((char *)&phi, sizeof(double));
      for (int j=0; j<headers[icomp].niatr; j++) 
	cout.write((char *)&(ivec[j]), sizeof(int));
      for (int j=0; j<headers[icomp].ndatr; j++) 
	cout.write((char *)&(dvec[j]), sizeof(double));
    }
    
    icomp++;
  }

  if (header.ndark) {
    vector<int>     ivec(max<int>(1, headers[icomp].niatr));
    vector<double>  dvec(max<int>(1, headers[icomp].ndatr));

    headers[icomp].write(&cout);

    for (int k=0; k<headers[icomp].nbod; k++) {

      mass = dark_particles[k].mass;
      for (int j=0; j<3; j++) {
	pos[j] = dark_particles[k].pos[j];
	vel[j] = dark_particles[k].vel[j];
      }

      dvec[0] = dark_particles[k].eps;
      phi     = dark_particles[k].phi;
      
      // Write phase space
      cout.write((char *)&mass, sizeof(double));
      for (int j=0; j<3; j++) cout.write((char *)&(pos[j]), sizeof(double));
      for (int j=0; j<3; j++) cout.write((char *)&(vel[j]), sizeof(double));
      cout.write((char *)&phi, sizeof(double));
      for (int j=0; j<headers[icomp].niatr; j++) 
	cout.write((char *)&(ivec[j]), sizeof(int));
      for (int j=0; j<headers[icomp].ndatr; j++) 
	cout.write((char *)&(dvec[j]), sizeof(double));
    }
    
    icomp++;
  }

  if (header.nstar) {
    vector<int>     ivec(max<int>(1, headers[icomp].niatr));
    vector<double>  dvec(max<int>(1, headers[icomp].ndatr));

    headers[icomp].write(&cout);

    for (int k=0; k<headers[icomp].nbod; k++) {

      mass = star_particles[k].mass;
      for (int j=0; j<3; j++) {
	pos[j] = star_particles[k].pos[j];
	vel[j] = star_particles[k].vel[j];
      }

      dvec[0] = star_particles[k].metals;
      dvec[1] = star_particles[k].tform;
      dvec[2] = star_particles[k].eps;
      phi     = star_particles[k].phi;

      // Write phase space
      cout.write((char *)&mass, sizeof(double));
      for (int j=0; j<3; j++) cout.write((char *)&(pos[j]), sizeof(double));
      for (int j=0; j<3; j++) cout.write((char *)&(vel[j]), sizeof(double));
      cout.write((char *)&phi, sizeof(double));
      for (int j=0; j<headers[icomp].niatr; j++) 
	cout.write((char *)&(ivec[j]), sizeof(int));
      for (int j=0; j<headers[icomp].ndatr; j++) 
	cout.write((char *)&(dvec[j]), sizeof(double));
    }
    
    icomp++;
  }

  return 0;
}
