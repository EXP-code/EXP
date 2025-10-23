/*
  Translate from ascii components into PSP using supplied header info

  MDWeinberg 05/15/04, 12/12/19
*/

#include <cstdlib>
#include <cstring>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <random>

using namespace std;

#include "header.H"
#include "cxxopts.H"		// Option parsing
#include "libvars.H"		// EXP library globals

int
main(int argc, char **argv)
{
  char *prog = argv[0];
  double time = 0.0;
  bool verbose = false;
  string outfile("new.psp");

  cxxopts::Options options("ascii2psp", "Construct a PSP file from ascii input files");

  options.add_options()
    ("h,help", "print this help message")
    ("v,verbose", "print verbose output messages")
    ("o,output", "output PSP file name",
     cxxopts::value<std::string>(outfile)->default_value("new.psp"))
    ("t,time", "desired time stamp",
     cxxopts::value<double>(time)->default_value("0.0"))
    ;

  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    if (myid==0) std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("help")) {
    std::cout << options.help() << std::endl;
    exit(-1);
  }

  if (vm.count("verbose")) {
    verbose = true;
  }

  ofstream out(outfile.c_str());
  if (!out) {
    cerr << "Error opening <" << outfile << "> for output\n";
    exit(-1);
  }

  const int lenbuf = 1024;
  char buf[lenbuf];
  vector<string> names;
  int N;

  cout << "Number of components: ";
  cin.getline(buf, lenbuf);
  istringstream sin(buf);
  sin >> N;
  
  for (int i=0; i<N; i++) {
    cout << "File name[" << i << "]: ";
    cin.getline(buf, lenbuf);
    if (cin)
      names.push_back(buf);
    else {
      cerr << "Error reading filename <" << buf << ">\n";
      exit(-1);
    }
  }

  std::vector<std::shared_ptr<ifstream>> in;
  std::vector<ComponentHeader> headers;

  int ntmp, ntot=0;
  ComponentHeader header;

  for (int i=0; i<N; i++) {
    auto in2 = std::make_shared<ifstream>(names[i].c_str());
    if (!*in2) {
      cerr << "Error opening file <" << names[i] << "> for input\n";
      exit(-1);
    }
    
    in2->getline(buf, lenbuf);
    if (!*in2) {
      cerr << "Error reading header from  <" << names[i] << ">\n";
      exit(-1);
    }

    istringstream ins(buf);

    ins >> header.nbod;
    ins >> header.niatr;
    ins >> header.ndatr;

    ntot += header.nbod;

    cout << "Component name: ";
    cin.getline(buf, lenbuf);
    string cname(buf);

    cout << "Component parameters: ";
    cin.getline(buf, lenbuf);
    string cparam(buf);

    cout << "Force name: ";
    cin.getline(buf, lenbuf);
    string fname(buf);

    cout << "Force parameters: ";
    cin.getline(buf, lenbuf);
    string fparam(buf);

    ostringstream outs;
    outs << cname << " : " << fname << " : " << cparam << " : " << fparam << '\0';
    strncpy(header.info.get(), outs.str().c_str(), header.ninfochar);

    in.push_back(in2);
    headers.push_back(header);
  }

  MasterHeader master;
  master.time = time;
  master.ntot = ntot;
  master.ncomp = N;
  
  // Write master header
  out.write((char *)&master, sizeof(MasterHeader));
  
  double mass, pos[3], vel[3], pot=0.0;

  for (int i=0; i<N; i++) {

    std::vector<int>     ivec(max<int>(1, headers[i].niatr));
    std::vector<double>  dvec(max<int>(1, headers[i].ndatr));

    headers[i].write(&out);

    for (int k=0; k<headers[i].nbod; k++) {

      in[i]->getline(buf, lenbuf);
      istringstream ins(buf);

      // Read phase space from file
      ins >> mass;
      for (int j=0; j<3; j++) ins >> pos[j];
      for (int j=0; j<3; j++) ins >> vel[j];

      for (int j=0; j<headers[i].niatr; j++) ins >> ivec[j];
      for (int j=0; j<headers[i].ndatr; j++) ins >> dvec[j];
      
      // Write phase space
      out.write((char *)&mass, sizeof(double));
      for (int j=0; j<3; j++) out.write((char *)&(pos[j]), sizeof(double));
      for (int j=0; j<3; j++) out.write((char *)&(vel[j]), sizeof(double));
      out.write((char *)&pot, sizeof(double));
      for (int j=0; j<headers[i].niatr; j++) 
	out.write((char *)&(ivec[j]), sizeof(int));
      for (int j=0; j<headers[i].ndatr; j++) 
	out.write((char *)&(dvec[j]), sizeof(double));
    }
  }

  return 0;
}
  
