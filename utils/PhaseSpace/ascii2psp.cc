/*
  Translate from ascii components into PSP using supplied header info

  MDWeinberg 05/15/04, 12/12/19
*/

using namespace std;

#include <cstdlib>
#include <cstring>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>

#include <header.H>
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
  cerr << prog << ": [-t time -v -h]\n\n";
  cerr << "    -t time         header time\n";
  cerr << "    -o outfile      output file\n";
  cerr << "    -h              print this help message\n";
  cerr << "    -v              verbose output\n\n";
  exit(0);
}


int
main(int argc, char **argv)
{
  char *prog = argv[0];
  double time = 0.0;
  bool verbose = false;
  string outfile("new.psp");

  // Parse command line

  while (1) {

    int c = getopt(argc, argv, "t:o:vh");

    if (c == -1) break;

    switch (c) {

    case 't':
      time = atof(optarg);
      break;

    case 'o':
      outfile.erase();
      outfile = optarg;
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

  vector<ifstream*> in;
  vector<ComponentHeader> headers;

  int ntmp, ntot=0;
  ComponentHeader header;

  for (int i=0; i<N; i++) {
    ifstream *in2 = new ifstream(names[i].c_str());
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

    vector<int>     ivec(max<int>(1, headers[i].niatr));
    vector<double>  dvec(max<int>(1, headers[i].ndatr));

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
  
