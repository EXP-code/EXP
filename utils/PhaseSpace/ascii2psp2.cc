/*
  Translate from ascii components into PSP using supplied header info.  New style.

  MDWeinberg 12/12/19
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <random>
#include <string>
#include <cstring>

#include <getopt.h>		// C-style option parsing

#include <yaml-cpp/yaml.h>	// YAML support
#include "libvars.H"		// EXP global library variables

#include "PSP.H"

//-------------
// Help message
//-------------

void Usage(char* prog) {
  cerr << prog << ": [-t time -v -h]\n\n";
  cerr << "    -t time         header time\n";
  cerr << "    -c config       YAML config file\n";
  cerr << "    -o outfile      output file\n";
  cerr << "    -I              add an index\n";
  cerr << "    -4              use float rather than double\n";
  cerr << "    -h              print this help message\n";
  cerr << "    -v              verbose output\n\n";
  exit(0);
}


int
main(int argc, char **argv)
{
  char *prog = argv[0];
  double time = 0.0;
  bool verbose = false, real4 = false, indexing = false;
  string outfile("new.psp"), config("new.config");

  // Parse command line

  while (1) {

    int c = getopt(argc, argv, "t:c:o:I4vh");

    if (c == -1) break;

    switch (c) {

    case 't':
      time = atof(optarg);
      break;

    case 'c':
      config.erase();
      config = optarg;
      break;

    case 'o':
      outfile.erase();
      outfile = optarg;
      break;

    case 'I':
      indexing = true;
      break;

    case '4':
      real4 = true;
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


  YAML::Node conf;

  std::ifstream fconf(config);

  try {
    conf = YAML::Load(fconf);
  }
  catch (YAML::Exception & error) {
    std::cout << "Error parsing component config file."
	      << std::endl
	      << error.what() << std::endl;
    exit(-1);
  }

  YAML::Node cconf  = conf["Components"];

  if (not cconf.IsSequence()) {
    std::cout << "Problem reading config file <" << config << ">"
	      << std::endl;
    exit(-1);
  }

  std::vector<std::ifstream>   in(cconf.size());
  std::vector<ComponentHeader> headers(cconf.size());

  const int lenbuf = 1024;
  char buf[lenbuf];
  int ntot = 0;

  for (int i=0; i<cconf.size(); i++) {

    auto bodies = cconf[i]["bodyfile"].as<std::string>();

    in[i].open(bodies);
    if (!in[i]) {
      std::cerr << "Error opening file <" << bodies << "> for input\n";
      exit(-1);
    }

    in[i].getline(buf, lenbuf);
    if (!in[i]) {
      cerr << "Error reading header from  <" << bodies << ">\n";
      exit(-1);
    }

    std::istringstream ins(buf);

    ins >> headers[i].nbod;
    ins >> headers[i].niatr;
    ins >> headers[i].ndatr;

    ntot += headers[i].nbod;

    std::ostringstream outs;
    outs << cconf[i] << std::endl;

    if (headers[i].ninfochar < outs.str().size()) {
      headers[i].ninfochar = outs.str().size();
      headers[i].info = std::make_shared<char>(headers[i].ninfochar+1);
    }

    // Copy to info string
    strncpy(headers[i].info.get(), outs.str().c_str(), headers[i].ninfochar);

    // DEBUGGING
    if (true) {
      std::cout << std::string(72, '-') << std::endl
		<< "Serialized YAML header looks like this:" << std::endl
		<< std::string(72, '-') << std::endl
		<< outs.str() << std::endl
		<< "Cur size=" << outs.str().size()
		<< " max size=" << headers[i].ninfochar << std::endl
		<< std::string(72, '-') << std::endl;
    }
  }

  MasterHeader master;

  master.time = time;
  master.ntot = ntot;
  master.ncomp = cconf.size();

  // Write master header
  //
  out.write((char *)&master, sizeof(MasterHeader));


  double mass, pos[3], vel[3], pot=0.0;

  for (int i=0; i<cconf.size(); i++) {

    const static unsigned long magic = 0xadbfabc0;
    unsigned long rsize;
    if (real4) rsize = sizeof(float);
    else       rsize = sizeof(double);
    unsigned long cmagic = magic + rsize;

    out.write((const char*)&cmagic, sizeof(unsigned long));

    if (not headers[i].write(&out)) {
      std::ostringstream sout;
      sout << "Error writing particle header at " << __FILE__ << ":" << __LINE__;
      throw std::runtime_error(sout.str());
    }

    std::vector<int>     ivec(max<int>(1, headers[i].niatr));
    std::vector<double>  dvec(max<int>(1, headers[i].ndatr));
    float                fv;

    for (int k=0; k<headers[i].nbod; k++) {

      in[i].getline(buf, lenbuf);
      istringstream ins(buf);

      // Read phase space from file
      ins >> mass;
      for (int j=0; j<3; j++) ins >> pos[j];
      for (int j=0; j<3; j++) ins >> vel[j];

      for (int j=0; j<headers[i].niatr; j++) ins >> ivec[j];
      for (int j=0; j<headers[i].ndatr; j++) ins >> dvec[j];

      // Write phase space
      //

      if (indexing) {      // Index?
	unsigned long indx = k + 1;
	out.write((char *)&indx, sizeof(unsigned long));
      }

      if (real4) {
	out.write((char *)&(fv=mass), sizeof(float));
	for (int j=0; j<3; j++) out.write((char *)&(fv=pos[j]), sizeof(float));
	for (int j=0; j<3; j++) out.write((char *)&(fv=vel[j]), sizeof(float));
	out.write((char *)&(fv=pot), sizeof(float));
	for (int j=0; j<headers[i].niatr; j++)
	  out.write((char *)&(ivec[j]), sizeof(int));
	for (int j=0; j<headers[i].ndatr; j++)
	  out.write((char *)&(fv=dvec[j]), sizeof(float));
      } else {
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
    // END: particle loop

    std::cout << "Wrote " << headers[i].nbod << " particles" << std::endl;

  }
  // END: Component loop

  return 0;
}
