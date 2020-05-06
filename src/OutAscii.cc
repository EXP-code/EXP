using namespace std;

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "expand.h"
#include <global.H>

#include <AxisymmetricBasis.H>
#include <OutAscii.H>

OutAscii::OutAscii(const YAML::Node& conf) : Output(conf)
{
  nint = 100;
  nintsub = std::numeric_limits<int>::max();
  nbeg = 0;
  name = "";
  accel = false;
  initialize();

  if (name.size()>0) {
				// Look for the fiducial component
    bool found = false;
    for (auto c : comp->components) {
      if ( !name.compare(c->name) ) {
	c0 = c;
	found = true;
      break;
      }
    }

    if (!found) {
      if (myid==0) 
	std::cerr << "Process " << myid << ": can't find desired component <"
		  << name << ">" << std::endl;
      MPI_Finalize();
      exit(35);
    }
  }
  else
    c0 = NULL;
}

void OutAscii::initialize()
{
  try {
    if (Output::conf["nint"])    nint  = Output::conf["nint"].as<int>();
    if (Output::conf["nintsub"])    nintsub  = Output::conf["nintsub"].as<int>();
    if (Output::conf["nbeg"])    nbeg  = Output::conf["nbeg"].as<int>();
    if (Output::conf["name"])    name  = Output::conf["name"].as<std::string>();
    if (Output::conf["accel"])   accel = Output::conf["name"].as<bool>();
    
    if (Output::conf["filename"])
      filename = Output::conf["filename"].as<std::string>();
    else {
      filename.erase();
      filename = outdir + "OUTASC." + runtag;
    }
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in OutAscii: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }
				// Determine last file

  if (restart && nbeg==0 && myid==0) {

    for (nbeg=0; nbeg<100000; nbeg++) {

				// Output name
      ostringstream fname;
      fname << filename << "." << setw(5) << setfill('0') << nbeg;

				// See if we can open file
      ifstream in(fname.str().c_str());

      if (!in) {
	cout << "OutAscii: will begin with nbeg=" << nbeg << endl;
	break;
      }
    }
  }
}


void OutAscii::Run(int n, int mstep, bool last)
{
  if (n % nint && !last) return;
  if (mstep % nintsub !=0) return;
  if (!c0) return;

  std::ofstream out;

  int nOK = 0;

  if (myid==0) {
				// Output name
    std::ostringstream fname;
    fname << filename << "." << setw(5) << setfill('0') << nbeg++;

				// Open file and write master header
    out.open(fname.str());

    if (out.fail()) {
      std::cerr << "OutAscii: can't open file <" << fname.str() 
		<< "> . . . quitting\n";
      nOK = 1;
    }
    
    if (nOK == 0) {
      out << "# Time=" << tnow << "\n";
      out << setw(10) << c0->NewTotal()
	  << setw(10) << c0->niattrib
	  << setw(10) << c0->ndattrib << "\n";
    }
  }

  MPI_Bcast(&nOK, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (nOK) {
    MPI_Finalize();
    exit(33);
  }

  c0->write_ascii(&out, accel);

  if (myid==0) {
    try {
      out.close();
    }
    catch (const ofstream::failure& e) {
      std::cout << "OutAscii: exception closing file: "
		<< e.what() << std::endl;
    }
  }

}

