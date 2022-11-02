using namespace std;

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "expand.H"
#include <global.H>

#include <AxisymmetricBasis.H>
#include <OutAscii.H>

const std::set<std::string> OutAscii::valid_keys = {
  "nint",
  "nintsub",
  "nbeg",
  "name",
  "accel",
  "filename"
};

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
  // Remove matched keys
  for (auto v : valid_keys) current_keys.erase(v);
  
  // Assign values from YAML
  //
  try {
    if (Output::conf["nint"])    nint     = Output::conf["nint"].as<int>();
#ifdef ALLOW_NINTSUB
    if (Output::conf["nintsub"]) nintsub  = Output::conf["nintsub"].as<int>();
    if (nintsub <= 0) nintsub = 1; // Sanity check
#else
    nintsub_warning("OutAscii");
#endif
    if (Output::conf["nbeg"])    nbeg     = Output::conf["nbeg"].as<int>();
    if (Output::conf["name"])    name     = Output::conf["name"].as<std::string>();
    if (Output::conf["accel"])   accel    = Output::conf["accel"].as<bool>();

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
  if (multistep>1 and mstep % nintsub !=0) return;
  if (!c0) return;

#ifdef HAVE_LIBCUDA
  if (use_cuda) {
    if (c0->force->cudaAware() and not comp->fetched[c0]) {
      comp->fetched[c0] = true;
      c0->CudaToParticles();
    }
  }
#endif

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
    
  }

  // Check that root has a good stream
  //
  MPI_Bcast(&nOK, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (nOK) {
    MPI_Finalize();
    exit(33);
  }

  // Update total particle number
  //
  c0->NewTotal();

  // Write the header
  //
  if (nOK == 0) {
      out << "# Time=" << tnow << "\n";
      out << setw(10) << c0->CurTotal()
	  << setw(10) << c0->niattrib
	  << setw(10) << c0->ndattrib << "\n";
  }
  
  // Dump the phase-space info into the file
  //
  c0->write_ascii(&out, accel);

  // Close file and done
  //
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

