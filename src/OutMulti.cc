#include <iostream>
#include <fstream>
#include <iomanip>

#include "expand.H"
#include "global.H"

#include "AxisymmetricBasis.H"
#include "OutMulti.H"

const std::set<std::string>
OutMulti::valid_keys = {
  "filename",
  "nint",
  "nintsub"
};

OutMulti::OutMulti(const YAML::Node& conf) : Output(conf)
{
  initialize();
}

void OutMulti::initialize()
{
  // Remove matched keys
  //
  for (auto v : valid_keys) current_keys.erase(v);
  
  // Assign values from YAML
  //
  try {
				// Get file name
    if (Output::conf["filename"])
      filename = Output::conf["filename"].as<std::string>();
    else {
      filename.erase();
      filename = outdir + "OUT.multi." + runtag;
    }

    if (Output::conf["nint"])
      nint = Output::conf["nint"].as<int>();
    else
      nint = 100;

    if (Output::conf["nintsub"]) {
      nintsub = Output::conf["nintsub"].as<int>();
      if (nintsub <= 0) nintsub = 1;
    } else
      nintsub = std::numeric_limits<int>::max();

  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in OutMulti: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    throw std::runtime_error("OutMulti::initialize: error parsing YAML");
  }

}


void OutMulti::Run(int n, int mstep, bool last)
{
  if (n % nint && !last && !dump_signal) return;
  if (restart  && n==0  && !dump_signal) return;
  if (multistep>1 && mstep % nintsub !=0 && !dump_signal) return;
  if (tnow <= prev) return;
  
  prev = tnow;			// Record current run time

  std::ofstream *out;
  
  vector<unsigned> counts(multistep+1, 0), histo(multistep+1, 0);

  for (auto c : comp->components) {

    PartMapItr it = c->Particles().begin();
    for (int q=0; q<c->Number(); q++) counts[c->Part((it++)->first)->level]++;
  }

  MPI_Reduce(&counts[0], &histo[0], multistep+1, MPI_UNSIGNED, MPI_SUM,
	     0, MPI_COMM_WORLD);

  int nOK = 0;

  if (myid==0) {
				// Open file and write master header
    ofstream out(filename.c_str(), ios::out | ios::app);

    if (!out) {
      std::cerr << "OutMulti: can't open file <" << filename.c_str() 
		<< "> . . . quitting" << std::endl;
    }
				// Dump the histogram
    for (int n=0; n<=multistep; n++)
      out << setw(18) << tnow << setw(5) << n 
	  << setw(8) << histo[n] << endl;
    out << endl;
  }

  MPI_Bcast(&nOK, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (nOK) {
    MPI_Finalize();
    exit(33);
  }

}

