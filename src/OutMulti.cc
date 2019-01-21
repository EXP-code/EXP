#include <iostream>
#include <fstream>
#include <iomanip>

#include "expand.h"
#include <global.H>

#include <AxisymmetricBasis.H>
#include <OutMulti.H>

OutMulti::OutMulti(const YAML::Node& conf) : Output(conf)
{
  initialize();
}

void OutMulti::initialize()
{
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
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in OutMulti: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }

}


void OutMulti::Run(int n, bool last)
{
  if (n % nint && !last && !dump_signal) return;
  if (restart  && n==0  && !dump_signal) return;

  ofstream *out;
  
  vector<unsigned> counts(multistep+1, 0), histo(multistep+1, 0);

  for (auto c : comp->components) {

    PartMapItr it = c->Particles().begin();
    for (int q=0; q<c->Number(); q++) counts[c->Part((it++)->first)->level]++;
  }

  MPI_Reduce(&counts[0], &histo[0], multistep+1, MPI_UNSIGNED, MPI_SUM,
	     0, MPI_COMM_WORLD);

  if (myid==0) {
				// Open file and write master header
    ofstream out(filename.c_str(), ios::out | ios::app);

    if (!out) {
      cerr << "OutMulti: can't open file <" << filename.c_str() 
	   << "> . . . quitting\n";
      MPI_Abort(MPI_COMM_WORLD, 33);
    }

				// Dump the histogram
    for (int n=0; n<=multistep; n++)
      out << setw(18) << tnow << setw(5) << n 
	  << setw(8) << histo[n] << endl;
    out << endl;
  }

  dump_signal = 0;
}

