#include <iostream>
#include <fstream>
#include <iomanip>

#include <expand.h>

#include <OutCoef.H>

OutCoef::OutCoef(const YAML::Node& conf) : Output(conf)
{
  nint = 10;
  filename = outdir + "outcoef." + runtag;
  tcomp = NULL;

  initialize();

  if (!tcomp) {
    if (myid==0) {
      cerr << "OutCoef: no component to trace\n";
      MPI_Abort(MPI_COMM_WORLD, 112);
    }
  }

  if (!(tcomp->force->HaveCoefDump())) {
    if (myid==0) {
      cerr << "OutCoef: no coefficients for this force\n";
    }
  }

}

void OutCoef::initialize()
{
  if (conf["filename"])     filename = conf["filename"].as<std::string>();
  if (conf["nint"])         nint     = conf["nint"].as<int>();
  if (conf["name"])
    {				// Search for desired component
      std::string tmp = conf["name"].as<std::string>();
      for (auto c : comp->components) {
	if (!(c->name.compare(tmp))) tcomp  = c;
      }
    }
}

void OutCoef::Run(int n, bool last)
{
  if (!(tcomp->force->HaveCoefDump())) return;
  if (n % nint != 0 && !last) return;

  MPI_Status status;

				// Open output file
  ofstream out;
  if (myid==0) {
    out.open(filename.c_str(), ios::out | ios::app);
    if (!out) {
      cout << "OutCoef: can't open file <" << filename << ">\n";
    }
  }
  
  tcomp->force->dump_coefs(out);

  out.close();
}
