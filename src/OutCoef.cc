#include <iostream>
#include <fstream>
#include <iomanip>

#include <expand.H>

#include <OutCoef.H>

const std::set<std::string> OutCoef::valid_keys = {
  "filename",
  "nint",
  "nintsub",
  "native",
  "name"
};

OutCoef::OutCoef(const YAML::Node& conf) : Output(conf)
{
  nint = 10;
  nintsub = std::numeric_limits<int>::max();
  filename = outdir + "outcoef." + runtag;
  native = false;
  tcomp = NULL;

  initialize();

  if (!tcomp) {
    if (myid==0) {
      std::cerr << "OutCoef: no component to trace\n";
    }
    MPI_Finalize();
    exit(112);
  }

  if (!(tcomp->force->HaveCoefDump())) {
    if (myid==0) {
      cerr << "OutCoef: no coefficients for this force\n";
    }
  }

}

void OutCoef::initialize()
{
  // Remove matched keys
  for (auto v : valid_keys) current_keys.erase(v);
  
  // Assign values from YAML
  //
  try {
    if (conf["filename"])     filename = conf["filename"].as<std::string>();
    if (conf["nint"])         nint     = conf["nint"].as<int>();
    if (conf["native"])       native   = conf["native"].as<bool>();
    if (conf["nintsub"]) {
      nintsub  = conf["nintsub"].as<int>();
      if (nintsub <= 0) nintsub = 1;
    }
    if (conf["name"])
      {				// Search for desired component
	std::string tmp = conf["name"].as<std::string>();
	for (auto c : comp->components) {
	  if (!(c->name.compare(tmp))) tcomp  = c;
	}
      }
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in OutCoef: "
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

void OutCoef::Run(int n, int mstep, bool last)
{
  // Coefficients are not available for output (pure playback mode)
  //
  if (tcomp->force->NoCoefs())         return;

  // No coefficient output defined for this force
  //
  if (!(tcomp->force->HaveCoefDump())) return;

  // Skip this master step
  //
  if (n % nint != 0 && !last) return;

  // Skip this sub step
  //
  if (multistep>1 and mstep % nintsub != 0) return;

  if (myid==0) {

    if (native) {

      // Open output file
      //
      std::ofstream out;

      // Now open the file on the root node only
      //
      out.open(filename.c_str(), ios::out | ios::app);
      if (!out) {
	cout << "OutCoef: can't open file <" << filename << ">\n";
      }
      
      // Make a bigger output buffer
      //
      const int bufsize = 16384;
      char mybuffer [bufsize];
      out.rdbuf()->pubsetbuf(mybuffer, bufsize);
      
      tcomp->force->dump_coefs(out);

    } else {
      tcomp->force->dump_coefs_h5(filename);
    }
  }

}
