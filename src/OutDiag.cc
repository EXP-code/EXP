#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "expand.H"

#include <AxisymmetricBasis.H>
#include <OutDiag.H>

const std::set<std::string> OutDiag::valid_keys = {
  "filename",
  "nint",
  "nintsub",
  "RMIN",
  "RMAX",
  "THETA",
  "PHI",
  "NUM"
};

OutDiag::OutDiag(const YAML::Node& conf) : Output(conf)
{
  if (myid) return;
				// Defaults
  RMIN = 1.0e-3;
  RMAX = 10.0;
  THETA = 0.5*M_PI;
  PHI = 1.0e-10;
  NUM = 100;

  for (auto c : comp->components) {
    if (c->force->geometry == PotAccel::sphere || 
	c->force->geometry == PotAccel::cylinder) 
      {
	lcomp.push_back(c);
      }
  }
  
  names.push_back("Rho");
  names.push_back("Pot");
  names.push_back("d(Pot)/dr)");
  names.push_back("d(Pot)/d cos(theta)");
  names.push_back("d(Pot)/d phi");

  initialize();
}

void OutDiag::initialize()
{
  // Check for unmatched keys
  //
  auto unmatched = YamlCheck(conf, valid_keys);
  if (unmatched.size())
    throw YamlConfigError("OutDiag", "parameter", unmatched, __FILE__, __LINE__);

  // Assign values from YAML
  //
  try {
				// Get file name
    if (Output::conf["filename"])
      filename = Output::conf["filename"].as<std::string>();
    else {
      filename.erase();
      filename = "ORBDIAG." + runtag;
    }

    if (Output::conf["nint"])
      nint = Output::conf["nint"].as<int>();
    else
      nint = 1;

    if (Output::conf["nintsub"]) {
      nintsub = Output::conf["nintsub"].as<int>();
      if (nintsub <= 0) nintsub = 1;
    } else
      nintsub = std::numeric_limits<int>::max();

    if (Output::conf["RMIN"])
      RMIN = Output::conf["RMIN"].as<double>();

    if (Output::conf["RMAX"])
      RMAX = Output::conf["RMAX"].as<double>();

    if (Output::conf["THETA"])
      THETA = Output::conf["THETA"].as<double>();
    
    if (Output::conf["PHI"])
      PHI = Output::conf["PHI"].as<double>();
    
    if (Output::conf["NUM"])
      NUM = Output::conf["NUM"].as<int>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in OutDiag: "
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


void OutDiag::header(ostream& out)
{
  int ncur = 0;

  out << "# " << ++ncur << ": " << "Radius\n";

  for (auto c : lcomp) {

    out << "# [" << c->id << "]\n";
    
    for (int i=0; i<5; i++) 
      out << "# " << setw(3) << ++ncur << ": " << names[i] << endl;
  }

  out << "#" << endl;
}


void OutDiag::Run(int n, int mstep, bool last)
{
  if (myid) return;
  if (n % nint && !last) return;
  if (multistep>1 and mstep % nintsub !=0) return;


  double r, dr, dens, dens0;
  double potl0, potl, potr, pott, potp;
    
  ostringstream outs;
  outs << outdir << filename.c_str() << "." << n;
  ofstream out(outs.str().c_str());
  if (!out) return;

  out.setf(ios::scientific);
  out.precision(6);

  header(out);

  // Determine potential and acceleration

  dr = (RMAX - RMIN)/(double)NUM;
  for (int i=0; i<=NUM; i++) {
    r = RMIN + dr*i;

    out << setw(15) << r;
    
    for (auto c : lcomp) {
      AxisymmetricBasis * q = (AxisymmetricBasis *)(c->force);
      q->determine_fields_at_point_sph(r, THETA, PHI,
				       &dens0, &potl0, 
				       &dens, &potl, &potr, &pott, &potp);
      out << setw(15) << dens
	  << setw(15) << potl
	  << setw(15) << potr
	  << setw(15) << pott
	  << setw(15) << potp;
    }
    out << endl;
  }

}

