#include "expand.H"

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

#include <numerical.H>
#include <orbit.H>
#include <massmodel.H>
#include <model3d.H>
#include <isothermal.H>
#include <hernquist.H>

#include <HaloBulge.H>

using namespace std;

const std::set<std::string>
HaloBulge::valid_keys = {
  "HMODEL",
  "INFILE",
  "MHALO",
  "RHALO",
  "RMODMIN",
  "RMOD",
  "RBCORE",
  "MBULGE",
  "RBULGE",
  "RBMODMIN",
  "RBMOD"
};


HaloBulge::HaloBulge(const YAML::Node& conf) : ExternalForce(conf)
{
				// Defaults
  HMODEL   = file;
  INFILE   = "w05";
  
  MHALO    = 1.0;
  RHALO    = 1.0;
  RMODMIN  = 1.0e-3;
  RMOD     = 20.0;
  
  RBCORE   = 1.0;
  MBULGE   = 1.0;
  RBULGE   = 1.0;
  RBMODMIN = 1.0e-3;
  RBMOD    = 20.0;

  initialize();

  //
  // Initialize model and orbit
  //

  switch (HMODEL) {
  case file:
    model = std::make_shared<SphericalModelTable>(INFILE);
    break;
  case isothermal:
    model = std::make_shared<IsothermalSphere>(1.0, RMODMIN, RMOD);
    break;
  case hernquist_model:
    model = std::make_shared<HernquistSphere>(1.0, RMODMIN, RMOD);
    break; 
  default:
    {
      std::ostringstream sout;
      sout << "No such HALO model type: " << (int)HMODEL << endl;
      throw GenericError(sout.str(), __FILE__, __LINE__, 1023, false);
    }
  }

  bmodel = std::make_shared<HernquistSphere>(RBCORE, RBMODMIN, RBMOD);
}


void * HaloBulge::determine_acceleration_and_potential_thread(void * arg)
{
  double r, potl, dpot, potlB, dpotB;

  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  PartMapItr it = cC->Particles().begin();
  unsigned long i;

  for (int q=nbeg; q<nend; q++) {

    i = it->first; it++;

    r = 0.0;
    for (int k=0; k<3; k++) r += cC->Pos(i, k)*cC->Pos(i, k);
    r = sqrt(r);

    model->get_pot_dpot(r/RHALO, potl, dpot);
    potl *= MHALO/RHALO;
    dpot *= MHALO/RHALO/RHALO;

    bmodel->get_pot_dpot(r/RBULGE, potlB, dpotB);
    potlB *= MBULGE/RBULGE;
    dpotB *= MBULGE/RBULGE/RBULGE;

    for (int k=0; k<3; k++) 
      cC->AddAcc(i, k, -(dpot + dpotB)*cC->Pos(i, k)/r );
    cC->AddPotExt(i,  potl + potlB );
    
  }

  return (NULL);
}


void HaloBulge::initialize()
{
  // Remove matched keys
  //
  for (auto v : valid_keys) current_keys.erase(v);
  
  // Assign values from YAML
  //
  try {
    if (conf["HMODEL"])     HMODEL   = conf["HMODEL"].as<int>();
    if (conf["INFILE"])     INFILE   = conf["INFILE"].as<std::string>();
    if (conf["MHALO"])      MHALO    = conf["MHALO"].as<double>();
    if (conf["RHALO"])      RHALO    = conf["RHALO"].as<double>();
    if (conf["RMODMIN"])    RMODMIN  = conf["RMODMIN"].as<double>();
    if (conf["RMOD"])       RMOD     = conf["RMOD"].as<double>();
    if (conf["RBCORE"])     RBCORE   = conf["RBCORE"].as<double>();
    if (conf["MBULGE"])     MBULGE   = conf["MBULGE"].as<double>();
    if (conf["RBULGE"])     RBULGE   = conf["RBULGE"].as<double>();
    if (conf["RBMODMIN"])   RBMODMIN = conf["RBMODMIN"].as<double>();
    if (conf["RBMOD"])      RBMOD    = conf["RBMOD"].as<double>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in HaloBulge: "
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
