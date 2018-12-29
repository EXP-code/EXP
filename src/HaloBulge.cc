using namespace std;

#include "expand.h"

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <Vector.h>

#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

#include <numerical.h>
#include <orbit.h>
#include <massmodel.h>
#include <model3d.h>
#include <isothermal.h>
#include <hernquist.h>

#include <HaloBulge.H>

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
    model = new SphericalModelTable(INFILE); // Halo model
    break;
  case isothermal:
    model = new IsothermalSphere(1.0, RMODMIN, RMOD); // Halo model
    break;
  case hernquist_model:
    model = new HernquistSphere(1.0, RMODMIN, RMOD); // Halo model
    break; 
  default:
    {
      std::ostringstream sout;
      sout << "No such HALO model type: " << (int)HMODEL << endl;
      throw GenericError(sout.str(), __FILE__, __LINE__);
    }
  }

  bmodel = new HernquistSphere(RBCORE, RBMODMIN, RBMOD);
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
}


void HaloBulge::initialize()
{
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
			   << error.what() << std::endl;
    MPI_Finalize();
    exit(-1);
  }
}
