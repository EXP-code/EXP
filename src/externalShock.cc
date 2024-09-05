
#include "expand.H"

#include <stdlib.h>
#include <string>
#include <numerical.H>

#include <externalShock.H>

externalShock::externalShock(const YAML::Node& conf) : ExternalForce(conf)
{
  E      = -0.5;
  K      = 1.0e-4;
  PER    = 0.25;
  AMPL   = 1.0;
  INFILE = "w05";

  initialize();

  model = std::make_shared<SphericalModelTable>(INFILE);
  t = std::make_shared<SphericalOrbit>(model, E, K);
}


externalShock::~externalShock()
{
  // None
}

void externalShock::initialize()
{
  try {
    if (conf["E"])      E      = conf["E"].as<double>();
    if (conf["K"])      K      = conf["K"].as<double>();
    if (conf["PER"])    PER    = conf["PER"].as<double>();
    if (conf["AMPL"])   AMPL   = conf["AMPL"].as<double>();
    if (conf["INFILE"]) INFILE = conf["INFILE"].as<std::string>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in externalShock: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    throw std::runtime_error("externalShock::initialize: error in parsing YAML");
  }
}


void * externalShock::determine_acceleration_and_potential_thread(void * arg)
{
  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  double w2, x, z;

  w2 = get_tidal_shock(tnow);

  PartMapItr it = cC->Particles().begin();
  unsigned long i;


  for (int q=0   ; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++)
    {
      i = (it++)->first;

      x = cC->Pos(i, 0);
      z = cC->Pos(i, 2);

      if (component->freeze(i)) continue;

      cC->AddAcc(i, 2, -w2*x );
      cC->AddPotExt(i, 0.5*w2*z*z );
    }

  return (NULL);
}


double externalShock::get_tidal_shock(double T)
{
  return AMPL * model->get_dpot2(t->get_angle(6, T*PER));
}


