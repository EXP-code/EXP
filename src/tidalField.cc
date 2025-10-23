
#include "expand.H"

#include "tidalField.H"

const std::set<std::string>
tidalField::valid_keys = {
  "hills_omega",
  "hills_p"
};

tidalField::tidalField(const YAML::Node& config) : ExternalForce(config)
{
  hills_omega = 0.5;
  hills_p = 0.5;
  
  initialize();
}

void tidalField::initialize()
{
  // Remove matched keys
  //
  for (auto v : valid_keys) current_keys.erase(v);
  
  // Assign values from YAML
  //
  try {
    if (conf["hills_omega"])    hills_omega        = conf["hills_omega"].as<double>();
    if (conf["hills_p"])        hills_p            = conf["hills_p"].as<double>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in tidalField: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    throw std::runtime_error("tidalField::initialize: error parsing YAML");
  }
}

void * tidalField::determine_acceleration_and_potential_thread(void * arg)
{
  int i;
  double s, c, w2, pp, pm, x, y, z;

  w2 = hills_omega*hills_omega;
  pm = 1.0-hills_p;
  pp = 1.0+hills_p;
	
  c = cos(2.0*hills_omega*tnow);
  s = sin(2.0*hills_omega*tnow);

  int nbeg, nend, id = *((int*)arg);

  for (auto cp : comp->components) {

    nbodies = cp->Number();
    nbeg = nbodies*id/nthrds;
    nend = nbodies*(id+1)/nthrds;

    PartMapItr it = cp->Particles().begin();
    unsigned long i;

    for (int q=0   ; q<nbeg; q++) it++;
    for (int q=nbeg; q<nend; q++) 
      {
	i = (it++)->first;

	if (cp->freeze(i)) continue;
	x = cp->Pos(i, 0);
	y = cp->Pos(i, 1);
	z = cp->Pos(i, 2);
	cp->AddAcc(i, 0, 0.5*w2*(pp*(c*x + s*y) - pm*x) );
	cp->AddAcc(i, 1, 0.5*w2*(pp*(s*x - c*y) - pm*y) );
	cp->AddAcc(i, 2, w2*z );
	cp->AddPotExt(i, 0.5*w2*z*z - 
		      0.25*w2*(pp*(c+s)*x*x + pp*(s-c)*y*y - pm*(x*x+y*y) ) );
      }
  }

  return (NULL);
}
