#include <limits>

#include "expand.H"

#include "gaussQ.H"
#include "TwoCenter.H"
#include "MixtureBasis.H"

const std::set<std::string> TwoCenter::valid_keys = {
  "nhisto",
  "basis"

};


TwoCenter::TwoCenter(Component* c0, const YAML::Node& conf) : PotAccel(c0, conf)
{
  nhisto = 0;
  inner  = vector<double>(3, 0);
  outer  = vector<double>(3, 0);

  // Get initialization info
  //
  initialize();

  // Force identity string
  //
  id = "TwoCenter (" + basis + ")";

  // Create dianostic output
  //
  if (nhisto) {
    dz = 1.0/nhisto;
    histo = vector<double>(nhisto, 0);
    ohisto = "histo_stc." + runtag;
  }
  
  // Generate two expansion grids
  //
  mix_in  = new MixtureBasis(*this, &inner, "EJ",
			     static_cast<mixFunc>(&TwoCenter::Cmixture));

  mix_out = new MixtureBasis(*this, &outer, "COM",
			     static_cast<mixFunc>(&TwoCenter::mixture));
  
  // Instantiate the force ("reflection" by hand)
  //
  if ( !basis.compare("bessel") ) {
    exp_in  = new Bessel(c0, conf, mix_in );
    exp_out = new Bessel(c0, conf, mix_out);
  }
  else if ( !basis.compare("sphereSL") ) {
    exp_in  = new Sphere(c0, conf, mix_in );
    exp_out = new Sphere(c0, conf, mix_out);
  }
  else if ( !basis.compare("cylinder") ) {
    exp_in  = new Cylinder(c0, conf, mix_in );
    exp_out = new Cylinder(c0, conf, mix_out);
  }
  else {
    ostringstream msg;
    msg << "The basis <" << id << "> cannot be used as a multicenter component";
    throw GenericError(msg.str(), __FILE__, __LINE__, 1032, false);
  }

  dof = exp_in->dof;
}


void TwoCenter::initialize()
{
  // Check for unmatched keys
  //
  auto unmatched = YamlCheck(conf, valid_keys);
  if (unmatched.size())
    throw YamlConfigError("TwoCenter", "parameter", unmatched, __FILE__, __LINE__, 1033, false);

  // Assign values from YAML
  //
  try {
    if (conf["nhisto"])         nhisto     = conf["nhisto"].as<int>();
    if (conf["basis"])          basis      = conf["basis"].as<string>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in TwoCenter: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    throw std::runtime_error("TwoCenter::initialze: error parsing YAML");
  }
}

TwoCenter::~TwoCenter(void)
{
  delete exp_in ;
  delete exp_out;
  delete mix_in ;
  delete mix_out;
}

void TwoCenter::determine_coefficients(Component *c) 
{
  // Assign the two (inner and outer) centers
  //
  for (int k=0; k<3; k++) {
    inner[k] = component->center[k];
    if (component->com_system) 
      outer[k] = component->com0[k];
    else 
      outer[k] = component->com[k];
  }

  // Compute the two expansions
  //
  exp_in ->determine_coefficients(c);
  exp_out->determine_coefficients(c);
}

void TwoCenter::determine_coefficients(void) 
{
  // Assign the two centers
  //
  for (int k=0; k<3; k++) {
    inner[k] = component->center[k];
    if (component->com_system) 
      outer[k] = component->com0[k];
    else
      outer[k] = component->com[k];
  }

  // Compute the two expansions
  //
  exp_in ->determine_coefficients();
  exp_out->determine_coefficients();
}

void TwoCenter::get_acceleration_and_potential(Component* curComp)
{
  cC = curComp;			// "Register" component
  nbodies = cC->Number();	// And retrieve number of bodies

				// Reset diagnostic distribution
  if (multistep==0 || mstep==0) reset_histo();
  
  bool use_external1 = use_external;

  for (int k=0; k<3; k++) {
    inner[k] = component->center[k];
    if (component->com_system) 
      outer[k] = component->com0[k];
    else 
      outer[k] = component->com[k];
  }

  if (use_external) exp_in -> SetExternal();
  exp_in->get_acceleration_and_potential(cC);
  exp_in->ClearExternal();

				// Reset external force flag
  use_external = use_external1;

  if (use_external) exp_out -> SetExternal();
  exp_out->get_acceleration_and_potential(cC);
  exp_out->ClearExternal();


  // Clear external potential flag
  use_external = false;

  // Write diagnostic file
  if (multistep==0 || mstep==0) write_histo();
}

void TwoCenter::accum_histo(double value)
{
  if (nhisto) {
    if (value<0.0 || value>1.0) {
      cerr << "TwoCenter::accum_histo: out of bounds value="
	   << value << endl;
    } else {
      unsigned indx = static_cast<unsigned>( floor(value/dz) );
      indx = min<unsigned>(indx, nhisto - 1);
      histo[indx] += 1.0;
    }
  }
}

void TwoCenter::reset_histo()
{
  if (nhisto) {
    if (multistep==0 || mstep==0)
      for (unsigned n=0; n<nhisto; n++) histo[n] = 0;
  }
} 

//
// GNUPLOT format output
//
void TwoCenter::write_histo()
{
  if (nhisto) {

    vector<double> histo0(nhisto);
    MPI_Reduce(&histo[0], &histo0[0], nhisto, MPI_DOUBLE, MPI_SUM, 0, 
	       MPI_COMM_WORLD);

    if (myid==0) {
      double sum = 0.0, val;
      for (unsigned n=0; n<nhisto; n++) sum += histo0[n];
      if (sum<=0.0) sum=1.0;

      ofstream out(ohisto.c_str(), ios::out | ios::app);
      out.precision(3);
      for (unsigned n=0; n<nhisto; n++)
	out << setw(16) << tnow 
	    << setw(12) << dz*(0.5+n)
	    << setw(12) << histo0[n]/sum
	    << endl;
      out << endl;
    }
  }
}
