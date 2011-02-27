#include <values.h>

#include "expand.h"

#include <gaussQ.h>
#include <TwoCenter.H>
#include <MixtureBasis.H>

TwoCenter::TwoCenter(string& line) : PotAccel(line)
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
    exp_in  = new Bessel(line, mix_in );
    exp_out = new Bessel(line, mix_out);
  }
  else if ( !basis.compare("c_brock") ) {
    exp_in  = new CBrock(line, mix_in );
    exp_out = new CBrock(line, mix_out);
  }
  else if ( !basis.compare("c_brock_disk") ) {
    exp_in  = new CBrockDisk(line, mix_in );
    exp_out = new CBrockDisk(line, mix_out);
  }
  else if ( !basis.compare("hernq") ) {
    exp_in  = new Hernquist(line, mix_in );
    exp_out = new Hernquist(line, mix_out);
  }
  else if ( !basis.compare("sphereSL") ) {
    exp_in  = new Sphere(line, mix_in );
    exp_out = new Sphere(line, mix_out);
  }
  else if ( !basis.compare("cylinder") ) {
    exp_in  = new Cylinder(line, mix_in );
    exp_out = new Cylinder(line, mix_out);
  }
  else {
    ostringstream msg;
    msg << "The basis <" << id << "> cannot be used as a multicenter component";
    bomb(msg.str());
  }

  exp_in ->RegisterComponent(component);
  exp_out->RegisterComponent(component);

  dof = exp_in->dof;
}


void TwoCenter::initialize()
{
  string val;
  if (get_value("nhisto", val))		nhisto = atoi(val.c_str());
  if (get_value("basis",  val))		basis  = val.c_str();
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
      outer[k] = component->comI[k];
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
      outer[k] = component->comI[k];
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
      outer[k] = component->comI[k];
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
