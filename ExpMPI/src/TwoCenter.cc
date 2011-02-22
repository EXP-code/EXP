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
  initialize();

				// Force identity string
  id = "TwoCenter (" + basis + ")";

				// Create dianostic output
  if (nhisto) {
    dz = 1.0/nhisto;
    histo = vector<double>(nhisto, 0);
    ohisto = "histo_stc." + runtag;
  }
				// Generate two expansion grids

  mix_ej  = new MixtureBasis(*this, &inner, "EJ",
			     static_cast<mixFunc>(&TwoCenter::Cmixture));

  mix_com = new MixtureBasis(*this, &outer, "COM",
			     static_cast<mixFunc>(&TwoCenter::mixture));
  
				// Instantiate the force ("reflection" by hand)
  if ( !basis.compare("bessel") ) {
    exp_ej  = new Bessel(line, mix_ej );
    exp_com = new Bessel(line, mix_com);
  }
  else if ( !basis.compare("c_brock") ) {
    exp_ej  = new CBrock(line, mix_ej );
    exp_com = new CBrock(line, mix_com);
  }
  else if ( !basis.compare("c_brock_disk") ) {
    exp_ej  = new CBrockDisk(line, mix_ej );
    exp_com = new CBrockDisk(line, mix_com);
  }
  else if ( !basis.compare("hernq") ) {
    exp_ej  = new Hernquist(line, mix_ej );
    exp_com = new Hernquist(line, mix_com);
  }
  else if ( !basis.compare("sphereSL") ) {
    exp_ej  = new Sphere(line, mix_ej );
    exp_com = new Sphere(line, mix_com);
  }
  else if ( !basis.compare("cylinder") ) {
    exp_ej  = new Cylinder(line, mix_ej );
    exp_com = new Cylinder(line, mix_com);
  }
  else {
    ostringstream msg;
    msg << "The basis <" << id << "> cannot be used as a multicenter component";
    bomb(msg.str());
  }

  exp_ej ->RegisterComponent(component);
  exp_com->RegisterComponent(component);

  dof = exp_ej->dof;
}


void TwoCenter::initialize()
{
  string val;
  if (get_value("nhisto", val))		nhisto = atoi(val.c_str());
  if (get_value("basis",  val))		basis  = val.c_str();
}

TwoCenter::~TwoCenter(void)
{
  delete exp_ej;
  delete exp_com;
  delete mix_ej;
  delete mix_com;
}

void TwoCenter::determine_coefficients(Component *c) 
{
  for (int k=0; k<3; k++) {
    inner[k] = component->center[k];
    if (component->com_system) 
      outer[k] = component->comI[k];
    else 
      outer[k] = component->com[k];
  }

  exp_ej->determine_coefficients(c);
  
  exp_com->determine_coefficients(c);
}

void TwoCenter::determine_coefficients(void) 
{
  for (int k=0; k<3; k++) {
    inner[k] = component->center[k];
    if (component->com_system) 
      outer[k] = component->comI[k];
    else
      outer[k] = component->com[k];
  }

  exp_ej->determine_coefficients();

  exp_com->determine_coefficients();
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

  if (use_external) exp_ej -> SetExternal();
  exp_ej->get_acceleration_and_potential(cC);
  exp_ej->ClearExternal();

				// Reset external force flag
  use_external = use_external1;

  if (use_external) exp_com -> SetExternal();
  exp_com->get_acceleration_and_potential(cC);
  exp_com->ClearExternal();


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
