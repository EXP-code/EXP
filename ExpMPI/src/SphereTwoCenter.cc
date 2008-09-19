#include <values.h>

#include "expand.h"

#include <gaussQ.h>
#include <SphereTwoCenter.H>
#include <MixtureSL.H>

SphereTwoCenter::SphereTwoCenter(string& line) : PotAccel(line)
{
  id = "SphereTwoCenter SL";
  nhisto = 0;
  center  = vector<double>(3, 0);

				// Get initialization info
  initialize();
				// Create dianostic output
  if (nhisto) {
    dz = 1.0/nhisto;
    histo = vector<double>(nhisto, 0);
    ohisto = "histo_stc." + runtag;
  }
				// Generate two expansion grids
  mix_ej  = new MixtureSL(*this, &center, "EJ",
			  static_cast<mixFunc>(&SphereTwoCenter::Cmixture));

  exp_ej  = new Sphere(line, mix_ej);
  exp_ej->RegisterComponent(component);

  mix_com = new MixtureSL(*this, &center, "COM",
			  static_cast<mixFunc>(&SphereTwoCenter::mixture));

  exp_com = new Sphere(line, mix_com);
  exp_com->RegisterComponent(component);
}


void SphereTwoCenter::initialize()
{
  string val;
  if (get_value("nhisto", val))		nhisto = atoi(val.c_str());
}

SphereTwoCenter::~SphereTwoCenter(void)
{
  delete exp_ej;
  delete exp_com;
  delete mix_ej;
  delete mix_com;
}

void SphereTwoCenter::determine_coefficients(Component *c) 
{
  for (int k=0; k<3; k++) center[k] = component->center[k];

  exp_ej->determine_coefficients(c);
  
  if (component->com_system)
    for (int k=0; k<3; k++) center[k] = 0.0;
  else
    for (int k=0; k<3; k++) center[k] = component->com[k];

  exp_com->determine_coefficients(c);
}

void SphereTwoCenter::determine_coefficients(void) 
{
  for (int k=0; k<3; k++) center[k] = component->center[k];

  exp_ej->determine_coefficients();

  if (component->com_system)
    for (int k=0; k<3; k++) center[k] = 0.0;
  else
    for (int k=0; k<3; k++) center[k] = component->com[k];

  exp_com->determine_coefficients();
}

void SphereTwoCenter::get_acceleration_and_potential(Component* curComp)
{
  cC = curComp;			// "Register" component
  nbodies = cC->Number();	// And retrieve number of bodies

				// Reset diagnostic distribution
  if (multistep==0 || mstep==0) reset_histo();
  
  bool use_external1 = use_external;

				// Set center to Component center
  for (int k=0; k<3; k++) center[k] = component->center[k];
  if (use_external) exp_ej -> SetExternal();
  exp_ej->get_acceleration_and_potential(cC);
  exp_ej->ClearExternal();

				// Reset external force flag
  use_external = use_external1;

				// Set center to Component center of mass
  if (component->com_system)
    for (int k=0; k<3; k++) center[k] = 0.0;
  else
    for (int k=0; k<3; k++) center[k] = component->com[k];

  if (use_external) exp_com -> SetExternal();
  exp_com->get_acceleration_and_potential(cC);
  exp_com->ClearExternal();


  // Clear external potential flag
  use_external = false;

  // Write diagnostic file
  if (multistep==0 || mstep==0) write_histo();
}

void SphereTwoCenter::accum_histo(double value)
{
  if (nhisto) {
    if (value<0.0 || value>1.0) {
      cerr << "SphereTwoCenter::accum_histo: out of bounds value="
	   << value << endl;
    } else {
      unsigned indx = static_cast<unsigned>( floor(value/dz) );
      histo[indx] += 1.0;
    }
  }
}

void SphereTwoCenter::reset_histo()
{
  if (nhisto) {
    if (multistep==0 || mstep==0)
      for (unsigned n=0; n<nhisto; n++) histo[n] = 0;
  }
} 

//
// GNUPLOT format output
//
void SphereTwoCenter::write_histo()
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
