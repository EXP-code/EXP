// This may look like C code, but it is really -*- C++ -*-

/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Follow unbound satellite orbit
 *
 *  Call sequence:
 *  -------------
 *
 *  Parameters:
 *  ----------
 *
 *
 *  Returns:
 *  -------
 *
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  revision MDW 01/10/93
 *  updated to use orbit classes
 *           MDW 07/15/94
 *
 *  major rewrite: incorpated in to SatelliteOrbit class 11/15/98
 *
 ****************************************************************************/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>

#include <kevin_complex.h>
#include <Vector.h>
#include <orbit.h>
#include <massmodel.h>

#include <isothermal.h>
#include <hernquist.h>
#include <model3d.h>
#include <interp.h>

#include <UnboundOrbit.H>

using namespace std;

#include <global.H>
				// External prototype for Euler matrix

Matrix return_euler_slater(double PHI, double THETA, double PSI, int BODY);

				// Default input parameters (storage)

static database_record init[] = {
  {"MODEL",	"int",		"0"},
  {"DIVERGE",	"int",		"0"},
  {"DIVEXPON",	"double",	"1.0"},
  {"RCORE",	"double",	"1.0"},
  {"E",		"double",	"0.0"},
  {"Rperi",	"double",	"0.1"},
  {"Redge",	"double",	"2.0"},
  {"deltaR",	"double",	"0.01"},
  {"RMODMIN",	"double",	"1.0e-3"},
  {"RMODMAX",	"double",	"20.0"},
  {"VROT",	"double",	"1.0"},
  {"rmin",	"double",	"-1.0"},
  {"rmax",	"double",	"-1.0"},
  {"PHIP",	"double",	"178.45"},
  {"THETA",	"double",	"114.89"},
  {"PSI",	"double",	"54.05"},
  {"INFILE",	"string",	"SLGridSph.model"},
  {"orbfile",	"bool",		"true"},
  {""		"",		""}
};

// ===================================================================
// Constructor
// ===================================================================

UnboundOrbit::UnboundOrbit(const string &conf)
{
  config = new ParamDatabase(init);
  config->parseFile(conf);

  m = 0;
  model = 0;

  switch (config->get<int>("MODEL")) {

  case file:
    m = new SphericalModelTable(config->get<string>("INFILE"), 
				config->get<int   >("DIVERGE"),
				config->get<double>("DIVEXPON"));
    model = m;
				// Assign filename to ID string
    Model3dNames[0] = config->get<string>("INFILE");
    break;

  case sing_isothermal:
    model = new SingIsothermalSphere(1.0, 
				     config->get<double>("RMODMIN"),
				     config->get<double>("RMODMAX"));
    break;

  case isothermal:
    model = new IsothermalSphere(config->get<double>("RCORE"),
				 config->get<double>("RMODMAX"),
				 config->get<double>("VROT"));
    break;

  case hernquist_model:
    model = new HernquistSphere(1.0, 
				config->get<double>("RMODMIN"), 
				config->get<double>("RMODMAX"));
    break; 

  default:
    cerr << "Illegal model: " << config->get<int>("MODEL") << '\n';
    exit(-1);
  }

  double rmin = config->get<double>("rmin");
  double rmax = config->get<double>("rmax");
  if (rmin < 0.0) rmin = model->get_min_radius();
  if (rmax < 0.0) rmax = model->get_max_radius();

// =================================
// Compute orbit
// =================================

  double Rperi = config->get<double>("Rperi");
  double E     = config->get<double>("E");
  double VTperi = sqrt(2.0*(E - model->get_pot(Rperi)));
  double J = Rperi*VTperi;

  R.push_back(Rperi);
  T.push_back(0.0);
  PHI.push_back(0.0);

  //
  // Trapezoidal increments
  //

  double rnext, rlast = Rperi;
  double tnext, tlast = 0.0;
  double phinext, philast = 0.0;
  double deltaR = config->get<double>("deltaR");

  //
  // First step
  //

  double Redge = config->get<double>("Redge");
  double denom = sqrt(2.0*(VTperi*VTperi/Rperi - model->get_dpot(Rperi)));

  rnext = rlast + deltaR;
  tnext = tlast + 2.0*sqrt(rnext - Rperi)/denom;
  phinext = philast + 2.0*sqrt(rnext - Rperi)/denom * J/(Rperi*Rperi);
  
  R.push_back(rnext);
  T.push_back(tnext);
  PHI.push_back(phinext);

  rlast = rnext;
  tlast = tnext;
  philast = phinext;

  while (R.back() < Redge) {
    rnext = rlast + deltaR;
    tnext = tlast + 0.5*(rnext - rlast)*
      (
       1.0/sqrt(2.0*(E - model->get_pot(rlast)) - J*J/(rlast*rlast)) +
       1.0/sqrt(2.0*(E - model->get_pot(rnext)) - J*J/(rnext*rnext)) 
       );
    phinext = philast + 0.5*(rnext - rlast)*
      (
       J/(rlast*rlast) /
       sqrt(2.0*(E - model->get_pot(rlast)) - J*J/(rlast*rlast)) +
       J/(rnext*rnext) /
       sqrt(2.0*(E - model->get_pot(rnext)) - J*J/(rnext*rnext)) 
       );

    rlast = rnext;
    tlast = tnext;
    philast = phinext;

    R.push_back(rnext);
    T.push_back(tnext);
    PHI.push_back(phinext);
  }
  
  ofstream out;

  if (config->get<bool>("orbfile")) {
    string orbfile = runtag + ".xyz";
    out.open(orbfile.c_str());
  }


  Three_Vector In, Out;
  double THETA   = config->get<double>("THETA")   * M_PI/180.0;
  double PSI     = config->get<double>("PSI")     * M_PI/180.0;
  double PHIP    = config->get<double>("PHIP")    * M_PI/180.0;
  Matrix Trans = return_euler_slater(PHIP, THETA, PSI, 1);
  
  for (unsigned i=R.size()-1; i>=1; i--) {
    In[1] = R[i]*cos(-PHI[i]);
    In[2] = R[i]*sin(-PHI[i]);
    In[3] = 0.0;
    Out  = Trans * In;

    Time.push_back(-T[i]);
    Xpos.push_back(Out[1]);
    Ypos.push_back(Out[2]);
    Zpos.push_back(Out[3]);
    
    if (out)
      out << setw(18) << -T[i]
	  << setw(18) << Out[1]
	  << setw(18) << Out[2]
	  << setw(18) << Out[3]
	  << endl;
  }

  for (unsigned i=0; i<R.size(); i++) {
    In[1] = R[i]*cos(PHI[i]);
    In[2] = R[i]*sin(PHI[i]);
    In[3] = 0.0;
    Out  = Trans * In;

    Time.push_back(T[i]);
    Xpos.push_back(Out[1]);
    Ypos.push_back(Out[2]);
    Zpos.push_back(Out[3]);
    
    if (out)
      out << setw(18) << T[i]
	  << setw(18) << Out[1]
	  << setw(18) << Out[2]
	  << setw(18) << Out[3]
	  << endl;
  }


  if (myid==0) {
    
    cout << "UnboundOrbit initiated with:" << endl
	 << setw(10) << "" << setw(10) << "model" 
	 << " = " << config->get<int   >("MODEL") << endl
	 << setw(10) << "" << setw(10) << "INFILE" 
	 << " = " << config->get<string>("INFILE") << endl
	 << setw(10) << "" << setw(10) << "E" 
	 << " = " << config->get<double>("E") << endl
	 << setw(10) << "" << setw(10) << "Rperi" 
	 << " = " << config->get<double>("Rperi") << endl
	 << setw(10) << "" << setw(10) << "Redge" 
	 << " = " << config->get<double>("Redge") << endl;
  }

}


// ===================================================================
// Destructior
// ===================================================================

UnboundOrbit::~UnboundOrbit(void)
{
  if (m)
    delete m;
  else 
    delete model;

  delete config;
}

Vector UnboundOrbit::get_satellite_orbit(double t)
{
  Vector ret(1, 3);
  t = max<double>(t, Time.front());
  t = min<double>(t, Time.back() );

  ret[1] = odd2(t, Time, Xpos);
  ret[2] = odd2(t, Time, Ypos);
  ret[3] = odd2(t, Time, Zpos);

  return ret;
}

void UnboundOrbit::get_satellite_orbit(double t, double *v)
{
  t = max<double>(t, Time.front());
  t = min<double>(t, Time.back() );

  v[0] = odd2(t, Time, Xpos);
  v[1] = odd2(t, Time, Ypos);
  v[2] = odd2(t, Time, Zpos);
}

