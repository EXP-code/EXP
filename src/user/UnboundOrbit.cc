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

#include <orbit.H>
#include <massmodel.H>

#include <isothermal.H>
#include <hernquist_model.H>
#include <model3d.H>
#include <interp.H>
#include <YamlCheck.H>
#include <EXPException.H>

#include <UnboundOrbit.H>

using namespace std;

#include <global.H>
				// External prototype for Euler matrix

Eigen::Matrix3d
return_euler_slater(double PHI, double THETA, double PSI, int BODY);

const std::set<std::string>
UnboundOrbit::valid_keys = {
  "MODEL",
  "DIVERGE",
  "DIVEXPON",
  "RCORE",
  "E",
  "Rperi",
  "Redge",
  "deltaR",
  "RMODMIN",
  "RMODMAX",
  "VROT",
  "rmin",
  "rmax",
  "PHIP",
  "THETA",
  "PSI",
  "INFILE",
  "orbfile"
};

// ===================================================================
// Constructor
// ===================================================================

UnboundOrbit::UnboundOrbit(const YAML::Node& conf)
{
  // Default parameters
  //
  int     MODEL          = 0;
  int     DIVERGE        = 0;
  double  DIVEXPON       = 1.0;
  double  RCORE          = 1.0;
  double  E              = 0.0;
  double  Rperi          = 0.1;
  double  Redge          = 2.0;
  double  deltaR         = 0.01;
  double  RMODMIN        = 1.0e-3;
  double  RMODMAX        = 20.0;
  double  VROT           = 1.0;
  double  rmin           = -1.0;
  double  rmax           = -1.0;
  double  PHIP           = 178.45;
  double  THETA          = 114.89;
  double  PSI            = 54.05;
  string  INFILE         = "SLGridSph.model";
  bool    orbfile        = true;

  // Check for unmatched keys
  //
  auto unmatched = YamlCheck(conf, valid_keys);
  if (unmatched.size())
    throw YamlConfigError("UnboundOrbit", "parameter", unmatched,
			  __FILE__, __LINE__);

  // Configured parameters
  //
  try {
    if (conf["MODEL"])      MODEL = conf["MODEL"].as<int>();
    if (conf["DIVERGE"])    DIVERGE = conf["DIVERGE"].as<int>();
    if (conf["DIVEXPON"])   DIVEXPON = conf["DIVEXPON"].as<double>();
    if (conf["RCORE"])      RCORE = conf["RCORE"].as<double>();
    if (conf["E"])          E = conf["E"].as<double>();
    if (conf["Rperi"])      Rperi = conf["Rperi"].as<double>();
    if (conf["Redge"])      Redge = conf["Redge"].as<double>();
    if (conf["deltaR"])     deltaR = conf["deltaR"].as<double>();
    if (conf["RMODMIN"])    RMODMIN = conf["RMODMIN"].as<double>();
    if (conf["RMODMAX"])    RMODMAX = conf["RMODMAX"].as<double>();
    if (conf["VROT"])       VROT = conf["VROT"].as<double>();
    if (conf["rmin"])       rmin = conf["rmin"].as<double>();
    if (conf["rmax"])       rmax = conf["rmax"].as<double>();
    if (conf["PHIP"])       PHIP = conf["PHIP"].as<double>();
    if (conf["THETA"])      THETA = conf["THETA"].as<double>();
    if (conf["PSI"])        PSI = conf["PSI"].as<double>();
    if (conf["INFILE"])     INFILE = conf["INFILE"].as<string>();
    if (conf["orbfile"])    orbfile = conf["orbfile"].as<bool>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in UnboundOrbit: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }


  m = 0;
  model = 0;

  switch (MODEL) {

  case file:
    m = std::make_shared<SphericalModelTable>(INFILE, DIVERGE, DIVEXPON);
    model = m;
				// Assign filename to ID string
    Model3dNames[0] = INFILE;
    break;

  case sing_isothermal:
    model = std::make_shared<SingIsothermalSphere>(1.0, 
						   RMODMIN,
						   RMODMAX);
    break;

  case isothermal:
    model = std::make_shared<IsothermalSphere>(RCORE,
					       RMODMAX,
					       VROT);
    break;

  case hernquist_model:
    model = std::make_shared<HernquistSphere>(1.0, 
					      RMODMIN, 
					      RMODMAX);
    break; 

  default:
    cerr << "Illegal model: " << MODEL << '\n';
    exit(-1);
  }

  if (rmin < 0.0) rmin = model->get_min_radius();
  if (rmax < 0.0) rmax = model->get_max_radius();

// =================================
// Compute orbit
// =================================

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

  //
  // First step
  //

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

  if (orbfile) {
    string orbfile = runtag + ".xyz";
    out.open(orbfile.c_str());
  }


  Eigen::Vector3d In, Out;
  THETA   = THETA   * M_PI/180.0;
  PSI     = PSI     * M_PI/180.0;
  PHIP    = PHIP    * M_PI/180.0;
  Eigen::Matrix3d Trans = return_euler_slater(PHIP, THETA, PSI, 1);
  
  for (unsigned i=R.size()-1; i>=1; i--) {
    In[0] = R[i]*cos(-PHI[i]);
    In[1] = R[i]*sin(-PHI[i]);
    In[2] = 0.0;
    Out  = Trans * In;

    Time.push_back(-T[i]);
    Xpos.push_back(Out[0]);
    Ypos.push_back(Out[1]);
    Zpos.push_back(Out[2]);
    
    if (out)
      out << setw(18) << -T[i]
	  << setw(18) << Out[0]
	  << setw(18) << Out[1]
	  << setw(18) << Out[2]
	  << endl;
  }

  for (unsigned i=0; i<R.size(); i++) {
    In[0] = R[i]*cos(PHI[i]);
    In[1] = R[i]*sin(PHI[i]);
    In[2] = 0.0;
    Out  = Trans * In;

    Time.push_back(T[i]);
    Xpos.push_back(Out[0]);
    Ypos.push_back(Out[1]);
    Zpos.push_back(Out[2]);
    
    if (out)
      out << setw(18) << T[i]
	  << setw(18) << Out[0]
	  << setw(18) << Out[1]
	  << setw(18) << Out[2]
	  << endl;
  }


  if (myid==0) {
    
    cout << "UnboundOrbit initialized with:" << endl
	 << setw(5) << "" << setw(10) << "THETA" << " = "
	 << THETA << endl
	 << setw(5) << "" << setw(10) << "PSI "  << " = " 
	 << PSI   << endl
	 << setw(5) << "" << setw(10) << "PHIP"  << " = " 
	 << PHIP << endl
	 << setw(5) << "" << setw(10) << "model" 
	 << " = " << MODEL << endl
	 << setw(5) << "" << setw(10) << "INFILE" 
	 << " = " << INFILE << endl
	 << setw(5) << "" << setw(10) << "E" 
	 << " = " << E << endl
	 << setw(5) << "" << setw(10) << "Rperi" 
	 << " = " << Rperi << endl
	 << setw(5) << "" << setw(10) << "Redge" 
	 << " = " << Redge << endl;
  }

}


// ===================================================================
// Destructior
// ===================================================================

UnboundOrbit::~UnboundOrbit(void)
{
  // Nothing
}

Eigen::Vector3d UnboundOrbit::get_satellite_orbit(double t)
{
  Eigen::Vector3d ret;
  t = max<double>(t, Time.front());
  t = min<double>(t, Time.back() );

  ret[0] = odd2(t, Time, Xpos);
  ret[1] = odd2(t, Time, Ypos);
  ret[2] = odd2(t, Time, Zpos);

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

