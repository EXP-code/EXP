
#ifndef _Initialize_H
#define _Initialize_H

#include <vector>

#include <boost/random.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/shared_ptr.hpp>

#include "Particle.H"

/*
  Create standard shock tube initial conditions.

  Tube is along x-axis, divided into two equal segments at x=Lx/2.

  x<Lx/2: density is D1 and temperature is T1
  x>Lx/2: density is D2 and temperature is T2

*/

/*void Initialize(std::vector<Particle>& p, double Mass,
		double D1, double T1, double D2, double T2,
		vector<double> &L, int nd=6);*/

/*
  Uniform temperature box
*/

void InitializeUniform(std::vector<Particle>& p, double Mass,
		       double T, vector<double> &L, int nd=6);
/*
  Uniform temperature box for shearing sheet equations of motion
*/

void InitializeShearingSheet(std::vector<Particle>& p, double Mass,
			     double T, vector<double> &L, 
			     double V0, double S0, int nd=6);

/*
  With a sound wave perturbation
*/

void InitializeUniformWave(std::vector<Particle>& p, double Mass,
			   double T, vector<double> &L, 
			   double Period, double Amp, int nd=6);

/*
  Uniform temperature box with shear in z direction 
  and initial sinusoidal perturbation
*/

void InitializeUniformShear(std::vector<Particle>& p, double Mass,
			    double T1, double D1, double T2, double D2,
			    vector<double> &L, double mach, double Phi,
			    double Period=1.0, double Amp=0.0,
			    int ni=0, int nd=6);


typedef boost::shared_ptr<boost::mt19937> gen_ptr;
typedef boost::shared_ptr<boost::uniform_real<> > uniform_ptr;
typedef boost::shared_ptr<boost::normal_distribution<> > normal_ptr;
typedef boost::variate_generator<boost::mt19937&, boost::uniform_real<> > unif_var;
typedef boost::shared_ptr<unif_var> unit_ptr;
typedef boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > norm_var;
typedef boost::shared_ptr<norm_var> norm_ptr;

extern unit_ptr Unit;
extern norm_ptr Norm;


#endif
