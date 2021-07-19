#ifndef _Models3d_H
#define _Models3d_H

enum Models3d {file, isothermal, sing_isothermal, low_sing_isothermal, 
	       hernquist_model, gen_polytrope, plummer};

static std::string Model3dNames[] = {"file", "Isothermal", "SingIsothermal",
				     "LowSingIsothermal", "Hernquist", 
				     "GeneralizedPolytrope", "PlummerSphere"};

#endif
