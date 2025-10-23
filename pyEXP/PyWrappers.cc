#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include "config_exp.h"

namespace py = pybind11;

extern void MSSAtoolkitClasses   (py::module &m);
extern void EDMDtoolkitClasses   (py::module &m);
extern void CoefficientClasses   (py::module &m);
extern void BasisFactoryClasses  (py::module &m);
extern void FieldGeneratorClasses(py::module &m);
extern void ParticleReaderClasses(py::module &m);
extern void UtilityClasses       (py::module &m);

PYBIND11_MODULE(pyEXP, m)
{
  m.doc() =
    "pyEXP\n"
    "=====\n"
    "Provides a collection of EXP tools for processing and analyzing\n"
    "simulation data using BFE techniques and MSSA.\n\n"
    "How to get started\n"
    "------------------\n"
    "The main documentation is the many docstrings embedded in the\n"
    "code and a set of examples provided with the EXP source.  We hope\n"
    "to provide a online reference guide in the future.\n\n"
    "We recommend beginning with the example Python scripts and IPython\n"
    "notebooks and adapting them to your own needs. You can explore\n"
    "the available classes and member functions using the usual Python\n"
    "``help'' function.  The classes are organized into seven submodules\n"
    "that are described briefly below.  Run 'help(pyEXP.xxxx) for each\n"
    "of the submodules below for more detailed usage info...\n\n"
    "The available submodules\n"
    "------------------------\n"
    "read\n"
    "     Read particle snapshots of various types.  Currently EXP,\n"
    "     Gadget, Tipsy, and Bonzai types are supported.\n"
    "basis\n"
    "     Create and apply specific biorthogonal bases to generate\n"
    "     coefficients from particle data and evaluate potential,\n"
    "     density, and force fields\n"
    "coefs\n"
    "     Classes for reading, passing, writing, converting, and\n"
    "     querying coefficient sets\n"
    "field\n"
    "     Create two- and three-dimension rectangular grids of fields\n"
    "     for visualization\n"
    "mssa\n"
    "     Tools to apply Multivariate Singular Spectrum Analysis (MSSA)\n"
    "     to the coefficients computed using the 'basis' classes\n"
    "edmd\n"
    "     Tools to apply the discrete Koopman operator analysis using\n"
    "     the extended Dynamical Mode Decomposition (EDMD) algorith to\n"
    "     approximate the Koopman operator.  The use of the Koopman class\n"
    "     echos that for mSSA.  As of this point, this is NOT a recommended\n"
    "     toolset.  If you have success with this, please post a message.\n"
    "util\n"
    "     Miscellaneous tools that support the others.  Currently this\n"
    "     include centering algorithms.  While EXP has native methods for\n"
    "     doing this, others will need to supply an estimated center\n\n"
    "Example workflow\n"
    "----------------\n"
    "To provide some context, suppose you want to read some snapshots,\n"
    "make some coefficients, and then analyze them with MSSA. The class\n"
    "construction would go something like this:\n"
    "  1. Create a reader instance for your simulation, call it 'reader'.\n"
    "  2. Create a basis designed to represent the a particular particle\n"
    "     type.  Star particles, perhaps, so let's call it 'disk'.\n"
    "  3. We then pass 'reader' to the createCoefficients member of 'disk'\n"
    "     to get coefficients for your snapshots, called 'coefs'\n"
    "  4. We might then want to explore dynamical patterns in these\n"
    "     coefficients by passing 'coefs' to 'expMSSA'.  'expMSSA' will\n"
    "     return principal signals as an updated coefficient object,\n"
    "     that we call 'newcoefs'\n"
    "  5. 'newcoefs' and 'disk' can then be passed to the FieldGenerator\n"
    "     to provide density, potential, force fields, etc. for the each\n"
    "     principal signal\n"
    "This is only one example of many possible uses.  There are many\n"
    "variants to this work flow, of course, and I expect that you will\n"
    "invent some interesting ones\n\n"
    "The source code has some sample Python scripts and notebooks for a\n"
    "quick start (check the pyEXP directory).\n\n"
    "Version information\n"
    "-----------------------\n"
    "Version:          "  PACKAGE_STRING "\n"
    "Repository URL:   "  PACKAGE_URL "\n"
    "GIT branch:       "  GIT_BRANCH "\n"
    "GIT commit:       "  GIT_COMMIT "\n"
    "Compile time:     "  COMPILE_TIME "\n\n"
    "History and provenance\n"
    "----------------------\n"
    "These EXP interface classes and the Python interface were written\n"
    "to fill a demand by collaborators for Python access to EXP.   EXP\n"
    "itself represents 20 years of experience applying BFE techniques to\n"
    "both simulation and perturbation theory.  Recent innovations include\n"
    "the application of MSSA to discover and characterize some key\n"
    "dynamical features in simulations that are hard to find `by eye'.\n\n"
    "Please send comments, suggestions, and particularly good cookies to:\n"
    "mdw@umass.edu (Martin Weinberg)\n\n";
  
  auto mod_coefs = m.def_submodule("coefs",
				   "Classes for reading, passing and "
				   "manipulating coefficient sets");

  auto mod_basis = m.def_submodule("basis",
				   "Create and apply specific biorthogonal "
				   "bases to generate coefficients from "
				   "particle data and evaluate potential, "
				   "density, and force fields");

  auto mod_field = m.def_submodule("field",
				   "Create two- and three-dimension rectangular "
				   "grids of fields for visualization");

  auto mod_read = m.def_submodule("read", "Read particle snapshots of various "
				  "types.  Currently EXP, Gadget, Tipsy, and "
				  "Bonzai types are supported.");

  auto mod_mssa = m.def_submodule("mssa", "Tools to apply Multivariate Singular "
				  "Spectrum Analysis (MSSA) to the coefficients "
				  "computed using the 'basis' classes");
  
  auto mod_edmd = m.def_submodule("edmd", "Tools to apply extended Dynamical Mode "
				  "Decomposition to the coefficients computed using "
				  "the 'basis' classes");
  
  auto mod_util = m.def_submodule("util", "Miscellaneous tools that support the "
				  "others.  Currently these contain several "
				  "centering algorithms.");

  CoefficientClasses(mod_coefs);
  BasisFactoryClasses(mod_basis);
  FieldGeneratorClasses(mod_field);
  ParticleReaderClasses(mod_read);
  MSSAtoolkitClasses(mod_mssa);
  EDMDtoolkitClasses(mod_edmd);
  UtilityClasses(mod_util);
}

