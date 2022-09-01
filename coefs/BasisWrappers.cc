#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <BasisFactory.H>

namespace py = pybind11;

void BasisFactoryClasses(py::module &m) {

  m.doc() = "BasisFactory class bindings. This module provides a factory\n"
    "class that will create biorthogonal bases from input configura-\n"
    "tion files.  These basis can then be used to compute coeffi-\n"
    "cients, provide field quantities such as forces and, together\n"
    "with the FieldGenerator, surfaces and fields for visualization.";

  using namespace Basis;

  class PyBasis : public Basis
  {
  protected:
    void all_eval(double r, double costh, double phi,
		  double& den0, double& den1,
		  double& pot0, double& pot1,
		  double& potr, double&pott, double& potp)override {
      PYBIND11_OVERRIDE_PURE(void, Basis, all_eval,
			     r, costh, phi, den0, den1, pot0, pot1,
			     potr, pott, potp);
    }
    
    void load_coefs(Coefs::CoefStrPtr coefs, double time) override {
      PYBIND11_OVERRIDE_PURE(void, Basis, load_coefs, coefs, time);
    }
    
    void set_coefs(Coefs::CoefStrPtr coefs) override {
      PYBIND11_OVERRIDE_PURE(void, Basis, set_coefs, coefs);
    }

  public:
    // Inherit the constructors
    using Basis::Basis;

    void getFields(double x, double y, double z,
		   double& tdens0, double& tpotl0, double& tdens, double& tpotl,
		   double& tpotx, double& tpoty, double& tpotz) override {
      PYBIND11_OVERRIDE_PURE(void, Basis, getFields,
			     x, y, z,
			     tdens0, tpotl0, tdens, tpotl,
			     tpotx, tpoty, tpotz
			     );
    }

    void accumulate(double x, double y, double z, double mass) override {
      PYBIND11_OVERRIDE_PURE(void, Basis, accumulate, x, y, z, mass);
    }

    void reset_coefs(void) override {
      PYBIND11_OVERRIDE_PURE(void, Basis, reset_coefs,);
    }

    void make_coefs(void) override {
      PYBIND11_OVERRIDE_PURE(void, Basis, make_coefs,);
    }

  };

  class PySphericalSL : public SphericalSL
  {
  protected:

    void all_eval(double r, double costh, double phi,
		  double& den0, double& den1,
		  double& pot0, double& pot1,
		  double& potr, double& pott, double& potp) override {
      PYBIND11_OVERRIDE(void, SphericalSL, all_eval,
			r, costh, phi, den0, den1, pot0, pot1,
			potr, pott, potp);
    }
    
    void load_coefs(Coefs::CoefStrPtr coefs, double time) override {
      PYBIND11_OVERRIDE(void, SphericalSL, load_coefs, coefs, time);
    }
    
    void set_coefs(Coefs::CoefStrPtr coefs) override {
      PYBIND11_OVERRIDE(void, SphericalSL, set_coefs, coefs);
    }

  public:

    // Inherit the constructors
    using SphericalSL::SphericalSL;
    void getFields(double x, double y, double z,
		   double& tdens0, double& tpotl0, double& tdens, double& tpotl,
		   double& tpotx, double& tpoty, double& tpotz) override {
      PYBIND11_OVERRIDE(void, SphericalSL, getFields,
			x, y, z,
			tdens0, tpotl0, tdens, tpotl,
			tpotx, tpoty, tpotz
			);
    }

    void accumulate(double x, double y, double z, double mass) override {
      PYBIND11_OVERRIDE(void, SphericalSL, accumulate, x, y, z, mass);
    }

    void reset_coefs(void) override {
      PYBIND11_OVERRIDE(void, SphericalSL, reset_coefs,);
    }

    void make_coefs(void) override {
      PYBIND11_OVERRIDE(void, SphericalSL, make_coefs,);
    }

  };

  class PyCylindricalSL : public CylindricalSL
  {
  protected:

    void all_eval(double r, double costh, double phi,
		  double& den0, double& den1,
		  double& pot0, double& pot1,
		  double& potr, double&pott, double& potp)override {
      PYBIND11_OVERRIDE(void, CylindricalSL, all_eval,
			r, costh, phi, den0, den1, pot0, pot1,
			potr, pott, potp);
    }
    
    void load_coefs(Coefs::CoefStrPtr coefs, double time) override {
      PYBIND11_OVERRIDE(void, CylindricalSL, load_coefs, coefs, time);
    }
    
    void set_coefs(Coefs::CoefStrPtr coefs) override {
      PYBIND11_OVERRIDE(void, CylindricalSL, set_coefs, coefs);
    }

  public:

    // Inherit the constructors
    using CylindricalSL::CylindricalSL;

    void getFields(double x, double y, double z,
		   double& tdens0, double& tpotl0, double& tdens, double& tpotl,
		   double& tpotx, double& tpoty, double& tpotz) override {
      PYBIND11_OVERRIDE(void, CylindricalSL, getFields,
			x, y, z,
			tdens0, tpotl0, tdens, tpotl,
			tpotx, tpoty, tpotz
			);
    }

    void accumulate(double x, double y, double z, double mass) override {
      PYBIND11_OVERRIDE(void, CylindricalSL, accumulate, x, y, z, mass);
    }

    void reset_coefs(void) override {
      PYBIND11_OVERRIDE(void, CylindricalSL, reset_coefs,);
    }

    void make_coefs(void) override {
      PYBIND11_OVERRIDE(void, CylindricalSL, make_coefs,);
    }

  };


  py::class_<Basis::Basis, std::shared_ptr<Basis::Basis>, PyBasis>(m, "Basis")
    .def(py::init<const std::string&>(),
	 "Initialize a biorthogonal basis from the configuration in the "
	 "provided YAML configuration", py::arg("YAMLstring"))
    .def("createCoefficients", &Basis::Basis::createCoefficients,
	 "Generate the coefficients from the supplied ParticleReader")
    .def("getFields",          &Basis::Basis::getFields,
	 "Return the field values for a cartesian position")
    .def("accumulate",         &Basis::Basis::accumulate,
	 "Add the contribution of a single particle to the coefficients")
    .def("getMass",            &Basis::Basis::getMass,
	 "Return the total mass of particles contributing the the current coefficient set")
    .def("reset_coefs",        &Basis::Basis::reset_coefs,
	 "Reset the coefficients to begin a generating a new set")
    .def("make_coefs",         &Basis::Basis::make_coefs,
	 "Create the coefficients after particle accumuluation is complete")
    .def("factory",            &Basis::Basis::factory_string,
	 "Generate a basis from a YAML configuration supplied as a string");

    py::class_<Basis::SphericalSL, std::shared_ptr<Basis::SphericalSL>, PySphericalSL, Basis::Basis>(m, "SphericalSL")
      .def(py::init<const std::string&>(), "Create a spherical Sturm-Liouville basis");

  py::class_<Basis::CylindricalSL, std::shared_ptr<Basis::CylindricalSL>, PyCylindricalSL, Basis::Basis>(m, "CylindricalSL")
    .def(py::init<const std::string&>(), "Create a cylindrical EOF basis");
}
