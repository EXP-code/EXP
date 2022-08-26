#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <BasisFactory.H>

namespace py = pybind11;

void BasisFactoryClasses(py::module &m) {

  m.doc() = "BasisFactory class bindings";

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


  py::class_<Basis::Basis, PyBasis>(m, "Basis")
    .def(py::init<const std::string&>())
    .def("createCoefficients", &Basis::Basis::createCoefficients)
    .def("getFields",          &Basis::Basis::getFields)
    .def("accumulate",         &Basis::Basis::accumulate)
    .def("getMass",            &Basis::Basis::getMass)
    .def("reset_coefs",        &Basis::Basis::reset_coefs)
    .def("make_coefs",         &Basis::Basis::make_coefs)
    .def("factory",            &Basis::Basis::factory_string);

  py::class_<Basis::SphericalSL, PySphericalSL>(m, "SphericalSL")
    .def(py::init<const std::string&>());

  py::class_<Basis::CylindricalSL, PyCylindricalSL>(m, "CylindricalSL")
    .def(py::init<const std::string&>());
}

