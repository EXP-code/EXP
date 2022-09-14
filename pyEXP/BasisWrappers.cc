#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <BasisFactory.H>

namespace py = pybind11;

void BasisFactoryClasses(py::module &m) {

  m.doc() = "BasisFactory class bindings\n\n"
    "This module provides a factory class that will create biorthogonal\n"
    "bases from input YAML configuration files.  Each basis can then be\n"
    "used to compute coefficients, provide field quantities such as\n"
    "forces and, together with the FieldGenerator, surfaces and fields for\n"
    "visualization.\n\n"
    "Two bases are currently implemented: SphericalSL, the Sturm-\n"
    "Liouiville spherical basis and the Cylindrical basis, which is\n"
    "created by computing empirical orthogonal functions over a densely\n"
    "sampled SphericalSL basis.  Each of these bases take a YAML\n"
    "configuration file as input. These are the same parameter lists used\n"
    "by EXP. See the EXP documentation and examples for more detail.\n\n"
    "The primary functions of these basis classes are:\n"
    "  1. To compute BFE coefficients from phase-space snapshots\n"
    "     using the ParticleReader class. See help(pyEXP.read).\n"
    "  2. To evaluate the fields from the basis and a coefficient\n"
    "     object. See help(pyEXP.coefs) and help(pyEXP.field).\n\n";

  using namespace Basis;

  class PyBasis : public Basis
  {
  protected:
    void all_eval(double r, double costh, double phi,
		  double& den0, double& den1,
		  double& pot0, double& pot1,
		  double& potr, double& pott, double& potp)override {
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

  class PyCylindrical : public Cylindrical
  {
  protected:

    void all_eval(double r, double costh, double phi,
		  double& den0, double& den1,
		  double& pot0, double& pot1,
		  double& potr, double& pott, double& potp)override {
      PYBIND11_OVERRIDE(void, Cylindrical, all_eval,
			r, costh, phi, den0, den1, pot0, pot1,
			potr, pott, potp);
    }
    
    void load_coefs(Coefs::CoefStrPtr coefs, double time) override {
      PYBIND11_OVERRIDE(void, Cylindrical, load_coefs, coefs, time);
    }
    
    void set_coefs(Coefs::CoefStrPtr coefs) override {
      PYBIND11_OVERRIDE(void, Cylindrical, set_coefs, coefs);
    }

  public:

    // Inherit the constructors
    using Cylindrical::Cylindrical;

    void getFields(double x, double y, double z,
		   double& tdens0, double& tpotl0, double& tdens, double& tpotl,
		   double& tpotx, double& tpoty, double& tpotz) override {
      PYBIND11_OVERRIDE(void, Cylindrical, getFields,
			x, y, z,
			tdens0, tpotl0, tdens, tpotl,
			tpotx, tpoty, tpotz
			);
    }

    void accumulate(double x, double y, double z, double mass) override {
      PYBIND11_OVERRIDE(void, Cylindrical, accumulate, x, y, z, mass);
    }

    void reset_coefs(void) override {
      PYBIND11_OVERRIDE(void, Cylindrical, reset_coefs,);
    }

    void make_coefs(void) override {
      PYBIND11_OVERRIDE(void, Cylindrical, make_coefs,);
    }

  };


  py::class_<Basis::Basis, std::shared_ptr<Basis::Basis>, PyBasis>(m, "Basis")
    .def(py::init<const std::string&>(),
	 "Initialize a biorthogonal basis from the configuration in the\n"
	 "provided YAML configuration", py::arg("YAMLstring"))
    .def("createFromReader", &Basis::Basis::createFromReader,
	 "Generate the coefficients from the supplied ParticleReader and\n"
	 "an optional expansion center location",
	 py::arg("reader"), 
	 py::arg("center") = std::vector<double>(3, 0.0))
    .def("createFromArray",
	 [](Basis::Basis& A, RowMatrixXd& data, double time,
	    std::vector<double> center)
	 {
	   return A.createFromArray(data, time, center);
	 },
	 "Generate the coefficients from an array with rows with the values\n"
	 "(mass, x, y, z, u, v, w), time, and an optional expansion center\n"
	 "location",
	 py::arg("data"), py::arg("time"),
	 py::arg("center") = std::vector<double>(3, 0.0))
    .def("getFields",
	 [](Basis::Basis& A, double x, double y, double z)
	 {
	   std::vector<double> ret(7);
	   A.getFields(x, y, z,
		       ret[0], ret[1], ret[2], ret[3],
		       ret[4], ret[5], ret[6]);
	   return ret;
	 },
	 "Return the density, potential, and forces for a cartesian position.\n"
	 "Field order is: dens0, potl0, dens, potl, fx, fy, fz. Dens0 and\n"
	 "potl0 are the fields evaluated for l=0 or m=0 and dens and potl\n"
	 "are evaluated for l>0 or m>0\n",
	 py::arg("x"), py::arg("y"), py::arg("z"))
    .def("accumulate",         &Basis::Basis::accumulate,
	 "Add the contribution of a single particle to the coefficients")
    .def("getMass",            &Basis::Basis::getMass,
	 "Return the total mass of particles contributing the the current\n"
	 "coefficient set")
    .def("reset_coefs",        &Basis::Basis::reset_coefs,
	 "Reset the coefficients to begin a generating a new set")
    .def("make_coefs",         &Basis::Basis::make_coefs,
	 "Create the coefficients after particle accumuluation is complete")
    .def("set_coefs",          &Basis::Basis::set_coefs,
	 "Install a new set of coefficients from a CoefStruct")
    .def("factory",            &Basis::Basis::factory_string,
	 "Generate a basis from a YAML configuration supplied as a string");

    py::class_<Basis::SphericalSL, std::shared_ptr<Basis::SphericalSL>, PySphericalSL, Basis::Basis>(m, "SphericalSL")
      .def(py::init<const std::string&>(), "Create a spherical Sturm-Liouville basis")
      .def("getBasis", &Basis::SphericalSL::getBasis,
	   "Evaluate the basis functions on a logarithmically spaced grid for"
	   "inspection",
	   py::arg("logxmin")=-3.0,
	   py::arg("logxmax")=0.5,
	   py::arg("numr")=400);

  py::class_<Basis::Cylindrical, std::shared_ptr<Basis::Cylindrical>, PyCylindrical, Basis::Basis>(m, "Cylindrical")
    .def(py::init<const std::string&>(), "Create a cylindrical EOF basis")
    .def("getBasis", &Basis::Cylindrical::getBasis,
	 "Evaluate the basis functions on a linearly spaced 2d-grid for"
	 "inspection",
	 py::arg("xmin")=0.0,
	 py::arg("xmax")=1.0,
	 py::arg("numr")=40,
	 py::arg("zmin")=-0.1,
	 py::arg("zmax")=0.1,
	 py::arg("numz")=40 );
}
