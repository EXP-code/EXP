#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <CoefContainer.H>

namespace py = pybind11;

void CoefContainerClasses(py::module &m) {

  m.doc() = "CoefContainer class bindings. This class stores, writes, and\n"
    "provides and interface to coefficients and table data for use\n"
    "by the other pyEXP classes.\n\n"
    "The base class, 'Coefs', provides a factory reader that will\n"
    "create one of the derived coefficient classes, SphCoef, CylCoef,\n"
    "or TblCoef, deducing the type from the input file.  The input\n"
    "files may be EXP native or HDF5 cofficient files.  The Basis\n"
    "factory, Basis::createCoefficients, will create set of coef-\n"
    "ficients from phase-space snapshots.  See Basis\n\n"
    "Once created, one may get a list of times, get the total gravi-\n"
    "tation power, and write a new HDF5 file.  Their main use is\n"
    "as a container object for MSSA (using expMSSA) and field\n"
    "visualization using the FieldGenerator class.  See\n"
    "help(pyEXP.mssa) and help(pyEXP.field) for more details.\n\n"
    "NB: the time list, returned with the Times() member function\n"
    "is truncated to 8 signficant figures so that it may be used\n"
    "as a dictionary key\n\n";

  using namespace Coefs;

  class PyCoefStruct : public CoefStruct
  {
  public:
    
    // Inherit the constructors
    using CoefStruct::CoefStruct;

    bool read(std::istream& in, bool exp_type) override {
      PYBIND11_OVERRIDE_PURE(bool, CoefStruct, read, in, exp_type);
    }
  };

  class PyCoefs : public Coefs
  {
  protected:
    void readNativeCoefs(const std::string& file) override {
      PYBIND11_OVERRIDE_PURE(void, Coefs, readNativeCoefs, file);
    }
    
    std::string getYAML() override {
      PYBIND11_OVERRIDE_PURE(std::string, Coefs, getYAML,);
    }
    
    void WriteH5Params(HighFive::File& file) override {
      PYBIND11_OVERRIDE_PURE(void, Coefs, WriteH5Params, file);
    }

    unsigned WriteH5Times(HighFive::Group& group, unsigned count) override {
      PYBIND11_OVERRIDE_PURE(unsigned, Coefs, WriteH5Times, group, count);
    }

  public:
    // Inherit the constructors
    using Coefs::Coefs;

    Eigen::MatrixXcd& operator()(double time) override {
      PYBIND11_OVERRIDE_PURE(Eigen::MatrixXcd&, Coefs, operator(), time);
    }

    std::shared_ptr<CoefStruct> getCoefStruct(double time) override {
      PYBIND11_OVERRIDE_PURE(std::shared_ptr<CoefStruct>, Coefs, getCoefStruct,
			     time);
    }

    void dump(int mmin, int mmax, int nmin, int nmax) override {
      PYBIND11_OVERRIDE_PURE(void, Coefs, dump,
			     mmin, mmax, nmin, nmax);
    }

    std::vector<double> Times() override {
      PYBIND11_OVERRIDE_PURE(std::vector<double>, Coefs, Times,);
    }

    void WriteH5Coefs(const std::string& prefix) override {
      PYBIND11_OVERRIDE_PURE(void, Coefs, WriteH5Coefs, prefix);
    }

    void ExtendH5Coefs(const std::string& prefix) override {
      PYBIND11_OVERRIDE_PURE(void, Coefs, ExtendH5Coefs, prefix);
    }

    Eigen::MatrixXd& Power() override {
      PYBIND11_OVERRIDE_PURE(Eigen::MatrixXd&, Coefs, Power,);
    }

    bool CompareStanzas(std::shared_ptr<Coefs> check) override {
      PYBIND11_OVERRIDE_PURE(bool, Coefs, CompareStanzas, check);
    }

    void clear() override {
      PYBIND11_OVERRIDE_PURE(void, Coefs, clear,);
    }

    void add(CoefStrPtr coef) override {
      PYBIND11_OVERRIDE_PURE(void, Coefs, add, coef);
    }

    std::vector<Key> makeKeys(Key k) override {
      PYBIND11_OVERRIDE_PURE(std::vector<Key>, Coefs, makeKeys, k);
    }
  };

  class PySphCoefs : public SphCoefs
  {
  protected:
    void readNativeCoefs(const std::string& file) override {
      PYBIND11_OVERRIDE(void, SphCoefs, readNativeCoefs, file);
    }
    
    std::string getYAML() override {
      PYBIND11_OVERRIDE(std::string, SphCoefs, getYAML,);
    }
    
    void WriteH5Params(HighFive::File& file) override {
      PYBIND11_OVERRIDE(void, SphCoefs, WriteH5Params, file);
    }

    unsigned WriteH5Times(HighFive::Group& group, unsigned count) override {
      PYBIND11_OVERRIDE(unsigned, SphCoefs, WriteH5Times, group, count);
    }

  public:
    // Inherit the constructors
    using SphCoefs::SphCoefs;

    Eigen::MatrixXcd& operator()(double time) override {
      PYBIND11_OVERRIDE(Eigen::MatrixXcd&, SphCoefs, operator(), time);
    }

    std::shared_ptr<CoefStruct> getCoefStruct(double time) override {
      PYBIND11_OVERRIDE(std::shared_ptr<CoefStruct>, SphCoefs, getCoefStruct,
			time);
    }

    void dump(int mmin, int mmax, int nmin, int nmax) override {
      PYBIND11_OVERRIDE(void, SphCoefs,dump,
			mmin, mmax, nmin, nmax);
    }

    std::vector<double> Times() override {
      PYBIND11_OVERRIDE(std::vector<double>, SphCoefs, Times,);
    }

    void WriteH5Coefs(const std::string& prefix) override {
      PYBIND11_OVERRIDE(void, SphCoefs, WriteH5Coefs, prefix);
    }

    void ExtendH5Coefs(const std::string& prefix) override {
      PYBIND11_OVERRIDE(void, SphCoefs, ExtendH5Coefs, prefix);
    }

    Eigen::MatrixXd& Power() override {
      PYBIND11_OVERRIDE(Eigen::MatrixXd&, SphCoefs, Power,);
    }

    bool CompareStanzas(std::shared_ptr<Coefs> check) override {
      PYBIND11_OVERRIDE(bool, SphCoefs, CompareStanzas, check);
    }

    void clear() override {
      PYBIND11_OVERRIDE(void, SphCoefs,	clear,);
    }

    void add(CoefStrPtr coef) override {
      PYBIND11_OVERRIDE(void, SphCoefs,	add, coef);
    }

    std::vector<Key> makeKeys(Key k) override {
      PYBIND11_OVERRIDE(std::vector<Key>, SphCoefs, makeKeys, k);
    }
  };

  class PyCylCoefs : public CylCoefs
  {
  protected:
    void readNativeCoefs(const std::string& file) override {
      PYBIND11_OVERRIDE(void, CylCoefs,	readNativeCoefs, file);
    }
    
    std::string getYAML() override {
      PYBIND11_OVERRIDE(std::string, CylCoefs, getYAML,);
    }
    
    void WriteH5Params(HighFive::File& file) override {
      PYBIND11_OVERRIDE(void, CylCoefs, WriteH5Params, file);
    }

    unsigned WriteH5Times(HighFive::Group& group, unsigned count) override {
      PYBIND11_OVERRIDE(unsigned, CylCoefs, WriteH5Times, group, count);
    }

  public:
    // Inherit the constructors
    using CylCoefs::CylCoefs;

    Eigen::MatrixXcd& operator()(double time) override {
      PYBIND11_OVERRIDE(Eigen::MatrixXcd&, CylCoefs, operator(), time);
    }

    std::shared_ptr<CoefStruct> getCoefStruct(double time) override {
      PYBIND11_OVERRIDE(std::shared_ptr<CoefStruct>, CylCoefs, getCoefStruct,
			time);
    }

    void dump(int mmin, int mmax, int nmin, int nmax) override {
      PYBIND11_OVERRIDE(void, CylCoefs, dump,
			mmin, mmax, nmin, nmax);
    }

    std::vector<double> Times() override {
      PYBIND11_OVERRIDE(std::vector<double>, CylCoefs, Times,);
    }

    void WriteH5Coefs(const std::string& prefix) override {
      PYBIND11_OVERRIDE(void, CylCoefs, WriteH5Coefs, prefix);
    }

    void ExtendH5Coefs(const std::string& prefix) override {
      PYBIND11_OVERRIDE(void, CylCoefs, ExtendH5Coefs, prefix);
    }

    Eigen::MatrixXd& Power() override {
      PYBIND11_OVERRIDE(Eigen::MatrixXd&, CylCoefs, Power,);
    }

    bool CompareStanzas(std::shared_ptr<Coefs> check) override {
      PYBIND11_OVERRIDE(bool, CylCoefs, CompareStanzas,	check);
    }

    void clear() override {
      PYBIND11_OVERRIDE(void, CylCoefs, clear,);
    }

    void add(CoefStrPtr coef) override {
      PYBIND11_OVERRIDE(void, CylCoefs,	add, coef);
    }

    std::vector<Key> makeKeys(Key k) override {
      PYBIND11_OVERRIDE(std::vector<Key>, CylCoefs, makeKeys, k);
    }
  };

  class PyTableData : public TableData
  {
  protected:
    void readNativeCoefs(const std::string& file) override {
      PYBIND11_OVERRIDE(void, TableData, readNativeCoefs, file);
    }
    
    std::string getYAML() override {
      PYBIND11_OVERRIDE(std::string, TableData, getYAML,);
    }
    
    void WriteH5Params(HighFive::File& file) override {
      PYBIND11_OVERRIDE(void, TableData, WriteH5Params, file);
    }

    unsigned WriteH5Times(HighFive::Group& group, unsigned count) override {
      PYBIND11_OVERRIDE(unsigned, TableData, WriteH5Times, group, count);
    }

  public:
    // Inherit the constructors
    using TableData::TableData;

    Eigen::MatrixXcd& operator()(double time) override {
      PYBIND11_OVERRIDE(Eigen::MatrixXcd&, TableData, operator(), time);
    }

    std::shared_ptr<CoefStruct> getCoefStruct(double time) override {
      PYBIND11_OVERRIDE(std::shared_ptr<CoefStruct>, TableData, getCoefStruct,
			time);
    }

    void dump(int mmin, int mmax, int nmin, int nmax) override {
      PYBIND11_OVERRIDE(void, TableData, dump,
			mmin, mmax, nmin, nmax);
    }

    std::vector<double> Times() override {
      PYBIND11_OVERRIDE(std::vector<double>, TableData, Times,);
    }

    void WriteH5Coefs(const std::string& prefix) override {
      PYBIND11_OVERRIDE(void, TableData, WriteH5Coefs, prefix);
    }

    void ExtendH5Coefs(const std::string& prefix) override {
      PYBIND11_OVERRIDE(void, TableData, ExtendH5Coefs, prefix);
    }

    Eigen::MatrixXd& Power() override {
      PYBIND11_OVERRIDE(Eigen::MatrixXd&, TableData, Power,);
    }

    bool CompareStanzas(std::shared_ptr<Coefs> check) override {
      PYBIND11_OVERRIDE(bool, TableData, CompareStanzas, check);
    }

    void clear() override {
      PYBIND11_OVERRIDE(void, TableData, clear,);
    }

    void add(CoefStrPtr coef) override {
      PYBIND11_OVERRIDE(void, TableData, add, coef);
    }

    std::vector<Key> makeKeys(Key k) override {
      PYBIND11_OVERRIDE(std::vector<Key>, TableData, makeKeys, k);
    }
  };

  py::class_<Coefs::CoefStruct, std::shared_ptr<Coefs::CoefStruct>, PyCoefStruct>(m, "CoefStruct")
    .def(py::init<>(), "Base class coefficient data structure object");

  py::class_<Coefs::SphStruct, std::shared_ptr<Coefs::SphStruct>, CoefStruct>(m, "SphStruct")
    .def(py::init<>(), "Spherical coefficient data structure object");

  py::class_<Coefs::CylStruct, std::shared_ptr<Coefs::CylStruct>, CoefStruct>(m, "CylStruct")
    .def(py::init<>(), "Cylindrical coefficient data structure object");

  py::class_<Coefs::Coefs, std::shared_ptr<Coefs::Coefs>, PyCoefs>(m, "Coefs")
    .def(py::init<std::string, bool>(),
	 "Base coefficient container class",
	 py::arg("type"), py::arg("verbose"))
    .def("operator()",     &Coefs::Coefs::operator(),
	 "Return the coefficient matrix for the desired time",
	 py::arg("time"))
    .def("add",            &Coefs::Coefs::add,
	 "Add a coefficient structure to the coefficient container",
	 py::arg("coef"))
    .def("getCoefStruct",  &Coefs::Coefs::getCoefStruct,
	 "Return the coefficient structure for the desired time",
	 py::arg("time"))
    .def("Times",          &Coefs::Coefs::Times,
	 "Return a list of times for coefficient sets current in the container")
    .def("WriteH5Coefs",   &Coefs::Coefs::WriteH5Coefs,
	 "Write the coefficients into an EXP HDF5 coefficieint file with the "
	 "given prefix name", py::arg("filename"))
    .def("ExtendH5Coefs",  &Coefs::Coefs::ExtendH5Coefs,
	 "Extend an existing EXP HDF5 coefficient file with given prefix "
	 "name using the coefficient in the container", py::arg("filename"))
    .def("Power",          &Coefs::Coefs::Power,
	 "Return a ndarray table of the full power for the top-level harmonic "
	 "index as function of time")
    .def("makeKeys",       &Coefs::Coefs::makeKeys,
	 "Return a vector/list of keys for an entire subspace of "
	 "subdimensional rank", py::arg("subkey"))
    .def("getGeometry",    &Coefs::Coefs::getGeometry,
	 "Return the coefficient geometry string")
    .def("getName",        &Coefs::Coefs::getName,
	 "Return the coefficient set nmenonic name")
    .def("setName",        &Coefs::Coefs::setName,
	 "Set or rename the coefficient set nmenonic name", py::arg("newname"))
    .def("CompareStanzas", &Coefs::Coefs::CompareStanzas,
	 "Check that the data in one Coefs set is identical to "
	 "that in another")
    .def_static("factory", &Coefs::Coefs::factory,
		"Deduce the type and read coefficients from a native or HDF5 file",
		py::arg("file"))
    .def_static("makecoefs", &Coefs::Coefs::makecoefs,
		"Create a new coefficient instance compatible with the "
		"supplied coefficient structure",
		py::arg("coef"), py::arg("name")="");

  py::class_<Coefs::SphCoefs, std::shared_ptr<Coefs::SphCoefs>, PySphCoefs, Coefs::Coefs>(m, "SphCoefs", "Container for spherical coefficients")
    .def(py::init<bool>());

  py::class_<Coefs::CylCoefs, std::shared_ptr<Coefs::CylCoefs>, PyCylCoefs, Coefs::Coefs>(m, "CylCoefs", "Container for cylindrical coefficients")
    .def(py::init<bool>());

  py::class_<Coefs::TableData, std::shared_ptr<Coefs::TableData>, PyTableData, Coefs::Coefs>(m, "TableData", "Container for simple data tables with multiple columns")
    .def(py::init<bool>());
}

