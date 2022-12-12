#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <CoefContainer.H>

namespace py = pybind11;

#include "TensorToArray.H"

void CoefficientClasses(py::module &m) {

  m.doc() = "Coefficient class bindings\n\n"
    "These classes store, write, and provide an interface to coefficients\n"
    "and table data for use by the other pyEXP classes.\n\n"
    "CoefStruct\n"
    "----------\n"
    "The CoefStruct class is low-level structure that stores the data\n"
    "and metadata specific to each geometry. These are spherical\n"
    "(SphStruct), cylindrical (CylStruct), and table data (TblStruct).\n"
    "EXP also knows about rectangular grids and slabs.  These may be\n"
    "added in a future release if there is a need.  Instances of these\n"
    "structures represent individual times points and are created,\n"
    "maintained, and interfaced by the Coefs class.  Access to the\n"
    "underlying data is provided to Python in case you need to change\n"
    "or rewrite the data for some reason.  We have also provided a\n"
    "create() member so that you can instaniate and load a coefficient\n"
    "structure using Python.  To do this, use the constructor to make\n"
    "a blank instance, assign the dimensions and use create() to create\n"
    "a data matrix of initially zero values.  The dimensions are \n",
    "(lmax, nmax) for SphStruct, (mmax,nmax) for a CylStruct, and\n"
    "(cols) for a TblStruct.\n\n"
    "Coefs\n"
    "-----\n"
    "The base class, 'Coefs', provides a factory reader that will\n"
    "create one of the derived coefficient classes, SphCoef, CylCoef,\n"
    "or TblCoef, deducing the type from the input file.  The input\n"
    "files may be EXP native or HDF5 cofficient files.  The Basis\n"
    "factory, Basis::createCoefficients, will create set of coef-\n"
    "ficients from phase-space snapshots.  See help(pyEXP.basis).\n"
    "Files which are not recognized as EXP coefficient files are\n"
    "assumed to be data files and are parsed by the TblCoefs class.\n"
    "The first column in data tables is interpreted as time and each\n"
    "successive column is interpreted as a new data field.\n\n"
    "Once created, you may get a list of times, get the total gravi-\n"
    "tation power from biothogonal basis coefficients, and write a new\n"
    "HDF5 file.  Their main use is as a container object for MSSA (using\n"
    "expMSSA) and field visualization using the FieldGenerator class.\n\n"
    "Updates\n"
    "-------\n"
    "The expMSSA class will update the contribution to the coefficients\n"
    "specified by key from each eigen component to the reconstructed\n"
    "series. Unspecified coefficients series will not be updated and\n"
    "their original data will be intact. For visualization, the series\n"
    "data in a Coefs object may be zeroed using the 'zerodata()' member\n"
    "function prior to an expMSSA update.  This allows one to include\n"
    "reconstructions that *only* include particular eigen components for\n"
    "the coefficients specified by key.  Then, one can visualize only the\n"
    "updated fields using 'FieldGenerator'. See help(pyEXP.mssa) and\n"
    "help(pyEXP.field) for more details.\n\n"
    "Dataset indexing\n"
    "----------------\n"
    "Coefficients and other auxilliary data from simulations are stored\n"
    "and retrieved by their time.  Internally, these are floating fixed-\n"
    "point values truncated to 8 signficant figures so that they may be\n"
    "used as dictionary/map keys.  The values in the time list, returned\n"
    "with the Times() member function contain the truncated values for\n"
    "reference.\n\n"
    "Object lifetime\n"
    "---------------\n"
    "As in native Python, the memory for created objects persists until\n"
    "it is no longer referenced.  For example, replacing the variable with\n"
    "a new set of coefficients will allow the memory to be deallocated if\n"
    "no other class instance holds a reference. Because coefficient sets\n"
    "can be large, the creation of a Coefs instance by the 'factory' will\n"
    "be passed to any other class that needs it.  The MSSA class, expMSSA,\n"
    "will hold a reference to the the Coefs object passed on creation, and\n"
    "it update the values of the coefficients on reconstruction, without\n"
    "copying. If you want to keep the initial set without change, we have\n"
    "provided a 'deepcopy()' member that provides a byte-by-byte copy.\n\n";

  using namespace CoefClasses;

  class PyCoefStruct : public CoefStruct
  {
  public:
    
    // Inherit the constructors
    using CoefStruct::CoefStruct;

    bool read(std::istream& in, bool exp_type, bool verbose) override {
      PYBIND11_OVERRIDE_PURE(bool, CoefStruct, read, in, exp_type, verbose);
    }

    void create() override {
      PYBIND11_OVERRIDE_PURE(void, CoefStruct, create,);
    }

    std::shared_ptr<CoefStruct> deepcopy() override {
      PYBIND11_OVERRIDE_PURE(std::shared_ptr<CoefStruct>, CoefStruct, deepcopy,);
    }

  };

  class PyCoefs : public Coefs
  {
  protected:
    void readNativeCoefs(const std::string& file, int stride, double tmin, double tmax) override {
      PYBIND11_OVERRIDE_PURE(void, Coefs, readNativeCoefs, file, stride, tmin, tmax);
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
    using CoefClasses::Coefs::Coefs;

    Eigen::MatrixXcd& operator()(double time) override {
      PYBIND11_OVERRIDE_PURE(Eigen::MatrixXcd&, Coefs, operator(), time);
    }

    void setMatrix(double time, const Eigen::MatrixXcd& mat) override {
      PYBIND11_OVERRIDE_PURE(void, Coefs, setMatrix, time, mat);
    }

    using ValueError = std::tuple<Eigen::MatrixXcd&, bool>;
    ValueError interpolate(double time) override {
      PYBIND11_OVERRIDE(ValueError, Coefs, interpolate, time);
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

    Eigen::MatrixXd& Power(int min, int max) override {
      PYBIND11_OVERRIDE_PURE(Eigen::MatrixXd&, Coefs, Power, min, max);
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

    std::shared_ptr<Coefs> deepcopy() override {
      PYBIND11_OVERRIDE_PURE(std::shared_ptr<Coefs>, Coefs, deepcopy,);
    }

    void zerodata() override {
      PYBIND11_OVERRIDE_PURE(void, Coefs, zerodata,);
    }
  };

  class PySphCoefs : public SphCoefs
  {
  protected:
    void readNativeCoefs(const std::string& file, int stride, double tmin, double tmax) override {
      PYBIND11_OVERRIDE(void, SphCoefs, readNativeCoefs, file, stride, tmin, tmax);
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

    void setMatrix(double time, const Eigen::MatrixXcd& mat) override {
      PYBIND11_OVERRIDE(void, SphCoefs, setMatrix, time, mat);
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

    Eigen::MatrixXd& Power(int min, int max) override {
      PYBIND11_OVERRIDE(Eigen::MatrixXd&, SphCoefs, Power, min, max);
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

    std::shared_ptr<Coefs> deepcopy() override {
      PYBIND11_OVERRIDE(std::shared_ptr<Coefs>, SphCoefs, deepcopy,);
    }

    void zerodata() override {
      PYBIND11_OVERRIDE(void, SphCoefs, zerodata,);
    }


  };

  class PyCylCoefs : public CylCoefs
  {
  protected:
    void readNativeCoefs(const std::string& file, int stride, double tmin, double tmax) override {
      PYBIND11_OVERRIDE(void, CylCoefs,	readNativeCoefs, file, stride, tmin, tmax);
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

    void setMatrix(double time, const Eigen::MatrixXcd& mat) override {
      PYBIND11_OVERRIDE(void, CylCoefs, setMatrix, time, mat);
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

    Eigen::MatrixXd& Power(int min, int max) override {
      PYBIND11_OVERRIDE(Eigen::MatrixXd&, CylCoefs, Power, min, max);
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

    std::shared_ptr<Coefs> deepcopy() override {
      PYBIND11_OVERRIDE(std::shared_ptr<Coefs>, CylCoefs, deepcopy,);
    }

    void zerodata() override {
      PYBIND11_OVERRIDE(void, CylCoefs, zerodata,);
    }

  };

  class PyTableData : public TableData
  {
  protected:
    void readNativeCoefs(const std::string& file, int stride, double tmin, double tmax) override {
      PYBIND11_OVERRIDE(void, TableData, readNativeCoefs, file, stride, tmin, tmax);
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

    void setMatrix(double time, const Eigen::MatrixXcd& mat) override {
      PYBIND11_OVERRIDE(void, TableData, setMatrix, time, mat);
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

    Eigen::MatrixXd& Power(int min, int max) override {
      PYBIND11_OVERRIDE(Eigen::MatrixXd&, TableData, Power, min, max);
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

    std::shared_ptr<Coefs> deepcopy() override {
      PYBIND11_OVERRIDE(std::shared_ptr<Coefs>, TableData, deepcopy,);
    }

    void zerodata() override {
      PYBIND11_OVERRIDE(void, TableData, zerodata,);
    }

  };

  py::class_<CoefClasses::CoefStruct, std::shared_ptr<CoefClasses::CoefStruct>, PyCoefStruct>(m, "CoefStruct")
    .def(py::init<>(), "Base class coefficient data structure object")
    .def("create",            &CoefStruct::create, 
	 "Initialize a coefficient zeroed structure from user supplied "
	 "dimensions")
    .def("deepcopy",          &CoefStruct::deepcopy,
	 "Make a new instance and copy all data into the new instance")
    .def_readonly("geometry", &CoefStruct::geom,  "The geometry type")
    .def_readwrite("time",    &CoefStruct::time,  "The data's time")
    .def("data",              &CoefStruct::data,
	 "Read-write access to the underlying data store. The Eigen complex\n"
	 "valued matrix is mapped to a numpy.ndarray of complex values.\n"
	 "Changes to the data array will automatically mapped back to the\n"
	 "C++ CoefStruct instance.");
  

  py::class_<CoefClasses::SphStruct, std::shared_ptr<CoefClasses::SphStruct>, CoefStruct>(m, "SphStruct")
    .def(py::init<>(), "Spherical coefficient data structure object");

  py::class_<CoefClasses::CylStruct, std::shared_ptr<CoefClasses::CylStruct>, CoefStruct>(m, "CylStruct")
    .def(py::init<>(), "Cylindrical coefficient data structure object");

  py::class_<CoefClasses::TblStruct, std::shared_ptr<CoefClasses::TblStruct>, CoefStruct>(m, "TblStruct")
    .def(py::init<>(), "Multicolumn table data structure object");

  py::class_<CoefClasses::Coefs, std::shared_ptr<CoefClasses::Coefs>, PyCoefs>(m, "Coefs")
    .def(py::init<std::string, bool>(),
	 "Base coefficient container class",
	 py::arg("type"), py::arg("verbose"))
    .def("__call__",       &CoefClasses::Coefs::operator(),
	 "Return the coefficient matrix for the desired time",
	 py::arg("time"))
    .def("setMatrix",      &CoefClasses::Coefs::setMatrix,
	 "Rewrite the coefficient matrix at the time provided. For those\n"
	 "tracking memory use, these will be copied,not passed by reference.",
	 py::arg("time"), py::arg("mat"))
    .def("add",            &CoefClasses::Coefs::add,
	 "Add a coefficient structure to the coefficient container",
	 py::arg("coef"))
    .def("getCoefStruct",  &CoefClasses::Coefs::getCoefStruct,
	 "Return the coefficient structure for the desired time",
	 py::arg("time"))
    .def("Times",          &CoefClasses::Coefs::Times,
	 "Return a list of times for coefficient sets current in the container")
    .def("WriteH5Coefs",   &CoefClasses::Coefs::WriteH5Coefs,
	 "Write the coefficients into an EXP HDF5 coefficieint file with the "
	 "given prefix name", py::arg("filename"))
    .def("ExtendH5Coefs",  &CoefClasses::Coefs::ExtendH5Coefs,
	 "Extend an existing EXP HDF5 coefficient file with given prefix "
	 "name using the coefficient in the container", py::arg("filename"))
    .def("Power",          &CoefClasses::Coefs::Power,
	 "Return a ndarray table of the full power for the top-level harmonic "
	 "index as function of time",
	 py::arg("min")=0, py::arg("max")=std::numeric_limits<int>::max())
    .def("makeKeys",       &CoefClasses::Coefs::makeKeys,
	 "Return a vector/list of keys for an entire subspace of "
	 "subdimensional rank", py::arg("subkey"))
    .def("getGeometry",    &CoefClasses::Coefs::getGeometry,
	 "Return the coefficient geometry string")
    .def("getName",        &CoefClasses::Coefs::getName,
	 "Return the coefficient set nmenonic name")
    .def("setName",        &CoefClasses::Coefs::setName,
	 "Set or rename the coefficient set nmenonic name", py::arg("newname"))
    .def("deepcopy",       &CoefClasses::Coefs::deepcopy,
	 "Return a byte-by-byte copy of the original data")
    .def("zerodata",       &CoefClasses::Coefs::zerodata,
	 "Zero all of the coefficient data, keeping sizes and metadata intact")
    .def("CompareStanzas", &CoefClasses::Coefs::CompareStanzas,
	 "Check that the data in one Coefs set is identical to "
	 "that in another")
    .def_static("factory", &CoefClasses::Coefs::factory,
		"Deduce the type and read coefficients from a native or HDF5 file",
		py::arg("file"), py::arg("stride")=1,
		py::arg("tmin")=-std::numeric_limits<double>::max(),
		py::arg("tmax")= std::numeric_limits<double>::max())
    .def_static("makecoefs", &CoefClasses::Coefs::makecoefs,
		"Create a new coefficient instance compatible with the "
		"supplied coefficient structure. You still need to call "
		"add() to add the coefficient structure to the container",
		py::arg("coef"), py::arg("name")="");

  py::class_<CoefClasses::SphCoefs, std::shared_ptr<CoefClasses::SphCoefs>, PySphCoefs, CoefClasses::Coefs>(m, "SphCoefs", "Container for spherical coefficients")
    .def(py::init<bool>())
    .def("getAllCoefs",
	 [](CoefClasses::SphCoefs& A)
	 {
	   auto M = A.getAllCoefs(); // Need a copy here
	   py::array_t<std::complex<double>> ret = make_ndarray<std::complex<double>>(M);
	   return ret;
	 },
	 "Return a 3d ndarray index by spherical index, radial index and time index.  The spherical "
	 "index serializes all pairs of (l, m). The index for (l, m) is: l*(l+1)/2 + m");

  py::class_<CoefClasses::CylCoefs, std::shared_ptr<CoefClasses::CylCoefs>, PyCylCoefs, CoefClasses::Coefs>(m, "CylCoefs", "Container for cylindrical coefficients")
    .def(py::init<bool>())
    .def("getAllCoefs",
	 [](CoefClasses::CylCoefs& A)
	 {
	   auto M = A.getAllCoefs(); // Need a copy here
	   py::array_t<std::complex<double>> ret = make_ndarray<std::complex<double>>(M);
	   return ret;
	 },
	 "Return a 3d ndarray index by azimuthal index, radial index and time index")
    .def("EvenOddPower",
	 [](CoefClasses::CylCoefs& A, int nodd, int min, int max)
	 {
	   return A.EvenOddPower(nodd, min, max);
	 },
	 "Get cylindrical coefficient power separated into vertically even and odd\n"
	 "contributions.  The default parameters will query the YAML config for the\n"
	 "value of ncylodd but this can be provided as an argument if it is\n"
	 "not explicit in your EXP::Cylinder configuration.",
	 py::arg("nodd")=-1, py::arg("min")=0,
	 py::arg("max") = std::numeric_limits<double>::max());


  py::class_<CoefClasses::TableData, std::shared_ptr<CoefClasses::TableData>, PyTableData, CoefClasses::Coefs>(m, "TableData", "Container for simple data tables with multiple columns")
    .def(py::init<bool>())
    .def("getAllCoefs",    &CoefClasses::TableData::getAllCoefs,
	 "Return a 2d ndarray index by column and time");
}

