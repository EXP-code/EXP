#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <CoefContainer.H>

namespace py = pybind11;

#include "TensorToArray.H"

void CoefficientClasses(py::module &m) {

  m.doc() =
    R"(
    Coefficient class bindings

    These classes store, write, and provide an interface to coefficients
    and table data for use by the other pyEXP classes.

    CoefStruct
    ----------
    The CoefStruct class is low-level structure that stores the data
    and metadata specific to each geometry. There are three groups of
    CoefStruct derived classes for biorthogonal basis coefficients,
    field data coeffients, and auxiliary table data. The biorthogonal
    classes are spherical (SphStruct), cylindrical (CylStruct), slab
    (SlabStruct), and cube (CubeStruct).  The field classes
    cylindrical (CylFldStruct), and spherical (SphFldStruct).  The
    table data is stored in TblStruct.

    Instances of these structures represent individual times points
    and are created, maintained, and interfaced by the Coefs class.
    Access to the underlying data is provided to Python in case you
    need to change or rewrite the data for some reason.  We have also
    provided a assign() member so that you can instaniate and load a
    coefficient structure using Python.  To do this, use the
    constructor to make a blank instance, assign the dimensions and
    use assign() to create a data matrix with the supplied matrix or
    array.  The dimen- sions are:
     1. (lmax, nmax) for SphStruct
     2. (mmax, nmax) for a CylStruct
     3. (nmaxx, nmaxy, nmaxz) for a SlabStruct
     4. (nmaxx, nmaxy, nmaxz) for a CubeStruct
     5. (nfld, lmax, nmax) for a SphFldStruct
     6. (nfld, mmax, nmax) for a CylFldStruct
     7. (cols) for a TblStruct.

    Coefs
    -----
    The base class, 'Coefs', provides a factory reader that will
    create one of the derived coefficient classes, SphCoefs, CylCoefs,
    SlabCoefs, CubeCoefs, TblCoefs, SphFldCoefs, and CylFldCoefs,
    deducing the type from the input file. The input files may be EXP
    native or HDF5 cofficient files.  Only biorthgonal basis
    coefficients have a native EXP type.  The Basis factory,
    Basis::createCoefficients, will create set of coefficients from
    phase-space snapshots.  See help(pyEXP.basis). Files which are not
    recognized as EXP coefficient files are assumed to be data files
    and are parsed by the TblCoefs class. The first column in data
    tables is interpreted as time and each successive column is
    interpreted as a new data field.

    Once created, you may get a list of times, get the total
    gravitational power from biothogonal basis coefficients and
    general power from the field coefficients, and write a new HDF5
    file.  Their primary use is as a container object for MSSA (using
    expMSSA) and field visualization using the FieldGenerator class.

    Updates
    -------
    The expMSSA class will update the contribution to the coefficients
    specified by key from each eigen component to the reconstructed
    series. Unspecified coefficients series will not be updated and
    their original data will be intact. For visualization, the series
    data in a Coefs object may be zeroed using the 'zerodata()' member
    function prior to an expMSSA update.  This allows one to include
    reconstructions that *only* include particular eigen components for
    the coefficients specified by key.  Then, one can visualize only the
    updated fields using 'FieldGenerator'. See help(pyEXP.mssa) and
    help(pyEXP.field) for more details.

    Dataset indexing
    ----------------
    Coefficients and other auxilliary data from simulations are stored
    and retrieved by their time field.  Internally, these are floating
    fixed-point values truncated to 8 signficant figures so that they
    may be used as dictionary/map keys.  The values in the time list,
    returned with the Times() member function contain the truncated
    values for reference.

    Object lifetime
    ---------------
    As in native Python, the memory for created objects persists until
    it is no longer referenced.  For example, replacing the variable
    with a new set of coefficients will allow the memory to be
    deallocated if no other class instance holds a reference. Because
    coefficient sets can be large, the creation of a Coefs instance by
    the 'factory' will be passed to any other class that needs it.
    The MSSA class, expMSSA, will hold a reference to the the Coefs
    object passed on creation, and it update the values of the
    coefficients on reconstruction, without copying. If you want to
    keep the initial set without change, we have provided a
    'deepcopy()' member that provides a byte-by-byte copy.

  )";

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

    Eigen::VectorXcd& getData(double time) override {
      PYBIND11_OVERRIDE_PURE(Eigen::VectorXcd&, Coefs, getData, time);
    }

    void setData(double time, Eigen::VectorXcd& array) override {
      PYBIND11_OVERRIDE_PURE(void, Coefs, setData, time, array);
    }

    std::shared_ptr<CoefStruct> getCoefStruct(double time) override {
      PYBIND11_OVERRIDE_PURE(std::shared_ptr<CoefStruct>, Coefs, getCoefStruct, time);
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

    Eigen::VectorXcd& getData(double time) override {
      PYBIND11_OVERRIDE(Eigen::VectorXcd&, SphCoefs, getData, time);
    }

    void setData(double time, Eigen::VectorXcd& array) override {
      PYBIND11_OVERRIDE(void, SphCoefs, setData, time, array);
    }

    std::shared_ptr<CoefStruct> getCoefStruct(double time) override {
      PYBIND11_OVERRIDE(std::shared_ptr<CoefStruct>, SphCoefs, getCoefStruct,
			time);
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

    Eigen::VectorXcd& getData(double time) override {
      PYBIND11_OVERRIDE(Eigen::VectorXcd&, CylCoefs, getData, time);
    }

    void setData(double time, Eigen::VectorXcd& array) override {
      PYBIND11_OVERRIDE(void, CylCoefs, setData, time, array);
    }

    std::shared_ptr<CoefStruct> getCoefStruct(double time) override {
      PYBIND11_OVERRIDE(std::shared_ptr<CoefStruct>, CylCoefs, getCoefStruct,
			time);
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

  class PySlabCoefs : public SlabCoefs
  {
  protected:
    void readNativeCoefs(const std::string& file, int stride, double tmin, double tmax) override {
      PYBIND11_OVERRIDE(void, SlabCoefs, readNativeCoefs, file, stride, tmin, tmax);
    }

    std::string getYAML() override {
      PYBIND11_OVERRIDE(std::string, SlabCoefs, getYAML,);
    }

    void WriteH5Params(HighFive::File& file) override {
      PYBIND11_OVERRIDE(void, SlabCoefs, WriteH5Params, file);
    }

    unsigned WriteH5Times(HighFive::Group& group, unsigned count) override {
      PYBIND11_OVERRIDE(unsigned, SlabCoefs, WriteH5Times, group, count);
    }

  public:
    // Inherit the constructors
    using SlabCoefs::SlabCoefs;

    Eigen::VectorXcd& getData(double time) override {
      PYBIND11_OVERRIDE(Eigen::VectorXcd&, SlabCoefs, getData, time);
    }

    void setData(double time, Eigen::VectorXcd& array) override {
      PYBIND11_OVERRIDE(void, SlabCoefs, setData, time, array);
    }

    std::shared_ptr<CoefStruct> getCoefStruct(double time) override {
      PYBIND11_OVERRIDE(std::shared_ptr<CoefStruct>, SlabCoefs, getCoefStruct,
			time);
    }

    std::vector<double> Times() override {
      PYBIND11_OVERRIDE(std::vector<double>, SlabCoefs, Times,);
    }

    void WriteH5Coefs(const std::string& prefix) override {
      PYBIND11_OVERRIDE(void, SlabCoefs, WriteH5Coefs, prefix);
    }

    void ExtendH5Coefs(const std::string& prefix) override {
      PYBIND11_OVERRIDE(void, SlabCoefs, ExtendH5Coefs, prefix);
    }

    Eigen::MatrixXd& Power(int min, int max) override {
      PYBIND11_OVERRIDE(Eigen::MatrixXd&, SlabCoefs, Power, min, max);
    }

    bool CompareStanzas(std::shared_ptr<Coefs> check) override {
      PYBIND11_OVERRIDE(bool, SlabCoefs, CompareStanzas, check);
    }

    void clear() override {
      PYBIND11_OVERRIDE(void, SlabCoefs,	clear,);
    }

    void add(CoefStrPtr coef) override {
      PYBIND11_OVERRIDE(void, SlabCoefs,	add, coef);
    }

    std::vector<Key> makeKeys(Key k) override {
      PYBIND11_OVERRIDE(std::vector<Key>, SlabCoefs, makeKeys, k);
    }

    std::shared_ptr<Coefs> deepcopy() override {
      PYBIND11_OVERRIDE(std::shared_ptr<Coefs>, SlabCoefs, deepcopy,);
    }

    void zerodata() override {
      PYBIND11_OVERRIDE(void, SlabCoefs, zerodata,);
    }


  };

  class PyCubeCoefs : public CubeCoefs
  {
  protected:
    void readNativeCoefs(const std::string& file, int stride, double tmin, double tmax) override {
      PYBIND11_OVERRIDE(void, CubeCoefs, readNativeCoefs, file, stride, tmin, tmax);
    }

    std::string getYAML() override {
      PYBIND11_OVERRIDE(std::string, CubeCoefs, getYAML,);
    }

    void WriteH5Params(HighFive::File& file) override {
      PYBIND11_OVERRIDE(void, CubeCoefs, WriteH5Params, file);
    }

    unsigned WriteH5Times(HighFive::Group& group, unsigned count) override {
      PYBIND11_OVERRIDE(unsigned, CubeCoefs, WriteH5Times, group, count);
    }

  public:
    // Inherit the constructors
    using CubeCoefs::CubeCoefs;

    Eigen::VectorXcd& getData(double time) override {
      PYBIND11_OVERRIDE(Eigen::VectorXcd&, CubeCoefs, getData, time);
    }

    void setData(double time, Eigen::VectorXcd& array) override {
      PYBIND11_OVERRIDE(void, CubeCoefs, setData, time, array);
    }

    std::shared_ptr<CoefStruct> getCoefStruct(double time) override {
      PYBIND11_OVERRIDE(std::shared_ptr<CoefStruct>, CubeCoefs, getCoefStruct,
			time);
    }

    std::vector<double> Times() override {
      PYBIND11_OVERRIDE(std::vector<double>, CubeCoefs, Times,);
    }

    void WriteH5Coefs(const std::string& prefix) override {
      PYBIND11_OVERRIDE(void, CubeCoefs, WriteH5Coefs, prefix);
    }

    void ExtendH5Coefs(const std::string& prefix) override {
      PYBIND11_OVERRIDE(void, CubeCoefs, ExtendH5Coefs, prefix);
    }

    Eigen::MatrixXd& Power(int min, int max) override {
      PYBIND11_OVERRIDE(Eigen::MatrixXd&, CubeCoefs, Power, min, max);
    }

    bool CompareStanzas(std::shared_ptr<Coefs> check) override {
      PYBIND11_OVERRIDE(bool, CubeCoefs, CompareStanzas, check);
    }

    void clear() override {
      PYBIND11_OVERRIDE(void, CubeCoefs,	clear,);
    }

    void add(CoefStrPtr coef) override {
      PYBIND11_OVERRIDE(void, CubeCoefs,	add, coef);
    }

    std::vector<Key> makeKeys(Key k) override {
      PYBIND11_OVERRIDE(std::vector<Key>, CubeCoefs, makeKeys, k);
    }

    std::shared_ptr<Coefs> deepcopy() override {
      PYBIND11_OVERRIDE(std::shared_ptr<Coefs>, CubeCoefs, deepcopy,);
    }

    void zerodata() override {
      PYBIND11_OVERRIDE(void, CubeCoefs, zerodata,);
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

    Eigen::VectorXcd& getData(double time) override {
      PYBIND11_OVERRIDE(Eigen::VectorXcd&, TableData, getData, time);
    }

    void setData(double time, Eigen::VectorXcd& array) override {
      PYBIND11_OVERRIDE(void, TableData, setData, time, array);
    }

    std::shared_ptr<CoefStruct> getCoefStruct(double time) override {
      PYBIND11_OVERRIDE(std::shared_ptr<CoefStruct>, TableData, getCoefStruct,
			time);
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

  class PyTrajectoryData : public TrajectoryData
  {
  protected:
    void readNativeCoefs(const std::string& file, int stride, double tmin, double tmax) override {
      PYBIND11_OVERRIDE(void, TrajectoryData, readNativeCoefs, file, stride, tmin, tmax);
    }

    std::string getYAML() override {
      PYBIND11_OVERRIDE(std::string, TrajectoryData, getYAML,);
    }

    void WriteH5Params(HighFive::File& file) override {
      PYBIND11_OVERRIDE(void, TrajectoryData, WriteH5Params, file);
    }

    unsigned WriteH5Times(HighFive::Group& group, unsigned count) override {
      PYBIND11_OVERRIDE(unsigned, TrajectoryData, WriteH5Times, group, count);
    }

  public:
    // Inherit the constructors
    using TrajectoryData::TrajectoryData;

    Eigen::VectorXcd& getData(double time) override {
      PYBIND11_OVERRIDE(Eigen::VectorXcd&, TrajectoryData, getData, time);
    }

    void setData(double time, Eigen::VectorXcd& array) override {
      PYBIND11_OVERRIDE(void, TrajectoryData, setData, time, array);
    }

    std::shared_ptr<CoefStruct> getCoefStruct(double time) override {
      PYBIND11_OVERRIDE(std::shared_ptr<CoefStruct>, TrajectoryData, getCoefStruct,
			time);
    }

    std::vector<double> Times() override {
      PYBIND11_OVERRIDE(std::vector<double>, TrajectoryData, Times,);
    }

    void WriteH5Coefs(const std::string& prefix) override {
      PYBIND11_OVERRIDE(void, TrajectoryData, WriteH5Coefs, prefix);
    }

    void ExtendH5Coefs(const std::string& prefix) override {
      PYBIND11_OVERRIDE(void, TrajectoryData, ExtendH5Coefs, prefix);
    }

    Eigen::MatrixXd& Power(int min, int max) override {
      PYBIND11_OVERRIDE(Eigen::MatrixXd&, TrajectoryData, Power, min, max);
    }

    bool CompareStanzas(std::shared_ptr<Coefs> check) override {
      PYBIND11_OVERRIDE(bool, TrajectoryData, CompareStanzas, check);
    }

    void clear() override {
      PYBIND11_OVERRIDE(void, TrajectoryData, clear,);
    }

    void add(CoefStrPtr coef) override {
      PYBIND11_OVERRIDE(void, TrajectoryData, add, coef);
    }

    std::vector<Key> makeKeys(Key k) override {
      PYBIND11_OVERRIDE(std::vector<Key>, TrajectoryData, makeKeys, k);
    }

    std::shared_ptr<Coefs> deepcopy() override {
      PYBIND11_OVERRIDE(std::shared_ptr<Coefs>, TrajectoryData, deepcopy,);
    }

    void zerodata() override {
      PYBIND11_OVERRIDE(void, TrajectoryData, zerodata,);
    }

  };

  py::class_<CoefClasses::CoefStruct, std::shared_ptr<CoefClasses::CoefStruct>, PyCoefStruct>
    (m, "CoefStruct")
    .def(py::init<>(),
	 R"(
         Base class coefficient data structure object

         Returns
         -------
         CoefStruct
         )")
    .def("create", &CoefStruct::create,
	 R"(
        Initialize a coefficient zeroed structure from user supplied dimensions

        Returns
        -------
        None
        )")
    .def("deepcopy", &CoefStruct::deepcopy,
        R"(
        Make a new instance and copy all data into the new instance

        Returns
        -------
        CoefStruct
            new CoefStruct instance with the copied data

        Notes
        -----
        This is useful if you would like to modify some coefficients
        while preserving your original coefficients.
        )")
    .def_readonly("geometry", &CoefStruct::geom,
		  R"(
                  str
                      geometry type
                  )")
    .def_readonly("time", &CoefStruct::time,
		   R"(
                   float
                       data's time stamp
                   )")
    .def_readonly("center", &CoefStruct::ctr,
    R"(
                float
                    data's center value
                )") 
    .def("getCoefTime", &CoefStruct::getTime,
        R"(
        Read-only access to the coefficient time

        Returns
        -------
        numpy.ndarray
            vector of center data

        See also
        --------
        setCoefTime : read-write access to the coefficient time
        )")
    .def("setCoefTime",
        static_cast<void (CoefStruct::*)(double&)>(&CoefStruct::setTime),
        py::arg("tval"),
        R"(
        Set the coefficient time

        Parameters
        ----------
        tval  : float
                time value

        Returns
        -------
        None

        Notes
        -----

        See also
        --------
        getCenter : read-only access to coefficient time
        )")
    .def("getCoefCenter", &CoefStruct::getCenter,
          R"(
          Read-only access to the center data
  
          Returns
          -------
          numpy.ndarray
              vector of center data
  
          See also
          --------
          setCenter : read-write access to the center data
          )")
    .def("setGravConstant", &CoefStruct::setGravConstant,
          py::arg("G"),
          R"(
          Set the gravitational constant
  
          Parameters
          ----------
          G  : float
             gravitational constant, default is 1.0
  
          Returns
          -------
          None
  
          Notes
          -----
          The gravitational constant is used for field evaluation for
          biorthogonal basis sets.  It will be set automatically when
          reading EXP coefficient files.
          )")
    .def("setCoefCenter",
          static_cast<void (CoefStruct::*)(std::vector<double>&)>(&CoefStruct::setCenter),
          py::arg("mat"),
          R"(
          Set the center vector
  
          Parameters
          ----------
          mat  : numpy.ndarray
                center vector
  
          Returns
          -------
          None
  
          Notes
          -----
  
          See also
          --------
          getCenter : read-only access to center data
          )")
    .def("getCoefs", &CoefStruct::getCoefs,
        R"(
        Read-only access to the underlying data store

        Returns
        -------
        numpy.ndarray
            complex-valued matrix as a flattened NumPy array of complex values

        See also
        --------
        setCoefs : read-write access to Coefs
        )")
  .def("setCoefs",		// Member function overload
        static_cast<void (CoefStruct::*)(Eigen::VectorXcd&)>(&CoefStruct::setCoefs),
        py::arg("mat"),
        R"(
        Set the coefficient matrix with the coefficient vector in the same form as returned
        by getCoefs

        Parameters
        ----------
        mat  : numpy.ndarray
             Flattened array of coefficients

        Returns
        -------
        None

        Notes
        -----
        The rank data array must match the rank of the CoefStruct.  Use getCoefs to create
        such an array with the correct rank.

        See also
        --------
        getCoefs : read-only access to Coefs
        )")
  .def("setCoefs",		// Member function overload
        static_cast<Eigen::Ref<Eigen::VectorXcd>(CoefStruct::*)()>(&CoefStruct::setCoefs),
        R"(
        Read-write access to the underlying data store

        Returns
        -------
        numpy.ndarray
            reference to a complex-valued matrix represented as a NumPy array of complex
            values

        Notes
        -----
        Changes made to the data array will be automatically mapped back to the C++
        CoefStruct instance.  You may use the setCoefs(array) call to set the data array
        directly.

        See also
        --------
        getCoefs : read-only access to Coefs
        )");


  py::class_<CoefClasses::SphStruct, std::shared_ptr<CoefClasses::SphStruct>, CoefStruct>
    (m, "SphStruct")
    .def(py::init<>(), "Spherical coefficient data structure object")
    .def("assign", &SphStruct::assign,
	      R"(
        Assign a coefficient matrix to CoefStruct.

        Parameters
        ----------
        mat  : numpy.ndarray
             Matrix of coefficients
        lmax : int
             angular order
        nmax : int
             radial order

        Returns
        -------
        None
        )");

  py::class_<CoefClasses::CylStruct, std::shared_ptr<CoefClasses::CylStruct>, CoefStruct>
    (m, "CylStruct")
    .def(py::init<>(), "Cylindrical coefficient data structure object")
    .def("assign", &CylStruct::assign,
	      R"(
        Assign a coefficient matrix to CoefStruct.

        Parameters
        ----------
        mat  : numpy.ndarray
             Matrix of coefficients
        mmax : int
             angular order
        nmax : int
             radial order

        Returns
        -------
        None
        )");

  py::class_<CoefClasses::SlabStruct, std::shared_ptr<CoefClasses::SlabStruct>, CoefStruct>
    (m, "SlabStruct")
    .def(py::init<>(), "Slab coefficient data structure object")
    .def("assign", &SlabStruct::assign,
	      R"(
        Assign a coefficient matrix to CoefStruct.

        Parameters
        ----------
        mat  : numpy.ndarray, complex
             complex-valued NumPy tensor of coefficient values

        Returns
        -------
        None

        Notes
        -----
        The dimensions are inferred from the 3-dimensional NumPy array
        (tensor)
        )");

  py::class_<CoefClasses::CubeStruct, std::shared_ptr<CoefClasses::CubeStruct>, CoefStruct>
    (m, "CubeStruct")
    .def(py::init<>(), "Cube coefficient data structure object")
    .def("assign", &CubeStruct::assign,
	      R"(
        Assign a coefficient matrix to CoefStruct.

        Parameters
        ----------
        mat  : numpy.ndarray, complex
             complex-valued NumPy tensor of coefficient values

        Returns
        -------
        None

        Notes
        -----
        The dimensions are inferred from the 3-dimensional NumPy array
        (tensor)
        )");

  py::class_<CoefClasses::TblStruct, std::shared_ptr<CoefClasses::TblStruct>, CoefStruct>
    (m, "TblStruct")
    .def(py::init<>(), "Multicolumn table data structure object")
    .def("assign", &TblStruct::assign,
	      R"(
        Assign a coefficient matrix to CoefStruct.

        Parameters
        ----------
        mat  : numpy.ndarray, complex
             complex-valued NumPy array of table values

        Returns
        -------
        None
        )");

  py::class_<CoefClasses::SphFldStruct, std::shared_ptr<CoefClasses::SphFldStruct>, CoefStruct>
    (m, "SphFldStruct")
    .def(py::init<>(), "Spherical field coefficient data structure object")
    .def("assign", &SphFldStruct::assign,
	R"(
        Assign a coefficient matrix to CoefStruct.

        Parameters
        ----------
        mat  : numpy.ndarray
             Flattened array of coefficients
        nfld : int
             number of data fields
        lmax : int
             angular order
        nmax : int
             radial order

        Returns
        -------
        None
        )");

  py::class_<CoefClasses::CylFldStruct, std::shared_ptr<CoefClasses::CylFldStruct>, CoefStruct>
    (m, "CylFldStruct")
    .def(py::init<>(), "Cylindrical field coefficient data structure object")
    .def("assign", &CylFldStruct::assign,
	      R"(
        Assign a flattened coefficient array to CylFldStruct.

        Parameters
        ----------
        mat  : numpy.ndarray
             Flattened array of coefficients
        nfld : int
             number of data fields
        mmax : int
             angular order
        nmax : int
             radial order

        Returns
        -------
        None
        )");

  py::class_<CoefClasses::Coefs, std::shared_ptr<CoefClasses::Coefs>, PyCoefs>
    (m, "Coefs")
    .def(py::init<std::string, bool>(),
         R"(
         Create a coefficient container

         Parameters
         ----------
         type : str
             type of coefficient container
         verbose : bool
             display verbose information.

         Returns
         -------
         Coefs instance

         Notes
         -----
         This container that holds, stores, and reads coefficient
         structures for a part or all of the snapshots in your
         simulation
         )",
         py::arg("type"),
         py::arg("verbose"))
    .def("__call__",
	 &CoefClasses::Coefs::getData,
         R"(
         Return the flattened coefficient structure for the desired time.

         Parameters
         ----------
         time : float
             the desired time

         Returns
         -------
         numpy.ndarray
             the flattened coefficient array at the requested time

         Notes
         -----
         This operator will return the 0-rank array if no coefficients
         are found at the requested time
         )",
         py::arg("time"))
    .def("setData",
         &CoefClasses::Coefs::setData,
         R"(
         Enter and/or rewrite the flattened coefficient array at the provided time

         Parameters
         ----------
         time : float
             snapshot time corresponding to the the coefficient matrix
             mat : numpy.ndarray
                 the new coefficient array.

         Returns
         -------
         None
         )",py::arg("time"), py::arg("array"))
    .def("add",
         &CoefClasses::Coefs::add,
         R"(
         Add a coefficient structure to the coefficient container

         Parameters
         ----------
         coef : CoefStruct
             coefficient structure to add

         Returns
         -------
         None

         Notes
         -----
         The time value is supplied by the CoefStruct

         See also
         --------
         CoefStruct : the coefficient structure for a single snaphot
         )",py::arg("coef"))
    .def("getCoefStruct",
         &CoefClasses::Coefs::getCoefStruct,
         R"(
         Return the coefficient structure for the desired time

         Parameters
         ----------
         time : float
             requested time

         Returns
         -------
         CoefStruct: coefficient structure

         Notes
         -----
         You will get a runtime error if the entry does not exist.
         )",py::arg("time"))
    .def("Times",
            &CoefClasses::Coefs::Times,
            R"(
            Return a list of times for coefficient sets currently in the container

            Returns
            -------
            list(float,...)
                list of times
            )")
    .def("setUnit",
            &CoefClasses::Coefs::setUnit,
            R"(
            Set the units for the coefficient structure.

            Parameters
            ----------
            name : str
               the name of physical quantity (G, Length, Mass, Time, etc)
            unit : str
               the unit string (scalar, mixed, kpc, Msun, Myr, km/s etc.).
               This field is optional and can be empty.
            value : float
	       the default value of the multiples of the unit

            Returns
            -------
            None
            )", py::arg("name"), py::arg("unit")="", py::arg("value")=1.0)
    .def("WriteH5Coefs",
	 [](CoefClasses::Coefs& self, const std::string& filename) {
	   if (self.getUnits().size()==1) {
	     std::cout << "Coefs::WriteH5Coefs: please set units for your coefficient set using the `setUnit()` member," << std::endl
		       << "                     one for each unit.  We suggest explicitly setting 'G', 'Length', 'Mass'," << std::endl
		       << "                     'Time', and optionally 'Velocity' before writing HDF5 coefficients" << std::endl;
	   }
	   self.WriteH5Coefs(filename);
	 },
	 R"(
            Write the coefficients into an EXP HDF5 coefficient file with the given prefix name.

            Parameters
            ----------
            filename : str
                the filename prefix.

            Returns
            -------
            None

            Notes
            -----
            This call will throw a runtime exception of the HDF5
            coefficient file already exists.  This is a safety
            feature.  If you'd like a new version of this file, delete
            the old before this call.
            )",
	 py::arg("filename"))
    .def("ExtendH5Coefs",
            &CoefClasses::Coefs::ExtendH5Coefs,
            R"(
            Extend an existing EXP HDF5 coefficient file with added data

            Parameters
            ----------
            filename : str 
                the filename prefix

            Returns
            -------
            None

            Notes
            -----
            You will get a runtime error if the H5 filename does not exist
            )",py::arg("filename"))
    .def("Power",
             &CoefClasses::Coefs::Power,
             R"(
             Table of the full power for the top-level harmonic index as a function of time

             Parameters
             ----------
             min : int, default=0
                 the minimum harmonic index
             max : int, default=inf
                 the maximum harmonic index

             Returns
             -------
             numpy.ndarray: 
                 table of coefficient power values
             )",py::arg("min")=0, py::arg("max")=std::numeric_limits<int>::max())
    .def("makeKeys",
         &CoefClasses::Coefs::makeKeys,
         R"(
         a vector/list of keys for an entire subspace of subdimensional rank

         Parameters
         ----------
         subkey : list(int,...)
             the subkey

         Returns
         -------
         list(list)
             list of keys with the provided subkey prefix

         Notes
         -----
         The subkey prefix is a list of ints with size smaller than
         the key.  The values will be used as leading key indices and
         all the matching trailing key values will be generated.
         )", py::arg("subkey"))
    .def("getGeometry",    &CoefClasses::Coefs::getGeometry,
         R"(
         The coefficient geometry string

         Returns
         -------
         str
             geometry name
         )")
    .def("getName",        &CoefClasses::Coefs::getName,
         R"(
         The coefficient set mnemonic name.

         Returns
         -------
         str
             mnemonic name
         )")
    .def("setName",        &CoefClasses::Coefs::setName,
         R"(
         Set or rename the coefficient set mnemonic name.

         Parameters
         ----------
         newname : str
             new mnemonic name

         Returns
         -------
         None
         )", py::arg("newname"))
    .def("deepcopy",       &CoefClasses::Coefs::deepcopy,
         R"(
         byte-by-byte copy of the original data

         Returns
         -------
         Coefs
             New copied instance of the Coefs object

         Notes
         -----
         Useful if you would like to change your coefficients or
         filter using mSSA but keep a copy of your original
         coefficient db for comparison

         See also
         --------
         pyEXP.mssa
         )")
    .def("zerodata",       &CoefClasses::Coefs::zerodata,
         R"(
         Zero all of the coefficient data, keeping sizes and metadata intact.

         Returns
         -------
         None
         )")
    .def("CompareStanzas", &CoefClasses::Coefs::CompareStanzas,
         R"(
         Check that the data in one Coefs set is identical to that in another

         Returns
         -------
         bool
             True if the data is identical, False otherwise.
         )")
    .def("getUnits", &CoefClasses::Coefs::getUnits,
         R"(
         Get the units of the coefficient data

         Returns
         -------
         list((str,str,float))
             list of
         )")
    .def_static("factory", &CoefClasses::Coefs::factory,
              R"(
              Deduce the type and read coefficients from a native or HDF5 file

              Parameters
              ----------
              file : str
                  the file path.
              stride : int, default=1
                  stride value
              tmin : float, default=-inf
                   minimum time value
              tmax : float, default=inf
                   maximum time value

            Returns
            -------
            Coefs
                the newly created Coefs object
            )",
            py::arg("file"), py::arg("stride")=1,
            py::arg("tmin")=-std::numeric_limits<double>::max(),
            py::arg("tmax")= std::numeric_limits<double>::max())
    .def_static("makecoefs", &CoefClasses::Coefs::makecoefs,
		R"(
                make a new coefficient container instance compatible

                Parameters
                ----------
                coef : CoefStruct
                    coefficient structure
                name : str, default=""
                     name of the coefficient instance

                Returns
                -------
                Coefs 
                    the newly created Coefs object

                Notes
                -----
                The type will be deduced from the supplied coefficient
                structure.  Subsequent additions of coefficients sets
                should use the addcoef() member

                See also
                --------
                addcoef : add coefficient structures to an existing coefficieint container
                )",
		py::arg("coef"), py::arg("name")="");


  py::class_<CoefClasses::SphCoefs, std::shared_ptr<CoefClasses::SphCoefs>, PySphCoefs, CoefClasses::Coefs>
    (m, "SphCoefs", "Container for spherical coefficients")
    .def(py::init<bool>(),
	 R"(
         Construct a null SphCoefs object

         Parameters
         ----------
         verbose : bool
             display verbose information.

         Returns
         -------
         SphCoefs instance
         )")
    .def("__call__",
	 &CoefClasses::SphCoefs::getMatrix,
         R"(
         Return the coefficient Matrix for the desired time.

         Parameters
         ----------
         time : float
             the desired time

         Returns
         -------
         numpy.ndarray
             the coefficient Matrix at the requested time

         Notes
         -----
         This operator will return the 0-rank matrix if no
         coefficients are found at the requested time
         )",
         py::arg("time"))
    .def("setMatrix",
	 &CoefClasses::SphCoefs::setMatrix,
         R"(
         Enter and/or rewrite the coefficient matrix at the provided time

         Parameters
         ----------
         time : float
             snapshot time corresponding to the the coefficient matrix
             mat : numpy.ndarray
                 the new coefficient array.

         Returns
         -------
         None
         )",
         py::arg("time"), py::arg("mat"))
    .def("getAllCoefs",
	 [](CoefClasses::SphCoefs& A)
	 {
	   Eigen::Tensor<std::complex<double>, 3> M = A.getAllCoefs(); // Need a copy here
	   py::array_t<std::complex<double>> ret = make_ndarray3<std::complex<double>>(M);
	   return ret;
	 },
	 R"(
        Provide a 3-dimensional ndarray indexed by spherical index, radial index,
        and time index

        Returns
        -------
        numpy.ndarray
            3-dimensional numpy array containing the spherical coefficients

        Notes
        -----
        The spherical index serializes all pairs of (l, m). The index
        for (l, m) is calculated as: l*(l+1)/2 + m.
        )");

  py::class_<CoefClasses::CylCoefs, std::shared_ptr<CoefClasses::CylCoefs>, PyCylCoefs, CoefClasses::Coefs>
    (m, "CylCoefs", "Container for cylindrical coefficients")
    .def(py::init<bool>(),
	 R"(
         Construct a null CylCoefs object

         Parameters
         ----------
         verbose : bool
             display verbose information.

         Returns
         -------
         CylCoefs instance
         )")
    .def("__call__",
	 &CoefClasses::CylCoefs::getMatrix,
         R"(
         Return the coefficient Matrix for the desired time.

         Parameters
         ----------
         time : float
             the desired time

         Returns
         -------
         numpy.ndarray
             the coefficient Matrix at the requested time

         Notes
         -----
         This operator will return the 0-rank matrix if no
         coefficients are found at the requested time
         )",
         py::arg("time"))
    .def("setMatrix",
	 &CoefClasses::CylCoefs::setMatrix,
         R"(
         Enter and/or rewrite the coefficient matrix at the provided time

         Parameters
         ----------
         time : float
             snapshot time corresponding to the the coefficient matrix
             mat : numpy.ndarray
                 the new coefficient array.

         Returns
         -------
         None
         )",
         py::arg("time"), py::arg("mat"))
    .def("getAllCoefs",
	 [](CoefClasses::CylCoefs& A)
	 {
	   Eigen::Tensor<std::complex<double>, 3> M = A.getAllCoefs(); // Need a copy here
	   py::array_t<std::complex<double>> ret = make_ndarray3<std::complex<double>>(M);
	   return ret;
	 },
	 R"(
         Provide a 3-dimensional ndarray indexed by azimuthal index, radial index, and time index

         Returns
         -------
         numpy.ndarray
             3-dimensional numpy array containing the cylindrical coefficients
         )")
    .def("EvenOddPower",
	 [](CoefClasses::CylCoefs& A, int nodd, int min, int max)
	 {
	   return A.EvenOddPower(nodd, min, max);
	 },
	 R"(
         Get cylindrical coefficient power separated into vertically even and odd contributions.

         Parameters
         ----------
         nodd : int, default=-1
             number of odd vertical modes to compute
         min : int, default=0
             minimum time index for power calculation
         max : int, default=inf
             maximum time index for power calculation

         Returns
         -------
         numpy.ndarray
             2-dimensional numpy array containing the even and odd
             power coefficients

         Notes
         -----
         The default parameters (nodd<0) will query the YAML config
         for the value of ncylodd, but this can be provided as an
         argument if it is not explicitly set in your EXP::Cylinder
         configuration. If in doubt, use the default.
         )",
	 py::arg("nodd")=-1, py::arg("min")=0,
	 py::arg("max")=std::numeric_limits<int>::max());


  py::class_<CoefClasses::SphFldCoefs, std::shared_ptr<CoefClasses::SphFldCoefs>, CoefClasses::Coefs>
    (m, "SphFldCoefs", "Container for spherical field coefficients")
    .def(py::init<bool>(),
	       R"(
         Construct a null SphFldCoefs object

         Parameters
         ----------
         verbose : bool
             display verbose information.

         Returns
         -------
         SphFldCoefs instance
         )")
    .def("__call__",
	 [](CoefClasses::SphFldCoefs& A, double time)
	 {
	   // Need a copy here
	   auto M = A.getMatrix(time);
	   return make_ndarray3<std::complex<double>>(M);
	 },
         R"(
         Return the coefficient tensor for the desired time.

         Parameters
         ----------
         time : float
             the desired time

         Returns
         -------
         numpy.ndarray
             the coefficient tensor at the requested time

         Notes
         -----
         This operator will return the 0-rank tensor if no
         coefficients are found at the requested time
         )",
         py::arg("time"))
    .def("setMatrix",
	 [](CoefClasses::SphFldCoefs& A, double time,
	    py::array_t<std::complex<double>> mat)
	 {
	   auto M = make_tensor3<std::complex<double>>(mat);
	   A.setMatrix(time, M);
	 },
         R"(
         Enter and/or rewrite the coefficient tensor at the provided time

         Parameters
         ----------
         time : float
             snapshot time corresponding to the the coefficient matrix
         mat : numpy.ndarray
             the new coefficient array.

         Returns
         -------
         None
         )",
         py::arg("time"), py::arg("mat"))
    .def("getAllCoefs",
	 [](CoefClasses::SphFldCoefs& A)
	 {
	   Eigen::Tensor<std::complex<double>, 4> M = A.getAllCoefs(); // Need a copy here
	   py::array_t<std::complex<double>> ret = make_ndarray4<std::complex<double>>(M);
	   return ret;
	 },
	 R"(
        Provide a 4-dimensional ndarray indexed by channel index, spherical index, radial index, and time index

        Returns
        -------
        numpy.ndarray
            4-dimensional numpy array containing the spherical coefficients

        Notes
        -----
        The spherical index serializes all pairs of (l, m) where l, m
        are the aximuthal indices. The index for (l, m) pair is
        calculated as: l*(l+1)/2 + m
        )");

  py::class_<CoefClasses::CylFldCoefs, std::shared_ptr<CoefClasses::CylFldCoefs>, CoefClasses::Coefs>
    (m, "CylFldCoefs", "Container for cylindrical field coefficients")
    .def(py::init<bool>(),
	 R"(
         Construct a null CylFldCoefs object

         Parameters
         ----------
         verbose : bool
             display verbose information.

         Returns
         -------
         CylFldCoefs instance
         )")
    .def("__call__",
	 [](CoefClasses::CylFldCoefs& A, double time)
	 {
	   auto M = A.getMatrix(time); // Need a copy here
	   return make_ndarray3<std::complex<double>>(M);
	 },
         R"(
         Return the coefficient tensor for the desired time.

         Parameters
         ----------
         time : float
             the desired time

         Returns
         -------
         numpy.ndarray
             the coefficient tensor at the requested time

         Notes
         -----
         This operator will return the 0-rank tensor if no
         coefficients are found at the requested time
         )",
         py::arg("time"))
    .def("setMatrix",
	 [](CoefClasses::CylFldCoefs& A, double time,
	    py::array_t<std::complex<double>> mat)
	 {
	   auto M = make_tensor3<std::complex<double>>(mat);
	   A.setMatrix(time, M);
	 },
         R"(
         Enter and/or rewrite the coefficient tensor at the provided time

         Parameters
         ----------
         time : float
             snapshot time corresponding to the the coefficient matrix
         mat : numpy.ndarray
             the new coefficient array.

         Returns
         -------
         None
         )",
         py::arg("time"), py::arg("mat"))
    .def("getAllCoefs",
	 [](CoefClasses::CylFldCoefs& A)
	 {
	   Eigen::Tensor<std::complex<double>, 4> M = A.getAllCoefs(); // Need a copy here
	   py::array_t<std::complex<double>> ret = make_ndarray4<std::complex<double>>(M);
	   return ret;
	 },
	 R"(
        Provide a 4-dimensional ndarray indexed by channel index, spherical index, radial index, and time index

        Returns
        -------
        numpy.ndarray
            4-dimensional numpy array containing the cylindrical coefficients
        )");


  py::class_<CoefClasses::SlabCoefs, std::shared_ptr<CoefClasses::SlabCoefs>, PySlabCoefs, CoefClasses::Coefs>
    (m, "SlabCoefs", "Container for cube coefficients")
    .def(py::init<bool>(),
	       R"(
         Construct a null SlabCoefs object

         Parameters
         ----------
         verbose : bool
             display verbose information.

         Returns
         -------
         SlabCoefs instance
         )")
    .def("__call__",
	 &CoefClasses::SlabCoefs::getTensor,
         R"(
         Return the coefficient tensor for the desired time.

         Parameters
         ----------
         time : float
             the desired time

         Returns
         -------
         numpy.ndarray
             the coefficient Matrix at the requested time

         Notes
         -----
         This operator will return the 0-rank tensor if no
         coefficients are found at the requested time
         )",
         py::arg("time"))
    .def("setTensor",
	 &CoefClasses::SlabCoefs::setTensor,
         R"(
         Enter and/or rewrite the coefficient tensor at the provided time

         Parameters
         ----------
         time : float
             snapshot time corresponding to the the coefficient tensor
             mat : numpy.ndarray
                 the new coefficient array.

         Returns
         -------
         None
         )",
         py::arg("time"), py::arg("tensor"))
    .def("getAllCoefs",
	 [](CoefClasses::SlabCoefs& A)
	 {
	   Eigen::Tensor<std::complex<double>, 4> M = A.getAllCoefs(); // Need a copy here
	   py::array_t<std::complex<double>> ret = make_ndarray4<std::complex<double>>(M);
	   return ret;
	 },
	 R"(
         Provide a 4-dimensional ndarray indexed by nx, ny, nz, and time indices.

         Returns
         -------
         numpy.ndarray
             4-dimensional numpy array containing the slab coefficients
         )")
    .def("PowerDim",
	 [](CoefClasses::SlabCoefs& A, std::string d, int min, int max)
	 {
	   return A.Power(d[0], min, max);
	 },
	 R"(
         Get power for the coefficient DB as a function of harmonic index for a
         given dimension.  This Power() member is equivalent to PowerDim('x').

         Parameters
         ----------
         d    : char
            dimension for power summary; one of 'x', 'y', or 'z'
         min  : int
            minimum index along requested dimension (default=0)
         max  : int
            maximum index along requested dimension (default=max int)

         Returns
         -------
         numpy.ndarray
             2-dimensional numpy array containing the power
         )", py::arg("d"), py::arg("min")=0,
	 py::arg("max")=std::numeric_limits<int>::max());


  py::class_<CoefClasses::CubeCoefs, std::shared_ptr<CoefClasses::CubeCoefs>, PyCubeCoefs, CoefClasses::Coefs>
    (m, "CubeCoefs", "Container for cube coefficients")
    .def(py::init<bool>(),
	       R"(
         Construct a null CubeCoefs object

         Parameters
         ----------
         verbose : bool
             display verbose information.

         Returns
         -------
         CubeCoefs instance
         )")
    .def("__call__",
	 &CoefClasses::CubeCoefs::getTensor,
         R"(
         Return the coefficient tensor for the desired time.

         Parameters
         ----------
         time : float
             the desired time

         Returns
         -------
         numpy.ndarray
             the coefficient Matrix at the requested time

         Notes
         -----
         This operator will return the 0-rank tensor if no
         coefficients are found at the requested time
         )",
         py::arg("time"))
    .def("setTensor",
	 &CoefClasses::CubeCoefs::setTensor,
         R"(
         Enter and/or rewrite the coefficient tensor at the provided time

         Parameters
         ----------
         time : float
             snapshot time corresponding to the the coefficient tensor
             mat : numpy.ndarray
                 the new coefficient array.

         Returns
         -------
         None
         )",
         py::arg("time"), py::arg("tensor"))
    .def("getAllCoefs",
	 [](CoefClasses::CubeCoefs& A)
	 {
	   Eigen::Tensor<std::complex<double>, 4> M = A.getAllCoefs(); // Need a copy here
	   py::array_t<std::complex<double>> ret = make_ndarray4<std::complex<double>>(M);
	   return ret;
	 },
	 R"(
         Provide a 4-dimensional ndarray indexed by nx, ny, nz, and time indices.

         Returns
         -------
         numpy.ndarray
             4-dimensional numpy array containing the cube coefficients
         )")
    .def("PowerDim",
	 [](CoefClasses::CubeCoefs& A, std::string d, int min, int max)
	 {
	   return A.Power(d[0], min, max);
	 },
	 R"(
         Get power for the coefficient DB as a function of harmonic index for a
         given dimension.  This Power() member is equivalent to PowerDim('x').

         Parameters
         ----------
         d    : char
            dimension for power summary; one of 'x', 'y', or 'z'
         min  : int
            minimum index along requested dimension (default=0)
         max  : int
            maximum index along requested dimension (default=max int)

         Returns
         -------
         numpy.ndarray
             2-dimensional numpy array containing the power
         )", py::arg("d"), py::arg("min")=0,
	 py::arg("max")=std::numeric_limits<int>::max());

  py::class_<CoefClasses::TableData, std::shared_ptr<CoefClasses::TableData>, PyTableData, CoefClasses::Coefs>
    (m, "TableData", "Container for simple data tables with multiple columns")
    .def(py::init<bool>(),
	 R"(
         Construct a null TableData object

         Parameters
         ----------
         verbose : bool
             display verbose information.

         Returns
         -------
         TableData instance
         )", py::arg("verbose")=true)
    .def(py::init<std::string&>(),
	 R"(
         Construct a TableData object from a data file

         Parameters
         ----------
         type : str
             ascii table data file

         Returns
         -------
         TableData instance
         )")
    .def(py::init<std::string&, bool>(),
	 R"(
         Construct a TableData object from a data file

         Parameters
         ----------
         type : str
             ascii table data file
         verbose : bool
             display verbose information.

         Returns
         -------
         TableData instance
         )", py::arg("filename"), py::arg("verbose")=true)
    .def(py::init<std::vector<double>&, std::vector<std::vector<double>>&, bool>(),
	 R"(
         Construct a TableData object from data arrays

         Parameters
         ----------
         time : ndarray
             time data
         array : ndarray
             data columns
         verbose : bool
             display verbose information.

         Returns
         -------
         TableData instance
         )", py::arg("time"), py::arg("array"), py::arg("verbose")=true)
    .def("getAllCoefs", &CoefClasses::TableData::getAllCoefs,
	 R"(
         Return a 2-dimensional ndarray indexed by column and time

         Returns
         -------
         numpy.ndarray
             2-dimensional numpy array containing the data table
         )");

  py::class_<CoefClasses::TrajectoryData, std::shared_ptr<CoefClasses::TrajectoryData>, PyTrajectoryData, CoefClasses::Coefs>
    (m, "TrajectoryData", "Container for trajectory/orbit data")
    .def(py::init<bool>(),
	 R"(
         Construct a null TrajectoryData object

         Parameters
         ----------
         verbose : bool
             display verbose information.

         Returns
         -------
         TrajectoryData instance
         )", py::arg("verbose")=true)
    .def(py::init<std::string&>(),
	 R"(
         Construct a TrajectoryData object from a data file

         Parameters
         ----------
         type : str
             ascii table data file

         Returns
         -------
         TrajectoryData instance
         )")
    .def(py::init<std::string&, bool>(),
	 R"(
         Construct a TrajectoryData object from a data file

         Parameters
         ----------
         type : str
             ascii table data file
         verbose : bool
             display verbose information.

         Returns
         -------
         TrajectoryData instance
         )", py::arg("filename"), py::arg("verbose")=true)
    .def(py::init<std::vector<double>&, std::vector<Eigen::MatrixXd>&, bool>(),
	 R"(
         Construct a TrajectoryData object from data arrays

         Parameters
         ----------
         time : ndarray
             time data
         array : ndarray
             data columns
         verbose : bool
             display verbose information.

         Returns
         -------
         TrajectoryData instance
         )", py::arg("time"), py::arg("array"), py::arg("verbose")=true)
    .def("getAllCoefs",
	 [](CoefClasses::TrajectoryData& A)
	 {
	   Eigen::Tensor<double, 3> M = A.getAllCoefs(); // Need a copy here
	   py::array_t<double> ret = make_ndarray3<double>(M);
	   return ret;
	 },

	 R"(
         Return a 3-dimensional ndarray indexed by column and time

         Returns
         -------
         numpy.ndarray
             2-dimensional numpy array containing the data table
         )");
}
