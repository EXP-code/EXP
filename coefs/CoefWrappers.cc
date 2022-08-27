#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <CoefFactory.H>

namespace py = pybind11;

void CoefFactoryClasses(py::module &m) {

  m.doc() = "CoefFactory class bindings";

  using namespace Coefs;

  class PyCoefs : public Coefs
  {
  protected:
    void readNativeCoefs(const std::string& file) override {
      PYBIND11_OVERRIDE_PURE(void, Coefs, readNativeCoefs, file);
    }
    
    std::string getYAML() override {
      PYBIND11_OVERRIDE_PURE(std::string, Coefs, getYaml,);
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

    void add(CoefStrPtr coef) override {
      PYBIND11_OVERRIDE_PURE(void, Coefs, add, coef);
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

    void add(CoefStrPtr coef) override {
      PYBIND11_OVERRIDE(void, SphCoefs,	add, coef);
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

    void add(CoefStrPtr coef) override {
      PYBIND11_OVERRIDE(void, CylCoefs,	add, coef);
    }

  };

  py::class_<Coefs::Coefs, std::shared_ptr<Coefs::Coefs>, PyCoefs>(m, "Coefs")
    .def(py::init<std::string, bool>())
    .def("operator()", &Coefs::Coefs::operator())
    .def("getCoefStruct", &Coefs::Coefs::getCoefStruct)
    .def("Times", &Coefs::Coefs::Times)
    .def("WriteH5Coefs", &Coefs::Coefs::WriteH5Coefs)
    .def("ExtendH5Coefs", &Coefs::Coefs::ExtendH5Coefs)
    .def("Power", &Coefs::Coefs::Power)
    .def("CompareStanzas", &Coefs::Coefs::CompareStanzas);

  py::class_<Coefs::SphCoefs, std::shared_ptr<Coefs::SphCoefs>, PySphCoefs>(m, "SphCoefs")
    .def(py::init<bool>())
    .def("operator()", &Coefs::SphCoefs::operator())
    .def("getCoefStruct", &Coefs::SphCoefs::getCoefStruct)
    .def("Times", &Coefs::SphCoefs::Times)
    .def("WriteH5Coefs", &Coefs::SphCoefs::WriteH5Coefs)
    .def("ExtendH5Coefs", &Coefs::SphCoefs::ExtendH5Coefs)
    .def("Power", &Coefs::SphCoefs::Power)
    .def("CompareStanzas", &Coefs::SphCoefs::CompareStanzas);

  py::class_<Coefs::CylCoefs, std::shared_ptr<Coefs::CylCoefs>, PyCylCoefs>(m, "CylCoefs")
    .def(py::init<bool>())
    .def("operator()", &Coefs::CylCoefs::operator())
    .def("getCoefStruct", &Coefs::CylCoefs::getCoefStruct)
    .def("Times", &Coefs::CylCoefs::Times)
    .def("WriteH5Coefs", &Coefs::CylCoefs::WriteH5Coefs)
    .def("ExtendH5Coefs", &Coefs::CylCoefs::ExtendH5Coefs)
    .def("Power", &Coefs::CylCoefs::Power)
    .def("CompareStanzas", &Coefs::CylCoefs::CompareStanzas);
}

