#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <ParticleReader.H>

namespace py = pybind11;

void ParticleReaderClasses(py::module &m) {

  m.doc() = "ParticleReader class bindings";

  using namespace PR;

  class PyParticleReader : public ParticleReader
  {
  public:

    // Inherit the constructors
    using ParticleReader::ParticleReader;

    void SelectType(const std::string& type) override {
      PYBIND11_OVERRIDE_PURE(void, ParticleReader, SelectType, type);
    }
    
    std::vector<std::string> GetTypes() override {
      PYBIND11_OVERRIDE_PURE(std::vector<std::string>, ParticleReader, GetTypes,);
    }
    
    unsigned long CurrentNumber() override {
      PYBIND11_OVERRIDE_PURE(unsigned long, ParticleReader, CurrentNumber,);
    }
    
    double CurrentTime() override {
      PYBIND11_OVERRIDE_PURE(double, ParticleReader, CurrentTime,);
    }
    
    const Particle* firstParticle() override {
      PYBIND11_OVERRIDE_PURE(const Particle*, ParticleReader, firstParticle,);
    }

    const Particle* nextParticle() override {
      PYBIND11_OVERRIDE_PURE(const Particle*, ParticleReader, nextParticle,);
    }

    void PrintSummary(std::ostream& out, bool stats, bool timeonly) override {
      PYBIND11_OVERRIDE_PURE(void, ParticleReader, PrintSummary,
			     out, stats, timeonly);
    }

  };

  class PyGadgetHDF5 : public GadgetHDF5
  {
  public:

    // Inherit the constructors
    using GadgetHDF5::GadgetHDF5;

    void SelectType(const std::string& type) override {
      PYBIND11_OVERRIDE(void, GadgetHDF5, SelectType, type);
    }
    
    std::vector<std::string> GetTypes() override {
      PYBIND11_OVERRIDE(std::vector<std::string>, GadgetHDF5, GetTypes,);
    }
    
    unsigned long CurrentNumber() override {
      PYBIND11_OVERRIDE(unsigned long, GadgetHDF5, CurrentNumber,);
    }
    
    double CurrentTime() override {
      PYBIND11_OVERRIDE(double, GadgetHDF5, CurrentTime,);
    }
    
    const Particle* firstParticle() override {
      PYBIND11_OVERRIDE(const Particle*, GadgetHDF5, firstParticle,);
    }

    const Particle* nextParticle() override {
      PYBIND11_OVERRIDE(const Particle*, GadgetHDF5, nextParticle,);
    }

    void PrintSummary(std::ostream& out, bool stats, bool timeonly) override {
      PYBIND11_OVERRIDE(void, GadgetHDF5, PrintSummary,
			out, stats, timeonly);
    }
    
  };

  class PyGadgetNative : public GadgetNative
  {
  public:

    // Inherit the constructors
    using GadgetNative::GadgetNative;

    void SelectType(const std::string& type) override {
      PYBIND11_OVERRIDE(void, GadgetNative, SelectType,	type);
    }
    
    std::vector<std::string> GetTypes() override {
      PYBIND11_OVERRIDE(std::vector<std::string>, GadgetNative, GetTypes,);
    }
    
    unsigned long CurrentNumber() override {
      PYBIND11_OVERRIDE(unsigned long, GadgetNative, CurrentNumber,);
    }
    
    double CurrentTime() override {
      PYBIND11_OVERRIDE(double, GadgetNative, CurrentTime,);
    }
    
    const Particle* firstParticle() override {
      PYBIND11_OVERRIDE(const Particle*, GadgetNative, firstParticle,);
    }

    const Particle* nextParticle() override {
      PYBIND11_OVERRIDE(const Particle*, GadgetNative, nextParticle,);
    }

    void PrintSummary(std::ostream& out, bool stats, bool timeonly) override {
      PYBIND11_OVERRIDE(void, GadgetNative, PrintSummary,
			out, stats, timeonly);
    }
    
  };

  class PyPSP : public PSP
  {
  public:

    // Inherit the constructors
    using PSP::PSP;

    void SelectType(const std::string& type) override {
      PYBIND11_OVERRIDE(void, PSP, SelectType, type);
    }

    PSPstanza* GetNamed(const std::string& name) override {
      PYBIND11_OVERRIDE(PSPstanza*, PSP, GetNamed, name);
    }
    
    PSPstanza* GetStanza() override {
      PYBIND11_OVERRIDE(PSPstanza*, PSP, GetStanza,);
    }
    
    PSPstanza* NextStanza() override {
      PYBIND11_OVERRIDE(PSPstanza*, PSP, NextStanza,);
    }
    
    const Particle* firstParticle() override {
      PYBIND11_OVERRIDE_PURE(const Particle*, ParticleReader, firstParticle,);
    }

    const Particle* nextParticle() override {
      PYBIND11_OVERRIDE_PURE(const Particle*, ParticleReader, nextParticle,);
    }

  };

  class PyPSPout : public PSPout
  {
  public:

    // Inherit the constructors
    using PSPout::PSPout;

    const Particle* firstParticle() override {
      PYBIND11_OVERRIDE(const Particle*, PSPout, firstParticle,);
    }

    const Particle* nextParticle() override {
      PYBIND11_OVERRIDE(const Particle*, PSPout, nextParticle,);
    }

  };

  class PyPSPspl : public PSPspl
  {
  public:

    // Inherit the constructors
    using PSPspl::PSPspl;

    const Particle* firstParticle() override {
      PYBIND11_OVERRIDE(const Particle*, PSPspl, firstParticle,);
    }

    const Particle* nextParticle() override {
      PYBIND11_OVERRIDE(const Particle*, PSPspl, nextParticle,);
    }
  };

  class PyTipsy : public Tipsy
  {
  public:

    // Inherit the constructors
    using Tipsy::Tipsy;

    void SelectType(const std::string& type) override {
      PYBIND11_OVERRIDE(void, Tipsy, SelectType, type);
    }
    
    std::vector<std::string> GetTypes() override {
      PYBIND11_OVERRIDE(std::vector<std::string>, Tipsy, GetTypes,);
    }
    
    unsigned long CurrentNumber() override {
      PYBIND11_OVERRIDE(unsigned long, Tipsy, CurrentNumber,);
    }
    
    double CurrentTime() override {
      PYBIND11_OVERRIDE(double, Tipsy, CurrentTime,);
    }
    
    const Particle* firstParticle() override {
      PYBIND11_OVERRIDE(const Particle*, Tipsy,	firstParticle,);
    }

    const Particle* nextParticle() override {
      PYBIND11_OVERRIDE(const Particle*, Tipsy,	nextParticle,);
    }

    void PrintSummary(std::ostream& out, bool stats, bool timeonly) override {
      PYBIND11_OVERRIDE(void, Tipsy, PrintSummary,
			out, stats, timeonly);
    }
    
  };


  py::class_<Particle>(m, "Particle")
    .def(py::init<>());

  py::class_<ParticleReader, PyParticleReader>(m, "ParticleReader")
    .def(py::init<>())
    .def("SelectType",      &ParticleReader::SelectType)
    .def("CurrentNumber",   &ParticleReader::CurrentNumber)
    .def("GetTypes",        &ParticleReader::GetTypes)
    .def("CurrentTime",     &ParticleReader::CurrentTime)
    .def("PrintSummary",    &ParticleReader::PrintSummary)
    .def("parseFileList",   &ParticleReader::parseFileList)
    .def("parseStringList", &ParticleReader::parseStringList)
    .def("createReader",    &ParticleReader::createReader);

  py::class_<GadgetHDF5, PyGadgetHDF5, ParticleReader>(m, "GadgetHDF5")
    .def(py::init<const std::vector<std::string>&, bool>());

  py::class_<GadgetNative, PyGadgetNative, ParticleReader>(m, "GadgetNative")
    .def(py::init<const std::vector<std::string>&, bool>());

  py::class_<PSP, PyPSP, ParticleReader>(m, "PSP")
    .def(py::init<bool>())
    .def("SelectType",      &PSP::SelectType)
    .def("PrintSummary",    &PSP::PrintSummary);

  py::class_<PSPout, PyPSPout, PSP>(m, "PSPout")
    .def(py::init<const std::vector<std::string>&, bool>())
    .def("CurrentNumber",   &PSPout::CurrentNumber)
    .def("GetTypes",        &PSPout::GetTypes)
    .def("CurrentTime",     &PSPout::CurrentTime);


  py::class_<PSPspl, PyPSPspl, PSP>(m, "PSPspl")
    .def(py::init<const std::vector<std::string>&, bool>())
    .def("CurrentNumber",   &PSPspl::CurrentNumber)
    .def("GetTypes",        &PSPspl::GetTypes)
    .def("CurrentTime",     &PSPspl::CurrentTime);

  py::class_<Tipsy, PyTipsy, ParticleReader> tipsy(m, "Tipsy");
  
  py::enum_<Tipsy::TipsyType>(tipsy, "TipsyType")
    .value("native", Tipsy::TipsyType::native)
    .value("xdr",    Tipsy::TipsyType::xdr   )
    .value("bonsai", Tipsy::TipsyType::bonsai)
    .export_values();

  tipsy.def(py::init<const std::string&, Tipsy::TipsyType, bool>())
    .def(py::init<const std::vector<std::string>&, Tipsy::TipsyType, bool>())
    .def("SelectType",      &Tipsy::SelectType)
    .def("CurrentNumber",   &Tipsy::CurrentNumber)
    .def("GetTypes",        &Tipsy::GetTypes)
    .def("CurrentTime",     &Tipsy::CurrentTime)
    .def("PrintSummary",    &Tipsy::PrintSummary);

}

