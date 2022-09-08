#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <ParticleReader.H>

namespace py = pybind11;

void ParticleReaderClasses(py::module &m) {

  m.doc() = "ParticleReader class bindings\n\n"
    "This collection of classes reads and converts your phase-space\n"
    "snapshots to iterable objects for generating basis coefficients.\n"
    "The particle fields themselves are not available in Python for\n"
    "inspection currently.\n\n"
    "The available particle readers are:\n"
    "  1. PSPout         The monolithic EXP phase-space snapshot format\n"
    "  2. PSPspl         Like PSPout, but split into multiple file chunks\n"
    "  3. GadgetNative   The original Gadget native format\n"
    "  4  GadgetHDF5     The newer HDF5 Gadget format\n"
    "  5. TipsyNative    The original Tipsy format\n"
    "  6. TipsyXDR       The original XDR Tipsy format\n"
    "  7. Bonsai         This is the Bonsai varient of Tipsy files\n\n"
    "We have a helper function, getReaders, to get a list to help you\n"
    "remember.  Try: pyEXP.read.ParticleReader.getReaders()\n\n"
    "Once the ParticleReader instance is created, you can select\n"
    "the type of particle you want to read; each reader has different\n"
    "mnemonics for these, as you know.  You can get the list of types\n"
    "using the GetTypes() member and select the desired type with the\n"
    "SelectType() member.  You can also get the current time and number\n"
    "of particles for the selected type using the CurrentTime() and\n"
    "CurrentNumber() members, respectively.  The main function of of\n"
    "this ParticleReader object in pyEXP is creating the coefficient\n"
    "expansion using Basis.createCoefficients.  See help(pyEXP.basis)\n\n";
			    

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


  py::class_<Particle> P(m, "Particle");

  P.def(py::init<>(), "The internal particle type");

  py::class_<ParticleReader, std::shared_ptr<ParticleReader>, PyParticleReader>
    pr(m, "ParticleReader");

  pr.def(py::init<>(), "The base class for particle reading");

  pr.def("SelectType",      &ParticleReader::SelectType,
	 "Select the particle type to read.  "
	 "Use GetTypes() to see the available types.");

  pr.def("CurrentNumber",   &ParticleReader::CurrentNumber,
	 "Return the current number of particles in the snapshot "
	 "for the selected type");
  
  pr.def("GetTypes",        &ParticleReader::GetTypes,
	 "View the available particle types");
  
  pr.def("CurrentTime",     &ParticleReader::CurrentTime,
	 "Return the time for the current snapshot");
  
  pr.def("PrintSummary",    &ParticleReader::PrintSummary,
	 "Print a summary of list of extents, center of mass, and "
	 "other global quantities for this snapshopt.  This requires "
	 "a read pass and may be time consuming");
  
  pr.def_static("parseFileList", &ParticleReader::parseFileList,
		py::doc("Read snapshot file names from a file and format into "
			"bunches for the reader using the provided delimiter "
			"string to separate the snapshot name from the "
			"processor index"),
		py::arg("file"), py::arg("delimiter")="");
  
  pr.def_static("parseStringList", &ParticleReader::parseStringList,
		py::doc("Format a list of snapshot file names into bunches "
			"for the reader"),
		py::arg("filelist"), py::arg("delimiter")="");

  pr.def_static("createReader", &ParticleReader::createReader,
	 py::doc("Create a particle reader from the provided type "
		 "string bunch list constructed by praseFileList or "
		 "parseStringList"),
	 py::arg("type"), py::arg("bunch"),
	 py::arg("myid")=0, py::arg("verbose")=false);

  pr.def_static("getReaders", []()
  {
    const std::vector<std::string> formats = {
      "PSPout", "PSPspl", "GadgetNative", "GadgetHDF5", "TipsyNative",
      "TipsyXDR", "Bonsai"};

    return formats;
    },
    py::doc("Returns the list of phase-space snapshot format types as a Python list"));

  py::class_<GadgetHDF5, std::shared_ptr<GadgetHDF5>, PyGadgetHDF5, ParticleReader>(m, "GadgetHDF5")
    .def(py::init<const std::vector<std::string>&, bool>(),
	 "Read Gadget HDF5 format snapshots");

  py::class_<GadgetNative, std::shared_ptr<GadgetNative>, PyGadgetNative, ParticleReader>(m, "GadgetNative")
    .def(py::init<const std::vector<std::string>&, bool>(), "Read Gadget native format snapshots");

  py::class_<PSP, std::shared_ptr<PSP>, PyPSP, ParticleReader>(m, "PSP")
    .def(py::init<bool>(), "Base class for PSP reader")
    .def("SelectType",      &PSP::SelectType)
    .def("PrintSummary",    &PSP::PrintSummary);

  py::class_<PSPout, std::shared_ptr<PSPout>, PyPSPout, PSP>(m, "PSPout")
    .def(py::init<const std::vector<std::string>&, bool>(),
	 "Reader for monolitic PSP format (single file written by root process)")
    .def("CurrentNumber",   &PSPout::CurrentNumber)
    .def("GetTypes",        &PSPout::GetTypes)
    .def("CurrentTime",     &PSPout::CurrentTime);
	 

  py::class_<PSPspl, std::shared_ptr<PSPspl>, PyPSPspl, PSP>(m, "PSPspl")
    .def(py::init<const std::vector<std::string>&, bool>(),
	 "Reader for split PSP files (multiple files written by each process)")
    .def("CurrentNumber",   &PSPspl::CurrentNumber)
    .def("GetTypes",        &PSPspl::GetTypes)
    .def("CurrentTime",     &PSPspl::CurrentTime);
	 
py::class_<Tipsy, std::shared_ptr<Tipsy>, PyTipsy, ParticleReader> tipsy(m, "Tipsy");
  
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

