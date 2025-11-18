#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "ParticleReader.H"

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
    "  3. PSPhdf5        The HDF5 version of the EXP phase-space format\n"
    "  4. GadgetNative   The original Gadget native format\n"
    "  5  GadgetHDF5     The newer HDF5 Gadget format\n"
    "  6. TipsyNative    The original Tipsy format\n"
    "  7. TipsyXDR       The original XDR Tipsy format\n"
    "  8. Bonsai         This is the Bonsai varient of Tipsy files\n\n"
    "We have a helper function, getReaders, to get a list to help you\n"
    "remember.  Try: pyEXP.read.ParticleReader.getReaders()\n\n"
    "Each reader can manage snapshots split into many files by parallel,\n"
    "per process writing.  The reader classes takes a list of lists on input.\n"
    "The outer list represents individual snapshots and the inner list is\n"
    "all files belonging to a single snapshot.  For unsplit snapshots, the\n"
    "inner list is a single file.  We provide two helper functions that will\n"
    "make the list of lists from the lexical sort and pattern match on a\n"
    "single list of filenames of the form 'prefix_xxxxx-yyyyy, where xxxxx\n"
    "is the integer index of the snapshot and yyyyy is the process index\n"
    "for a single snapshot.  The value for 'prefix' is arbitrary.  The\n"
    "delimiter between xxxxx and yyyyy may be user specified.  See routines\n"
    "'parseFileList()' and 'parseStringList()'.\n\n"
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

  class PyPSPhdf5 : public PSPhdf5
  {
  public:

    // Inherit the constructors
    using PSPhdf5::PSPhdf5;

    void SelectType(const std::string& type) override {
      PYBIND11_OVERRIDE(void, PSPhdf5, SelectType, type);
    }
    
    std::vector<std::string> GetTypes() override {
      PYBIND11_OVERRIDE(std::vector<std::string>, PSPhdf5, GetTypes,);
    }
    
    unsigned long CurrentNumber() override {
      PYBIND11_OVERRIDE(unsigned long, PSPhdf5, CurrentNumber,);
    }
    
    double CurrentTime() override {
      PYBIND11_OVERRIDE(double, PSPhdf5, CurrentTime,);
    }
    
    const Particle* firstParticle() override {
      PYBIND11_OVERRIDE(const Particle*, PSPhdf5, firstParticle,);
    }

    const Particle* nextParticle() override {
      PYBIND11_OVERRIDE(const Particle*, PSPhdf5, nextParticle,);
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

  };


  py::class_<Particle> P(m, "Particle");

  P.def(py::init<>(), "The internal particle type");

  py::class_<ParticleReader, std::shared_ptr<ParticleReader>, PyParticleReader>
    pr(m, "ParticleReader");

  pr.def(py::init<>(), "The base class for particle reading");

  pr.def("SelectType",      &ParticleReader::SelectType,
	 R"(
         Select the particle type to read.  

         Parameters
         ----------
         type : str
             particle type to select

         Returns
         -------
         None

         Notes
         -----
         Use GetTypes() to see the available types
         )");

  pr.def("CurrentNumber",   &ParticleReader::CurrentNumber,
	 R"(
         Number of particles in the snapshot for the selected type

         Returns
         -------
         int
             number of particles
         )");
  
  pr.def("GetTypes",        &ParticleReader::GetTypes,
	 R"(
         View the available particle types

         Returns
         -------
         list(str)
             Available types to select

         See also
         --------
         SelectType
         )");
  
  pr.def("CurrentTime",     &ParticleReader::CurrentTime,
	 R"(
	 Return the time for the current snapshot

	 Returns
         -------
         float
             the current time
	 )");
  
  pr.def("PrintSummary",
	 [](ParticleReader& A, bool stats, bool timeonly)
	 { A.PrintSummary(std::cout, stats, timeonly); },
	 R"(
         Summarize global phase-space features

         Print a summary of list of extents, center of mass, and
	 other global quantities for this snapshopt.  This requires
	 a read pass and may be time consuming.

         Parameters
         ----------
         stats : bool, default=True
             compute ensemble properties of positions and velocities
         timeonly : bool, default=False
             report current time only
         )",
	 py::arg("stats")=true, py::arg("timeonly")=false);
  
  pr.def_static("parseFileList", &ParticleReader::parseFileList,
		py::doc(R"(
                        Group files into times and segments for reader

                        Read snapshot file names from a file and format into bunches for the 
                        reader using the provided delimiter string to separate the snapshot name 
                        from the processor index

                        Parameters
                        ----------
                        file : str
                            file containing the file list
                        delimiter : str, default=" "
                            string that delimits filename fields

                        Returns
                        -------
                        list(list(str))
                            List of bunches for each snapshot. Bunches are lists of phase-space 
                            partitions for each time.
                        )"),
		py::arg("file"), py::arg("delimiter")=" ");
  
  pr.def_static("parseStringList", &ParticleReader::parseStringList,
		py::doc(R"(
			Format a list of snapshot file names into bunches for the reader

                        As in parseFileList but for a vector of file name strings

                        Parameters
                        ----------
                        file : list(str)
                            list of file names
                        delimiter : str, default=" "
                            string that delimits filename fields

                        Returns
                        -------
                        list(list(str))
                            List of bunches for each snapshot. Bunches are lists of phase-space 
                            partitions for each time.
                        )"),		
		py::arg("filelist"), py::arg("delimiter")=" ");

  pr.def_static("createReader", &ParticleReader::createReader,
		py::doc(R"(
                        Create a particle reader from the provided type

		        Uses the string bunch list constructed by praseFileList or parseStringList

                        Parameters
                        ----------
                        type : str
                            component type
                        bunch : list(str)
                            list of file segments to process
                        myid : int, default=0
                            MPI processor id (use 0 if not using MPI)
                        verbose : bool, default=False
                            verbose, diagnostic output

                        Returns
                        -------
                        ParticleReader
                        )"),
		py::arg("type"), py::arg("bunch"),
		py::arg("myid")=0, py::arg("verbose")=false);

  pr.def_static("getReaders", []()
  {
    const std::vector<std::string> formats = {
      "PSPout", "PSPspl", "PSPhdf5", "GadgetNative", "GadgetHDF5",
      "TipsyNative", "TipsyXDR", "Bonsai"};

    return formats;
    },
    py::doc(R"(
            Returns the list of phase-space snapshot format types as a Python list

            Returns
            -------
            list(str)
                List of valid format strings
            )")
    );

  py::class_<GadgetHDF5, std::shared_ptr<GadgetHDF5>, PyGadgetHDF5, ParticleReader>(m, "GadgetHDF5")
    .def(py::init<const std::vector<std::string>&, bool>(),
	 R"(
         Read Gadget HDF5 format snapshots

         Parameters
         ----------
         files : list(str)
             List of files with phase-space segments comprising a single snapshot
         verbose : bool, default=False
             Verbose, diagnostic output

         Returns
         -------
         ParticleReader
         )", py::arg("files"), py::arg("verbose")=false);

  py::class_<GadgetNative, std::shared_ptr<GadgetNative>, PyGadgetNative, ParticleReader>(m, "GadgetNative")
    .def(py::init<const std::vector<std::string>&, bool>(),
	 R"(
         Read Gadget native format snapshots

         Parameters
         ----------
         files : list(str)
             List of files with phase-space segments comprising a single snapshot
         verbose : bool, default=False
             Verbose, diagnostic output

         Returns
         -------
         ParticleReader
         )", py::arg("files"), py::arg("verbose")=false);

  py::class_<PSP, std::shared_ptr<PSP>, PyPSP, ParticleReader>(m, "PSP")
    .def(py::init<bool>(), "Base class for PSP reader")
    .def("SelectType",      &PSP::SelectType)
    .def("PrintSummary",
	 [](PSP& A, bool stats, bool timeonly)
	 { A.PrintSummary(std::cout, stats, timeonly); },
	 "Print a summary of list of extents, center of mass, and "
	 "other global quantities for this snapshopt.  This requires "
	 "a read pass and may be time consuming",
	 py::arg("stats")=true, py::arg("timeonly")=false);


  py::class_<PSPout, std::shared_ptr<PSPout>, PyPSPout, PSP>(m, "PSPout")
    .def(py::init<const std::vector<std::string>&, bool>(),
	 R"(
         Read PSP ascii monolithic format snapshot files

         Parameters
         ----------
         files : list(str)
             List of files with phase-space segments comprising a single snapshot
         verbose : bool, default=False
             Verbose, diagnostic output

         Returns
         -------
         ParticleReader
         )", py::arg("files"), py::arg("verbose")=false)
    .def("CurrentNumber",   &PSPout::CurrentNumber)
    .def("GetTypes",        &PSPout::GetTypes)
    .def("CurrentTime",     &PSPout::CurrentTime);
	 

  py::class_<PSPspl, std::shared_ptr<PSPspl>, PyPSPspl, PSP>(m, "PSPspl")
    .def(py::init<const std::vector<std::string>&, bool>(),
	 R"(
	 Reader for split PSP files (multiple files written by each process)

         Parameters
         ----------
         files : list(str)
             List of files with phase-space segments comprising a single snapshot
         verbose : bool, default=False
             Verbose, diagnostic output

         Returns
         -------
         ParticleReader
         )", py::arg("files"), py::arg("verbose")=false)
    .def("CurrentNumber",   &PSPspl::CurrentNumber)
    .def("GetTypes",        &PSPspl::GetTypes)
    .def("CurrentTime",     &PSPspl::CurrentTime);
	 
  py::class_<PSPhdf5, std::shared_ptr<PSPhdf5>, PyPSPhdf5, PSP>(m, "PSPhdf5")
    .def(py::init<const std::vector<std::string>&, bool>(),
	 R"(
	 Reader for HDF5 PSP files (both single and multiple files)

         Parameters
         ----------
         files : list(str)
             List of files with phase-space segments comprising a single snapshot
         verbose : bool, default=False
             Verbose, diagnostic output

         Returns
         -------
         ParticleReader
         )", py::arg("files"), py::arg("verbose")=false)
    .def("SelectType",      &PSPhdf5::SelectType)
    .def("NumFiles",        &PSPhdf5::NumFiles)
    .def("CurrentNumber",   &PSPhdf5::CurrentNumber)
    .def("GetTypes",        &PSPhdf5::GetTypes)
    .def("CurrentTime",     &PSPhdf5::CurrentTime);
	 
  py::class_<Tipsy, std::shared_ptr<Tipsy>, PyTipsy, ParticleReader> tipsy(m, "Tipsy");
  
  py::enum_<Tipsy::TipsyType>(tipsy, "TipsyType")
    .value("native", Tipsy::TipsyType::native)
    .value("xdr",    Tipsy::TipsyType::xdr   )
    .value("bonsai", Tipsy::TipsyType::bonsai)
    .export_values();

  tipsy.def(py::init<const std::string&, Tipsy::TipsyType, bool>(),
	    R"(
            Read Tipsy format snapshots

            Parameters
            ----------
            file : str
                   The Tipsy snapshot file
            type : TipsyType
                   The Tipsy file type (native, xdr, bonsai)
            verbose : bool, default=False
                   Verbose, diagnostic output

            Returns
            -------
            ParticleReader

            )", py::arg("file"), py::arg("type"), py::arg("verbose")=false)
    .def(py::init<const std::vector<std::string>&, Tipsy::TipsyType, bool>())
    .def("SelectType",      &Tipsy::SelectType)
    .def("CurrentNumber",   &Tipsy::CurrentNumber)
    .def("GetTypes",        &Tipsy::GetTypes)
    .def("CurrentTime",     &Tipsy::CurrentTime);

}
