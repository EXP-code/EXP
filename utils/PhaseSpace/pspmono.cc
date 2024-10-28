/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Compute the monopole (spherical) model from the input PSP files
 *
 *
 *  Call sequence:
 *  -------------
 *
 *  Parameters:
 *  ----------
 *
 *
 *  Returns:
 *  -------
 *
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  MDW 01/28/22
 *
 ***************************************************************************/

#include <cmath>
#include <cstdlib>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>

#include <ParticleReader.H>
#include <MakeModel.H>
#include <localmpi.H>
#include <cxxopts.H>		// Option parsing
#include <libvars.H>		// EXP library globals
#include <interp.H>

int
main(int argc, char **argv)
{
  // MPI initialization
  //
  local_init_mpi(argc, argv);

  // Parameter assignment
  //
  double       RMIN;
  double       RMAX;
  int          NREPORT;
  int          RNUM;

  std::string  COMP;
  std::string  CURDIR;
  std::string  OUTFILE;
  std::string  fileType;

  std::vector<std::string>  INFILE, ORIENTFILE;

  std::vector<double> p0;

  const char* desc = 
    "=======================================================\n"		\
    "Compute spherical mass model from input files          \n"		\
    "=======================================================\n"		;

  cxxopts::Options options(argv[0], desc);

  options.add_options()
    ("h,help", "This help message")
    ("F,filetype", "input file type (one of: PSPout, PSPspl, GadgetNative, GadgetHDF5)",
     cxxopts::value<std::string>(fileType)->default_value("PSPout"))
    ("RMIN", "Minimum model radius",
     cxxopts::value<double>(RMIN)->default_value("0.0"))
    ("RMAX", "Maximum model radius",
     cxxopts::value<double>(RMAX)->default_value("1.0"))
    ("NREPORT", "Interval for reporting processing progress",
     cxxopts::value<int>(NREPORT)->default_value("0"))
    ("RNUM", "Size of radial grid for model",
     cxxopts::value<int>(RNUM)->default_value("1000"))
    ("INFILE", "Phase-space file(s)",
     cxxopts::value<std::vector<std::string>>(INFILE))
    ("CURDIR", "Alternative directory",
     cxxopts::value<string>(CURDIR)->default_value(""))
    ("COMP", "Compute monopole model for this component name",
     cxxopts::value<std::string>(COMP)->default_value("stars"))
    ("ORIENTFILE", "EXP generated orient file for center selection (ignored if null)",
     cxxopts::value<std::vector<std::string>>(ORIENTFILE))
    ("OUTFILE", "Output model file",
     cxxopts::value<std::string>(OUTFILE)->default_value("model.file"))
    ("CENTER", "Phase-space center",
     cxxopts::value<std::vector<double>>(p0))
    ;

  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("help")) {
    if (myid==0) std::cout << options.help() << std::endl;
    MPI_Finalize();
    return 0;
  }

  if (vm.count("CENTER")) {
    if (myid==0) {
      std::cout << "Using center:";
      for (auto v : p0) std::cout << " " << v;
      std::cout << std::endl;
    }
  } else {
    p0 = {0.0, 0.0, 0.0};
  }

  bool LOGR;
  if (RMIN>1.0e-08) {
    LOGR = true;
  } else {
    LOGR = false;
  }

  // Create the output model
  //
  auto cmodel = std::make_shared<MakeModel>(RNUM, RMIN, RMAX, LOGR);

  // Per file weight
  //
  double weight = 1.0/INFILE.size();

  // Loop through files
  // 
  for (size_t n=0; n<INFILE.size(); n++) {

    PR::PRptr psp;
    double time;

    try {
      psp = PR::ParticleReader::createReader(fileType, {INFILE[n]}, myid, true);

      psp->SelectType(COMP);

      time = psp->CurrentTime();

      if (myid==0) {
	std::cout << "File: " << INFILE[n] << std::endl;
	std::cout << "Found dump at time: " << time << std::endl;
      }
    }
    catch (const std::runtime_error& error) {
      if (myid==0)
	std::cerr << "pspmono: error opening snapshot in file <"
		  << INFILE[n] << ">" << std::endl
		  << "pspmono: " << error.what() << std::endl;
      MPI_Finalize();
      exit(-1);
    }

    
    // ===================================================================
    // Use orient file?
    // ===================================================================
  
    if (ORIENTFILE.size()) {

      std::ifstream in(ORIENTFILE[n]);
      if (!in) {
	if (myid==0)
	  std::cerr << "Couldn't open desired orient file: " << ORIENTFILE[n]
		    << std::endl;
	MPI_Finalize();
	exit(-1);
      }

      double val;
      static int buf_size = 1024;
      auto buf = std::make_unique<char[]>(buf_size);
      
      std::vector<double> or_time, or_c[3];

      while (!in.eof()) {
	in.getline(buf.get(), buf_size);
	if (in.eof()) break;
	
	istringstream ins(buf.get());
	ins >> val;
	or_time.push_back(val);
	for (int k=0; k<5; k++) ins >> val;
	for (int k=0; k<3; k++) {
	  ins >> val;
	  or_c[k].push_back(val);
	}
	for (int k=0; k<3; k++) {
	  ins >> val;
	  *(or_c[k].end()-1) += val;
	}
      }
    
      if (time < or_time.front())
	for (int k=0; k<3; k++) p0[k] = or_c[k].front();
      else if (time > or_time.back())
	for (int k=0; k<3; k++) p0[k] = or_c[k].back();
      else
	for (int k=0; k<3; k++) p0[k] = odd2(time, or_time, or_c[k]);
    }

    //============================================================
    // Build model
    //============================================================
    
    int N = 0;

    for (auto part=psp->firstParticle(); part!=0; part=psp->nextParticle()) {

      if (N++ % numprocs == myid) {
	
	double r = 0.0; 
	for (int k=0; k<3; k++)
	  r += (part->pos[k] - p0[k])*(part->pos[k] - p0[k]);
	r = sqrt(r);
	
	cmodel->AddPoint(r, part->mass*weight);
      }
      
      if (myid==0 and NREPORT) {
	if (!((N+1)%NREPORT)) std::cout << "\rProcessed: " 
					<< std::setw(10) << N+1 << std::flush;
      }
    }
  }

  auto hmodel = cmodel->Compute();

  if (myid==0) {
    cmodel->WriteModel(OUTFILE);
  }
  
  MPI_Finalize();
  
  return 0;
}
