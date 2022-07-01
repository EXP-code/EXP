/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Use a particle reader to read a body file and print ascii
 *  components for EXP input
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
 *  MDW 07/01/22
 *
 ***************************************************************************/

#include <cmath>
#include <cstdlib>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>

#include <ParticleReader.H>	// Read n-body snapshots
#include <cxxopts.H>		// Option parsing
#include <libvars.H>		// EXP library globals

int
main(int argc, char **argv)
{
  // Parameter assignment
  //
  int          NREPORT;

  std::string  COMP;
  std::string  CURDIR;
  std::string  OUTFILE;
  std::string  fileType;
  double       mscale, lscale, vscale;

  std::vector<std::string>  INFILE;
  
  const char* desc = 
    "=======================================================\n"		\
    "Compute EXP ascii files from n-body snapshots          \n"		\
    "=======================================================\n"		;

  cxxopts::Options options(argv[0], desc);

  options.add_options()
    ("h,help", "This help message")
    ("F,filetype", "input file type (one of: PSPout, PSPspl, GadgetNative, GadgetHDF5)",
     cxxopts::value<std::string>(fileType)->default_value("PSPout"))
    ("n,NREPORT", "Interval for reporting processing progress",
     cxxopts::value<int>(NREPORT)->default_value("0"))
    ("i,INFILE", "Phase-space file(s)",
     cxxopts::value<std::vector<std::string>>(INFILE))
    ("d,CURDIR", "Alternative directory",
     cxxopts::value<string>(CURDIR)->default_value(""))
    ("c,COMP", "Compute wake for this component name",
     cxxopts::value<std::string>(COMP)->default_value("stars"))
    ("o,OUTFILE", "Output model file",
     cxxopts::value<std::string>(OUTFILE)->default_value("model.file"))
    ("m,mscale", "Mass scale factor",
     cxxopts::value<double>(mscale)->default_value("1.0"))
    ("l,lscale", "Length scale factor",
     cxxopts::value<double>(lscale)->default_value("1.0"))
    ("v,vscale", "Velocity scale factor",
     cxxopts::value<double>(lscale)->default_value("1.0"))
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
    return 0;
  }

  PR::PRptr snap;
  double time;

  std::ofstream out(OUTFILE);
  if (!out) {
    std::cerr << "part2ascii: error opening output file <"
	      << OUTFILE << ">" << std::endl;
    exit(-1);
  }

  try {
    snap = PR::ParticleReader::createReader(fileType, INFILE, myid, true);
    
    snap->SelectType(COMP);

    time = snap->CurrentTime();

    std::cout << "File(s) ";
    for (auto s : INFILE) std::cout << s << " ";
    std::cout << std::endl;
    std::cout << "Found dump at time: " << time << std::endl;
  }
  catch (const std::runtime_error& error) {
    std::cerr << "part2ascii: error opening snapshot in files ";
    for (auto s : INFILE) std::cerr << s << " ";
    std::cerr << std::endl
	      << "part2ascii: " << error.what() << std::endl;
    exit(-1);
  }
  
  out << std::setw(15) << snap->CurrentNumber()
      << std::setw(10) << 0
      << std::setw(10) << 0
      << std::endl;

  int N=0;

  for (auto part=snap->firstParticle(); part!=0; part=snap->nextParticle(), N++) {
    out << std::setw(18) << part->indx
	<< std::setw(18) << part->mass*mscale;
    for (int k=0; k<3; k++) out << std::setw(18) << part->pos[k]*lscale;
    for (int k=0; k<3; k++) out << std::setw(18) << part->vel[k]*vscale;
    out << std::endl;
    
    if (NREPORT) {
      if (!((N+1)%NREPORT)) std::cout << "\rProcessed: " 
				      << std::setw(10) << N+1 << std::flush;
    }
  }

  return 0;
}
