#include <filesystem>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>
#include <array>

#include <cxxopts.H>
#include <libvars.H>

#include <ParticleReader.H>

int main(int argc, char **argv)
{
  const double pi = 3.14159265358979323846;
  std::string ascii, type, files, delim, comp;
  double rmax = 0.05;
  int nbins = 80;

  // Parse Command line
  //
  const std::string overview = "A quick test of particle reading through comparison with an ascii bods file\n";

  cxxopts::Options options(argv[0], overview);

  options.add_options()
    ("h,help", "produce this help message")
    ("v, verbose", "provide details of running pipeline")
    ("files", "snapshot file list",
     cxxopts::value<std::string>(files)->default_value("file.list"))
    ("delim", "grouping delimiter",
     cxxopts::value<std::string>(delim))
    ("bods", "ascii bods file",
     cxxopts::value<std::string>(ascii)->default_value("bods"))
    ("t,type", "snapshot type",
     cxxopts::value<std::string>(type)->default_value("PSPspl"))
    ("c,comp", "The component name",
     cxxopts::value<std::string>(comp)->default_value("dark"))
     ;
  
  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("help")) {
    std::cout << options.help() << std::endl;
    return 1;
  }

  bool verbose = vm.count("verbose");

  // First open ascii file and read positions
  //
  std::map<unsigned, std::array<double, 4>> mpos;

  try {
    // If the file exists, open it
    //
    if (std::filesystem::exists(ascii)) {
      std::ifstream bods(ascii);
      
      unsigned indx;
      std::string line;
      std::getline(bods, line);	// Header line; discard

      while (std::getline(bods, line)) {
	std::istringstream sin(line);
	sin >> indx;
	for (int k=0; k<4; k++) sin >> mpos[indx][k];
      }
    }
  }
  catch (const std::runtime_error& error) {
    std::cout << "testread: found a problem reading ascii file" << std::endl
	      << error.what() << std::endl;
    return(1);
  }

  try {

    unsigned total = 0;
    double maxdif = 0.0;
    std::vector<double> histo(nbins, 0.0);
    double delta = rmax/nbins;

    // Loop through list of snapshots
    //
    for (auto batch : PR::ParticleReader::parseFileList(files, delim)) {
      
      // Create a particle reader object
      //
      auto reader = PR::ParticleReader::createReader(type, batch);
      reader->SelectType(comp);

      for (auto p=reader->firstParticle(); p!=0; p=reader->nextParticle()) {

	// Find the entry in the ascii-read db
	auto it = mpos.find(p->indx);
	if (it == mpos.end()) {
	  std::cout << "Index=" << p->indx << " not found" << std::endl;
	} else {
	  if (fabs(p->mass - it->second[0]) > 1.0e-6) {
	    std::cout << "Index=" << p->indx << " weird mass=" << p->mass
		      << std::endl;
	  }
	  
	  double dif = 0.0;
	  for (int k=0; k<3; k++) {
	    double d = it->second[k+1] - p->pos[k];
	    dif += d*d;
	  }
	  dif = sqrt(dif);
	  
	  if (dif>maxdif) maxdif = dif;

	  if (fabs(dif) > 1.0e-6) {
	    std::cout << "Index=" << p->indx << " weird pos dif=" << dif
		      << std::endl;
	  }
	  total++;

	  double r = sqrt(p->pos[0]*p->pos[0] + p->pos[1]*p->pos[1]);
	  int indx = floor(r/delta);
	  if (indx >= 0 and indx < nbins) histo[indx] += p->mass;
	}
      }
    }

    std::cout << "Read " << total << " particles with max dif="
	      << maxdif << std::endl;

    std::ofstream out("histo.out");
    for (int n=0; n<nbins; n++) {
      out << std::setw(18) << delta*(0.5 + n)
	  << std::setw(18) << histo[n] / (pi*delta*delta*(2.0*n + 1.0))
	  << std::endl;
    }

  }
  catch (const std::runtime_error& error) {
    std::cout << "testread: found a problem reading snap files" << std::endl
	      << error.what() << std::endl;
    return(1);
  }

  // Done
  //
  return(0);
}
