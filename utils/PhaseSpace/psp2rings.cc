/*
  Separate a psp structure and make a kinematic Fourier coefficients
  series in rings

  MDWeinberg 10/18/19
*/

#include <functional>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <random>
#include <string>
#include <list>
#include <map>

#include <header.H>
#include <PSP.H>
#include <FileUtils.H>

#include <Progress.H>
#include <cxxopts.H>
#include <libvars.H>		// EXP library globals


//! Coefficient file header
struct RingCoefHeader
{
  double time;
  int nrings;
  int mnum;
};

class RingCoefs
{
 private:

  double time;
  int nrad, mnum;

public:

  //! Coefficient data
  std::map<int, std::array<std::vector<double>, 3>> cos_c, sin_c;

  //! Constructor
  RingCoefs(double time, int mmin, int mmax, int nrad) :
    time(time), nrad(nrad)
  {
    for (int m=mmin; m<=mmax; m++) {
      for (size_t k=0; k<3; k++) {
	cos_c[m][k].resize(nrad, 0);
	if (m) sin_c[m][k].resize(nrad, 0);
      }
    }
  }

  void write(std::ostream& out);
};

typedef std::shared_ptr<RingCoefs> RingCoefPtr;

void RingCoefs::write(std::ostream& out)
{
  mnum = cos_c.size();

  RingCoefHeader header;
  header.time   = time;
  header.nrings = nrad;
  header.mnum   = mnum;
  
  out.write((const char *)&header, sizeof(RingCoefHeader));

  for (auto d : cos_c) {
    out.write((const char *)&d.first, sizeof(int));
    out.write((const char *)d.second[0].data(), sizeof(double)*nrad);
    out.write((const char *)d.second[1].data(), sizeof(double)*nrad);
    out.write((const char *)d.second[2].data(), sizeof(double)*nrad);
    if (d.first) {
      out.write((const char *)sin_c[d.first][0].data(), sizeof(double)*nrad);
      out.write((const char *)sin_c[d.first][1].data(), sizeof(double)*nrad);
      out.write((const char *)sin_c[d.first][2].data(), sizeof(double)*nrad);
    }
  }
}

int
main(int ac, char **av)
{
  char *prog = av[0];
  bool verbose = false;
  std::string cname, tname, new_dir, suffix, work_dir, runtag;
  int axis, numb, comp, mmin, mmax, ibeg, iend;
  double pmin, pmax;

  // Parse command line
  //
  cxxopts::Options options(prog, "Separate a psp structure and make a kinematic Fourier coefficients series in rings.\n");

  options.add_options()
   ("h,help", "produce help message")
   ("v,verbose", "verbose output")
   ("i,beg", "initial snapshot index",
     cxxopts::value<int>(ibeg)->default_value("0"))
   ("e,end", "final snapshot index",
     cxxopts::value<int>(iend)->default_value(std::to_string(std::numeric_limits<int>::max())))
   ("m,mmin", "minimum Fourier component in bin",
     cxxopts::value<int>(mmin)->default_value("1"))
   ("M,mmax", "maximum Fourier component in bin",
     cxxopts::value<int>(mmax)->default_value("4"))
   ("r,rmin", "minimum bin radius",
     cxxopts::value<double>(pmin)->default_value("0.0"))
   ("R,rmax", "maximum bin radius",
     cxxopts::value<double>(pmax)->default_value("0.04"))
   ("b,bins", "number of bins",
     cxxopts::value<int>(numb)->default_value("40"))
   ("c,name", "component name",
     cxxopts::value<std::string>(cname)->default_value("comp"))
   ("d,dir", "rewrite directory location for SPL files",
     cxxopts::value<std::string>(new_dir)->default_value("./"))
   ("w,work", "working directory for output file",
     cxxopts::value<std::string>(work_dir)->default_value("."))
   ("t,type", "PSP output type (OUT or SPL)",
     cxxopts::value<std::string>(tname)->default_value("OUT"))
   ("T,runtag", "Runtag id",
     cxxopts::value<std::string>(runtag)->default_value("run0"))
   ("s,suffix", "Output file suffix",
     cxxopts::value<std::string>(suffix)->default_value("ring_coefs"))
    ;

  cxxopts::ParseResult vm;

  try {
    vm = options.parse(ac, av);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("help")) {
    std::cout << options.help() << std::endl;
    std::cout << "Example: " << std::endl;
    std::cout << "\t" << av[0]
	      << " --output=out.bod" << std::endl;
    return 1;
  }

  if (vm.count("verbose")) {
    verbose = true;
  }

  int n;
  for (n=ibeg; n<=iend; n++) {

    std::ostringstream fname;
    fname << tname << "." << runtag << "."
	  << std::setw(5) << std::setfill('0') << n;

    std::string file = fname.str();

    if (!FileExists(file)) {
      std::cerr << "Error opening file <" << file << "> for input"
		<< std::endl;
      break;
    }
  }
  iend = n-1;
  std::cerr << "Assuming last file has index <" << iend << ">"
	    << std::endl;

  std::string outcoefs = work_dir + "/" + runtag + "." + suffix;
  std::ofstream out(outcoefs);
  if (!out) {
    std::cerr << "Error opening file <" << outcoefs << "> for output"
	      << std::endl;
    exit(-1);
  }

				// Write ring bin data
  double dp = (pmax - pmin)/numb;
  out.write((const char *)&numb, sizeof(int));
  for (int p=0; p<numb; p++) {
    double b = pmin + dp*p;
    out.write((const char *)&b, sizeof(double));
  }
  for (int p=0; p<numb; p++) {
    double b = pmin + dp*(p+1);
    out.write((const char *)&b, sizeof(double));
  }

  progress::progress_display progress(iend - ibeg + 1);

  for (int n=ibeg; n<=iend; n++) {

    std::ostringstream fname;
    fname << tname << "." << runtag << "."
	  << std::setw(5) << std::setfill('0') << n;

    std::string file = fname.str();

    if (!FileExists(file)) {
      cerr << "Error opening file <" << file << "> for input\n";
      break;
    }

    if (verbose) cerr << "Using filename: " << file << endl;
    else ++progress;

				// Parse the PSP file
				// ------------------
    PSPptr psp;
    if (file.find("SPL") != std::string::npos)
      psp = std::make_shared<PSPspl>(file, new_dir);
    else
      psp = std::make_shared<PSPout>(file);


				// Now write a summary
				// -------------------
    if (verbose) {
      
      psp->PrintSummary(cerr);
    
      std::cerr << "\nPSP file <" << file << "> has time <" 
	   << psp->CurrentTime() << ">\n";
    }

				// Dump ascii for each component
				// -----------------------------
    std::vector<double> pos(3), vel(3);

				// Make the arrays
				// ---------------

    RingCoefs rings(psp->CurrentTime(), mmin, mmax, numb);

    std::vector<float> value(numb, 0), bmass(numb, 0);
    std::array<std::map<int, std::vector<float>>, 3> vel_c, vel_s;

    PSPstanza *stanza;
    SParticle* part;

    for (stanza=psp->GetStanza(); stanza!=0; stanza=psp->NextStanza()) {
    
      if (stanza->name != cname) continue;

      for (part=psp->GetParticle(); part!=0; part=psp->NextParticle()) {
	
	// Cylindrical radius
	//
	double val = 0.0;
	for (int k=0; k<2; k++) val += part->pos(k) * part->pos(k);
	val = sqrt(val);

	int iv = static_cast<int>( floor( (val - pmin)/dp ) );

	if (iv < 0 || iv >= numb) continue;

	double mass = part->mass();
	
	bmass[iv] += mass;

	// Make cylindrical velocity bins
	//
	double phi  = std::atan2(part->pos(1), part->pos(0));
	double cosp = std::cos(phi);
	double sinp = std::sin(phi);

	// uvec vr:  cos(phi), sin(phi)
	double vr = cosp*part->vel(0) + sinp*part->vel(1);

	// uvec vt: -sin(phi), cos(phi)
	double vt = -sinp*part->vel(0) + cosp*part->vel(1);
	
	// vertical
	double vz = part->vel(2);

	// Add to grid
	for (int m=mmin; m<=mmax; m++) {
	  double cosm = std::cos(phi*m);
	  rings.cos_c[m][0][iv] += mass*vr*cosm;
	  rings.cos_c[m][1][iv] += mass*vt*cosm;
	  rings.cos_c[m][2][iv] += mass*vz*cosm;
	  if (m) {
	    double sinm = std::sin(phi*m);
	    rings.sin_c[m][0][iv] += mass*vr*sinm;
	    rings.sin_c[m][1][iv] += mass*vt*sinm;
	    rings.sin_c[m][2][iv] += mass*vz*sinm;
	  }
	}
      }

    }

    // Normalize
    //
    for (int m=mmin; m<=mmax; m++) {
      for (int b=0; b<numb; b++) {
	if (bmass[b]>0.0) {
	  for (int k=0; k<3; k++) {
	    rings.cos_c[m][k][b] /= bmass[b];
	    if (m) rings.sin_c[m][k][b] /= bmass[b];
	  }
	}
      }
    }

    rings.write(out);
  }
  std::cout << std::endl;

  return 0;
}
