/*
  Separate a psp structure and make a histogram

  MDWeinberg 03/15/10, 11/24/19
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <string>
#include <list>

#include <Progress.H>
#include <StringTok.H>
#include <cxxopts.H>
#include <libvars.H>
#include <header.H>
#include <PSP.H>


int
main(int argc, char **argv)
{
  char *prog = argv[0];
  double time=1e20;
  bool mweight = false;
  bool areal   = false;
  bool verbose = false;
  std::string cname("comp"), new_dir("./"), file;
  double xmin = 0.0, xmax = 1.0;
  double ymin = 0.0, ymax = 1.0;
  double zmin = -100.0, zmax = 100.0;
  int numx = 40;
  int numy = 40;
  int comp = 9;

  //--------------------
  // Parse command line
  //--------------------

  cxxopts::Options options(argv[0], "Compute histograms from PSP phase-space field quantities\n");

  options.add_options()
    ("h,help", "print this help message")
    ("t,time", "use PSP closest to given time",
     cxxopts::value<double>(time)->default_value("1.0e20"))
    ("x,xmin", "minimum x component",
     cxxopts::value<double>(xmin)->default_value("0.0"))
    ("X,xmax", "maximum x component",
     cxxopts::value<double>(xmax)->default_value("1.0"))
    ("y,ymin", "minimum y component",
     cxxopts::value<double>(ymin)->default_value("0.0"))
    ("Y,ymax", "maximum y component",
     cxxopts::value<double>(ymax)->default_value("1.0"))
    ("z,zmin", "minimum z component",
     cxxopts::value<double>(zmin)->default_value("0.0"))
    ("Z,zmax", "maximum z component",
     cxxopts::value<double>(zmax)->default_value("1.0"))
    ("1,numx", "number of bins in x direction",
     cxxopts::value<int>(numx)->default_value("40"))
    ("2,numy", "number of bins in y direction",
     cxxopts::value<int>(numy)->default_value("40"))
    ("o,comp", "component name",
     cxxopts::value<std::string>(cname)->default_value("comp"))
    ("i,index", "attribute index",
     cxxopts::value<int>(comp)->default_value("0"))
    ("d,dir", "output date location directory",
     cxxopts::value<std::string>()->default_value("./"))
    ("f,file", "input PSP file",
     cxxopts::value<std::string>(file)->default_value("out.psp"))
    ("m,mweight", "use mass-weighted values")
    ("a,areal",   "areal average")
    ("v,verbose", "verbose output")
    ;

  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    if (myid==0) std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("help")) {
    std::cout << options.help() << std::endl;
    exit(0);
  }

  if (vm.count("mweight")) mweight = true;
  if (vm.count("areal")  ) areal   = true;
  if (vm.count("verbose")) verbose = true;

  if (verbose) cerr << "Using filename: " << file << endl;


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
    
    cerr << "\nPSP file named <" << file << "> has time <" 
	 << psp->CurrentTime() << ">\n";
  }

				// Dump ascii for each component
				// -----------------------------
  
  double rtmp, mass, fac, val, dx=(xmax - xmin)/numx, dy=(ymax - ymin)/numy;
  double mtot = 0.0;
  vector<double> pos(3), vel(3);
  int itmp, ix, iy;

				// Make the array
				// --------------

  std::vector< std::vector<float> > value(numy), bmass(numy);
  for (int j=0; j<numy; j++) {
    value[j] = vector<float>(numx, 0);
    bmass[j] = vector<float>(numx, 0);
  }

  for (auto stanza=psp->GetStanza(); stanza!=0; stanza=psp->NextStanza()) {
    
    if (stanza->name != cname) continue;

    progress::progress_display progress(stanza->comp.nbod, std::cerr);

    for (auto part=psp->GetParticle(); part!=0; part=psp->NextParticle()) {

      if (part->pos(0)>=xmin && part->pos(0)<=xmax &&
	  part->pos(1)>=ymin && part->pos(1)<=ymax &&
	  part->pos(2)>=zmin && part->pos(2)<=zmax) {

	ix = static_cast<int>( floor( (part->pos(0) - xmin)/dx ) );
	iy = static_cast<int>( floor( (part->pos(1) - ymin)/dy ) );
      
	mtot += part->mass();

	if (mweight) {
	  bmass[iy][ix] += part->mass();
	  fac = mass;
	} else {
	  bmass[iy][ix] += 1.0;
	  fac = 1.0;
	}
	
	if (comp == 0)
	  value[iy][ix] += fac*part->mass();
	else if (comp <= 3)
	  value[iy][ix] += fac*part->pos(comp-1);
	else if (comp <= 6)
	  value[iy][ix] += fac*part->vel(comp-4);
	else if (comp <= 7 + part->niatr())
	  value[iy][ix] += fac*part->iatr(comp-7);
	else
	  value[iy][ix] += fac*part->datr(comp-7-part->niatr());
      }

      ++progress;
    }
  }
  
  //
  // Output
  //
  cout << std::setw(16) << numx << std::setw(16) << numy << std::endl;

  float f, m=0.0, rhomin=1e20, rhomax=0.0;

  for (int j=0; j<numy; j++) {
    for (int i=0; i<numx; i++) {
      m += bmass[j][i];
      if (areal)  {
	f = value[j][i]/dx*dy;
      } else {
	if (bmass[j][i] > 0.0) f = value[j][i]/bmass[j][i];
	else f = 0.0;
      }

      std::cout << setw(16) << xmin + dx*(0.5+i)
		<< setw(16) << ymin + dy*(0.5+j)
		<< setw(16) << f
		<< setw(16) << bmass[j][i]
		<< setw(16) << value[j][i]
		<< std::endl;
      
      rhomin = min<float>(f, rhomin);
      rhomax = max<float>(f, rhomax);
    }
  }
    
  cerr << "Total mass=" << mtot
       << ", minD=" << rhomin 
       << ", maxD=" << rhomax << endl;

  return 0;
}
  
