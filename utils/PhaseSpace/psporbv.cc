/*
  Look for radial mode

  This is an example of using PSP classes to look for a particular
  feature based on a background model

  MDWeinberg 01/08/2022
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <memory>
#include <string>
#include <cmath>
#include <list>

#include <Progress.H>
#include <cxxopts.H>
#include <libvars.H>
#include <header.H>

#include <massmodel.H>
#include <localmpi.H>
#include <orbit.H>
#include <PSP.H>


int
main(int argc, char **argv)
{
  char *prog = argv[0];
  bool verbose = false;
  std::string cname("comp"), new_dir("./"), modfile("SLGridSph.model");
  std::string psfile, outfile;
  double rmin = 0.0, rmax = 1.0, maxkap = 0.5, KTOL, frac, vcut;

  // MPI initialization
  //
  local_init_mpi(argc, argv);

  // Parse command line
  //
  cxxopts::Options options(argv[0], "Compute orbits that satisfy particular constraints in orbital phase and (Energy, Kappa)\n");

  options.add_options()
    ("h,help", "print this help message")
    ("r,rmin", "minimum apocentric radius",
     cxxopts::value<double>(rmin)->default_value("0.0"))
    ("R,rmax", "maximum apocentric radius",
     cxxopts::value<double>(rmax)->default_value("1.0"))
    ("F,frac", "maximum fraction of apocentric radius",
     cxxopts::value<double>(frac)->default_value("0.5"))
    ("K,maxK", "maximum kappa value",
     cxxopts::value<double>(maxkap)->default_value("0.5"))
    ("V,vcut", "maximum ratio of vrad to vtan",
     cxxopts::value<double>(vcut)->default_value("-2.0"))
    ("t,tol", "minimum kappa value",
     cxxopts::value<double>(KTOL)->default_value("0.005"))
    ("c,comp", "component name",
     cxxopts::value<std::string>(cname)->default_value("comp"))
    ("d,dir", "output date location directory",
     cxxopts::value<std::string>()->default_value("./"))
    ("f,PSP", "input PSP file",
     cxxopts::value<std::string>(psfile)->default_value("out.psp"))
    ("o,outfile", "output data file",
     cxxopts::value<std::string>(outfile)->default_value("orbits.dat"))
    ("m,model", "SphericalModel file",
     cxxopts::value<std::string>(modfile)->default_value("SLGridSph.mod"))
    ("v,verbose", "verbose output")
    ;

  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    MPI_Finalize();
    exit(-1);
  }

  if (vm.count("help")) {
    if (myid==0) std::cout << options.help() << std::endl;
    MPI_Finalize();
    exit(0);
  }

  if (vm.count("verbose")) verbose = true;

  if (verbose) {
    if (myid==0) std::cerr << "Using filename: " << psfile << std::endl;
  }

				// Parse the PSP file
				// ------------------
  PSPptr psp;
  if (psfile.find("SPL") != std::string::npos)
    psp = std::make_shared<PSPspl>(psfile, new_dir);
  else
    psp = std::make_shared<PSPout>(psfile);


				// Now write a summary
				// -------------------
  if (verbose) {

    psp->PrintSummary(cerr);
    
    if (myid==0)
      std::cerr << std::endl << "PSP file named <" << psfile << "> has time <" 
		<< psp->CurrentTime() << ">" << std::endl;
  }
  
				// Create the orbit instance
				// -------------------------
  SphModTblPtr model = std::make_shared<SphericalModelTable>(modfile);
  SphericalOrbit orb(model);
  
  double Emin = model->get_pot(model->get_min_radius());
  double Emax = model->get_pot(model->get_max_radius());

				// Dump ascii for each component
				// -----------------------------
  std::ostringstream sout;
  sout << outfile << "." << myid;

  std::ofstream out(sout.str());
  const int wid = 18;

  if (out) {
    std::vector<std::string> labels =
      { "radius", "energy", "kappa", "I_r", "I_p",
	"L_x", "L_y", "L_z", "Omega_1", "Omega_2", "r/r_apo", "r/r_peri",
	"w1", "v_rad/v_tan", "phi", "theta", "v_rad", "v_tan", "index" };

    // Separating header
    //
    auto sep = [wid, labels, out](int v) {
		 for (int i=0; i<labels.size(); i++) {
		   std::ostringstream sout;
		   if (i==0) out << "#"; else out << "+";
		   if (v==0) sout << std::setw(wid-1) << std::setfill('-') << '-';
		   if (v==1) sout << std::setw(wid-2) << std::right << labels[i] << ' ';
		   if (v==2) { std::ostringstream sout2; sout2 << "[" << i+1 << "] ";
		     sout << std::setw(wid-1) << std::right << sout2.str();
		   }
		   out << sout.str();
		 }
		 out << std::endl;
	       };
    
    out << sep(0) << sep(1) << sep(2) << sep(0);
  } else {
    std::cerr << "Error opening <" << sout.str() << std::endl;
  }

  unsigned good = 0, all = 0;

  std::shared_ptr<progress::progress_display> progress;

  for (auto stanza=psp->GetStanza(); stanza!=0; stanza=psp->NextStanza()) {
    
    if (stanza->name != cname) continue;

    // Make the progress bar for root process only
    if (myid==0)
      progress = std::make_shared<progress::progress_display>
	(stanza->comp.nbod/numprocs, std::cerr);
    
    for (auto part=psp->GetParticle(); part!=0; part=psp->NextParticle()) {
      
      if (all++ % numprocs != myid) continue;

      double r2 = 0.0, v2 = 0.0;
      for (int k=0; k<3; k++) {
	r2 += part->pos(k)*part->pos(k);
	v2 += part->vel(k)*part->vel(k);
      }
      double r = sqrt(r2);

      double Jx = part->pos(1)*part->vel(2) - part->pos(2)*part->vel(1);
      double Jy = part->pos(2)*part->vel(0) - part->pos(0)*part->vel(2);
      double Jz = part->pos(0)*part->vel(1) - part->pos(1)*part->vel(0);

      double J  = sqrt(Jx*Jx + Jy*Jy + Jz*Jz);
      double E  = 0.5*v2 + model->get_pot(r);

      if (E > Emin and E < Emax) {

	orb.new_orbit(E, 0.5);
	double kappa = J/orb.Jmax();

	if (kappa > KTOL and kappa < 1.0 - KTOL and kappa < maxkap) {

	  orb.new_orbit(E, kappa);
	  
	  double apo = orb.apo(), peri = orb.peri();

	  if (apo > rmin and apo < rmax) {

	    if (r/apo > frac) {
	
	      double theta = acos(part->pos(2)/r);
	      double phi   = atan2(part->pos(1), part->pos(0));

	      double vrad  =
		sin(theta)*cos(phi)*part->vel(0) +
		sin(theta)*sin(phi)*part->vel(1) +
		cos(theta)*part->vel(2) ;
	      
	      double vtan = sqrt(v2 - vrad*vrad);
	      double omg0 = orb.get_freq(0);
	      double omg1 = orb.get_freq(1);
	      double w1   = orb.get_angle(r);
	      
	      if (fabs(vrad/vtan) < vcut) {

		out << std::setw(wid) << r
		    << std::setw(wid) << E
		    << std::setw(wid) << kappa
		    << std::setw(wid) << orb.get_action(1)
		    << std::setw(wid) << orb.get_action(2)
		    << std::setw(wid) << Jx
		    << std::setw(wid) << Jy
		    << std::setw(wid) << Jz
		    << std::setw(wid) << omg0
		    << std::setw(wid) << omg1
		    << std::setw(wid) << r/apo
		    << std::setw(wid) << r/peri
		    << std::setw(wid) << w1
		    << std::setw(wid) << vrad/vtan
		    << std::setw(wid) << phi
		    << std::setw(wid) << theta - 0.5*M_PI
		    << std::setw(wid) << vrad
		    << std::setw(wid) << vtan
		    << std::setw(wid) << part->indx()
		    << std::endl;
		
		good++;
	      }
	    }
	  }
	}
      }
	      
      if (myid==0) ++(*progress);
    }
  }
  
  if (myid==0) {
    MPI_Reduce(MPI_IN_PLACE, &good, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
    std::cout << "Found: " << good << "/" << all << std::endl;
  } else {
    MPI_Reduce(&good, 0, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  }

  MPI_Finalize();

  return 0;
}
  
