/*
  Compute simple statistics from psp dump

  MDWeinberg 06/10/02
*/

#include <unistd.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <tuple>
#include <list>

using namespace std;

#include "Species.H"

#include "StringTok.H"
#include "cxxopts.H"		// Option parsing
#include "libvars.H"		// EXP library globals
#include "header.H"
#include "PSP.H"

#include "atomic_constants.H"

int
main(int ac, char **av)
{
  char *prog = av[0];
  double Lunit, Munit, Tunit;
  int sindx, eindx, icons, econs;
  std::string cname;
  bool verbose = false, trace = false;

  // Parse command line

  cxxopts::Options options(prog, "Compute simple statistics from psp dump");

  options.add_options()
   ("h,help", "produce help message")
   ("v,verbose", "verbose output")
   ("trace", "trace element method")
   ("OUT", "assume that PSP files are in original format")
   ("SPL", "assume that PSP files are in split format")
   ("s,species", "position of species index",
     cxxopts::value<int>(sindx)->default_value("-1"))
   ("e,electrons", "position of electron index",
     cxxopts::value<int>(eindx)->default_value("10"))
   ("I,consI", "position of ion conservation (-1 to ignore)",
     cxxopts::value<int>(icons)->default_value("-1"))
   ("E,consE", "position of electron conservation (-1 to ignore)",
     cxxopts::value<int>(econs)->default_value("-1"))
   ("c,name", "component name",
     cxxopts::value<std::string>(cname)->default_value("gas"))
   ("L,Lunit", "physical length unit in pc",
     cxxopts::value<double>(Lunit)->default_value("1.0"))
   ("M,Munit", "physical mass unit in solar masses",
     cxxopts::value<double>(Munit)->default_value("0.1"))
   ("T,Tunit", "physical time unit in years",
     cxxopts::value<double>(Tunit)->default_value("1.0e+05"))
   ("f,files", "input files",
     cxxopts::value< std::vector<std::string> >())
    ;


  double mu0 = 1.0/(0.76/atomic_mass[1] + 0.24/atomic_mass[2]);

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
	      << " -c gas -s 0 -f OUT.run.00001" << std::endl;
    return 1;
  }

  if (vm.count("verbose")) {
    verbose = true;
  }

  if (vm.count("trace")) {
    trace = true;
  }


  // Convert system units to cgs
  //
  Lunit *= pc;
  Munit *= msun;
  Tunit *= year;

				// Velocity and energy system units
  double Vunit = Lunit/Tunit;
  double Eunit = Munit*Vunit*Vunit;

  // Get file arguments
  //
  std::vector<std::string> files = vm["files"].as< std::vector<std::string> >();

  bool first = true;
  
  for (auto file : files ) {
    
    std::ifstream in(file);
    if (in) {
      std::cerr << "Error opening file <" << file << "> for input" << std::endl;
      exit(-1);
    }
    in.close();

    if (verbose) cerr << "Using filename: " << file << endl;


				// Parse the PSP file
				// ------------------
    PSPptr psp;
    if (vm.count("SPL")) psp = std::make_shared<PSPspl>(file);
    else                 psp = std::make_shared<PSPout>(file);


				// Now write a summary
				// -------------------
    if (verbose) {

      psp->PrintSummary(std::cerr);
    
      std::cerr << std::endl << "Phase space has time <" 
		<< psp->CurrentTime() << ">" << std::endl;
    }
    
				// Will contain array for each gas species
				// ---------------------------------------
    typedef std::tuple<double, double, double, double,
    //                 ^       ^       ^       ^
    //                 |       |       |       |
    // 0: mass --------+       |       |       |
    // 1: Ion KE --------------+       |       |
    // 2: Electron KE -----------------+       |
    // 3: True number -------------------------+
    //
		       double, double, unsigned> shtup;
    //                 ^       ^       ^
    //                 |       |       |
    // 4: Ion cons E---+       |       |
    // 5: Electron cons E------+       |
    // 6: Superparticle count ---------+
    //

    typedef std::map<speciesKey, shtup> shist;
    const shtup tupzero(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0);

    PSPstanza *stanza;
    SParticle* part;
    double rtmp;

    for (stanza=psp->GetStanza(); stanza!=0; stanza=psp->NextStanza()) {
      
      if (stanza->name != cname) continue;


				// Setup stats for each component
				// -----------------------------

      double com1[3] = {0.0, 0.0, 0.0};
      double cov1[3] = {0.0, 0.0, 0.0};
      double dsp1[3] = {0.0, 0.0, 0.0};
      double covE[3] = {0.0, 0.0, 0.0};
      double dspE[3] = {0.0, 0.0, 0.0};
      double ang1[3] = {0.0, 0.0, 0.0};
      double KE1     = 0.0;
      double PE1     = 0.0;
      double EE1     = 0.0;
      double Iv2     = 0.0;
      double Ev2     = 0.0;
      double mass1   = 0.0;

      shist  hist1;		// For gas only

				// Phase space stuff
				// -----------------
      double ms;
      double pos[3];
      double vel[3];
      double mom[3];
      double pot;

				// Print the header

      cout << "Comp name: " << stanza->name << endl << endl
	   << "     Bodies\t\t"
	   << setw(15) << stanza->comp.nbod 
	   << setw(10) << stanza->comp.niatr 
	   << setw(10) << stanza->comp.ndatr 
	   << endl;


      for (part=psp->GetParticle(); part!=0; part=psp->NextParticle()) {

	mom[0] = part->pos(1)*part->vel(2) - part->pos(2)*part->vel(1);
	mom[1] = part->pos(2)*part->vel(0) - part->pos(0)*part->vel(2);
	mom[2] = part->pos(0)*part->vel(1) - part->pos(1)*part->vel(0);

				// Accumulate statistics
	double ms = part->mass();
	mass1 += ms;
	for (int i=0; i<3; i++) com1[i] += ms*part->pos(i);
	for (int i=0; i<3; i++) cov1[i] += ms*part->vel(i);
	for (int i=0; i<3; i++) dsp1[i] += ms*part->vel(i)*part->vel(i);
	for (int i=0; i<3; i++) ang1[i] += ms*mom[i];
	rtmp = 0.0;
	for (int i=0; i<3; i++) rtmp += part->vel(i)*part->vel(i);
	KE1 += 0.5*ms*rtmp;
	PE1 += 0.5*ms*part->phi();
	Iv2 += ms*rtmp;

	if (sindx >= 0) {
	  KeyConvert kc(part->iatr(sindx));
	  speciesKey k = kc.getKey();

	  if (hist1.find(k) == hist1.end()) hist1[k] = tupzero;

	  double ke = 0.0;
	  for (int i=0; i<3; i++) {
	    double t = part->datr(eindx+i);
	    covE[i] += ms*t;
	    dspE[i] += ms*t*t;
	    ke += t*t;
	  }
	  if (trace)
	    EE1 += ke * 0.5*ms*atomic_mass[0]/mu0;
	  else
	    EE1 += ke * 0.5*ms*atomic_mass[0]/atomic_mass[k.first];
	  Ev2 += ms * ke;

				// Mass
	  std::get<0>(hist1[k]) += ms;
				// Ion KE
	  std::get<1>(hist1[k]) += 0.5 * ms * rtmp * Eunit;
				// Electron KE
	  if (trace)
	    std::get<2>(hist1[k]) += 0.5 * ke * ms * atomic_mass[0]/mu0 * Eunit;
	  else
	    std::get<2>(hist1[k]) += 0.5 * ke * ms * atomic_mass[0]/atomic_mass[k.first] * Eunit;
				// True particle number
	  if (trace)
	    std::get<3>(hist1[k]) += ms * Munit/(mu0*amu);
	  else
	    std::get<3>(hist1[k]) += ms * Munit/(atomic_mass[k.first]*amu);
				// Ion energy conservation
	  if (icons>=0) std::get<4>(hist1[k]) += part->datr(icons) * Eunit;
				// Electron energy conservation
	  if (econs>=0) std::get<5>(hist1[k]) += part->datr(econs) * Eunit;
				// Superparticle count
	  std::get<6>(hist1[k]) ++;
	}

      }
    
      cout  << "     COM\t\t";
      for (int i=0; i<3; i++) cout << setw(15) << com1[i]/mass1;
      cout << endl;
      cout  << "     COVi\t\t";
      for (int i=0; i<3; i++) cout << setw(15) << cov1[i]/mass1;
      cout << endl;
      cout  << "     DSPi\t\t";
      for (int i=0; i<3; i++) cout << setw(15) << std::sqrt(dsp1[i]/mass1 - cov1[i]*cov1[i]/mass1/mass1);
      cout << endl;
      cout  << "     COVe\t\t";
      for (int i=0; i<3; i++) cout << setw(15) << covE[i]/mass1;
      cout << endl;
      cout  << "     DSPe\t\t";
      for (int i=0; i<3; i++) cout << setw(15) << std::sqrt(dspE[i]/mass1 - covE[i]*covE[i]/mass1/mass1);
      cout << endl;
      cout  << "     Ang mom\t\t";
      for (int i=0; i<3; i++) cout << setw(15) << ang1[i];
      cout << endl << endl;
      
      // Get the system time from the PSP header
      //
      double systime = psp->CurrentTime();

      cout << "     System time\t\t"        << systime       << std::endl
	   << "     Kinetic energy\t\t"     << KE1           << std::endl
	   << "     Potential energy\t\t"   << PE1           << std::endl
	   << "     Virial ratio\t\t"       << -2.0*KE1/PE1  << std::endl
	   << "     Electron energy\t\t"    << EE1           << std::endl
	   << "     Particle mass\t\t"      << mass1         << std::endl
	   << "     Mean ion vel\t\t"       << sqrt(Iv2/mass1) << std::endl
	   << "     Mean elec vel\t\t"      << sqrt(Ev2/mass1) << std::endl
	   << "     Mean ion-elec vel\t\t"  << sqrt((Iv2+Ev2)/mass1) << std::endl
	   << std::endl;

      if (sindx >= 0) {
				// Header
	cout   << std::right
	       << std::setw(16) << "Species" << " : "
	       << std::setw(16) << "Mass"
	       << std::setw(16) << "Ion temp"
	       << std::setw(16) << "Elec temp";

	if (icons>=0)
	  cout << std::setw(16) << "dE(Ion)/KE";

	if (econs>=0)
	  cout << std::setw(16) << "dE(Elec)/KE";

	cout   << std::setw(10) << "Count"
	       << endl;
				// Species loop
	for (auto v : hist1) {
	  speciesKey k = v.first;
	  std::ostringstream sout;
	  sout   << "<" << k.first << "," << k.second << ">";
	  cout   << std::right << std::setw(16) << sout.str()
		 << " : "
		 << std::setw(16) << std::get<0>(v.second)
		 << std::setw(16) << std::get<1>(v.second)/(1.5*std::get<3>(v.second)*boltz)
		 << std::setw(16) << std::get<2>(v.second)/(1.5*std::get<3>(v.second)*boltz);

	  if (icons>=0)
	    cout << std::setw(16) << std::get<4>(v.second)/std::get<1>(v.second);

	  if (econs>=0)
	    cout << std::setw(16) << std::get<5>(v.second)/std::get<2>(v.second);
	  
	  cout   << std::setw(10) << std::get<6>(v.second)
	         << std::endl;
	}
      }
    }
    std::cout << std::endl;
  }
  
  return 0;
}
  
