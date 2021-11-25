/*
  Compute simple statistics from psp dump

  MDWeinberg 07/24/19
*/

using namespace std;

#include <unistd.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <list>

#include <Species.H>
#include <yaml-cpp/yaml.h>
#include <cxxopts.H>

#include <header.H>
#include <PSP.H>

#include "atomic_constants.H"

int
main(int ac, char **av)
{
  char *prog = av[0];
  double Lunit, Munit, Tunit;
  int sindx, eindx, icons, econs;
  std::string cname, species, new_dir("./");
  bool verbose = false;

  // Parse command line
  //
  cxxopts::Options options(prog, "Compute simple statistics from psp dump\n");

  options.add_options()
   ("h,help", "produce help message")
   ("v,verbose", "verbose output")
   ("species", "position of species index",
     cxxopts::value<std::string>(species)->default_value("species.yml"))
   ("e,electrons", "position of electron index",
     cxxopts::value<int>(eindx)->default_value("10"))
   ("I,consI", "position of ion conservation (-1 to ignore)",
     cxxopts::value<int>(icons)->default_value("-1"))
   ("E,consE", "position of electron conservation (-1 to ignore)",
     cxxopts::value<int>(econs)->default_value("-1"))
   ("c,name", "component name",
     cxxopts::value<std::string>(cname)->default_value("gas"))
   ("d,dir", "rewrite data directory for SPL files",
     cxxopts::value<std::string>(new_dir))
   ("L,Lunit", "physical length unit in pc",
     cxxopts::value<double>(Lunit)->default_value("1.0"))
   ("M,Munit", "physical mass unit in solar masses",
     cxxopts::value<double>(Munit)->default_value("1.0"))
   ("T,Tunit", "physical time unit in years",
     cxxopts::value<double>(Tunit)->default_value("1.0e+03"))
   ("f,files", "input files",
     cxxopts::value< std::vector<std::string> >())
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
	      << " -c gas -s 0 -f OUT.run.00001" << std::endl;
    return 1;
  }

  if (vm.count("verbose")) {
    verbose = true;
  }

  // Convert system units to cgs
  //
  Lunit *= pc;
  Munit *= msun;
  Tunit *= year;

				// Velocity and energy system units
  double Vunit = Lunit/Tunit;
  double Eunit = Munit*Vunit*Vunit;

				// Get species file

  YAML::Node config;
  try {
    config = YAML::LoadFile(species);
  }
  catch (YAML::Exception & error) {
    std::cerr << "pspstatT: error opening species file <"
	      << species << "> . . . quitting" << std::endl;
    exit(1);
  }

  if (not config["species_map"]) {
    std::cerr << "pspstatT: no species map in config"
	      << " . . . quitting"
	      << std::endl;
    exit(-1);
  }
  
  std::string type;

  if (config["species_map"]["type"]) {
    type = config["species_map"]["type"].as<std::string>();
  } else {
    std::cerr << "pspstatT: no <type> key found . . . "
	      << "quitting" << std::endl;
    exit(-2);
  }    

  if (type.compare("trace") != 0) {
    std::cerr << "pspstatT: expected <trace> and found <" << type << "> . . . "
	      << "quitting" << std::endl;
    exit(-3);
  }

  
  int use_cons = -1;

  if (config["species_map"]["cons"]) {
    use_cons = config["species_map"]["cons"].as<int>();
  } else {
    std::cerr << "pspstatT: no <cons> key found . . . "
	      << "quitting" << std::endl;
    exit(-4);
  }

  int use_elec = -1;

  if (config["species_map"]["elec"]) {
    use_elec = config["species_map"]["elec"].as<int>();
  } else {
    std::cerr << "pspstatT: no <elec> key found . . . "
	      << "quitting" << std::endl;
    exit(-5);
  }

  std::map<speciesKey, int> SpList;
  std::set<unsigned short>  ZList;

  if (config["species_map"]["elements"]) {
    YAML::Node elems = config["species_map"]["elements"];

    for (YAML::const_iterator it=elems.begin(); it!=elems.end(); ++it) {
      auto dd = it->as<std::vector<int>>();
      speciesKey key(dd[0], dd[1]);
      SpList[key] = dd[2];
      ZList.insert(key.first);
    }
  } else {
    std::cerr << "pspstatT: no <elements> key found . . . "
	      << "quitting" << std::endl;
    exit(-6);
  }


  // Get file arguments
  //
  std::vector<std::string> files = vm["files"].as< std::vector<std::string> >();

  bool first = true;
  
  for (auto file : files ) {

    if (verbose) cerr << "Using filename: " << file << endl;

				// Parse the PSP file
				// ------------------
    std::shared_ptr<PSP> psp;
    if (file.find("SPL") != std::string::npos)
      psp = std::make_shared<PSPout>(file);
    else
      psp = std::make_shared<PSPspl>(file, new_dir);

				// Now write a summary
				// -------------------
    if (verbose) psp->PrintSummary(cerr);

				// Will contain array for each gas species
				// ---------------------------------------
    double mass=0.0, KEi=0.0, KEe=0.0, Icons=0.0, Econs=0.0, NumbI=0.0, NumbE=0.0;
    unsigned int N=0;

				// Quantiles for printing rank
				// ---------------------------
    const std::vector<double> pval = {0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.99};

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

				// Phase space stuff
				// -----------------
      double ms;
      double pos[3];
      double vel[3];
      double mom[3];
      double pot;

				// Print the header

      std::cout << "Comp name: " << stanza->name << endl << endl
		<< "     Bodies\t\t" << std::left
		<< setw(15) << stanza->comp.nbod 
		<< setw(15) << stanza->comp.niatr 
		<< setw(15) << stanza->comp.ndatr 
		<< endl;

      std::vector<double> evel, ivel;

      for (part=psp->GetParticle(); part!=0; part=psp->NextParticle()) {

	//
	// Accumulate statistics
	//

				// Molectural weight and electron
				// fraction
	double Mu=0.0, Eta=0.0;
	for (auto v : SpList) {
	  speciesKey k = v.first;
	  double   one = part->datr(v.second)/atomic_mass[k.first];
	  Mu  += one;
	  Eta += one * (k.second - 1);
	}
	Eta /= Mu;
	Mu   = 1.0/Mu;

	double ms = part->mass();
	mass1 += ms;

	// Angular momentum
	//
	mom[0] = part->pos(1)*part->vel(2) - part->pos(2)*part->vel(1);
	mom[1] = part->pos(2)*part->vel(0) - part->pos(0)*part->vel(2);
	mom[2] = part->pos(0)*part->vel(1) - part->pos(1)*part->vel(0);

	for (int i=0; i<3; i++) com1[i] += ms*part->pos(i);
	for (int i=0; i<3; i++) cov1[i] += ms*part->vel(i);
	for (int i=0; i<3; i++) dsp1[i] += ms*part->vel(i)*part->vel(i);
	for (int i=0; i<3; i++) ang1[i] += ms*mom[i];

	rtmp = 0.0;
	for (int i=0; i<3; i++) rtmp += part->vel(i)*part->vel(i);
	KE1 += 0.5*ms*rtmp;
	PE1 += 0.5*ms*part->phi();
	Iv2 += ms*rtmp;

	double ke = 0.0;
	for (int i=0; i<3; i++) {
	  double t = part->datr(use_elec+i);
	  covE[i] += ms*t;
	  dspE[i] += ms*t*t;
	  ke += t*t;
	}

	EE1 += ke * 0.5*ms*atomic_mass[0]*Eta/Mu;
	Ev2 += ms * ke;

	ivel.push_back(0.5*Mu*amu*rtmp*Vunit*Vunit/eV);
	evel.push_back(0.5*me*ke*Vunit*Vunit/eV);

	// Mass
	mass += ms;

	// Ion KE
	KEi  += 0.5 * ms * rtmp * Eunit;

	// Electron KE
	KEe  += 0.5 * ke * ms * atomic_mass[0] * Eta/Mu * Eunit;

	// True ion number
	NumbI += ms * Munit/(Mu*amu);

	// True ion number
	NumbE += ms * Munit/(Mu*amu) * Eta;

	// Ion energy conservation
	if (icons>=0) Icons += part->datr(icons) * Eunit;
	
	// Electron energy conservation
	if (econs>=0) Econs += part->datr(econs) * Eunit;

	// Superparticle count
	N++;
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

      cout << "     Mass\t\t"               << mass          << std::endl
	   << "     KE(ion)\t\t"            << KEi           << std::endl
	   << "     KE(elec)\t\t"           << KEe           << std::endl
	   << "     T(ion)\t\t"             << KEi/(1.5*NumbI*boltz) << std::endl
	   << "     T(elec)\t\t"            << KEe/(1.5*NumbE*boltz) << std::endl
	   << "     N(ion)\t\t"             << NumbI         << std::endl
	   << "     N(elec)\t\t"            << NumbE         << std::endl;
      if (icons>=0) {
	cout << "     dE(ion)\t\t"          << Icons         << std::endl;
	cout << "     dE/E(ion)\t\t"        << Icons/KEi     << std::endl;
      }
      if (econs>=0) {
	cout << "     dE(elec)\t\t"         << Econs         << std::endl;
	cout << "     dE/E(elec)\t\t"       << Econs/KEe     << std::endl;
      }
      cout   << "     Total energy\t"       << KEi + KEe + Econs + Icons << std::endl;
      cout   << "     N(spart)\t\t"         << N             << std::endl;

      std::cout << std::endl
		<< "------------------------" << std::endl
		<< "Ranked energies" << std::endl
		<< "------------------------" << std::endl << std::left
		<< std::setw(8) << "Rank" << std::setw(18) << "Ion (eV)"
		<< std::setw(18) << "Elec (eV)"
		<< std::endl
		<< std::setw(8) << "----" << std::setw(18) << "-----------"
		<< std::setw(18) << "-----------"
		<< std::endl;
      
      std::sort(ivel.begin(), ivel.end());
      std::sort(evel.begin(), evel.end());

      for (int i=0; i<5; i++) {
	std::cout << std::setw( 8) << i
		  << std::setw(18) << ivel[i]
		  << std::setw(18) << evel[i] << std::endl;
      }

      for (auto v : pval) {
	int I = floor(v*evel.size());
	std::cout << std::setw( 8) << I
		  << std::setw(18) << ivel[I]
		  << std::setw(18) << evel[I] << std::endl;
      }

      for (int i=evel.size()-5; i<evel.size(); i++) {
	std::cout << std::setw( 8) << i
		  << std::setw(18) << ivel[i]
		  << std::setw(18) << evel[i] << std::endl;
      }

    }

    std::cout << std::endl;
  }
  
  return 0;
}
  
