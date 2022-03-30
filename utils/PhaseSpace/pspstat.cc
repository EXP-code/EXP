/*
  Compute simple statistics from psp dump

  MDWeinberg 06/10/02, 11/24/19
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <memory>
#include <limits>
#include <list>

#include <Eigen/Eigen>

#include <StringTok.H>
#include <libvars.H>		// EXP library globals
#include <cxxopts.H>		// Option parsing
#include <header.H>
#include <PSP.H>

int
main(int argc, char **argv)
{
  char *prog = argv[0];
  bool verbose = false;
  std::string new_dir;

  // Parse command line

  cxxopts::Options options("ascii2psp", "Construct a PSP file from ascii input files");

  std::vector<std::string> pos_names = {"file"};

  options.parse_positional(pos_names.begin(), pos_names.end());

  options.add_options()
    ("h,help", "print this help message")
    ("v,verbose", "print verbose output messages")
    ("d,dir", "directory containing PSP files",
     cxxopts::value<std::string>(new_dir)->default_value("./"))
    ("f,file", "PSP file",
     cxxopts::value<std::string>(pos_names[0]))
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
    exit(-1);
  }

  if (vm.count("verbose")) {
    verbose = true;
  }



  std::ifstream in;
  std::string filename = vm["file"].as<std::string>();
  
				// Parse the PSP file
				// ------------------
  PSPptr psp;
  if (filename.find("SPL") != std::string::npos)
    psp = std::make_shared<PSPspl>(filename, new_dir);
  else
    psp = std::make_shared<PSPout>(filename);
  

				// Now write a summary
				// -------------------
  if (verbose) {
    psp->PrintSummary(cerr);
    
    std::cout << std::endl << "PSP file <" << filename << "> has time <" 
	      << psp->CurrentTime() << ">" << std::endl;
  }

				// Setup stats for all components
				// ------------------------------

  Eigen::Vector3d com = Eigen::Vector3d::Zero();
  Eigen::Vector3d cov = Eigen::Vector3d::Zero();
  Eigen::Vector3d ang = Eigen::Vector3d::Zero();

  double KE     = 0.0;
  double PE     = 0.0;
  double mass   = 0.0;
  int   totbod  = 0;
  
  PSPstanza *stanza;
  SParticle* part;
  double rtmp;

  for (stanza=psp->GetStanza(); stanza!=0; stanza=psp->NextStanza()) {


				// Setup stats for each component
				// -----------------------------

    Eigen::Matrix3d moments2 = Eigen::Matrix3d::Zero();
    Eigen::Vector3d com1     = Eigen::Vector3d::Zero();
    Eigen::Vector3d cov1     = Eigen::Vector3d::Zero();
    Eigen::Vector3d ang1     = Eigen::Vector3d::Zero();
    Eigen::Vector3d pmn1, pmx1;

    // Set to maximum
    pmn1 << std::numeric_limits<double>::max(),
      std::numeric_limits<double>::max(),
      std::numeric_limits<double>::max();

    // Set to minimum
    pmx1 << -std::numeric_limits<double>::max(),
      -std::numeric_limits<double>::max(),
      -std::numeric_limits<double>::max();

    double KE1     = 0.0;
    double PE1     = 0.0;
    double mass1   = 0.0;

				// Phase space stuff
				// -----------------
    double ms;
    double pos[3];
    double vel[3];
    double mom[3];
    double pot;
				// Print the header

    cout << "Comp name: " << stanza->name << endl
	 << "     Bodies:\t\t"
	 << setw(15) << stanza->comp.nbod 
	 << setw(10) << stanza->comp.niatr 
	 << setw(10) << stanza->comp.ndatr 
	 << endl;

    totbod += stanza->comp.nbod;

    for (part=psp->GetParticle(); part!=0; part=psp->NextParticle()) {

      mom[0] = part->pos(1)*part->vel(2) - part->pos(2)*part->vel(1);
      mom[1] = part->pos(2)*part->vel(0) - part->pos(0)*part->vel(2);
      mom[2] = part->pos(0)*part->vel(1) - part->pos(1)*part->vel(0);

				// Accumulate statistics
      double ms = part->mass();
      mass1 += ms;
      for (int i=0; i<3; i++) com1[i] += ms*part->pos(i);
      for (int i=0; i<3; i++) cov1[i] += ms*part->vel(i);
      for (int i=0; i<3; i++) ang1[i] += ms*mom[i];
      rtmp = 0.0;
      for (int i=0; i<3; i++) rtmp += part->vel(i)*part->vel(i);
      KE1 += 0.5*ms*rtmp;
      PE1 += 0.5*ms*part->phi();

      for (int i=0; i<3; i++) {
	pmn1[i] = std::min<double>(pmn1[i], part->pos(i));
	pmx1[i] = std::max<double>(pmx1[i], part->pos(i));
      }

      for (int i=0; i<3; i++) {
	for (int j=0; j<3; j++) {
	  moments2(i, j) += ms*part->pos(i)*part->pos(j);
	}
      }

    }
    
    if (mass1>0.0) {

      cout  << "     MIN:\t\t";
      for (int i=0; i<3; i++) cout << setw(15) << pmn1[i];
      cout << endl;
      cout  << "     MAX:\t\t";
      for (int i=0; i<3; i++) cout << setw(15) << pmx1[i];
      cout << endl;
      cout  << "     COM:\t\t";
      for (int i=0; i<3; i++) cout << setw(15) << com1[i]/mass1;
      cout << endl;
      cout  << "     COV:\t\t";
      for (int i=0; i<3; i++) cout << setw(15) << cov1[i]/mass1;
      cout << endl;
      cout  << "     Ang mom:\t\t";
      for (int i=0; i<3; i++) cout << setw(15) << ang1[i];
      cout << endl;
      cout  << "     Stats:\t\tKE=" << KE1 << " PE=" << PE1 << " -2T/W=" << -2.0*KE1/PE1
	    << " Mass=" << mass1 << endl;

      moments2 /= mass1;
      
      for (int i=0; i<3; i++) {
	for (int j=0; j<3; j++) {
	  moments2(i, j) += -com1[i]*com1[j]/(mass1*mass1);
	}
      }

      cout << "2nd moment matrix:" << endl << moments2 << endl;
      Eigen::JacobiSVD<Eigen::Matrix3d> svd(moments2, Eigen::ComputeThinU | Eigen::ComputeThinV);
      cout << "Singular values:" << endl << svd.singularValues() << endl;
      cout << "Left vectors:" << endl << svd.matrixU() << endl;
      cout << "Right vectors:" << endl << svd.matrixV() << endl;
    }

    mass += mass1;
    for (int i=0; i<3; i++) com[i] += com1[i];
    for (int i=0; i<3; i++) cov[i] += cov1[i];
    for (int i=0; i<3; i++) ang[i] += ang1[i];
    KE += KE1;
    PE += PE1;
  }
  
  cout << endl << "Total:" << endl
       << "     Bodies:\t\t"
       << setw(15) << totbod << endl; 
  cout  << "     COM:\t\t";
  for (int i=0; i<3; i++) cout << setw(15) << com[i]/mass;
  cout << endl;
  cout  << "     COV:\t\t";
  for (int i=0; i<3; i++) cout << setw(15) << cov[i]/mass;
  cout << endl;
  cout  << "     Ang mom:\t\t";
  for (int i=0; i<3; i++) cout << setw(15) << ang[i];
  cout << endl;
  cout  << "     Stats:\t\tKE=" << KE << " PE=" << PE << " -2T/W=" << -2.0*KE/PE
	<< " Mass=" << mass << endl;
  
  return 0;
}
