#include <stdlib.h>
#include <sys/types.h>
#include <getopt.h>
#include <time.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <Timer.h>

#include "expand.h"
#include "UserTreeDSMC.H"
#include "CollideLTE.H"

using namespace std;

//
// Physical units
//

static double pc = 3.086e18;		// cm
static double a0 = 2.0*0.054e-7;	// cm (2xBohr radius)
static double boltz = 1.381e-16;	// cgs
static double year = 365.25*24*3600;	// seconds
static double mp = 1.67e-24;		// g
static double msun = 1.989e33;		// g

double UserTreeDSMC::Lunit = 3.0e5*pc;
double UserTreeDSMC::Tunit = 5e8*year;
double UserTreeDSMC::Vunit = Lunit/Tunit;
double UserTreeDSMC::Munit = 1.0e12*msun/mp;
double UserTreeDSMC::Eunit = Munit*Vunit*Vunit;

UserTreeDSMC::UserTreeDSMC(string& line) : ExternalForce(line)
{
				// Default parameter values
  ncell = 64;
  cnum = 0;
  diamfac = 1.0;
  taufrac = 0.2;
  boxsize = 1.0;
  comp_name = "gas disk";
  nsteps = -1;

				// Initialize using input parameters
  initialize();

				// Look for the fiducial component
  bool found = false;
  list<Component*>::iterator cc;
  Component *c;
  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    c = *cc;
    if ( !comp_name.compare(c->name) ) {
      c0 = c;
      found = true;
      break;
    }
  }

  Vunit = Lunit/Tunit;
  Eunit = Munit*Vunit*Vunit;

				// Diameter*Bohr radius in Lunits
  double diam = diamfac*a0/(Lunit);

  pHOT::s[0] = pHOT::s[1] = pHOT::s[2] = boxsize;

  pCell::bucket = ncell;

  volume = pHOT::s[0] * pHOT::s[1] * pHOT::s[2];


  Collide::CNUM = cnum;
  collide = new CollideLTE(diam, nthrds);

  //
  // Timers: set precision to microseconds
  //
  
  driftTime.Microseconds();
  bcTime.Microseconds();
  partnTime.Microseconds();
  treeTime.Microseconds();
  collideTime.Microseconds();
  densityTime.Microseconds();
  gatherTime.Microseconds();

}

UserTreeDSMC::~UserTreeDSMC()
{
  delete collide;
}

void UserTreeDSMC::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine TreeDSMC initialized, "
       << "Lunit=" << Lunit << ", Tunit=" << Tunit << ", Munit=" << Munit
       << ", cnum=" << cnum << ", diamfac=" << diamfac 
       << ", taufrac=" << taufrac << ", compname=" << comp_name;
  if (nsteps>0) cout << ", with diagnostic output";
  cout << endl;

  print_divider();
}

void UserTreeDSMC::initialize()
{
  string val;

  if (get_value("Lunit", val))		Lunit = atof(val.c_str());
  if (get_value("Tunit", val))		Tunit = atof(val.c_str());
  if (get_value("Munit", val))		Munit = atof(val.c_str());
  if (get_value("cnum", val))		cnum = atoi(val.c_str());
  if (get_value("diamfac", val))	diamfac = atof(val.c_str());
  if (get_value("taufrac", val))	taufrac = atof(val.c_str());
  if (get_value("nsteps", val))		nsteps = atoi(val.c_str());
  if (get_value("compname", val))	comp_name = val;
}


void UserTreeDSMC::determine_acceleration_and_potential(void)
{
  static bool firstime = true;

  //
  // Make the cells
  //

  if (firstime) {
    tree.sendParticles(c0->Particles());
    tree.makeTree();
    firstime = false;
  }

  // DEBUG
  // tree.densCheck();
  
  //
  // Compute time step
  //
  double minvol = tree.minVol();
  double medianvol = tree.medianVol();
  double minsize = pow(minvol, 0.3333333);
  double mediansize = pow(medianvol, 0.3333333);
  double tau = taufrac*mediansize/sqrt(max<double>(T1, T2));
  double ElostTot = 0.0;

  MPI_Bcast(&tau, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (myid==0) 
    cout << setw(10) << "Length: " << setw(16) << mediansize << endl
	 << setw(10) << "Tau: "    << setw(16) << tau        << endl;

  TimeElapsed driftSoFar, bcSoFar, partnSoFar;
  TimeElapsed treeSoFar, collideSoFar, densitySoFar, gatherSoFar;


  //
  // Test repartition
  //
  partnTime.start();
  tree.Repartition();
  partnSoFar = partnTime.stop();
  
  //
  // Sort the particles into cells
  //
  treeTime.start();
  tree.makeTree();
  treeSoFar = treeTime.stop();
  
  //
  // Evaluate collisions among the particles
  //
  unsigned col = collide->collide(tree, Fn, tau);
    
  collideSoFar = collideTime.stop();

				// Debug
  // collide->Debug(tau*(n+1));

  //
  // Periodically display the current progress
  //
  //
  if(nsteps>1 && n%nsteps == 0 ) {
    
    static unsigned tmpc = 0;
    unsigned medianNumb = collide->medianNumber();
    unsigned medianColl = collide->medianColl();

    ostringstream sout;
    sout << "tmp.out." << tmpc++;
    string filen = sout.str();
    tree.testFrontier(filen);

    pHOT_iterator ic(tree);   // Begin frontier iteration to add up KE

    double KEtot1=0.0, KEtot=0.0, totE, dspE;
    unsigned ncells = tree.Number();
    for (unsigned j=0; j<ncells; j++) {
      ic.nextCell();
      ic.KE(totE, dspE);
      KEtot1 += totE;
    }
    double Elost1=collide->Elost(), Elost=0.0;

    MPI_Reduce(&KEtot1, &KEtot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Elost1, &Elost, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (myid==0) {
      ostringstream sout;
      sout << runtag << ".DSMC_log";
      ofstream mout(sout.str().c_str(), ios::app);

      mout << "Finished " 
	   << setw(digit1) << n << " of " << nsteps << " steps; " 
	   << setw(digit2) << collide->total() << " coll; "
	   << setw(digit2) << collide->errors() << " coll errs; "
	   << medianNumb << " num/cell; "
	   << medianColl << " coll/cell; "
	   << tree.Number() << " cells";
      mout << endl << endl;
	
      ElostTot += Elost;

      mout << "Energy (system):" << endl
	   << "       Lost=" << Elost << endl
	   << " Total loss=" << ElostTot << endl
	   << "   Total KE=" << KEtot << endl
	   << "      Ratio=" << Elost/KEtot << endl << endl
	   << "Timing (secs):" << endl
	   << "      drift=" << driftSoFar()*1.0e-6 << endl
	   << "         bc=" << bcSoFar()*1.0e-6 << endl
	   << "  partition=" << partnSoFar()*1.0e-6 << endl
	   << "       tree=" << treeSoFar()*1.0e-6 << endl
	   << "    collide=" << collideSoFar()*1.0e-6 << endl
	   << "    density=" << densitySoFar()*1.0e-6 << endl
	   << "     gather=" << gatherSoFar()*1.0e-6 << endl
	   << endl;

      driftTime.reset();
      bcTime.reset();
      partnTime.reset();
      treeTime.reset();
      collideTime.reset();
      densityTime.reset();
    }

  }

}


extern "C" {
  ExternalForce *makerTreeDSMC(string& line)
  {
    return new UserTreeDSMC(line);
  }
}

class proxytreedsmc { 
public:
  proxytreedsmc()
  {
    factory["usertreedsmc"] = makerTreeDSMC;
  }
};

proxytreedsmc p;
