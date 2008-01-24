#include <sys/timeb.h>
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

extern "C"
void *
dsmc_thread_call(void *atp)
{
  thrd_pass_dsmc *tp = (thrd_pass_dsmc *)atp;
  UserTreeDSMC *p = (UserTreeDSMC *)tp->p;
  p->timestep_thread((void*)&tp->arg);
  return NULL;
}

void UserTreeDSMC::dsmc_thread_fork(pHOT* tree)
{
  int errcode;
  void *retval;
  
  td = new thrd_pass_dsmc [nthrds];
  t = new pthread_t [nthrds];

  if (!td) {
    cerr << "Process " << myid 
         << ": dsmc_thread_fork: error allocating memory for thread counters\n";
    exit(18);
  }
  if (!t) {
    cerr << "Process " << myid
         << ": dsmc_thread_fork: error allocating memory for thread\n";
    exit(18);
  }

                                // Make the <nthrds> threads
  for (int i=0; i<nthrds; i++) {
    td[i].p = this;
    td[i].arg.tree = tree;
    td[i].arg.id = i;

    errcode =  pthread_create(&t[i], 0, dsmc_thread_call, &td[i]);
    if (errcode) {
      cerr << "Process " << myid;
      cerr << " UserTreeDSMC: cannot make thread " << i
	   << ", errcode=" << errcode << endl;
      exit(19);
    }
  }
    
                                // Collapse the threads
  for (int i=0; i<nthrds; i++) {
    if ((errcode=pthread_join(t[i], &retval))) {
      cerr << "Process " << myid;
      cerr << " UserTreeDSMC: thread join " << i
           << " failed, errcode=" << errcode << endl;
      exit(20);
    }
  }
  
  delete [] td;
  delete [] t;
}


void * UserTreeDSMC::timestep_thread(void * arg)
{
  pHOT *tree = (pHOT*)((dsmc_pass_arguments*)arg)->tree;
  int id = (int)((dsmc_pass_arguments*)arg)->id;

  // Loop over cells, cell time-of-flight time
  // for each particle

  unsigned ncells = tree->Number(), nbods;
  double L, vel;

  pHOT_iterator pt(*tree);

  for (unsigned j=0; j<ncells; j++) {
    nbods = pt.nextCell();
    if ( (int)(j%nthrds) == id ) {
      L = pow(pt.Volume(), 0.3333333333333);
      for (unsigned i=0; i<nbods; i++) {
	vel = 1.0e-40;
	for (unsigned k=0; k<3; k++)
	  vel = max<double>(vel, fabs(pt.Body(i)->vel[k]));
	pt.Body(i)->dtreq = L/vel;
      }
    }
  }

  return (NULL);
}


//
// Physical units
//

static double pc = 3.086e18;		// cm
static double a0 = 2.0*0.054e-7;	// cm (2xBohr radius)
// static double boltz = 1.381e-16;	// cgs
// static double year = 365.25*24*3600;	// seconds
static double mp = 1.67e-24;		// g
static double msun = 1.989e33;		// g

double UserTreeDSMC::Lunit = 3.0e5*pc;
double UserTreeDSMC::Munit = 1.0e12*msun;
double UserTreeDSMC::Tunit = sqrt(Lunit*Lunit*Lunit/(Munit*6.673e-08));
double UserTreeDSMC::Vunit = Lunit/Tunit;
double UserTreeDSMC::Eunit = Munit*Vunit*Vunit;

UserTreeDSMC::UserTreeDSMC(string& line) : ExternalForce(line)
{
  id = "TreeDSMC";		// ID string

				// Default parameter values
  ncell = 64;
  cnum = 0;
  diamfac = 1.0;
  boxsize = 1.0;
  comp_name = "gas disk";
  nsteps = -1;
  use_temp = -1;
  use_dens = -1;
  frontier = false;
  mfpstat = false;
  use_multi = false;
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
  diam = diamfac*a0/(Lunit);
				// Number of protons per mass unit
  collfrac = Munit/mp;

  pHOT::sides[0] = pHOT::sides[1] = pHOT::sides[2] = 2.0*boxsize;

  pCell::bucket = ncell;

  volume = pHOT::sides[0] * pHOT::sides[1] * pHOT::sides[2];


  Collide::CNUM = cnum;
  collide = new CollideLTE(diam, nthrds);
  collide->set_temp_dens(use_temp, use_dens);
  ElostTot = 0.0;

  //
  // Timers: set precision to microseconds
  //
  
  partnTime.Microseconds();
  treeTime.Microseconds();
  collideTime.Microseconds();

  if (mfpstat) {
    quant.push_back(0.05);
    quant.push_back(0.1);
    quant.push_back(0.2);
    quant.push_back(0.5);
    quant.push_back(0.8);
    quant.push_back(0.9);
    quant.push_back(0.95);
  }       

  userinfo();
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
       << ", cnum=" << cnum << ", diamfac=" << diamfac << ", diam=" << diam
       << ", compname=" << comp_name;
  if (nsteps>0) cout << ", with diagnostic output";
  if (use_temp>=0) cout << ", temp at pos=" << use_temp;
  if (use_dens>=0) cout << ", dens at pos=" << use_dens;
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
  if (get_value("nsteps", val))		nsteps = atoi(val.c_str());
  if (get_value("compname", val))	comp_name = val;
  if (get_value("use_temp", val))	use_temp = atoi(val.c_str());
  if (get_value("use_dens", val))	use_dens = atoi(val.c_str());
  if (get_value("frontier", val))	frontier = atoi(val.c_str()) ? true : false;
  if (get_value("mfpstat", val))	mfpstat = atoi(val.c_str()) ? true : false;
  if (get_value("use_multi", val))	use_multi = atoi(val.c_str()) ? true : false;
}


void UserTreeDSMC::determine_acceleration_and_potential(void)
{
  static bool firstime = true;
  static unsigned nrep = 0;

  //
  // Make the cells
  //

  if (firstime) {
    c0->Tree()->Repartition();
    MPI_Barrier(MPI_COMM_WORLD);
    c0->Tree()->makeTree();
    c0->Tree()->checkBounds(2.0, "AFTER makeTree (first time)");

    stepnum = 0;
    curtime = tnow;

    firstime = false;

  } else {

    if (tnow-curtime > 1.0e-12) {
      stepnum++;
      curtime = tnow;
    }

  }

  // DEBUG
  // tree.densCheck();
  
  //
  // Compute time step
  //
  double minvol = c0->Tree()->minVol();
  // double medianvol = c0->Tree()->medianVol();
  // double minsize = pow(minvol, 0.3333333);
  // double mediansize = pow(medianvol, 0.3333333);
  double tau = dtime*mintvl[mlevel]/Mstep;

  MPI_Bcast(&tau, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  TimeElapsed partnSoFar, treeSoFar, collideSoFar;

  //
  // Test repartition
  //
  partnTime.start();
  {
    ostringstream sout;
    sout << "before Repartition [" << nrep << "], " 
	 << __FILE__ << ": " << __LINE__;
    c0->Tree()->checkBounds(2.0, sout.str().c_str());
  }

  c0->Tree()->Repartition();
  {
    ostringstream sout;
    sout << "after Repartition [" << nrep << "], " 
	 << __FILE__ << ": " << __LINE__;
    c0->Tree()->checkBounds(2.0, sout.str().c_str());
  }
  MPI_Barrier(MPI_COMM_WORLD);
  partnSoFar = partnTime.stop();
  
  //
  // Sort the particles into cells
  //
  treeTime.start();
  c0->Tree()->makeTree();
  {
    ostringstream sout;
    sout << "after makeTree [" << nrep << "], " 
	 << __FILE__ << ": " << __LINE__;
    c0->Tree()->checkBounds(2.0, sout.str().c_str());
  }
				// Time step request
  if (use_multi) dsmc_thread_fork(c0->Tree());
  treeSoFar = treeTime.stop();
  
  //
  // Evaluate collisions among the particles
  //

  collideTime.start();

  {
    ostringstream sout;
    sout << "before Collide [" << nrep << "], " 
	 << __FILE__ << ": " << __LINE__;
    c0->Tree()->checkBounds(2.0, sout.str().c_str());
  }

  // unsigned col = 
  collide->collide(*c0->Tree(), collfrac, tau);
    
  {
    ostringstream sout;
    sout << "after Collide [" << nrep << "], " 
	 << __FILE__ << ": " << __LINE__;
    c0->Tree()->checkBounds(2.0, sout.str().c_str());
  }

  collideSoFar = collideTime.stop();

				// Debug
  // collide->Debug(tau*(n+1));

  //
  // Periodically display the current progress
  //
  //
  if(nsteps>1 && stepnum%nsteps == 0 ) {
    
    collide->Debug(tnow);

    unsigned medianNumb = collide->medianNumber();
    unsigned medianColl = collide->medianColl();

    if (frontier) {
      ostringstream sout;
      sout << runtag << ".DSMC.frontier." << stepnum;
      string filen = sout.str();
      c0->Tree()->testFrontier(filen);
    }

    if (mfpstat) {

      collide->mfpsizeQuantile(quant, mfp_, ts_);

      if (myid==0) {

	ostringstream sout;
	sout << runtag << ".DSMC.mfpstat." << stepnum;
	string filen = sout.str();
	ofstream out(filen.c_str());
	if (out) {
	  
	  out << setw(18) << "Quantiles: " 
	      << setw(18) << "MFP/size"
	      << setw(18) << "Flight/size"
	      << endl;
	  for (unsigned nq=0; nq<quant.size(); nq++)
	    out << setw(18) << quant[nq] << ": " 
		<< setw(18) << mfp_[nq] << setw(18) << ts_[nq] << endl;
	}
      }
    }

				// Begin frontier iteration to add up KE
    pHOT_iterator ic(*c0->Tree());

    double KEtot1=collide->Etotal(), KEtot=0.0;
    double Elost1=collide->Elost(),  Elost=0.0;

    MPI_Reduce(&KEtot1, &KEtot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Elost1, &Elost, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (myid==0) {
      ostringstream sout;
      sout << runtag << ".DSMC_log";
      ofstream mout(sout.str().c_str(), ios::app);
      mout << "Summary:" << endl << left
	   << setw(6) << " " << setw(20) << stepnum << "steps" << endl
	   << setw(6) << " " << setw(20) << collide->total() << "collisions" << endl
	   << setw(6) << " " << setw(20) << collide->errors() << "collision errors" << endl
	   << setw(6) << " " << setw(20) << medianNumb << "number/cell" << endl
	   << setw(6) << " " << setw(20) << medianColl << "collision/cell" << endl
	   << setw(6) << " " << setw(20) << c0->Tree()->TotalNumber() << "cells" << endl
	   << endl;
	
      ElostTot += Elost;

      mout << "Energy (system):" << endl
	   << "       Lost=" << Elost << endl
	   << " Total loss=" << ElostTot << endl
	   << "   Total KE=" << KEtot << endl;

      if (KEtot<=0.0)
	mout << "      Ratio= XXXX" << endl << endl;
      else
	mout << "      Ratio=" << Elost/KEtot << endl << endl;

      mout << "Timing (secs):" << endl
	   << "  partition=" << partnSoFar()*1.0e-6 << endl
	   << "       tree=" << treeSoFar()*1.0e-6 << endl
	   << "    collide=" << collideSoFar()*1.0e-6 << endl
	   << endl;

      collide->tsdiag(mout);

      partnTime.reset();
      treeTime.reset();
      collideTime.reset();
    }

  }
				// Remake level lists because particles
				// will (usually) have been exchanged 
				// between nodes
  c0->reset_level_lists();

  MPI_Barrier(MPI_COMM_WORLD);
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
