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
  
  if (nthrds==1) {
    thrd_pass_dsmc td;

    td.p = this;
    td.arg.tree = tree;
    td.arg.id = 0;

    dsmc_thread_call(&td);
    
    return;
  }

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
  double L, vel, DT;

  pHOT_iterator pt(*tree);

  for (unsigned j=0; j<ncells; j++) {
    nbods = pt.nextCell();
    if ( (int)(j%nthrds) == id ) {
      L = pow(pt.Volume(), 0.3333333333333);
      for (unsigned i=0; i<nbods; i++) {
				// Time of flight criterion
	vel = 1.0e-40;
	for (unsigned k=0; k<3; k++)
	  vel = max<double>(vel, fabs(pt.Body(i)->vel[k]));
	DT = L/vel;
				// Cooling criterion
	int sz = pt.Body(i)->dattrib.size();
	if (use_delt>=0 && use_delt<sz)
	  DT = min<double>(DT, coolfrac*pt.Body(i)->dattrib[use_delt]);
	
	pt.Body(i)->dtreq = DT;
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
  epsm = -1.0;
  diamfac = 1.0;
  boxsize = 1.0;
  boxratio = 1.0;
  comp_name = "gas disk";
  nsteps = -1;
  use_temp = -1;
  use_dens = -1;
  use_delt = -1;
  coolfrac = 0.1;
  frontier = false;
  mfpstat = false;
  use_multi = false;
  use_pullin = false;
  cba = true;
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
  pHOT::offst[0] = pHOT::offst[1] = pHOT::offst[2] = boxsize;

				// For cylindrical disk
  pHOT::sides[2] = boxratio;
  pHOT::offst[2] = boxratio;

  pCell::bucket = ncell;

  volume = pHOT::sides[0] * pHOT::sides[1] * pHOT::sides[2];


  Collide::CBA = cba;
  Collide::PULLIN = use_pullin;
  Collide::CNUM = cnum;
  Collide::EPSMratio = epsm;
  collide = new CollideLTE(diam, nthrds);
  collide->set_temp_dens(use_temp, use_dens);
  collide->set_timestep(use_delt);
  ElostTotCollide = ElostTotEPSM = 0.0;

  //
  // Timers: set precision to microseconds
  //
  
  partnTime.Microseconds();
  treeTime.Microseconds();
  collideTime.Microseconds();

  //
  // Quantiles for distribution diagnostic
  //
  quant.push_back(0.01);
  quant.push_back(0.05);
  quant.push_back(0.1);
  quant.push_back(0.2);
  quant.push_back(0.5);
  quant.push_back(0.8);
  quant.push_back(0.9);
  quant.push_back(0.95);
  quant.push_back(0.99);

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
       << ", epsm=" << epsm << ", boxsize=" << boxsize << ", boxratio=" << boxratio
       << ", compname=" << comp_name;
  if (nsteps>0) cout << ", with diagnostic output";
  if (use_temp>=0) cout << ", temp at pos=" << use_temp;
  if (use_dens>=0) cout << ", dens at pos=" << use_dens;
  if (use_pullin)  cout << ", Pullin algorithm enabled";
  if (use_multi) {
    cout << ", multistep enabled";
    if (use_delt>=0) 
      cout << ", time step at pos=" << use_delt << ", coolfrac=" << coolfrac;
  }
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
  if (get_value("epsm", val))		epsm = atof(val.c_str());
  if (get_value("diamfac", val))	diamfac = atof(val.c_str());
  if (get_value("boxsize", val))	boxsize = atof(val.c_str());
  if (get_value("boxratio", val))	boxratio = atof(val.c_str());
  if (get_value("coolfrac", val))	coolfrac = atof(val.c_str());
  if (get_value("nsteps", val))		nsteps = atoi(val.c_str());
  if (get_value("compname", val))	comp_name = val;
  if (get_value("use_temp", val))	use_temp = atoi(val.c_str());
  if (get_value("use_dens", val))	use_dens = atoi(val.c_str());
  if (get_value("use_delt", val))	use_delt = atoi(val.c_str());
  if (get_value("frontier", val))	frontier = atoi(val.c_str()) ? true : false;
  if (get_value("mfpstat", val))	mfpstat = atoi(val.c_str()) ? true : false;
  if (get_value("use_multi", val))	use_multi = atoi(val.c_str()) ? true : false;
  if (get_value("use_pullin", val))	use_pullin = atoi(val.c_str()) ? true : false;
  if (get_value("cba", val))		cba = atoi(val.c_str()) ? true : false;
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
  // double minvol = c0->Tree()->minVol();
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

  //
  // Periodically display the current progress
  //
  //
  if(nsteps>1 && stepnum%nsteps == 0 ) {
    
				// Uncomment for debug
    // collide->Debug(tnow);

    unsigned medianNumb = collide->medianNumber();
    collide->collQuantile(quant, coll_);

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
	      << setw(18) << "Collide/size"
	      << endl;
	  for (unsigned nq=0; nq<quant.size(); nq++)
	    out << setw(18) << quant[nq] << ": " 
		<< setw(18) << mfp_[nq] 
		<< setw(18) << ts_[nq] 
		<< setw(18) << coll_[nq] << endl;
	}
      }
    }

				// Begin frontier iteration to add up KE
    // pHOT_iterator ic(*c0->Tree());

    double KEtot1=collide->Etotal(), KEtot=0.0;
    double Elost1, Elost2, ElostC=0.0, ElostE=0.0;

    collide->Elost(&Elost1, &Elost2);

    MPI_Reduce(&KEtot1, &KEtot,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Elost1, &ElostC, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Elost2, &ElostE, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    unsigned cellBods = c0->Tree()->checkNumber();

    if (myid==0) {

      unsigned coll_total = collide->total()/nsteps;
      unsigned coll_error = collide->errors()/nsteps;
      unsigned epsm_total = collide->EPSMtotal()/nsteps;
      unsigned epsm_cells = collide->EPSMcells()/nsteps;

      ostringstream sout;
      sout << runtag << ".DSMC_log";
      ofstream mout(sout.str().c_str(), ios::app);

      unsigned prec = mout.precision(2);
      mout << "Summary:" << endl << left << "--------" << endl
	   << setw(6) << " " << setw(20) << tnow << "current time" << endl
	   << setw(6) << " " << setw(20) << stepnum << "step number" << endl
	   << setw(6) << " " << setw(20) << coll_total << "collisions" << endl
	   << setw(6) << " " << setw(20) << coll_error << "collision errors (" 
	   << 100.0*coll_error/(1.0e-08+coll_total) << "%)" << endl;
      if (epsm>0) mout << setw(6) << " " << setw(20) << epsm_total << "EPSM particles (" 
		       << 100.0*epsm_total/c0->nbodies_tot << "%)" << endl;
      mout << setw(6) << " " << setw(20) << medianNumb << "number/cell" << endl
	   << setw(6) << " " << setw(20) << c0->Tree()->TotalNumber() << "total # cells" << endl;
      if (epsm>0) mout << setw(6) << " " << setw(20) << epsm_cells << "EPSM cells (" 
		       << 100.0*epsm_cells/c0->Tree()->TotalNumber() << "%)" << endl;
      mout << setw(6) << " " << setw(20) << coll_[0] << "collision/cell @ 1%" << endl
	   << setw(6) << " " << setw(20) << coll_[1] << "collision/cell @ 5%" << endl
	   << setw(6) << " " << setw(20) << coll_[4] << "collision/cell @ 50%" << endl
	   << setw(6) << " " << setw(20) << coll_[7] << "collision/cell @ 95%" << endl
	   << setw(6) << " " << setw(20) << coll_[8] << "collision/cell @ 99%" << endl
	   << setw(6) << " " << setw(20) << c0->Tree()->CellCount(0.05) 
	   << "occupation @ 5%" << endl
	   << setw(6) << " " << setw(20) << c0->Tree()->CellCount(0.50) 
	   << "occupation @ 50%" << endl
	   << setw(6) << " " << setw(20) << c0->Tree()->CellCount(0.95) 
	   << "occupation @ 95%" << endl
	   << setw(6) << " " << setw(20) << cellBods
	   << "total number in cells" << endl
	   << endl;
      mout.precision(prec);
	
      ElostTotCollide += ElostC;
      ElostTotEPSM    += ElostE;

      mout << "Energy (system):" << endl << " Lost collide =" << ElostC << endl;
      if (epsm>0) mout << "    Lost EPSM =" << ElostE << endl;
      mout << "   Total loss =" << ElostTotCollide+ElostTotEPSM << endl;
      if (epsm>0) mout << "   Total EPSM =" << ElostTotEPSM << endl;
      mout << "     Total KE =" << KEtot << endl;

      if (KEtot<=0.0)
	mout << "      Ratio= XXXX" << endl << endl;
      else
	mout << "      Ratio=" << (ElostC+ElostE)/KEtot << endl << endl;

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

				// Debugging disk
				//
  // triggered_cell_body_dump(0.01, 0.002);

  MPI_Barrier(MPI_COMM_WORLD);
}

void UserTreeDSMC::triggered_cell_body_dump(double time, double radius)
{
  static bool done = false;
  if (tnow<time) return;
  if (done) return;

  unsigned cnt=0;
  vector<double> p(3);
  pHOT_iterator c(*c0->Tree());

  for (unsigned n=0; n<c0->Tree()->Number(); n++) {
    unsigned tnum = c.nextCell();
    double r2 = 0.0;

    c.Cell()->MeanPos(p);
    for (unsigned k=0; k<3; k++) r2 += p[k]*p[k];

    if (r2 < radius*radius) {
      ostringstream ostr;
      ostr << runtag << ".testcell." << myid << "." << cnt++;
      ofstream out(ostr.str().c_str());

      for (unsigned j=0; j<tnum; j++) {
	for (unsigned k=0; k<3; k++) 
	  out << setw(18) << c0->Tree()->Body(c.Cell()->bods[j])->pos[k];
	for (unsigned k=0; k<3; k++) 
	  out << setw(18) << c0->Tree()->Body(c.Cell()->bods[j])->vel[k];
	out << endl;
      }
    }
  }

  done = true;
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
