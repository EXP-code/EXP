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
  ncell = 7;
  Ncell = 64;
  cnum = 0;
  epsm = -1.0;
  diamfac = 1.0;
  boxsize = 1.0;
  boxratio = 1.0;
  jitter = 0.0;
  comp_name = "gas disk";
  nsteps = -1;
  use_temp = -1;
  use_dens = -1;
  use_delt = -1;
  use_exes = -1;
  coolfrac = 0.1;
  frontier = false;
  mfpstat = false;
  dryrun = false;
  nocool = false;
  use_multi = false;
  use_pullin = false;
  ntc = true;
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

  if (!found) {
    cerr << "UserTreeDSMC: process " << myid 
	 << ": can't find fiducial component <" << comp_name << ">" << endl;
    MPI_Abort(MPI_COMM_WORLD, 35);
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
  pHOT::sides[2] *= boxratio;
  pHOT::offst[2] *= boxratio;

				// Jitter factor of the origin offset
				// 
  for (unsigned k=0; k<3; k++) pHOT::jittr[k] = jitter*pHOT::offst[k];

  pCell::bucket = ncell;
  pCell::Bucket = Ncell;

  volume = pHOT::sides[0] * pHOT::sides[1] * pHOT::sides[2];


  //
  // Sanity check on excess attribute if excess calculation is
  // desired
  //
  if (use_exes>=0) {

    int ok1 = 1, ok;

    map<unsigned long, Particle>::iterator p = c0->Particles().begin();
    map<unsigned long, Particle>::iterator pend = c0->Particles().end();
    for (; p!=pend; p++) {
      if (use_exes >= static_cast<int>(p->second.dattrib.size())) {
	ok1 = 0;
	break;
      }
    }

    MPI_Allreduce(&ok1, &ok, 1, MPI_INT, MPI_PROD, MPI_COMM_WORLD);

    // Turn off excess calculation
    // if particles have incompatible attributes
    if (ok==0) {
      if (myid==0) {
	cout << "UserTreeDSMC: excess calculation requested but some" << endl
	     << "particles have incompatible float attribute counts." << endl
	     << "Attribute #" << use_exes << ". Continuing without excess."
	     << endl;
      }
      use_exes = -1;
    }
  }

  //
  // Set collision parameters
  //
  Collide::NTC = ntc;
  Collide::CBA = cba;
  Collide::PULLIN = use_pullin;
  Collide::CNUM = cnum;
  Collide::EPSMratio = epsm;
  Collide::DRYRUN = dryrun;
  Collide::NOCOOL = nocool;
  Collide::MFPDIAG = mfpstat;
				// Create the collision instance
  collide = new CollideLTE(diam, nthrds);
  collide->set_temp_dens(use_temp, use_dens);
  collide->set_timestep(use_delt);
  collide->set_excess(use_exes);
  ElostTotCollide = ElostTotEPSM = 0.0;

  //
  // Timers: set precision to microseconds
  //
  
  partnTime.Microseconds();
  tree1Time.Microseconds();
  tree2Time.Microseconds();
  tstepTime.Microseconds();
  llistTime.Microseconds();
  collideTime.Microseconds();

  //
  // Quantiles for distribution diagnstic
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
       << ", epsm=" << epsm << ", boxsize=" << boxsize 
       << ", ncell=" << ncell << ", Ncell=" << Ncell
       << ", boxratio=" << boxratio << ", jitter=" << jitter 
       << ", compname=" << comp_name;
  if (nsteps>0) cout << ", with diagnostic output";
  if (use_temp>=0) cout << ", temp at pos=" << use_temp;
  if (use_dens>=0) cout << ", dens at pos=" << use_dens;
  if (use_exes>=0) cout << ", excess at pos=" << use_exes;
  if (use_pullin)  cout << ", Pullin algorithm enabled";
  if (dryrun)      cout << ", collisions disabled";
  if (nocool)      cout << ", cooling disabled";
  if (ntc)         cout << ", using NTC";
  else             cout << ", NTC disabled";
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
  if (get_value("jitter", val))		jitter = atof(val.c_str());
  if (get_value("coolfrac", val))	coolfrac = atof(val.c_str());
  if (get_value("nsteps", val))		nsteps = atoi(val.c_str());
  if (get_value("ncell", val))		ncell = atoi(val.c_str());
  if (get_value("Ncell", val))		Ncell = atoi(val.c_str());
  if (get_value("compname", val))	comp_name = val;
  if (get_value("use_temp", val))	use_temp = atoi(val.c_str());
  if (get_value("use_dens", val))	use_dens = atoi(val.c_str());
  if (get_value("use_delt", val))	use_delt = atoi(val.c_str());
  if (get_value("use_exes", val))	use_exes = atoi(val.c_str());
  if (get_value("frontier", val))	frontier = atoi(val.c_str()) ? true : false;
  if (get_value("mfpstat", val))	mfpstat = atoi(val.c_str()) ? true : false;
  if (get_value("dryrun", val))		dryrun = atoi(val.c_str()) ? true : false;
  if (get_value("nocool", val))		nocool = atoi(val.c_str()) ? true : false;
  if (get_value("use_multi", val))	use_multi = atoi(val.c_str()) ? true : false;
  if (get_value("use_pullin", val))	use_pullin = atoi(val.c_str()) ? true : false;
  if (get_value("cba", val))		cba = atoi(val.c_str()) ? true : false;
  if (get_value("ntc", val))		ntc = atoi(val.c_str()) ? true : false;
}


void UserTreeDSMC::determine_acceleration_and_potential(void)
{
  static bool firstime = true;
  static unsigned nrep = 0;

  //
  // Only compute DSMC when passed the fiducial component
  //

  if (cC != c0) return;

  //
  // Make the cells
  //

  if (firstime) {
    c0->Tree()->Repartition(); nrep++;
    c0->Tree()->makeTree();
    c0->Tree()->makeCellLevelList();
#ifdef DEBUG
    c0->Tree()->checkBounds(2.0, "AFTER makeTree (first time)");
#endif

#ifdef RECTIFICATION
    c0->Tree()->Rectify();	// This is only a test!
#endif

    stepnum = 0;
    curtime = tnow;
    llevel = 0;
    clevel = 0;

    firstime = false;

#ifdef DEBUG
    cout << "Computed partition and tree [firstime]" << endl;
#endif
  } else {

    if (tnow-curtime < 1.0e-14) {
      if (myid==0) {
	cout << "UserTreeDSMC: attempt to redo step at T=" << tnow << endl;
      }
      return; 			// Don't do this time step again!
    }

    stepnum++;
    curtime = tnow;
    clevel = llevel;
    llevel = mlevel;
  }

#ifdef DEBUG
  c0->Tree()->densCheck();
#endif
  
#ifdef DEBUG
  if (c0->Tree()->checkParticles()) {
    cout << "After init only: Particle check ok [" << clevel << "]" << endl;
  } else {
    cout << "After init only: Particle check FAILED [" << clevel << "]" << endl;
  }
#endif

  //
  // Only run diagnostis every nsteps
  //
  bool diagstep = (nsteps>0 && stepnum%nsteps == 0);

  //
  // Compute time step
  //
  // double tau = dtime*mintvl[mlevel]/Mstep;

  // Now, DSMC is computed on the smallest step, every step
  double tau = dtime*mintvl[multistep]/Mstep;


  MPI_Bcast(&tau, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  TimeElapsed partnSoFar, tree1SoFar, tree2SoFar, tstepSoFar, collideSoFar;

  //
  // Test repartition
  //
  partnTime.start();
  if (mlevel==0) {
    if (0) {
      ostringstream sout;
      sout << "before Repartition [" << nrep << "], " 
	   << __FILE__ << ": " << __LINE__;
      c0->Tree()->checkBounds(2.0, sout.str().c_str());
    }
    
#ifdef RECTIFICATION
    c0->Tree()->Rectify();	// This is only a test!
#endif

    c0->Tree()->Repartition(); nrep++;
    if (0) {
      ostringstream sout;
      sout << "after Repartition [" << nrep << "], " 
	   << __FILE__ << ": " << __LINE__;
      c0->Tree()->checkBounds(2.0, sout.str().c_str());
    }

#ifdef DEBUG
    cout << "Computed partition and tree [" << mlevel << "]" << endl;
#endif
  }
  partnSoFar = partnTime.stop();
  
  //
  // Sort the particles into cells
  //
  if (mlevel==0) {
    tree1Time.start();
    c0->Tree()->makeTree();
    c0->Tree()->makeCellLevelList();
#ifdef DEBUG
    cout << "Made tree and level list [" << mlevel << "]" << endl;
    if (c0->Tree()->checkParticles()) {
      cout << "Particle check on new tree ok [" << mlevel << "]" << endl;
    } else {
      cout << "Particle check on new tree FAILED [" << mlevel << "]" << endl;
    }
#endif
    tree1SoFar = tree1Time.stop();
  } else {
    tree2Time.start();
#ifdef DEBUG
    cout << "About to adjust tree [" << clevel << "]" << endl;
#endif
    c0->Tree()->adjustTree(clevel);
#ifdef DEBUG
    cout << "About to adjust level list [" << clevel << "]" << endl;
#endif
    c0->Tree()->adjustCellLevelList(clevel);
#ifdef DEBUG
    cout << "Adjusted tree and level list [" << clevel << "]" << endl;
#endif
    tree2SoFar = tree2Time.stop();
  }

  if (0) {
    cout << "Process " << myid << ": sanity check" << endl;;
    int cnt=0;
    for (unsigned M=0; M<=multistep; M++) {
      if (c0->Tree()->clevels[M].size()) {
	for (set<pCell*>::iterator it=c0->Tree()->clevels[M].begin();
	     it!=c0->Tree()->clevels[M].end(); it++) {
	  cout << hex << (*it) << "/" << &((*it)->bods) << dec 
	       << ": M=" << M << ", bods=" << (*it)->bods.size() << endl;
	  if (cnt++>10) break;
	}
      }
      if (cnt>10) break;
    }
  }

  if (0) {
    ostringstream sout;
    sout << "after makeTree [" << nrep << "], " 
	 << __FILE__ << ": " << __LINE__;
    c0->Tree()->checkBounds(2.0, sout.str().c_str());
  }

  //
  // Evaluate collisions among the particles
  //

  collideTime.start();

  if (0) {
    ostringstream sout;
    sout << "before Collide [" << nrep << "], " 
	 << __FILE__ << ": " << __LINE__;
    c0->Tree()->checkBounds(2.0, sout.str().c_str());
  }

  // unsigned col = 
  collide->collide(*c0->Tree(), collfrac, tau, mlevel, diagstep);
    
  if (0) {
    ostringstream sout;
    sout << "after Collide [" << nrep << "], " 
	 << __FILE__ << ": " << __LINE__;
    c0->Tree()->checkBounds(2.0, sout.str().c_str());
  }

  collideSoFar = collideTime.stop();


				// Time step request
  tstepTime.start();
  if (use_multi) collide->compute_timestep(c0->Tree(), coolfrac);
  tstepSoFar = tstepTime.stop();
  
  //
  // Periodically display the current progress
  //
  //
  if (diagstep) {
				// Uncomment for debug
    // collide->Debug(tnow);

    unsigned medianNumb = collide->medianNumber();
    collide->collQuantile(quant, coll_);
    collide->mfpsizeQuantile(quant, mfp_, ts_, nsel_, rate_);
      
    if (frontier) {
      ostringstream sout;
      sout << runtag << ".DSMC.frontier." << stepnum;
      string filen = sout.str();
      c0->Tree()->testFrontier(filen);
    }

    if (mfpstat && myid==0) {

      ostringstream sout;
      sout << runtag << ".DSMC.mfpstat." << stepnum;
      string filen = sout.str();
      ofstream out(filen.c_str());
      if (out) {
	  
	out << setw(18) << "Quantiles: " 
	    << setw(18) << "MFP/size"
	    << setw(18) << "Flight/size"
	    << setw(18) << "Collions/cell"
	    << setw(18) << "Number/Nsel"
	    << setw(18) << "Energy ratio"
	    << endl;
	for (unsigned nq=0; nq<quant.size(); nq++)
	  out << setw(18) << quant[nq] << ": " 
	      << setw(18) << mfp_[nq] 
	      << setw(18) << ts_[nq] 
	      << setw(18) << coll_[nq] 
	      << setw(18) << nsel_[nq] 
	      << setw(18) << rate_[nq] 
	      << endl;
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
    unsigned oobBods  = c0->Tree()->oobNumber();

    if (myid==0) {

      unsigned sell_total = collide->select();
      unsigned coll_total = collide->total();
      unsigned coll_error = collide->errors();
      unsigned epsm_total = collide->EPSMtotal();
      unsigned epsm_cells = collide->EPSMcells();

      vector<double> disp;
      collide->dispersion(disp);
      double dmean = (disp[0]+disp[1]+disp[2])/3.0;

      ostringstream sout;
      sout << runtag << ".DSMC_log";
      ofstream mout(sout.str().c_str(), ios::app);

      mout << "Summary:" << endl << left << "--------" << endl << scientific
	   << setw(6) << " " << setw(20) << tnow << "current time" << endl
	   << setw(6) << " " << setw(20) << stepnum << "step number" << endl
	   << setw(6) << " " << setw(20) << sell_total << "targets" << endl
	   << setw(6) << " " << setw(20) << coll_total << "collisions" << endl
	   << setw(6) << " " << setw(20) << coll_error << "collision errors (" 
	   << setprecision(2) << fixed 
	   << 100.0*coll_error/(1.0e-08+coll_total) << "%)" << endl
	   << setw(6) << " " << setw(20) << oobBods << "out-of-bounds" << endl
	   << endl;

      collide->colldeTime(mout);


      if (epsm>0) mout << setw(6) << " " << setw(20) << epsm_total 
		       << "EPSM particles ("
		       << 100.0*epsm_total/c0->nbodies_tot << "%)" 
		       << scientific << endl;
      mout << setw(6) << " " << setw(20) << medianNumb << "number/cell" << endl
	   << setw(6) << " " << setw(20) << c0->Tree()->TotalNumber() 
	   << "total # cells" << endl;

      if (epsm>0) mout << setw(6) << " " << setw(20) << epsm_cells 
		       << "EPSM cells (" << setprecision(2) << fixed 
		       << 100.0*epsm_cells/c0->Tree()->TotalNumber() 
		       << "%)" << scientific << endl;

      if (mfpstat)
	mout << setw(6) << " " << setw(20) << nsel_[1] << "body/collision @ 5%" << endl
	     << setw(6) << " " << setw(20) << nsel_[4] << "body/collision @ 50%" << endl
	     << setw(6) << " " << setw(20) << nsel_[7] << "body/collision @ 95%" << endl;

      mout << fixed << setprecision(0)
	   << setw(6) << " " << setw(20) << coll_[1] << "collision/cell @ 5%" << endl
	   << setw(6) << " " << setw(20) << coll_[4] << "collision/cell @ 50%" << endl
	   << setw(6) << " " << setw(20) << coll_[7] << "collision/cell @ 95%" << endl
	   << setw(6) << " " << setw(20) << c0->Tree()->CellCount(0.05) 
	   << "occupation @ 5%" << endl
	   << setw(6) << " " << setw(20) << c0->Tree()->CellCount(0.50) 
	   << "occupation @ 50%" << endl
	   << setw(6) << " " << setw(20) << c0->Tree()->CellCount(0.95) 
	   << "occupation @ 95%" << endl
	   << setw(6) << " " << setw(20) << cellBods
	   << "total number in cells" << endl
	   << endl << setprecision(4);
	
      ElostTotCollide += ElostC;
      ElostTotEPSM    += ElostE;

      mout << scientific << "Energy (system):" << endl 
	   << " Lost collide =" << ElostC << endl;
      if (epsm>0) mout << "    Lost EPSM =" << ElostE << endl;
      mout << "   Total loss =" << ElostTotCollide+ElostTotEPSM << endl;
      if (epsm>0) mout << "   Total EPSM =" << ElostTotEPSM << endl;
      mout << "     Total KE =" << KEtot << endl;
      if (KEtot<=0.0) mout << "         Ratio= XXXX" << endl;
      else mout << "    Ratio lost=" << (ElostC+ElostE)/KEtot << endl;
      mout << "     3-D disp =" << disp[0] << ", " << disp[1] 
	   << ", " << disp[2] << endl;
      if (dmean>0.0) {
	mout << "   Disp ratio =" << disp[0]/dmean << ", " 
	     << disp[1]/dmean << ", " << disp[2]/dmean << endl << endl;
      }
	
      double keymake, xchange, convert, overlap;
      c0->Tree()->adjustTiming(keymake, xchange, convert, overlap);      

      mout << "Timing (secs) at mlevel=" << mlevel << ":" << endl
	   << "  partition=" << partnSoFar()*1.0e-6 << endl
	   << "  make tree=" << tree1SoFar()*1.0e-6 << endl
	   << "adjust tree=" << tree2SoFar()*1.0e-6 << endl
	   << "      *** keymake=" << keymake << endl
	   << "      *** xchange=" << xchange << endl
	   << "      *** convert=" << convert << endl
	   << "      *** overlap=" << overlap << endl
	   << "  timesteps=" << tstepSoFar()*1.0e-6 << endl
	   << "  step list=" << llistTime.getTime().getRealTime()*1.0e-6 
	   << endl
	   << "    collide=" << collideSoFar()*1.0e-6 << endl
	   << endl;

      collide->tsdiag(mout);

      partnTime.reset();
      tree1Time.reset();
      tree2Time.reset();
      tstepTime.reset();
      llistTime.reset();
      collideTime.reset();
    }

  }

#ifdef DEBUG
  if (c0->Tree()->checkParticles()) {
    cout << "Before level list: Particle check ok [" << clevel << "]" << endl;
  } else {
    cout << "Before level list: Particle check FAILED [" << clevel << "]" << endl;
  }
#endif

				// Remake level lists because particles
				// will (usually) have been exchanged 
				// between nodes
  llistTime.start();
  c0->reset_level_lists();
  llistTime.stop();

#ifdef DEBUG
  if (c0->Tree()->checkParticles()) {
    cout << "After level list: Particle check ok [" << clevel << "]" << endl;
  } else {
    cout << "After level list: Particle check FAILED [" << clevel << "]" << endl;
  }
#endif

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
    c.nextCell();
    double r2 = 0.0;

    c.Cell()->MeanPos(p);
    for (unsigned k=0; k<3; k++) r2 += p[k]*p[k];

    if (r2 < radius*radius) {
      ostringstream ostr;
      ostr << runtag << ".testcell." << myid << "." << cnt++;
      ofstream out(ostr.str().c_str());

      for (set<unsigned>::iterator j=c.Cell()->bods.begin();
	   j!=c.Cell()->bods.end(); j++) {
	for (unsigned k=0; k<3; k++) 
	  out << setw(18) << c0->Tree()->Body(*j)->pos[k];
	for (unsigned k=0; k<3; k++) 
	  out << setw(18) << c0->Tree()->Body(*j)->vel[k];
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
