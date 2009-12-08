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

#include <expand.h>
#include <ExternalCollection.H>
#include <UserTreeDSMC.H>
#include <CollideLTE.H>

#ifdef USE_GPTL
#include <gptl.h>
#endif

using namespace std;

//
// Physical units
//

static double pc = 3.086e18;		// cm
static double a0 = 2.0*0.054e-7;	// cm (2xBohr radius)
static double boltz = 1.381e-16;	// cgs
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
  ncell      = 7;		// 
  Ncell      = 64;
  cnum       = 0;
  madj       = 512;		// No tree pruning by default
  wght       = 1;		// Cell time partitioning by default
  epsm       = -1.0;
  diamfac    = 1.0;
  boxsize    = 1.0;
  boxratio   = 1.0;
  jitter     = 0.0;
  comp_name  = "gas disk";
  nsteps     = -1;
  msteps     = -1;
  use_temp   = -1;
  use_dens   = -1;
  use_delt   = -1;
  use_exes   = -1;
  coolfrac   = 0.1;
  frontier   = false;
  tsdiag     = false;
  tspow      = 4;
  mfpstat    = false;
  cbadiag    = false;
  dryrun     = false;
  nocool     = false;
  use_multi  = false;
  use_pullin = false;
  ntc        = true;
  cba        = true;
  tube       = false;
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

  if (tube) {
    pHOT::sides[0] = pHOT::sides[1] = pHOT::sides[2] = boxsize;
    pHOT::offst[0] = pHOT::offst[1] = pHOT::offst[2] = 0.0;
    pHOT::sides[2] *= boxratio;
  } else {

    pHOT::sides[0] = pHOT::sides[1] = pHOT::sides[2] = 2.0*boxsize;
    pHOT::offst[0] = pHOT::offst[1] = pHOT::offst[2] = boxsize;

				// For rectangular slab
    pHOT::sides[2] *= boxratio;
    pHOT::offst[2] *= boxratio;
  }

				// Jitter factor of the dimensions
				// 
  for (unsigned k=0; k<3; k++) pHOT::jittr[k] = jitter*pHOT::sides[k];

  pCell::bucket = ncell;
  pCell::Bucket = Ncell;

  volume = pHOT::sides[0] * pHOT::sides[1] * pHOT::sides[2];


  //
  // Sanity check on excess attribute if excess calculation is
  // desired
  //
  if (use_exes>=0) {

    int ok1 = 1, ok;

    PartMapItr p = c0->Particles().begin();
    PartMapItr pend = c0->Particles().end();
    for (; p!=pend; p++) {
      if (use_exes >= static_cast<int>(p->second.dattrib.size())) {
	ok1 = 0;
	break;
      }
    }

    MPI_Allreduce(&ok1, &ok, 1, MPI_INT, MPI_PROD, MPI_COMM_WORLD);

    // Turn off excess calculation
    // if particles have incompatible attributes
    //
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
  Collide::NTC     = ntc;
  Collide::CBA     = cba;
  Collide::CBADIAG = cbadiag;
  Collide::PULLIN  = use_pullin;
  Collide::CNUM    = cnum;
  Collide::EPSMratio = epsm;
  Collide::DRYRUN  = dryrun;
  Collide::NOCOOL  = nocool;
  Collide::TSDIAG  = tsdiag;
  Collide::TSPOW   = tspow;
  Collide::MFPDIAG = mfpstat;
				// Create the collision instance
  collide = new CollideLTE(this, diam, nthrds);
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
  quant.push_back(0.0);		// 0
  quant.push_back(0.01);	// 1
  quant.push_back(0.05);	// 2
  quant.push_back(0.1);		// 3
  quant.push_back(0.2);		// 4
  quant.push_back(0.5);		// 5
  quant.push_back(0.8);		// 6
  quant.push_back(0.9);		// 7
  quant.push_back(0.95);	// 8
  quant.push_back(0.99);	// 9
  quant.push_back(1.0);		// 10

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
       << ", madj=" << madj << ", epsm=" << epsm << ", boxsize=" << boxsize 
       << ", ncell=" << ncell << ", Ncell=" << Ncell << ", wght=" << wght
       << ", boxratio=" << boxratio << ", jitter=" << jitter 
       << ", compname=" << comp_name;
  if (msteps>=0) 
    cout << ", with diagnostic output at levels <= " << msteps;
  else if (nsteps>0) 
    cout << ", with diagnostic output every " << nsteps << " steps";
  if (use_temp>=0) cout << ", temp at pos=" << use_temp;
  if (use_dens>=0) cout << ", dens at pos=" << use_dens;
  if (use_exes>=0) cout << ", excess at pos=" << use_exes;
  if (use_pullin)  cout << ", Pullin algorithm enabled";
  if (dryrun)      cout << ", collisions disabled";
  if (nocool)      cout << ", cooling disabled";
  if (epsm>0.0)    cout << ", using EPSM";
  else             cout << ", EPSM disabled";
  if (ntc)         cout << ", using NTC";
  else             cout << ", NTC disabled";
  if (cba)         cout << ", using CBA";
  else             cout << ", CBA disabled";
  if (cba && cbadiag)     
                   cout << " with diagnostics";
  if (tube)         cout << ", using TUBE mode";
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
  if (get_value("madj", val))		madj = atoi(val.c_str());
  if (get_value("wght", val))		wght = atoi(val.c_str());
  if (get_value("epsm", val))		epsm = atof(val.c_str());
  if (get_value("diamfac", val))	diamfac = atof(val.c_str());
  if (get_value("boxsize", val))	boxsize = atof(val.c_str());
  if (get_value("boxratio", val))	boxratio = atof(val.c_str());
  if (get_value("jitter", val))		jitter = atof(val.c_str());
  if (get_value("coolfrac", val))	coolfrac = atof(val.c_str());
  if (get_value("nsteps", val))		nsteps = atoi(val.c_str());
  if (get_value("msteps", val))		msteps = atoi(val.c_str());
  if (get_value("ncell", val))		ncell = atoi(val.c_str());
  if (get_value("Ncell", val))		Ncell = atoi(val.c_str());
  if (get_value("compname", val))	comp_name = val;
  if (get_value("use_temp", val))	use_temp = atoi(val.c_str());
  if (get_value("use_dens", val))	use_dens = atoi(val.c_str());
  if (get_value("use_delt", val))	use_delt = atoi(val.c_str());
  if (get_value("use_exes", val))	use_exes = atoi(val.c_str());
  if (get_value("frontier", val))	frontier = atoi(val.c_str()) ? true : false;
  if (get_value("tspow", val))		tspow = atoi(val.c_str());
  if (get_value("tsdiag", val))		tsdiag = atoi(val.c_str()) ? true : false;
  if (get_value("mfpstat", val))	mfpstat = atoi(val.c_str()) ? true : false;
  if (get_value("cbadiag", val))	cbadiag = atoi(val.c_str()) ? true : false;
  if (get_value("dryrun", val))		dryrun = atoi(val.c_str()) ? true : false;
  if (get_value("nocool", val))		nocool = atoi(val.c_str()) ? true : false;
  if (get_value("use_multi", val))	use_multi = atoi(val.c_str()) ? true : false;
  if (get_value("use_pullin", val))	use_pullin = atoi(val.c_str()) ? true : false;
  if (get_value("cba", val))		cba = atoi(val.c_str()) ? true : false;
  if (get_value("ntc", val))		ntc = atoi(val.c_str()) ? true : false;
  if (get_value("tube", val))		tube = atoi(val.c_str()) ? true : false;
}


void UserTreeDSMC::determine_acceleration_and_potential(void)
{
#ifdef USE_GPTL
  GPTLstart("UserTreeDSMC::determine_acceleration_and_potential");
#endif

  static bool firstime = true;
  static unsigned nrep = 0;

  //
  // Only compute DSMC when passed the fiducial component
  //

  if (cC != c0) {
#ifdef USE_GPTL
    GPTLstop("UserTreeDSMC::determine_acceleration_and_potential");
#endif
    return;
  }

  //
  // Make the cells
  //

  if (firstime) {
    c0->Tree()->setWeights(wght ? true : false);
    c0->Tree()->Repartition(0); nrep++;
    c0->Tree()->makeTree();
    c0->Tree()->makeCellLevelList();
#ifdef DEBUG
    c0->Tree()->checkBounds(2.0, "AFTER makeTree (first time)");
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
#ifdef USE_GPTL
      GPTLstop("UserTreeDSMC::determine_acceleration_and_potential");
#endif
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
  // Only run diagnostics every nsteps
  //
  bool diagstep = (nsteps>0 && stepnum%nsteps == 0);

  //
  // Diagnostics run at levels <= msteps (takes prececdence over nsteps)
  //
  if (msteps>=0) 
    diagstep = (mlevel <= static_cast<unsigned>(msteps)) ? true : false;

  //
  // Compute time step
  //
  // Now, DSMC is computed on the smallest step, every step
  double tau = dtime*mintvl[multistep]/Mstep;


  MPI_Bcast(&tau, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  TimeElapsed partnSoFar, tree1SoFar, tree2SoFar, tstepSoFar, collideSoFar;

  //
  // Sort the particles into cells
  //
  if (mlevel<=madj) {

#ifdef USE_GPTL
    GPTLstart("UserTreeDSMC::pHOT_1");
    GPTLstart("UserTreeDSMC::waiting");
    MPI_Barrier(MPI_COMM_WORLD);
    GPTLstop("UserTreeDSMC::waiting");
    GPTLstart("UserTreeDSMC::repart");
#endif

    partnTime.start();
    c0->Tree()->Repartition(mlevel); nrep++;
    partnSoFar = partnTime.stop();

    tree1Time.start();
#ifdef USE_GPTL
    GPTLstop("UserTreeDSMC::repart");
    GPTLstart("UserTreeDSMC::makeTree");
#endif
    c0->Tree()->makeTree();
#ifdef USE_GPTL
    GPTLstop("UserTreeDSMC::makeTree");
    GPTLstart("UserTreeDSMC::makeCLL");
#endif
    c0->Tree()->makeCellLevelList();
#ifdef USE_GPTL
    GPTLstop("UserTreeDSMC::makeCLL");
    GPTLstart("UserTreeDSMC::pcheck");
#endif
#ifdef DEBUG
    cout << "Made partition, tree and level list [" << mlevel << "]" << endl;
    if (c0->Tree()->checkParticles()) {
      cout << "Particle check on new tree ok [" << mlevel << "]" << endl;
    } else {
      cout << "Particle check on new tree FAILED [" << mlevel << "]" << endl;
    }
#endif
    tree1SoFar = tree1Time.stop();
    
#ifdef USE_GPTL
    GPTLstop("UserTreeDSMC::pcheck");
    GPTLstop("UserTreeDSMC::pHOT_1");
#endif

  } else {

#ifdef USE_GPTL
    GPTLstart("UserTreeDSMC::pHOT_2");
    GPTLstart("UserTreeDSMC::adjustTree");
#endif

    tree2Time.start();
#ifdef DEBUG
    cout << "About to adjust tree [" << clevel << "]" << endl;
#endif
    c0->Tree()->adjustTree(clevel);
#ifdef DEBUG
    cout << "About to adjust level list [" << clevel << "]" << endl;
#endif
#ifdef USE_GPTL
    GPTLstop("UserTreeDSMC::adjustTree");
    GPTLstart("UserTreeDSMC::adjustCLL");
#endif
    c0->Tree()->adjustCellLevelList(clevel);
#ifdef DEBUG
    cout << "Adjusted tree and level list [" << clevel << "]" << endl;
#endif
    tree2SoFar = tree2Time.stop();
    
#ifdef USE_GPTL
    GPTLstop("UserTreeDSMC::adjustCLL");
    GPTLstop("UserTreeDSMC::pHOT_2");
#endif

  }

  if (0) {
    cout << "Process " << myid << ": sanity check" << endl;;
    for (unsigned M=0; M<=multistep; M++) {
      int cnt=0;
      if (c0->Tree()->clevels[M].size()) {
	for (set<pCell*>::iterator it=c0->Tree()->clevels[M].begin();
	     it!=c0->Tree()->clevels[M].end(); it++) {
	  cout << hex << (*it) << "/" << &((*it)->bods) << dec 
	       << ": M=" << M << ", bods=" << (*it)->bods.size() << endl;
	  if (cnt++>10) break;
	}
	if (cnt>10) break;
      }
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

#ifdef USE_GPTL
  GPTLstart("UserTreeDSMC::collide");
#endif

  collide->collide(*c0->Tree(), collfrac, tau, mlevel, diagstep);
    
#ifdef USE_GPTL
  GPTLstop("UserTreeDSMC::collide");
#endif

  if (0) {
    ostringstream sout;
    sout << "after Collide [" << nrep << "], " 
	 << __FILE__ << ": " << __LINE__;
    c0->Tree()->checkBounds(2.0, sout.str().c_str());
  }

  collideSoFar = collideTime.stop();


				// Time step request
#ifdef USE_GPTL
  GPTLstart("UserTreeDSMC::collide_timestep");
#endif


  tstepTime.start();
  if (use_multi) collide->compute_timestep(c0->Tree(), coolfrac);
  tstepSoFar = tstepTime.stop();
  
#ifdef USE_GPTL
  GPTLstop("UserTreeDSMC::collide_timestep");
#endif

  //
  // Periodically display the current progress
  //
  //
  if (diagstep) {
#ifdef USE_GPTL
    GPTLstart("UserTreeDSMC::collide_diag");
#endif

				// Uncomment for debug
    collide->Debug(tnow);

    unsigned medianNumb = collide->medianNumber();
    unsigned collnum=0, coolnum=0;
    if (mfpstat) {
      collide->collQuantile(quant, coll_);
      collide->mfpsizeQuantile(quant, mfp_, ts_, nsel_, cool_, rate_,
			       collnum, coolnum);
    }

    double ExesCOLL, ExesEPSM;
    if (use_exes>=0) collide->energyExcess(ExesCOLL, ExesEPSM);
      
    if (frontier) {
      ostringstream sout;
      sout << outdir << runtag << ".DSMC_frontier";
      string filen = sout.str();
      c0->Tree()->testFrontier(filen);
    }

    vector<unsigned> ncells, bodies;
    c0->Tree()->countFrontier(ncells, bodies);

    if (mfpstat && myid==0) {

				// Generate the file name
      ostringstream sout;
      sout << outdir << runtag << ".DSMC_mfpstat";
      string filen = sout.str();

				// Check for existence
      ifstream in(filen.c_str());
      if (in.fail()) {
				// Write a new file
	ofstream out(filen.c_str());
	if (out) {
	  out << left
	      << setw(14) << "# Time"
	      << setw(14) << "Quantiles" 
	      << setw(14) << "Bodies"
	      << setw(14) << "MFP/size"
	      << setw(14) << "Flight/size"
	      << setw(14) << "Collions/cell"
	      << setw(14) << "Nsel/Number"
	      << setw(14) << "Energy ratio"
	      << setw(14) << "Excess ratio"
	      << endl;
	}
      }
      in.close();

				// Open old file to write a stanza
				// 
      ofstream out(filen.c_str(), ios::app);
      if (out) {
	for (unsigned nq=0; nq<quant.size(); nq++) {
	  out << setw(14) << tnow
	      << setw(14) << quant[nq]
	      << setw(14) << collnum
	      << setw(14) << mfp_[nq] 
	      << setw(14) << ts_[nq] 
	      << setw(14) << coll_[nq] 
	      << setw(14) << nsel_[nq] 
	      << setw(14) << cool_[nq] 
	      << setw(14) << rate_[nq] 
	      << endl;
	}
	out << endl;
      }
    }

				// Overall statistics
				// 
    double KEtotl=collide->Etotal(), KEtot=0.0;
    double Mtotal=collide->Mtotal(), Mtotl=0.0;
    double Elost1, Elost2, ElostC=0.0, ElostE=0.0;

    collide->Elost(&Elost1, &Elost2);

    MPI_Reduce(&KEtotl, &KEtot,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Mtotal, &Mtotl,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Elost1, &ElostC, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Elost2, &ElostE, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

				// Computing mass-weighted temperature
    const double f_H = 0.76;
    double mm = f_H*mp + (1.0-f_H)*4.0*mp;
    double meanT = 0.0;
    if (Mtotl>0.0) meanT = 2.0*KEtot/Mtotl*Eunit/3.0 * mm/Munit/boltz;

    unsigned cellBods = c0->Tree()->checkNumber();
    unsigned oobBods  = c0->Tree()->oobNumber();

    double Mass;
    unsigned Counts;
    c0->Tree()->totalMass(Counts, Mass);

    double zerorate = c0->Tree()->checkAdjust();

				// Check frontier for mass at or below 
				// current level
    double cmass1=0.0, cmass=0.0;
    pHOT_iterator pit(*(c0->Tree()));

    while (pit.nextCell()) {
      pCell *cc = pit.Cell();
      if (cc->maxplev >= mlevel && cc->count > 1) cmass1 += cc->state[0];
    }

    MPI_Reduce(&cmass1, &cmass, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
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
      sout << outdir << runtag << ".DSMC_log";
      ofstream mout(sout.str().c_str(), ios::app);

      mout << "Summary:" << endl << left << "--------" << endl << scientific
	   << setw(6) << " " << setw(20) << tnow       << "current time" << endl
	   << setw(6) << " " << setw(20) << mlevel     << "current level" << endl
	   << setw(6) << " " << setw(20) << Counts     << "total counts" << endl
	   << setw(6) << " " << setw(20) << Mass       << "total mass" << endl
	   << setw(6) << " " << setw(20) << meanT      << "mass-weighted temperature" << endl
	   << setw(6) << " " << setw(20) << Mtotl      << "accumulated mass" << endl
	   << setw(6) << " " << setw(20) << cmass      << "mass at this level" << endl
	   << setw(6) << " " << setw(20) << zerorate   << "zero partition" << endl
	   << setw(6) << " " << setw(20) << stepnum    << "step number" << endl
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

      if (mfpstat) {
	mout << setw(6) << " " << setw(20) << "--------" << "--------------------" << endl
	     << setw(6) << " " << setw(20) << nsel_[0] << "collision/body @  0%" << endl
	     << setw(6) << " " << setw(20) << nsel_[2] << "collision/body @  5%" << endl
	     << setw(6) << " " << setw(20) << nsel_[5] << "collision/body @ 50%" << endl
	     << setw(6) << " " << setw(20) << nsel_[8] << "collision/body @ 95%" << endl
	     << setw(6) << " " << setw(20) << nsel_[10] << "collision/body @100%" << endl;

	mout << fixed << setprecision(0)
	     << setw(6) << " " << setw(20) << "--------" << "--------------------" << endl

	     << setw(6) << " " << setw(20) << coll_[0] << "collision/cell @  0%" << endl
	     << setw(6) << " " << setw(20) << coll_[2] << "collision/cell @  5%" << endl
	     << setw(6) << " " << setw(20) << coll_[5] << "collision/cell @ 50%" << endl
	     << setw(6) << " " << setw(20) << coll_[8] << "collision/cell @ 95%" << endl
	     << setw(6) << " " << setw(20) << coll_[10] << "collision/cell @100%" << endl
	     << setw(6) << " " << setw(20) << "--------" << "--------------------" << endl

	     << setw(6) << " " << setw(20) << c0->Tree()->CellCount(0.0) 
	     << "occupation @  0%" << endl
	     << setw(6) << " " << setw(20) << c0->Tree()->CellCount(0.05) 
	     << "occupation @  5%" << endl
	     << setw(6) << " " << setw(20) << c0->Tree()->CellCount(0.50) 
	     << "occupation @ 50%" << endl
	     << setw(6) << " " << setw(20) << c0->Tree()->CellCount(0.95) 
	     << "occupation @ 95%" << endl
	     << setw(6) << " " << setw(20) << c0->Tree()->CellCount(1.0) 
	     << "occupation @100%" << endl
	     << setw(6) << " " << setw(20) << "--------" << "--------------------" << endl
	     << setw(6) << " " << setw(20) << cellBods
	     << "total number in cells" << endl
	     << endl << setprecision(4);
      }
	
      ElostTotCollide += ElostC;
      ElostTotEPSM    += ElostE;

      mout << scientific << "Energy (system):" << endl 
	   << " Lost collide =" << ElostC << endl;
      if (epsm>0) mout << "    Lost EPSM =" << ElostE << endl;
      mout << "   Total loss =" << ElostTotCollide+ElostTotEPSM << endl;
      if (epsm>0) mout << "   Total EPSM =" << ElostTotEPSM << endl;
      mout << "     Total KE =" << KEtot << endl;
      if (use_exes>=0) {
	mout << "  COLL excess =" << ExesCOLL << endl;
	if (epsm>0) mout << "  EPSM excess =" << ExesEPSM << endl;
      }
      if (KEtot<=0.0) mout << "         Ratio= XXXX" << endl;
      else mout << "    Ratio lost=" << (ElostC+ElostE)/KEtot << endl;
      mout << "     3-D disp =" << disp[0] << ", " << disp[1] 
	   << ", " << disp[2] << endl;
      if (dmean>0.0) {
	mout << "   Disp ratio =" << disp[0]/dmean << ", " 
	     << disp[1]/dmean << ", " << disp[2]/dmean << endl << endl;
      }
	
      unsigned sumcells=0, sumbodies=0;
      mout << endl;
      mout << "-----------------------------------------------------" << endl;
      mout << "-----Cell/body diagnostics---------------------------" << endl;
      mout << "-----------------------------------------------------" << endl;
      mout << right << setw(8) << "Level" << setw(15) << "Scale(x)"
	   << setw(10) << "Cells" << setw(10) << "Bodies" << endl;
      mout << "-----------------------------------------------------" << endl;
      for (unsigned n=0; n<ncells.size(); n++) {
	mout << setw(8) << n << setw(15) << pHOT::sides[0]/(1<<n)
	     << setw(10) << ncells[n] << setw(10) << bodies[n]
	     << endl;
	sumcells  += ncells[n];
	sumbodies += bodies[n];
      }
      mout << "-----------------------------------------------------" << endl;
      mout << setw(8) << "TOTALS" 
	   << setw(15) << "**********"
	   << setw(10) << sumcells << setw(10) << sumbodies << endl;
      mout << "-----------------------------------------------------" << endl;
      mout << left << endl;
      
      double keymake, xchange, prepare, convert, overlap, update, scatter;
      double repartn, tadjust, keycall, keycomp, keybods, keywait;
      unsigned numbods;
      c0->Tree()->adjustTiming(keymake, xchange, prepare, 
			       convert, overlap, update, 
			       scatter, repartn, tadjust,
			       keycall, keycomp, keybods,
			       keywait, numbods);

      mout << "Timing (secs) at mlevel=" << mlevel << ":" << endl
	   << "  partition=" << partnSoFar()*1.0e-6 << endl
	   << "  make tree=" << tree1SoFar()*1.0e-6 << endl
	   << "adjust tree=" << tree2SoFar()*1.0e-6 << endl
	   << "      *** keymake=" << keymake << endl
	   << "      *** keycall=" << keycall << endl
	   << "      *** keycomp=" << keycomp << endl
	   << "      *** keybods=" << keybods << endl
	   << "      *** keywait=" << keywait << endl
	   << "      *** xchange=" << xchange << endl
	   << "      *** prepare=" << prepare << endl
	   << "      *** convert=" << convert << endl
	   << "      *** overlap=" << overlap << endl
	   << "      *** cupdate=" << update  << endl
	   << "      *** scatter=" << scatter << endl
	   << "      *** repartn=" << repartn << endl
	   << "      *** tadjust=" << tadjust << endl
	   << "      *** numbods=" << numbods << endl
	   << "  timesteps=" << tstepSoFar()*1.0e-6 << endl
	   << "  step list=" << llistTime.getTime().getRealTime()*1.0e-6 
	   << endl
	   << "    collide=" << collideSoFar()*1.0e-6 << endl
	   << endl;

      collide->tsdiag(mout);

      {
	sout.str("");
	sout << outdir << runtag << ".DSMC_log.0";
	ofstream nout(sout.str().c_str(), ios::app);

	nout << "Timing (secs) at mlevel=" << mlevel << " and T=" << tnow << endl
	     << "      *** keymake=" << keymake << endl
	     << "      *** keycall=" << keycall << endl
	     << "      *** keycomp=" << keycomp << endl
	     << "      *** keybods=" << keybods << endl
	     << "      *** keywait=" << keywait << endl
	     << "      *** xchange=" << xchange << endl
	     << "      *** prepare=" << prepare << endl
	     << "      *** convert=" << convert << endl
	     << "      *** overlap=" << overlap << endl
	     << "      *** cupdate=" << update  << endl
	     << "      *** scatter=" << scatter << endl
	     << "      *** repartn=" << repartn << endl
	     << "      *** tadjust=" << tadjust << endl
	     << "      *** numbods=" << numbods << endl
	     << endl;
      }
      
      partnTime.reset();
      tree1Time.reset();
      tree2Time.reset();
      tstepTime.reset();
      llistTime.reset();
      collideTime.reset();

    } else {
      
      ostringstream sout;
      sout << outdir << runtag << ".DSMC_log." << myid;
      ofstream mout(sout.str().c_str(), ios::app);

      double keymake, xchange, prepare, convert, overlap, update, scatter;
      double repartn, tadjust, keycall, keycomp, keybods, keywait;
      unsigned numbods;
      c0->Tree()->adjustTiming(keymake, xchange, prepare, 
			       convert, overlap, update, 
			       scatter, repartn, tadjust,
			       keycall, keycomp, keybods,
			       keywait, numbods);

      mout << "Timing (secs) at mlevel=" << mlevel << " and T=" << tnow << endl
	   << "      *** keymake=" << keymake << endl
	   << "      *** keycall=" << keycall << endl
	   << "      *** keycomp=" << keycomp << endl
	   << "      *** keybods=" << keybods << endl
	   << "      *** keywait=" << keywait << endl
	   << "      *** xchange=" << xchange << endl
	   << "      *** prepare=" << prepare << endl
	   << "      *** convert=" << convert << endl
	   << "      *** overlap=" << overlap << endl
	   << "      *** cupdate=" << update  << endl
	   << "      *** scatter=" << scatter << endl
	   << "      *** repartn=" << repartn << endl
	   << "      *** tadjust=" << tadjust << endl
	   << "      *** numbods=" << numbods << endl
	   << endl;
    }

#ifdef USE_GPTL
    GPTLstop("UserTreeDSMC::collide_diag");
#endif

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

#ifdef USE_GPTL
  GPTLstop("UserTreeDSMC::determine_acceleration_and_potential");
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
      ostr << outdir << runtag << ".testcell." << myid << "." << cnt++;
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
