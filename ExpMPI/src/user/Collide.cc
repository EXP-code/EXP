#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>

using namespace std;

#include "Timer.h"
#include "global.H"
#include "pHOT.H"
#include "UserTreeDSMC.H"
#include "Collide.H"

				// Use the original Pullin velocity 
				// selection algorithm
bool Collide::PULLIN = false;
				// Print out sorted cell parameters
bool Collide::SORTED = false;
				// Print out T-rho plane for cells 
				// with mass weighting
bool Collide::PHASE = false;
				// Extra debugging output
bool Collide::EXTRA = false;
				// Turn off collisions for testing
bool Collide::DRYRUN = false;
				// Turn off cooling for testing
bool Collide::NOCOOL = false;
				// Ensemble-based excess cooling
bool Collide::ENSEXES = true;
				// Time step diagnostics
bool Collide::TSDIAG = false;
				// CBA length scale diagnostics
bool Collide::CBADIAG = false;
				// Mean free path diagnostics
bool Collide::MFPDIAG = false;
				// Sample based on maximum (true) or estimate
bool Collide::NTC = false;	// from variance (false);

				// Temperature floor in EPSM
double Collide::TFLOOR = 1000.0;

int Collide::TSPOW = 4;		// Power of two interval for KE/cool histogram

				// Proton mass (g)
const double mp = 1.67262158e-24;
				// Boltzmann constant (cgs)
const double boltz = 1.3810e-16;

extern "C"
void *
collide_thread_call(void *atp)
{
  thrd_pass_Collide *tp = (thrd_pass_Collide *)atp;
  Collide *p = (Collide *)tp->p;
  p -> collide_thread((void*)&tp->arg);
  return NULL;
}

void Collide::collide_thread_fork(pHOT* tree, double Fn, double tau)
{
  int errcode;
  void *retval;
  
  if (nthrds==1) {
    thrd_pass_Collide td;
    
    td.p = this;
    td.arg.tree = tree;
    td.arg.fn = Fn;
    td.arg.tau = tau;
    td.arg.id = 0;

    collide_thread_call(&td);

    return;
  }

  td = new thrd_pass_Collide [nthrds];
  t = new pthread_t [nthrds];

  if (!td) {
    cerr << "Process " << myid 
         << ": collide_thread_fork: error allocating memory for thread counters\n";
    exit(18);
  }
  if (!t) {
    cerr << "Process " << myid
         << ": collide_thread_fork: error allocating memory for thread\n";
    exit(18);
  }

                                // Make the <nthrds> threads
  for (int i=0; i<nthrds; i++) {
    td[i].p = this;
    td[i].arg.tree = tree;
    td[i].arg.fn = Fn;
    td[i].arg.tau = tau;
    td[i].arg.id = i;

    errcode =  pthread_create(&t[i], 0, collide_thread_call, &td[i]);
    if (errcode) {
      cerr << "Process " << myid;
      cerr << " collide: cannot make thread " << i
	   << ", errcode=" << errcode << endl;
      exit(19);
    }
  }
    
  waitTime.start();

                                // Collapse the threads
  for (int i=0; i<nthrds; i++) {
    if ((errcode=pthread_join(t[i], &retval))) {
      cerr << "Process " << myid;
      cerr << " collide: thread join " << i
           << " failed, errcode=" << errcode << endl;
      exit(20);
    }
    if (i==0) {
      waitSoFar = waitTime.stop();
      joinTime.start();
    }
  }
  
  joinSoFar = joinTime.stop();

  delete [] td;
  delete [] t;
}


int Collide::CNUM = 0;
bool Collide::CBA = true;
double Collide::EPSMratio = -1.0;

Collide::Collide(double diameter, int nth)
{
  nthrds = nth;

  colcntT = vector< vector<unsigned> > (nthrds);
  numcntT = vector< vector<unsigned> > (nthrds);
  tdispT  = vector< vector<double> >   (nthrds);
  error1T = vector<unsigned> (nthrds, 0);
  sel1T   = vector<unsigned> (nthrds, 0);
  col1T   = vector<unsigned> (nthrds, 0);
  epsm1T  = vector<unsigned> (nthrds, 0);
  Nepsm1T = vector<unsigned> (nthrds, 0);
  tmassT  = vector<double>   (nthrds, 0);
  decolT  = vector<double>   (nthrds, 0);
  decelT  = vector<double>   (nthrds, 0);
  exesCT  = vector<double>   (nthrds, 0);
  exesET  = vector<double>   (nthrds, 0);

  if (MFPDIAG) {
    tsratT  = vector< vector<double> > (nthrds);
    keratT  = vector< vector<double> > (nthrds);
    deratT  = vector< vector<double> > (nthrds);
    tdensT  = vector< vector<double> > (nthrds);
    tvolcT  = vector< vector<double> > (nthrds);
    ttempT  = vector< vector<double> > (nthrds);
    tdeltT  = vector< vector<double> > (nthrds);
    tselnT  = vector< vector<double> > (nthrds);
    tphaseT = vector< vector<Precord> > (nthrds);
    tmfpstT = vector< vector<Precord> > (nthrds);
  }

  cellist = vector< vector<pCell*> > (nthrds);

  diam0 = diam = diameter;
  seltot = 0;			// Count estimated collision targets
  coltot = 0;			// Count total collisions
  errtot = 0;			// Count errors in inelastic computation
  epsmcells = 0;		// Count cells in EPSM regime
  epsmtot = 0;			// Count particles in EPSM regime

				// Default cooling rate 
				// (if not set by derived class)
  coolrate = vector<double>(nthrds, 0.0);

				// EPSM diagnostics
  lostSoFar_EPSM = vector<double>(nthrds, 0.0);

  diagTime.Microseconds();
  snglTime.Microseconds();
  forkTime.Microseconds();
  waitTime.Microseconds();
  joinTime.Microseconds();

  stepcount = 0;
  bodycount = 0;

  initTime = vector<Timer>(nthrds);
  collTime = vector<Timer>(nthrds);
  elasTime = vector<Timer>(nthrds);
  stat1Time = vector<Timer>(nthrds);
  stat2Time = vector<Timer>(nthrds);
  stat3Time = vector<Timer>(nthrds);
  coolTime = vector<Timer>(nthrds);
  initSoFar = vector<TimeElapsed>(nthrds);
  collSoFar = vector<TimeElapsed>(nthrds);
  elasSoFar = vector<TimeElapsed>(nthrds);
  stat1SoFar = vector<TimeElapsed>(nthrds);
  stat2SoFar = vector<TimeElapsed>(nthrds);
  stat3SoFar = vector<TimeElapsed>(nthrds);
  coolSoFar = vector<TimeElapsed>(nthrds);
  collCnt = vector<int>(nthrds, 0);
  for (int n=0; n<nthrds; n++) {
    initTime[n].Microseconds();
    collTime[n].Microseconds();
    elasTime[n].Microseconds();
    stat1Time[n].Microseconds();
    stat2Time[n].Microseconds();
    stat3Time[n].Microseconds();
    coolTime[n].Microseconds();
  }
  
  if (TSDIAG) {
    tdiag  = vector<unsigned>(numdiag, 0);
    tdiag1 = vector<unsigned>(numdiag, 0);
    tdiag0 = vector<unsigned>(numdiag, 0);
    tdiagT = vector< vector<unsigned> > (nthrds);

    Eover  = vector<double>(numdiag, 0);
    Eover1 = vector<double>(numdiag, 0);
    Eover0 = vector<double>(numdiag, 0);
    EoverT = vector< vector<double> > (nthrds);
  }

  tcool  = vector<unsigned>(numdiag, 0);
  tcool1 = vector<unsigned>(numdiag, 0);
  tcool0 = vector<unsigned>(numdiag, 0);
  tcoolT = vector< vector<unsigned> > (nthrds);

  if (CBADIAG) {
    Cover  = vector<double>(numdiag, 0);
    Cover1 = vector<double>(numdiag, 0);
    Cover0 = vector<double>(numdiag, 0);
    CoverT = vector< vector<double> > (nthrds);
  }

  for (int n=0; n<nthrds; n++) {
    if (TSDIAG) {
      tdiagT[n] = vector<unsigned>(numdiag, 0);
      EoverT[n] = vector<double>(numdiag, 0);
    }
    if (CBADIAG) {
      CoverT[n] = vector<double>(numdiag, 0);
    }
    tcoolT[n] = vector<unsigned>(numdiag, 0);
    tdispT[n] = vector<double>(3, 0);
  }

  disptot = vector<double>(3, 0);
  masstot = 0.0;

  use_temp = -1;
  use_dens = -1;
  use_delt = -1;
  use_exes = -1;

  gen = new ACG(11+myid);
  unit = new Uniform(0.0, 1.0, gen);
  norm = new Normal(0.0, 1.0, gen);

  if (MFPDIAG) {
    prec = vector<Precord>(nthrds);
    for (int n=0; n<nthrds; n++)
      prec[n].second = vector<double>(Nmfp, 0);
  }

}

Collide::~Collide()
{
  delete gen;
  delete unit;
  delete norm;
}

void Collide::debug_list(pHOT& tree)
{
  unsigned ncells = tree.Number();
  pHOT_iterator c(tree);
  for (int cid=0; cid<numprocs; cid++) {
    if (myid == cid) {
      ostringstream sout;
      sout << "==== Collide " << myid << " ncells=" << ncells;
      cout << setw(70) << setfill('=') << left << sout.str() 
	   << endl << setfill(' ');

      for (int n=0; n<nthrds; n++) {
	int nbeg = ncells*(n  )/nthrds;
	int nend = ncells*(n+1)/nthrds;
	for (int j=nbeg; j<nend; j++) {
	  int tnum = c.nextCell();
	  cout << setw(8)  << j
	       << setw(12) << c.Cell()
	       << setw(12) << cellist[n][j-nbeg]
	       << setw(12) << c.Cell()->bods.size()
	       << setw(12) << tnum << endl;
	}
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}


unsigned Collide::collide(pHOT& tree, double Fn, double tau, int mlevel,
			  bool diag)
{
  snglTime.start();
				// Initialize diagnostic counters
				// 
  if (diag) pre_collide_diag();
				// Make cellist
				// 
  for (int n=0; n<nthrds; n++) cellist[n].clear();
  ncells = 0;
  set<pCell*>::iterator ic;
  for (unsigned M=mlevel; M<=multistep; M++) {
    if (tree.clevels[M].size())	// Don't queue null cells
      for (ic=tree.clevels[M].begin(); ic!=tree.clevels[M].end(); ic++) {
	cellist[(ncells++)%nthrds].push_back(*ic);
	bodycount += (*ic)->bods.size();
      }
  }
  stepcount++;
      
#ifdef DEBUG
  debug_list(tree);
#endif
  snglTime.stop();

  forkTime.start();
  if (0) {
    ostringstream sout;
    sout << "before fork, " << __FILE__ << ": " << __LINE__;
    tree.checkBounds(2.0, sout.str().c_str());
  }
  collide_thread_fork(&tree, Fn, tau);
  if (0) {
    ostringstream sout;
    sout << "after fork, " << __FILE__ << ": " << __LINE__;
    tree.checkBounds(2.0, sout.str().c_str());
  }
  forkSoFar = forkTime.stop();

  snglTime.start();
  unsigned col = 0;
  if (diag) col = post_collide_diag();
  snglSoFar = snglTime.stop();

  return( col );
}

void Collide::dispersion(vector<double>& disp)
{
  disp = disptot;
  if (masstot>0.0) {
    for (unsigned k=0; k<3; k++) disp[k] /= masstot;
  }
  for (unsigned k=0; k<3; k++) disptot[k] = 0.0;
  masstot = 0.0;
}


void * Collide::collide_thread(void * arg)
{
  pHOT *tree = (pHOT*)((thrd_pass_arguments*)arg)->tree;
  double Fn = (double)((thrd_pass_arguments*)arg)->fn;
  double tau = (double)((thrd_pass_arguments*)arg)->tau;
  int id = (int)((thrd_pass_arguments*)arg)->id;

				// Work vectors
  vector<double> vcm(3), vrel(3), crel(3);
  pCell *c;

  // Loop over cells, processing collisions in each cell
  //
  for (unsigned j=0; j<cellist[id].size(); j++ ) {

    // Number of particles in this cell
    //
    c = cellist[id][j];
    unsigned number = c->bods.size();
    numcntT[id].push_back(number);

    // Nbody list
    //
    vector<unsigned> bodx;
    set<unsigned>::iterator ib = c->bods.begin();
    for (ib=c->bods.begin(); ib!=c->bods.end(); ib++) bodx.push_back(*ib);

    // Skip cells with only one particle
    //
    if( number < 2 ) {
      colcntT[id].push_back(0);
      continue;  // Skip to the next cell
    }

    stat1Time[id].start();

    // Energy lost in this cell
    //
    decolT[id] = 0.0;
    decelT[id] = 0.0;

    // Compute 1.5 times the mean relative velocity in each MACRO cell
    //
    pCell *samp = c->sample;
    double crm=samp->CRMavg(), crmax=0.0;

    if (!NTC || crm<0.0) {
      crm = 0.0;
      if (samp->state[0]>0.0) {
	for (unsigned k=0; k<3; k++) 
	  crm += (samp->state[1+k] - 
		  samp->state[4+k]*samp->state[4+k]/samp->state[0])/samp->state[0];
      }
      crm = sqrt(fabs(2.0*crm));
      if (NTC) crm *= 1.5;
    }
    
    stat1SoFar[id] = stat1Time[id].stop();
    stat2Time[id].start();

    // KE in the cell
    //
    double kedsp=0.0;
    if (MFPDIAG) {
      if (c->state[0]>0.0) {
	for (unsigned k=0; k<3; k++) 
	  kedsp += 
	    0.5*(c->state[1+k] - c->state[4+k]*c->state[4+k]/c->state[0]);
      }
    }
    
    // Volume in the cell
    //
    double volc = c->Volume();

    // Mass in the cell
    //
    double mass = c->Mass();

    if (mass <= 0.0) continue;

    // Fiducial cross section
    //
    diam = diam0;
    double cross  = M_PI*diam*diam;

    // Determine cross section based on fixed number of collisions
    //
    if (CNUM) {
      cross = 2.0*CNUM*volc/(Fn*mass*tau*crm*number*(number-1));
      diam = sqrt(cross/M_PI);
    }

    double diamCBA = sqrt(Fn*mass)*diam;

    if (MFPDIAG) {

      // Diagnostic: MFP to linear cell size ratio 
      //
      tsratT[id].push_back(crm*tau/pow(volc,0.33333333));
      tdensT[id].push_back(number/volc);
      tvolcT[id].push_back(volc);
    
      double posx, posy, posz;
      c->MeanPos(posx, posy, posz);
      
      // MFP = 1/(n*cross_section)
      // MFP/side = MFP/vol^(1/3) = vol^(2/3)/(number*cross_section)
      
      prec[id].first = pow(volc, 0.66666667)/(Fn*mass*cross*number);
      prec[id].second[0] = sqrt(posx*posx + posy*posy);
      prec[id].second[1] = posz;
      prec[id].second[2] = sqrt(posx*posx+posy*posy*+posz*posz);
      prec[id].second[3] = mass/volc;
      prec[id].second[4] = volc;
      
      tmfpstT[id].push_back(prec[id]);
    }
      
    // Determine number of candidate collision pairs
    // to be selected in this cell
    //
    double coeff  = 0.5*number*(number-1)*Fn*mass/volc*cross*tau;
    double select = coeff*crm;

    if (TSDIAG) {		// Diagnose time step in this cell
      double vmass;
      vector<double> V1, V2;
      c->Vel(vmass, V1, V2);
      double scale = c->Scale();
      double taudiag = 1.0e40;
      for (int k=0; k<3; k++) {	// Time of flight
	taudiag = min<double>
	  (pHOT::sides[k]*scale/(fabs(V1[k]/vmass)+sqrt(V2[k]/vmass)+1.0e-40), 
	   taudiag);
      }
      
      int indx = (int)floor(log(taudiag/tau)/log(4.0) + 5);
      if (indx<0 ) indx = 0;
      if (indx>10) indx = 10;
      tdiagT[id][indx]++;
    }
				// Number per selection ratio
    if (MFPDIAG)
      tselnT[id].push_back(select/number);

    stat2SoFar[id] = stat2Time[id].stop();
    

    // double length = pow(c->Volume(), 0.3333333);

				// Number of pairs to be selected
    unsigned nsel = (int)floor(select+0.5);
    
    initTime[id].start();
    initialize_cell(c, crm, tau, select, id);
    collCnt[id]++;
    initSoFar[id] = initTime[id].stop();

				// No collisions, primarily for testing . . .
    if (DRYRUN) continue;

    collTime[id].start();
				// If more than EPSMratio collisions per
				// particle, assume equipartition
    if (static_cast<double>(number)/select < EPSMratio) {

      EPSM(tree, c, id);

    } else {

      unsigned colc = 0;

      // Loop over total number of candidate collision pairs
      //
      for (unsigned i=0; i<nsel; i++ ) {

	// Pick two particles at random out of this cell
	//
	unsigned k1 = min<int>((int)floor((*unit)()*number), number-1);
	unsigned k2 = ((int)floor((*unit)()*(number-1)) + k1 + 1) % number;
	Particle* p1 = tree->Body(bodx[k1]); // First particle
	Particle* p2 = tree->Body(bodx[k2]); // Second particle
	
	// Calculate pair's relative speed (pre-collision)
	//
	double cr = 0.0;
	for (int k=0; k<3; k++) {
	  crel[k] = p1->vel[k] - p2->vel[k];
	  cr += crel[k]*crel[k];
	}
	cr = sqrt(cr);
	if (NTC) crmax = max<double>(crmax, cr);
	
	// Accept or reject candidate pair according to relative speed
	//
	bool ok = false;
	if (NTC)
	  ok = ( cr/crm > (*unit)() );
	else
	  ok = true;

	if (ok) {
	  elasTime[id].start();

	  // If pair accepted, select post-collision velocities
	  //
	  colc++;			// Collision counter

	  // Do inelastic stuff
	  //
	  error1T[id] += inelastic(tree, p1, p2, &cr, id);
	  // May update relative velocity to reflect
	  // excitation of internal degrees of freedom
	  
	  // Center of mass velocity
	  //
	  double tmass = p1->mass + p2->mass;
	  for(unsigned k=0; k<3; k++)
	    vcm[k] = (p1->mass*p1->vel[k] + p2->mass*p2->vel[k]) / tmass;


	  double cos_th = 1.0 - 2.0*(*unit)();       // Cosine and sine of
	  double sin_th = sqrt(1.0 - cos_th*cos_th); // collision angle theta
	  double phi = 2.0*M_PI*(*unit)();           // Collision angle phi

	  vrel[0] = cr*cos_th;             // Compute post-collision
	  vrel[1] = cr*sin_th*cos(phi);    // relative velocity
	  vrel[2] = cr*sin_th*sin(phi);

	  // Update post-collision velocities
	  // 
	  for(unsigned k=0; k<3; k++ ) {
	    p1->vel[k] = vcm[k] + p2->mass/tmass*vrel[k];
	    p2->vel[k] = vcm[k] - p1->mass/tmass*vrel[k];
	  }

	  if (CBA) {

	    // Calculate pair's relative speed (post-collision)
	    //
	    cr = 0.0;
	    for (int k=0; k<3; k++) {
	      crel[k] = p1->vel[k] - p2->vel[k] - crel[k];
	      cr += crel[k]*crel[k];
	    }
	    cr = sqrt(cr);
	    
	    // Displacement
	    //
	    if (cr>0.0) {
	      double displ;
	      for (int k=0; k<3; k++) {
		displ = crel[k]*diamCBA/cr;
		p1->pos[k] += displ;
		p2->pos[k] -= displ;
	      }

	      if (CBADIAG) {

		double rat = fabs(displ)/pow(volc,0.33333333);
		int indx = (int)floor(log(rat)/log(4.0) + 5);

		if (indx<0 ) indx = 0;
		if (indx>10) indx = 10;

		CoverT[id][indx] += tmass;
	      }
	    }
	  }
	  
	  elasSoFar[id] = elasTime[id].stop();
	} // Loop over pairs

      }

      if (NTC) samp->CRMadd(crmax);

      // Count collisions
      //
      colcntT[id].push_back(colc);
      sel1T[id] += nsel;
      col1T[id] += colc;

    }
    collSoFar[id] = collTime[id].stop();

    // Compute dispersion diagnostics
    //
    stat3Time[id].start();

    double tmass = 0.0;
    vector<double> velm(3, 0.0), velm2(3, 0.0);
    for (unsigned j=0; j<number; j++) {
      Particle* p = tree->Body(bodx[j]);
      for (unsigned k=0; k<3; k++) {
	velm[k]  += p->mass*p->vel[k];
	velm2[k] += p->mass*p->vel[k]*p->vel[k];
      }
      tmass += p->mass;
    }

    if (tmass>0.0) {
      for (unsigned k=0; k<3; k++) {
	velm[k] /= tmass;
	velm2[k] = velm2[k] - velm[k]*velm[k]*tmass;
	if (velm2[k]>0.0) {
	  tdispT[id][k] += velm2[k];
	  tmassT[id]    += tmass;
	}
      }
    }
    
    // Energy lost from this cell compared to target
    //
    if (coolrate[id]>0.0) {
      if (mass>0.0) {
	double dE, excessT = decolT[id] + decelT[id];

				// Diagnostic
	if (MFPDIAG && kedsp>0.0) {
	  keratT[id].push_back((kedsp   - coolrate[id])/kedsp);
	  deratT[id].push_back((excessT - coolrate[id])/kedsp);
	}
	
	if (use_exes>=0) {	// Spread excess energy into the cell

	  if (ENSEXES) 		// All the energy is spread
	    dE = (excessT    - coolrate[id])/mass;
	  else			// Only the EPSM excess is spread
	    dE = (decelT[id] - coolrate[id])/mass;

	  for (unsigned j=0; j<number; j++) {
	    Particle* p = tree->Body(bodx[j]);
	    p->dattrib[use_exes] += dE*p->mass;
	  }
	}
      }
    }
    exesCT[id] += decolT[id];
    exesET[id] += decelT[id];

    stat3SoFar[id] = stat3Time[id].stop();

  } // Loop over cells

  return (NULL);
}


unsigned Collide::medianNumber() 
{
  MPI_Status s;

  if (myid==0) {
    unsigned num;
    for (int n=1; n<numprocs; n++) {
      MPI_Recv(&num, 1, MPI_UNSIGNED, n, 39, MPI_COMM_WORLD, &s);
      vector<unsigned> tmp(num);
      MPI_Recv(&tmp[0], num, MPI_UNSIGNED, n, 40, MPI_COMM_WORLD, &s);
      numcnt.insert(numcnt.end(), tmp.begin(), tmp.end());
    }

    std::sort(numcnt.begin(), numcnt.end()); 

    if (EXTRA) {
      ofstream out("tmp.numcnt");
      for (unsigned j=0; j<numcnt.size(); j++)
	out << setw(8) << j << setw(18) << numcnt[j] << endl;
    }

    return numcnt[numcnt.size()/2]; 

  } else {
    unsigned num = numcnt.size();
    MPI_Send(&num, 1, MPI_UNSIGNED, 0, 39, MPI_COMM_WORLD);
    MPI_Send(&numcnt[0], num, MPI_UNSIGNED, 0, 40, MPI_COMM_WORLD);

    return 0;
  }
}

unsigned Collide::medianColl() 
{ 
  MPI_Status s;

  if (myid==0) {
    unsigned num;
    for (int n=1; n<numprocs; n++) {
      MPI_Recv(&num, 1, MPI_UNSIGNED, n, 39, MPI_COMM_WORLD, &s);
      vector<unsigned> tmp(num);
      MPI_Recv(&tmp[0], num, MPI_UNSIGNED, n, 40, MPI_COMM_WORLD, &s);
      colcnt.insert(colcnt.end(), tmp.begin(), tmp.end());
    }

    std::sort(colcnt.begin(), colcnt.end()); 

    if (EXTRA) {
      ostringstream ostr;
      ostr << runtag << ".colcnt";
      ofstream out(ostr.str().c_str());
      for (unsigned j=0; j<colcnt.size(); j++)
	out << setw(8) << j << setw(18) << colcnt[j] << endl;
    }

    return colcnt[colcnt.size()/2]; 

  } else {
    unsigned num = colcnt.size();
    MPI_Send(&num, 1, MPI_UNSIGNED, 0, 39, MPI_COMM_WORLD);
    MPI_Send(&colcnt[0], num, MPI_UNSIGNED, 0, 40, MPI_COMM_WORLD);
    return 0;
  }

}

void Collide::collQuantile(vector<double>& quantiles, vector<double>& coll_)
{
  MPI_Status s;

  if (myid==0) {
    unsigned num;
    for (int n=1; n<numprocs; n++) {
      MPI_Recv(&num, 1, MPI_UNSIGNED, n, 39, MPI_COMM_WORLD, &s);
      vector<unsigned> tmp(num);
      MPI_Recv(&tmp[0], num, MPI_UNSIGNED, n, 40, MPI_COMM_WORLD, &s);
      colcnt.insert(colcnt.end(), tmp.begin(), tmp.end());
    }

    std::sort(colcnt.begin(), colcnt.end()); 

    coll_ = vector<double>(quantiles.size());
    for (unsigned j=0; j<quantiles.size(); j++)
      coll_[j] = colcnt[(unsigned)floor(quantiles[j]*colcnt.size())];

    ostringstream ostr;
    ostr << runtag << ".coll_counts";
    ifstream in(ostr.str().c_str());
    in.close();
    if (in.fail()) {
      ofstream out(ostr.str().c_str());
      out << left
	  << setw(14) << "# Time" 
	  << setw(10) << "Quantiles"
	  << setw(10) << "Counts"
	  << endl;
    }

    ofstream out(ostr.str().c_str(), ios::app);
    out << setw(14) << tnow
	<< setw(10) << 0.0
	<< setw(10) << colcnt.front() << endl;
    for (unsigned j=0; j<quantiles.size(); j++) {
      out << setw(14) << tnow
	  << setw(10) << quantiles[j] 
	  << setw(10) << coll_[j] << endl;
    }
    out << setw(14) << tnow
	<< setw(10) << 1.0
	<< setw(10) << colcnt.back() << endl << endl;
    
  } else {
    unsigned num = colcnt.size();
    MPI_Send(&num, 1, MPI_UNSIGNED, 0, 39, MPI_COMM_WORLD);
    MPI_Send(&colcnt[0], num, MPI_UNSIGNED, 0, 40, MPI_COMM_WORLD);
  }
}

void Collide::mfpsizeQuantile(vector<double>& quantiles, 
			      vector<double>& mfp_, 
			      vector<double>& ts_,
			      vector<double>& coll_,
			      vector<double>& cool_,
			      vector<double>& rate_,
			      unsigned &collnum, unsigned &coolnum) 
{
  if (!MFPDIAG) return;

  MPI_Status s;

  if (myid==0) {
    unsigned nmb, num;
    for (int n=1; n<numprocs; n++) {
      MPI_Recv(&nmb, 1, MPI_UNSIGNED, n, 38, MPI_COMM_WORLD, &s);
      MPI_Recv(&num, 1, MPI_UNSIGNED, n, 39, MPI_COMM_WORLD, &s);
      vector<double> tmb(nmb), tmp(num);
      MPI_Recv(&tmp[0], num, MPI_DOUBLE, n, 40, MPI_COMM_WORLD, &s);
      tsrat.insert(tsrat.end(), tmp.begin(), tmp.end());
      MPI_Recv(&tmp[0], num, MPI_DOUBLE, n, 41, MPI_COMM_WORLD, &s);
      tdens.insert(tdens.end(), tmp.begin(), tmp.end());
      MPI_Recv(&tmp[0], num, MPI_DOUBLE, n, 42, MPI_COMM_WORLD, &s);
      tvolc.insert(tvolc.end(), tmp.begin(), tmp.end());
      MPI_Recv(&tmp[0], num, MPI_DOUBLE, n, 43, MPI_COMM_WORLD, &s);
      ttemp.insert(ttemp.end(), tmp.begin(), tmp.end());
      MPI_Recv(&tmp[0], num, MPI_DOUBLE, n, 44, MPI_COMM_WORLD, &s);
      tdelt.insert(tdelt.end(), tmp.begin(), tmp.end());
      MPI_Recv(&tmp[0], num, MPI_DOUBLE, n, 45, MPI_COMM_WORLD, &s);
      tseln.insert(tseln.end(), tmp.begin(), tmp.end());
      MPI_Recv(&tmb[0], nmb, MPI_DOUBLE, n, 46, MPI_COMM_WORLD, &s);
      kerat.insert(kerat.end(), tmb.begin(), tmb.end());
      MPI_Recv(&tmb[0], nmb, MPI_DOUBLE, n, 47, MPI_COMM_WORLD, &s);
      derat.insert(derat.end(), tmb.begin(), tmb.end());

      vector<Precord> tmp2(num);

      MPI_Recv(&tmp[0], num, MPI_DOUBLE, n, 48, MPI_COMM_WORLD, &s);
      for (unsigned k=0; k<num; k++) {
				// Load density
	tmp2[k].first = tmp[k];
				// Initialize record
	tmp2[k].second = vector<double>(Nphase, 0);
      }
      for (unsigned l=0; l<Nphase; l++) {
	MPI_Recv(&tmp[0], num, MPI_DOUBLE, n, 49+l, MPI_COMM_WORLD, &s);
	for (unsigned k=0; k<num; k++) tmp2[k].second[l] = tmp[k];
      }
      tphase.insert(tphase.end(), tmp2.begin(), tmp2.end());

      MPI_Recv(&tmp[0], num, MPI_DOUBLE, n, 49+Nphase, MPI_COMM_WORLD, &s);
      for (unsigned k=0; k<num; k++) {
				// Load mfp
	tmp2[k].first = tmp[k];
				// Initialize record
	tmp2[k].second = vector<double>(Nmfp, 0);
      }
      for (unsigned l=0; l<Nmfp; l++) {
	MPI_Recv(&tmp[0], num, MPI_DOUBLE, n, 50+Nphase+l, MPI_COMM_WORLD, &s);
	for (unsigned k=0; k<num; k++) tmp2[k].second[l] = tmp[k];
      }
      tmfpst.insert(tmfpst.end(), tmp2.begin(), tmp2.end());
    }

    std::sort(tsrat.begin(),  tsrat.end()); 
    std::sort(tdens.begin(),  tdens.end()); 
    std::sort(tvolc.begin(),  tvolc.end()); 
    std::sort(ttemp.begin(),  ttemp.end()); 
    std::sort(tdelt.begin(),  tdelt.end()); 
    std::sort(tseln.begin(),  tseln.end()); 
    std::sort(tphase.begin(), tphase.end());
    std::sort(tmfpst.begin(), tmfpst.end());
    std::sort(kerat.begin(),  kerat.end());
    std::sort(derat.begin(),  derat.end());

    collnum = tseln.size();
    coolnum = kerat.size();

    mfp_  = vector<double>(quantiles.size());
    ts_   = vector<double>(quantiles.size());
    coll_ = vector<double>(quantiles.size());
    cool_ = vector<double>(quantiles.size());
    rate_ = vector<double>(quantiles.size());
    for (unsigned j=0; j<quantiles.size(); j++) {
      if (tmfpst.size())
	mfp_[j]  = tmfpst[(unsigned)floor(quantiles[j]*tmfpst.size())].first;
      else
	mfp_[j] = 0;
      if (tsrat.size())
	ts_[j]   = tsrat [(unsigned)floor(quantiles[j]*tsrat.size()) ];
      else
	ts_[j]   = 0;
      if (tseln.size())
	coll_[j] = tseln [(unsigned)floor(quantiles[j]*tseln.size()) ];
      else
	coll_[j] = 0;
      if (kerat.size())
	cool_[j] = kerat [(unsigned)floor(quantiles[j]*kerat.size()) ];
      else
	cool_[j] = 0;
      if (derat.size())
	rate_[j] = derat [(unsigned)floor(quantiles[j]*derat.size()) ];
      else
	rate_[j] = 0;
    }

    if (SORTED) {
      ostringstream ostr;
      ostr << runtag << ".collide." << this_step;
      ofstream out(ostr.str().c_str());
      out << left << setw(8) << "# N" // Header
	  << setw(18) << "| MFP/L"
	  << setw(18) << "| Cyl radius (MFP)"
	  << setw(18) << "| Vertical (MFP)"
	  << setw(18) << "| Sph radius (MFP)"
	  << setw(18) << "| Density(MFP)"
	  << setw(18) << "| Volume(MFP)"
	  << setw(18) << "| TOF/TS"
	  << setw(18) << "| Density"
	  << setw(18) << "| Cell vol"
	  << setw(18) << "| Cell temp"
	  << setw(18) << "| Cool/part"
	  << setw(18) << "| Number/Nsel"
	  << endl;
      out << "# " << setw(6) << 1;
      for (unsigned k=2; k<13; k++) out << "| " << setw(16) << k;
      out << endl;
      for (unsigned j=0; j<tmfpst.size(); j++)
	out << setw(8) << j 
	    << setw(18) << tmfpst[j].first
	    << setw(18) << tmfpst[j].second[0]
	    << setw(18) << tmfpst[j].second[1]
	    << setw(18) << tmfpst[j].second[2]
	    << setw(18) << tmfpst[j].second[3]
	    << setw(18) << tmfpst[j].second[4]
	    << setw(18) << tsrat[j] 
	    << setw(18) << tdens[j] 
	    << setw(18) << tvolc[j] 
	    << setw(18) << ttemp[j] 
	    << setw(18) << tdelt[j] 
	    << setw(18) << tseln[j] 
	    << endl;
    }


    if (PHASE) {
      ostringstream ostr;
      ostr << runtag << ".phase." << this_step;
      ofstream out(ostr.str().c_str());
      out << left << setw(8) << "# N" // Header
	  << setw(18) << "| Density"
	  << setw(18) << "| Temp"
	  << setw(18) << "| Number"
	  << setw(18) << "| Mass"
	  << setw(18) << "| Volume"
	  << endl;
      out << "# " << setw(6) << 1;
      for (unsigned k=2; k<7; k++) out << "| " << setw(16) << k;
      out << endl;
      for (unsigned j=0; j<tphase.size(); j++) {
	out << setw(8) << j << setw(18) << tphase[j].first;
	for (unsigned k=0; k<Nphase; k++) 
	  out << setw(18) << tphase[j].second[k];
	out << endl;
      }
    }
    
  } else {
    unsigned num = tmfpst.size();
    unsigned nmb = derat.size();
    MPI_Send(&nmb, 1, MPI_UNSIGNED, 0, 38, MPI_COMM_WORLD);
    MPI_Send(&num, 1, MPI_UNSIGNED, 0, 39, MPI_COMM_WORLD);
    MPI_Send(&tsrat[0],  num, MPI_DOUBLE, 0, 40, MPI_COMM_WORLD);
    MPI_Send(&tdens[0],  num, MPI_DOUBLE, 0, 41, MPI_COMM_WORLD);
    MPI_Send(&tvolc[0],  num, MPI_DOUBLE, 0, 42, MPI_COMM_WORLD);
    MPI_Send(&ttemp[0],  num, MPI_DOUBLE, 0, 43, MPI_COMM_WORLD);
    MPI_Send(&tdelt[0],  num, MPI_DOUBLE, 0, 44, MPI_COMM_WORLD);
    MPI_Send(&tseln[0],  num, MPI_DOUBLE, 0, 45, MPI_COMM_WORLD);
    MPI_Send(&kerat[0],  nmb, MPI_DOUBLE, 0, 46, MPI_COMM_WORLD);
    MPI_Send(&derat[0],  nmb, MPI_DOUBLE, 0, 47, MPI_COMM_WORLD);

    vector<double> tmp(num);

    for (unsigned k=0; k<num; k++) tmp[k] = tphase[k].first;
    MPI_Send(&tmp[0],  num, MPI_DOUBLE, 0, 48, MPI_COMM_WORLD);
    for (unsigned l=0; l<Nphase; l++) {
      for (unsigned k=0; k<num; k++) tmp[k] = tphase[k].second[l];
      MPI_Send(&tmp[0],  num, MPI_DOUBLE, 0, 49+l, MPI_COMM_WORLD);
    }

    for (unsigned k=0; k<num; k++) tmp[k] = tmfpst[k].first;
    MPI_Send(&tmp[0],  num, MPI_DOUBLE, 0, 49+Nphase, MPI_COMM_WORLD);
    for (unsigned l=0; l<Nmfp; l++) {
      for (unsigned k=0; k<num; k++) tmp[k] = tmfpst[k].second[l];
      MPI_Send(&tmp[0],  num, MPI_DOUBLE, 0, 50+Nphase+l, MPI_COMM_WORLD);
    }
  }
}

void Collide::EPSM(pHOT* tree, pCell* cell, int id)
{
  if (cell->bods.size()<2) return;

				// Compute mean and variance in each dimension
				// 
  vector<double> mvel(3, 0.0), disp(3, 0.0);
  double mass = 0.0;
  double Exes = 0.0;
  set<unsigned>::iterator ib;
  unsigned nbods = cell->bods.size();
  vector<unsigned> bodx;
  for (ib=cell->bods.begin(); ib!=cell->bods.end(); ib++) {
    Particle* p = tree->Body(*ib);
    bodx.push_back(*ib);
    if (p->mass<=0.0 || isnan(p->mass)) {
      cout << "[crazy mass]";
    }
    for (unsigned k=0; k<3; k++) {
      mvel[k] += p->mass*p->vel[k];
      disp[k] += p->mass*p->vel[k]*p->vel[k];
      if (fabs(p->pos[k])>1.0e6 || isnan(p->pos[k])) {
	cout << "[crazy pos]";
      }
      if (fabs(p->vel[k])>1.0e6 || isnan(p->vel[k])) {
	cout << "[crazy vel]";
      }
    }
    mass += p->mass;
				// Compute the total 
				// undercooled (-) or overcooled (+) 
				// energy.  That is Exes must be added
				// to the internal energy before cooling
				// at this step.
    if (use_exes>=0) {
      Exes += p->dattrib[use_exes];
      p->dattrib[use_exes] = 0;
    }
  }

				// Can't do anything if the gas has no mass
  if (mass<=0.0) return;

  double Einternal = 0.0, Enew;
  for (unsigned k=0; k<3; k++) {
    mvel[k] /= mass;
				// Disp is variance here
    disp[k] = (disp[k] - mvel[k]*mvel[k]*mass)/mass;

				// Crazy value?
    if (disp[k]<0.0) disp[k] = 0.0;

				// Total kinetic energy in COV frame
    Einternal += 0.5*mass*disp[k];
  }
				// Can't collide if with no internal energy
  if (Einternal<=0.0) return;

				// Correct 1d vel. disp. after cooling
				// 
  double Emin = 1.5*boltz*TFLOOR * mass/mp * 
    UserTreeDSMC::Munit/UserTreeDSMC::Eunit;

				// Einternal+Exes is the amount that
				// is *really* available for cooling
				// after excess or deficit is included
				//
				// Exes will be - if more energy still
				// needs to be removed and + if too much
				// energy was removed by cooling last step
				// 
  if (Einternal + Exes - Emin > coolrate[id]) {
    Enew = Einternal + Exes - coolrate[id];
  } else {
    Enew = min<double>(Emin, Einternal);

    decelT[id] += Einternal - Enew + Exes - coolrate[id];

    if (TSDIAG) {
      if (coolrate[id]-Exes>0.0) {
	/*
	int indx = (int)floor(log(Einternal/(coolrate[id]-Exes)) /
			      (log(2.0)*TSPOW) + 5);
	*/
	int indx = (int)floor(log(Einternal/coolrate[id]) /
			      (log(2.0)*TSPOW) + 5);
	if (indx<0 ) indx = 0;
	if (indx>10) indx = 10;

	EoverT[id][indx] += mass;
      }
    }
  }
				// Compute the mean 1d vel.disp. from the
				// new internal energy value
  double mdisp = sqrt(Enew/mass/3.0);

				// Sanity check
				// 
  if (mdisp<=0.0 || isnan(mdisp) || isinf(mdisp)) {
    cout << "Process " << myid  << " id " << id 
	 << ": crazy values, mdisp=" << mdisp << " Enew=" << Enew
	 << " Eint=" << Einternal << " nbods=" << nbods << endl;
    return;
  }
				// Realize new velocities for all particles
				// 
  if (PULLIN) {
    double R=0.0, T=0.0;	// [Shuts up the compile-time warnings]
    const double sqrt3 = sqrt(3.0);

    if (nbods==2) {
      Particle* p1 = tree->Body(bodx[0]);
      Particle* p2 = tree->Body(bodx[1]);
      for (unsigned k=0; k<3; k++) {
	R = (*unit)();
	if ((*unit)()>0.5)
	  p1->vel[k] = mvel[k] + mdisp;
	else 
	  p1->vel[k] = mvel[k] - mdisp;
	p2->vel[k] = 2.0*mvel[k] - p1->vel[k];
      }

    } else if (nbods==3) {
      Particle* p1 = tree->Body(bodx[0]);
      Particle* p2 = tree->Body(bodx[1]);
      Particle* p3 = tree->Body(bodx[2]);
      double v2, v3;
      for (unsigned k=0; k<3; k++) {
	T = 2.0*M_PI*(*unit)();
	v2 = M_SQRT2*mdisp*cos(T);
	v3 = M_SQRT2*mdisp*sin(T);
	p1->vel[k] = mvel[k] - M_SQRT2*v2/sqrt3;
	p2->vel[k] = p1->vel[k] + (sqrt3*v2 - v3)/M_SQRT2;
	p3->vel[k] = p2->vel[k] + M_SQRT2*v3;
      }
    } else if (nbods==4) {
      Particle* p1 = tree->Body(bodx[0]);
      Particle* p2 = tree->Body(bodx[1]);
      Particle* p3 = tree->Body(bodx[2]);
      Particle* p4 = tree->Body(bodx[3]);
      double v2, v3, e2, e4, v4;
      for (unsigned k=0; k<3; k++) {
	R = (*unit)();
	e2 = mdisp*mdisp*(1.0 - R*R);
	T = 2.0*M_PI*(*unit)();
	v2 = sqrt(2.0*e2)*cos(T);
	v3 = sqrt(2.0*e2)*sin(T);
	p1->vel[k] = mvel[k] - sqrt3*v2/2.0;
	p2->vel[k] = p1->vel[k] + (2.0*v2 - M_SQRT2*v3)/sqrt3;
	e4 = mdisp*mdisp*R*R;
	if ((*unit)()>0.5) v4 =  sqrt(2.0*e4);
	else               v4 = -sqrt(2.0*e4);
	p3->vel[k] = p2->vel[k] + (sqrt3*v3 - v4)/M_SQRT2;
	p4->vel[k] = p3->vel[k] + M_SQRT2*v4;
      }

    } else {

      Particle *Pm1, *P00, *Pp1;
      vector<double> Tk, v(nbods), e(nbods);
      int kmax, dim, jj;
      bool Even = (nbods/2*2 == nbods);

      for (int k=0; k<3; k++) {
	if (Even) { 
				// Even
	  kmax = nbods;
	  dim = kmax/2-1;
	  Tk = vector<double>(dim);
	  for (int m=0; m<dim; m++) 
	    Tk[m] = pow((*unit)(), 1.0/(kmax/2 - m - 1.5));
	} else {			
				// Odd
	  kmax = nbods-1;
	  dim = kmax/2-1;
	  Tk = vector<double>(dim);
	  for (int m=0; m<dim; m++) 
	    Tk[m] = pow((*unit)(), 1.0/(kmax/2 - m - 1.0));
	}
      
	e[1] = mdisp*mdisp*(1.0 - Tk[0]);
	T = 2.0*M_PI*(*unit)();
	v[1] = sqrt(2.0*e[1])*cos(T);
	v[2] = sqrt(2.0*e[1])*sin(T);

	P00 = tree->Body(bodx[0]);
	Pp1 = tree->Body(bodx[1]);

	P00->vel[k] = mvel[k] - sqrt(nbods-1)*v[1]/sqrt(nbods);
	Pp1->vel[k] = P00->vel[k] + (sqrt(nbods)*v[1] - sqrt(nbods-2)*v[2])/sqrt(nbods-1);

	double prod = 1.0;
	for (int j=4; j<kmax-1; j+=2) {
	  jj = j-1;

	  Pm1 = tree->Body(bodx[jj-2]);
	  P00 = tree->Body(bodx[jj-1]);
	  Pp1 = tree->Body(bodx[jj  ]);

	  prod *= Tk[j/2-2];
	  e[jj] = mdisp*mdisp*(1.0 - Tk[j/2-1])*prod;
	  T = 2.0*M_PI*(*unit)();
	  v[jj]   = sqrt(2.0*e[jj])*cos(T);
	  v[jj+1] = sqrt(2.0*e[jj])*sin(T);
	  
	  P00->vel[k] = Pm1->vel[k] + 
	    (sqrt(3.0+nbods-j)*v[jj-1] - sqrt(1.0+nbods-j)*v[jj]  )/sqrt(2.0+nbods-j);
	  Pp1->vel[k] = P00->vel[k] +
	    (sqrt(2.0+nbods-j)*v[jj  ] - sqrt(    nbods-j)*v[jj+1])/sqrt(1.0+nbods-j);
	}

	prod *= Tk[kmax/2-2];
	e[kmax-1] = mdisp*mdisp*prod;

	if (Even) {
	  if ((*unit)()>0.5) v[nbods-1] =  sqrt(2.0*e[kmax-1]);
	  else               v[nbods-1] = -sqrt(2.0*e[kmax-1]);
	} else {
	  T = 2.0*M_PI*(*unit)();
	  v[nbods-2] = sqrt(2.0*e[kmax-1])*cos(T);
	  v[nbods-1] = sqrt(2.0*e[kmax-1])*sin(T);

	  Pm1 = tree->Body(bodx[nbods-4]);
	  P00 = tree->Body(bodx[nbods-3]);

	  P00->vel[k] = Pm1->vel[k] + (2.0*v[nbods-3] - M_SQRT2*v[nbods-2])/sqrt3;
	}

	Pm1 = tree->Body(bodx[nbods-3]);
	P00 = tree->Body(bodx[nbods-2]);
	Pp1 = tree->Body(bodx[nbods-1]);

	P00->vel[k] = Pm1->vel[k] + (sqrt3*v[nbods-2] - v[nbods-1])/M_SQRT2;
	Pp1->vel[k] = P00->vel[k] + M_SQRT2*v[nbods-1];
      }
    }

    // End Pullin algorithm

  } else {

				// Realize a distribution with internal
				// dispersion only
    vector<double> Tmvel(3, 0.0);
    vector<double> Tdisp(3, 0.0);
    for (unsigned j=0; j<nbods; j++) {
      Particle* p = tree->Body(bodx[j]);
      for (unsigned k=0; k<3; k++) {
	p->vel[k] = mdisp*(*norm)();
	Tmvel[k] += p->mass*p->vel[k];
	Tdisp[k] += p->mass*p->vel[k]*p->vel[k];
	if (fabs(p->vel[k])>1e6 || isnan(p->vel[k])) {
	  cout << "[Collide crazy vel indx=" << p->indx 
	       << ", key=" << hex << p->key << dec << "]";
	}
      }
    }
				// Compute mean and variance
				// 
    double Tmdisp = 0.0;
    for (unsigned k=0; k<3; k++) {
      Tmvel[k] /= mass;
      Tdisp[k] = (Tdisp[k] - Tmvel[k]*Tmvel[k]*mass)/mass;
      Tmdisp += Tdisp[k];
    }
    Tmdisp = sqrt(Tmdisp/3.0);

				// Sanity check
				// 
    if (Tmdisp<=0.0 || isnan(Tmdisp) || isinf(Tmdisp)) {
      cout << "Process " << myid  << " id " << id 
	   << ": crazy values, Tmdisp=" << Tmdisp << " mdisp=" << mdisp 
	   << " nbods=" << nbods << endl;
      return;
    }
				// Enforce energy and momentum conservation
				// 
    for (unsigned j=0; j<nbods; j++) {
      Particle* p = tree->Body(bodx[j]);
      for (unsigned k=0; k<3; k++)
	p->vel[k] = mvel[k] + (p->vel[k]-Tmvel[k])*mdisp/Tmdisp;
    }

  }

				// Record diagnostics
				// 
  lostSoFar_EPSM[id] += Einternal - Enew;
  epsm1T[id] += nbods;
  Nepsm1T[id]++;
}

void Collide::list_sizes()
{
  string sname = runtag + ".collide_storage";
  for (int n=0; n<numprocs; n++) {
    if (myid==n) {
      ofstream out(sname.c_str(), ios::app);
      if (out) {
	out << setw(18) << tnow
	    << setw(6)  << myid;
	list_sizes_proc(&out);
	out << endl;
	if (myid==numprocs-1) out << endl;
	out.close();
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}


void Collide::list_sizes_proc(ostream* out)
{
  *out << setw(12) << numcnt.size()
       << setw(12) << colcnt.size()
       << setw(12) << tsrat.size()
       << setw(12) << derat.size()
       << setw(12) << tdens.size()
       << setw(12) << tvolc.size()
       << setw(12) << ttemp.size()
       << setw(12) << tdelt.size()
       << setw(12) << tseln.size()
       << setw(12) << tphase.size()
       << setw(12) << (tphaseT.size() ? tphaseT[0].size() : (size_t)0)
       << setw(12) << tmfpst.size()
       << setw(12) << (tmfpstT.size() ? tmfpstT[0].size() : (size_t)0)
       << setw(12) << (numcntT.size() ? numcntT[0].size() : (size_t)0)
       << setw(12) << (colcntT.size() ? colcntT[0].size() : (size_t)0)
       << setw(12) << error1T.size()
       << setw(12) << col1T.size()
       << setw(12) << epsm1T.size()
       << setw(12) << Nepsm1T.size()
       << setw(12) << KEtotT.size()
       << setw(12) << KElostT.size()
       << setw(12) << tmassT.size()
       << setw(12) << decelT.size()
       << setw(12) << (mfpratT.size() ? mfpratT[0].size() : (size_t)0)
       << setw(12) << (tsratT.size() ? tsratT[0].size() : (size_t)0)
       << setw(12) << (tdensT.size() ? tdensT[0].size() : (size_t)0)
       << setw(12) << (tvolcT.size() ? tvolcT[0].size() : (size_t)0)
       << setw(12) << (ttempT.size() ? ttempT[0].size() : (size_t)0)
       << setw(12) << (tselnT.size() ? tselnT[0].size() : (size_t)0)
       << setw(12) << (deratT.size() ? deratT[0].size() : (size_t)0)
       << setw(12) << (tdeltT.size() ? tdeltT[0].size() : (size_t)0)
       << setw(12) << (tdispT.size() ? tdispT[0].size() : (size_t)0)
       << setw(12) << tdiag.size()
       << setw(12) << tdiag1.size()
       << setw(12) << tdiag0.size()
       << setw(12) << tcool.size()
       << setw(12) << tcool1.size()
       << setw(12) << tcool0.size()
       << setw(12) << (cellist.size() ? cellist[0].size() : (size_t)0)
       << setw(12) << disptot.size()
       << setw(12) << lostSoFar_EPSM.size();
}


void Collide::colldeTime(ostream& out) 
{ 
  out << "-----------------------------------------------------" << endl;
  out << "-----Collide timing----------------------------------" << endl;
  out << "-----------------------------------------------------" << endl;
  out << left
      << setw(18) << "Thread time"  << threadTime() << endl
      << setw(18) << "Joined time"  << joinedTime() << endl
      << setw(18) << "Waiting time" << waitngTime() << endl
      << setw(18) << "Join time"    << joinngTime() << endl
      << setw(18) << "Diag time"    << dgnoscTime() << endl
      << setw(18) << "Body count"   << bodycount    
      << " [" << bodycount/stepcount << "]" << endl
      << setw(18) << "Step count"   << stepcount    << endl
      << endl;
  
#ifdef DEBUG
  out << left << setw(4) << "#"
      << setw(10) << "Init"
      << setw(10) << "Collide"
      << setw(10) << "Stat1"
      << setw(10) << "Stat2"
      << setw(10) << "Stat3"
      << setw(10) << "Cool"
      << setw(10) << "Inelastic"
      << setw(10) << "Cell count"
      << endl;
#endif
  
  for (int n=0; n<nthrds; n++) {
#ifdef DEBUG
    out << setw(4) << n 
	<< setw(10) << initSoFar[n]()*1.0e-6 
	<< setw(10) << collSoFar[n]()*1.0e-6 
	<< setw(10) << stat1SoFar[n]()*1.0e-6 
	<< setw(10) << stat2SoFar[n]()*1.0e-6 
	<< setw(10) << stat3SoFar[n]()*1.0e-6 
	<< setw(10) << coolSoFar[n]()*1.0e-6 
	<< setw(10) << elasSoFar[n]()*1.0e-6 
	<< setw(10) << collCnt[n] << endl
	<< setw(4) << "*" << setprecision(4)
	<< setw(10) << initSoFar[n]()*1.0e-6/stepcount
	<< setw(10) << collSoFar[n]()*1.0e-6/stepcount
	<< setw(10) << stat1SoFar[n]()*1.0e-6/stepcount
	<< setw(10) << stat2SoFar[n]()*1.0e-6/stepcount
	<< setw(10) << stat3SoFar[n]()*1.0e-6/stepcount
	<< setw(10) << coolSoFar[n]()*1.0e-6/stepcount
	<< setw(10) << elasSoFar[n]()*1.0e-6/stepcount
	<< setw(10) << collCnt[n]/stepcount << endl
	<< endl << setprecision(2);
#endif
    initTime[n].reset(); 
    collTime[n].reset(); 
    stat1Time[n].reset();
    stat2Time[n].reset();
    stat3Time[n].reset();
    coolTime[n].reset(); 
    elasTime[n].reset(); 
    coolTime[n].reset(); 
    collCnt[n] = 0;
  }

  stepcount = 0;
  bodycount = 0;

  out << "-----------------------------------------------------" << endl;
}

void Collide::tsdiag(ostream& out) 
{
  if (!TSDIAG) return;

  if (tdiag.size()==numdiag) {

    out << "-----------------------------------------------------" << endl;
    out << "-----Time step diagnostics---------------------------" << endl;
    out << "-----------------------------------------------------" << endl;
    out << right << setw(8) << "2^n" << setw(15) << "TS ratio"
	<< setw(15) << "Size/Vel";
    if (use_delt>=0) out << setw(15) << "Kinetic/Cool";
    out << endl << setprecision(3);
    out << "-----------------------------------------------------" << endl;
    for (unsigned k=0; k<numdiag; k++) {
      double rat = pow(4.0, -5.0+k);
      out << setw(8)  << -10+2*static_cast<int>(k)
	  << ((rat<1.0e-02 && rat>0.0) ? scientific : fixed)
	  << setw(15) << rat
	  << setw(15) << tdiag[k];
      if (use_delt>=0) out << setw(15) << tcool[k];
      out << endl;
      tdiag[k] = 0;
      if (use_delt>=0) tcool[k] = 0;
    }
    out << "-----------------------------------------------------" << endl;
    out << left;
  }


  if (Eover.size()==numdiag) {

    double emass = 0.0;
    for (unsigned k=0; k<numdiag; k++) emass += Eover[k];
    if (emass>0.0) {
      out << "-----Cooling rate diagnostics------------------------" << endl;
      out << "-----------------------------------------------------" << endl;
      out << right << setw(8) << "2^n" << setw(15) << "Ratio"
	  << setw(15) << "KE/Cool(%)" << endl;
      out << "-----------------------------------------------------" << endl;
      
      for (unsigned k=0; k<numdiag; k++) {
	double rat = pow(pow(2.0, TSPOW), -5.0+k);
	double val = Eover[k]*100.0/emass;
	out << setw(8)  << TSPOW*(-5 + static_cast<int>(k))
	    << ((rat<1.0e-02 && rat>0.0) ? scientific : fixed)
	    << setw(15) << rat
	    << ((val<1.0e-02 && val>0.0) ? scientific : fixed)
	    << setw(15) << val << endl;
	Eover[k] = 0;
      }
      out << "-----------------------------------------------------" << endl;
    }
    out << left;
  }

  if (Cover.size()==numdiag) {

    double cmass = 0.0;
    for (unsigned k=0; k<numdiag; k++) cmass += Cover[k];
    if (cmass>0.0) {
      out << "-----CBA scale diagnostics--------------------------" << endl;
      out << "-----------------------------------------------------" << endl;
      out << right << setw(8) << "2^n" << setw(15) << "Ratio"
	  << setw(15) << "Diam/Side(%)" << endl;
      out << "-----------------------------------------------------" << endl;
      
      for (unsigned k=0; k<numdiag; k++) {
	double rat = pow(4.0, -5.0+k);
	double val = Cover[k]*100.0/cmass;
	out << setw(8)  << -10+2*static_cast<int>(k)
	    << ((rat<1.0e-02 && rat>0.0) ? scientific : fixed)
	    << setw(15) << rat
	    << ((val<1.0e-02 && val>0.0) ? scientific : fixed)
	    << setw(15) << val << endl;
	Cover[k] = 0;
      }
      out << "-----------------------------------------------------" << endl;
    }
    out << left;
  }

}


// For timestep computation

extern "C"
void *
tstep_thread_call(void *atp)
{
  thrd_pass_tstep *tp = (thrd_pass_tstep *)atp;
  Collide *p = (Collide *)tp->p;
  p->timestep_thread((void*)&tp->arg);
  return NULL;
}

void Collide::compute_timestep(pHOT* tree, double coolfrac)
{
  int errcode;
  void *retval;
  
  if (nthrds==1) {
    thrd_pass_tstep td;

    td.p = this;
    td.arg.tree = tree;
    td.arg.coolfrac = coolfrac;
    td.arg.id = 0;

    tstep_thread_call(&td);
    
    return;
  }

  tdT = new thrd_pass_tstep [nthrds];
  t = new pthread_t [nthrds];

  if (!tdT) {
    cerr << "Process " << myid 
         << ": Collide::tstep_thread_call: error allocating memory for thread counters\n";
    exit(18);
  }
  if (!t) {
    cerr << "Process " << myid
         << ": Collide::tstep_thread_call: error allocating memory for thread\n";
    exit(18);
  }

                                // Make the <nthrds> threads
  for (int i=0; i<nthrds; i++) {
    tdT[i].p = this;
    tdT[i].arg.tree = tree;
    tdT[i].arg.coolfrac = coolfrac;
    tdT[i].arg.id = i;

    errcode =  pthread_create(&t[i], 0, tstep_thread_call, &tdT[i]);
    if (errcode) {
      cerr << "Process " << myid;
      cerr << " Collide: cannot make thread " << i
	   << ", errcode=" << errcode << endl;
      exit(19);
    }
  }
    
                                // Collapse the threads
  for (int i=0; i<nthrds; i++) {
    if ((errcode=pthread_join(t[i], &retval))) {
      cerr << "Process " << myid;
      cerr << " Collide::tstep_thread_call: thread join " << i
           << " failed, errcode=" << errcode << endl;
      exit(20);
    }
  }
  
  delete [] tdT;
  delete [] t;
}


void * Collide::timestep_thread(void * arg)
{
  pHOT* tree = (pHOT* )((tstep_pass_arguments*)arg)->tree;
  double coolfrac = (double)((tstep_pass_arguments*)arg)->coolfrac;
  int id = (int)((tstep_pass_arguments*)arg)->id;

  // Loop over cells, cell time-of-flight time
  // for each particle

  pCell *c;
  Particle *p;
  unsigned nbods;
  double L, DT, mscale;

  // Loop over cells, processing collisions in each cell
  //
  for (unsigned j=0; j<cellist[id].size(); j++ ) {

    // Number of particles in this cell
    //
    c = cellist[id][j];
    nbods = c->bods.size();
    L = c->Scale();
    
    for (set<unsigned>::iterator i=c->bods.begin(); i!=c->bods.end(); i++) {
				// Current particle
      p = tree->Body(*i);
				// Compute time of flight criterion
      DT = 1.0e40;
      mscale = 1.0e40;
      for (unsigned k=0; k<3; k++) {
	DT = min<double>
	  (pHOT::sides[k]*L/(fabs(p->vel[k])+1.0e-40), DT);
	mscale = min<double>(pHOT::sides[k]*L, mscale);
      }
				// Size scale for multistep timestep calc.
      p->scale = mscale;
				// Compute cooling criterion timestep
      double v = p->dattrib[use_delt];
      if (use_delt>=0 && v>0.0)	DT = min<double>(DT, v);
      
      p->dtreq = coolfrac*DT;
    }
  }
  
  return (NULL);
}

void Collide::energyExcess(double& ExesColl, double& ExesEPSM)
{
				// Sum up from all the threads
				// 
  for (int n=1; n<nthrds; n++) {
    exesCT[0] += exesCT[n];
    exesET[0] += exesET[n];
  }
				// Sum reduce result to root node
				// 
  MPI_Reduce(&exesCT[0], &ExesColl, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&exesET[0], &ExesEPSM, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

				// Zero out the thread accumulators
				// 
  for (int n=0; n<nthrds; n++) exesCT[0] = exesET[n] = 0.0;
}


void Collide::pre_collide_diag()
{
				// Clean thread variables
  diagTime.start();
  for (int n=0; n<nthrds; n++) {
    error1T[n] = 0;
    sel1T[n] = 0;
    col1T[n] = 0;
    epsm1T[n] = 0;
    Nepsm1T[n] = 0;
    tmassT[n] = 0;
    decolT[n] = 0;
    decelT[n] = 0;

    // For computing cell occupation #
    colcntT[n].clear();	// and collision counts
    numcntT[n].clear();

    // For computing MFP to cell size ratio 
    // and drift ratio
    if (MFPDIAG) {
      tsratT[n].clear();
      keratT[n].clear();
      deratT[n].clear();
      tdensT[n].clear();
      tvolcT[n].clear();
      ttempT[n].clear();
      tdeltT[n].clear();
      tselnT[n].clear();
      tphaseT[n].clear();
      tmfpstT[n].clear();
    }
    
    for (unsigned k=0; k<numdiag; k++) {
      if (TSDIAG) {
	tdiagT[n][k] = 0;
	EoverT[n][k] = 0;
      }
      if (CBADIAG)     CoverT[n][k] = 0;
      if (use_delt>=0) tcoolT[n][k] = 0;
    }
      
    for (unsigned k=0; k<3; k++) tdispT[n][k] = 0;
  }
  if (TSDIAG) {
    for (unsigned k=0; k<numdiag; k++) tdiag1[k] = tdiag0[k] = 0;
    for (unsigned k=0; k<numdiag; k++) Eover1[k] = Eover0[k] = 0;
  }
  if (CBADIAG) {
    for (unsigned k=0; k<numdiag; k++) Cover1[k] = Cover0[k] = 0;
  }
  if (use_delt>=0) 
    for (unsigned k=0; k<numdiag; k++) tcool1[k] = tcool0[k] = 0;
  
  diagTime.stop();
  
#ifdef DEBUG
  list_sizes();
#endif
}


unsigned Collide::post_collide_diag()
{
  unsigned sel=0, col=0;

  diagTime.start();
				// Diagnostics
  unsigned error1=0, error=0;

  unsigned sel1=0, col1=0;	// Count number of collisions
  unsigned epsm1=0, epsm=0, Nepsm1=0, Nepsm=0;

				// Dispersion test
  double mass1 = 0, mass0 = 0;
  vector<double> disp1(3, 0), disp0(3, 0);

  numcnt.clear();
  colcnt.clear();

  for (int n=0; n<nthrds; n++) {
    error1 += error1T[n];
    sel1   += sel1T[n];
    col1   += col1T[n];
    epsm1  += epsm1T[n];
    Nepsm1 += Nepsm1T[n];
    numcnt.insert(numcnt.end(), numcntT[n].begin(), numcntT[n].end());
    colcnt.insert(colcnt.end(), colcntT[n].begin(), colcntT[n].end());
    if (TSDIAG) {
      for (unsigned k=0; k<numdiag; k++) tdiag1[k] += tdiagT[n][k];
      for (unsigned k=0; k<numdiag; k++) Eover1[k] += EoverT[n][k];
    }
    if (CBADIAG)
      for (unsigned k=0; k<numdiag; k++) Cover1[k] += CoverT[n][k];
    if (use_delt>=0) 
      for (unsigned k=0; k<numdiag; k++) tcool1[k] += tcoolT[n][k];
  }

				// For computing MFP to cell size ratio 
				// and drift ratio (diagnostic only)
  if (MFPDIAG) {
    tsrat.clear();
    kerat.clear();
    derat.clear();
    tdens.clear();
    tvolc.clear();
    ttemp.clear();
    tdelt.clear();
    tseln.clear();
    tphase.clear();
    tmfpst.clear();

    for (int n=0; n<nthrds; n++) {
      tsrat. insert(tsrat.end(),   tsratT[n].begin(),  tsratT[n].end());
      kerat. insert(kerat.end(),   keratT[n].begin(),  keratT[n].end());
      derat. insert(derat.end(),   deratT[n].begin(),  deratT[n].end());
      tdens. insert(tdens.end(),   tdensT[n].begin(),  tdensT[n].end());
      tvolc. insert(tvolc.end(),   tvolcT[n].begin(),  tvolcT[n].end());
      ttemp. insert(ttemp.end(),   ttempT[n].begin(),  ttempT[n].end());
      tdelt. insert(tdelt.end(),   tdeltT[n].begin(),  tdeltT[n].end());
      tseln. insert(tseln.end(),   tselnT[n].begin(),  tselnT[n].end());
      tphase.insert(tphase.end(), tphaseT[n].begin(), tphaseT[n].end());
      tmfpst.insert(tmfpst.end(), tmfpstT[n].begin(), tmfpstT[n].end());
    }
  }

  for (int n=0; n<nthrds; n++) {
    for (unsigned k=0; k<3; k++) disp1[k] += tdispT[n][k];
    mass1 += tmassT[n];
  }

  MPI_Reduce(&sel1, &sel, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&col1, &col, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&epsm1, &epsm, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&Nepsm1, &Nepsm, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&error1, &error, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&ncells, &numtot, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  if (TSDIAG) {
    MPI_Reduce(&tdiag1[0], &tdiag0[0], numdiag, MPI_UNSIGNED, MPI_SUM, 0, 
	       MPI_COMM_WORLD);
    MPI_Reduce(&Eover1[0], &Eover0[0], numdiag, MPI_DOUBLE,   MPI_SUM, 0, 
	       MPI_COMM_WORLD);
  }
  if (CBADIAG) {
    MPI_Reduce(&Cover1[0], &Cover0[0], numdiag, MPI_DOUBLE,   MPI_SUM, 0, 
	       MPI_COMM_WORLD);
  }
  if (use_delt>=0)
    MPI_Reduce(&tcool1[0], &tcool0[0], numdiag, MPI_UNSIGNED, MPI_SUM, 0, 
	       MPI_COMM_WORLD);
  MPI_Reduce(&disp1[0], &disp0[0], 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&mass1, &mass0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  seltot    += sel;
  coltot    += col;
  epsmtot   += epsm;
  epsmcells += Nepsm;
  errtot    += error;
  if (TSDIAG) {
    for (unsigned k=0; k<numdiag; k++) {
      tdiag[k] += tdiag0[k];
      Eover[k] += Eover0[k];
    }
  }
  if (CBADIAG) {
    for (unsigned k=0; k<numdiag; k++)
      Cover[k] += Cover0[k];
  }
  if (use_delt>=0)
    for (unsigned k=0; k<numdiag; k++) tcool[k] += tcool0[k];
  for (unsigned k=0; k<3; k++) disptot[k] += disp0[k];
  masstot   += mass0;
  diagSoFar = diagTime.stop();

  return col;
}

