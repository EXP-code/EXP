#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>

using namespace std;

#include "Timer.h"
#include "pHOT.H"
#include "Collide.H"

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

Collide::Collide(double diameter, int nth)
{
  nthrds = nth;

  colcntT = vector< vector<unsigned> > (nthrds);
  numcntT = vector< vector<unsigned> > (nthrds);
  error1T = vector<unsigned> (nthrds, 0);
  col1T = vector<unsigned> (nthrds, 0);
  KEtotT = vector<double> (nthrds, 0);
  KElostT = vector<double> (nthrds, 0);

  mfpratT = vector< vector<double> > (nthrds);
  tsratT  = vector< vector<double> > (nthrds);
  tdensT  = vector< vector<double> > (nthrds);
  tvolcT  = vector< vector<double> > (nthrds);

  cellist = vector< vector<pCell*> > (nthrds);

  diam0 = diam = diameter;
  coltot = 0;			// Count total collisions
  errtot = 0;			// Count errors in inelastic computation


  snglTime.Microseconds();
  forkTime.Microseconds();
  waitTime.Microseconds();
  joinTime.Microseconds();

  collTime = vector<Timer>(nthrds);
  collSoFar = vector<TimeElapsed>(nthrds);
  collCnt = vector<int>(nthrds, 0);
  for (int n=0; n<nthrds; n++) collTime[n].Microseconds();
  
  tdiag  = vector<unsigned>(numdiag, 0);
  tdiag1 = vector<unsigned>(numdiag, 0);
  tdiag0 = vector<unsigned>(numdiag, 0);
  tdiagT = vector< vector<unsigned> > (nthrds);
  for (int n=0; n<nthrds; n++) 
    tdiagT[n] = vector<unsigned>(numdiag, 0);

  use_temp = -1;
  use_dens = -1;

  gen = new ACG(11+myid);
  unit = new Uniform(0.0, 1.0, gen);
}

Collide::~Collide()
{
  delete gen;
  delete unit;
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


unsigned Collide::collide(pHOT& tree, double Fn, double tau)
{
  snglTime.start();

				// Clean thread variables
  for (int n=0; n<nthrds; n++) {
    error1T[n] = 0;
    col1T[n] = 0;
				// Diagnostics
    KEtotT[n] = 0.0;
    KElostT[n] = 0.0;
				// For computing cell occupation #
    colcntT[n].clear();		// and collision counts
    numcntT[n].clear();
				// For computing MFP to cell size ratio 
				// and drift ratio
    mfpratT[n].clear();
    tsratT[n].clear();
    tdensT[n].clear();
    tvolcT[n].clear();

    for (unsigned k=0; k<numdiag; k++) tdiagT[n][k] = 0;
  }
  for (unsigned k=0; k<numdiag; k++) tdiag1[k] = tdiag0[k] = 0;


				// Make cellist
  unsigned ncells = tree.Number();
  pHOT_iterator c(tree);
  for (int n=0; n<nthrds; n++) {
    cellist[n].clear();
    int nbeg = ncells*(n  )/nthrds;
    int nend = ncells*(n+1)/nthrds;
    for (int j=nbeg; j<nend; j++) {
      c.nextCell();
      cellist[n].push_back(c.Cell());
    }
  }
      
#ifdef DEBUG
  debug_list(tree);
#endif

  snglTime.stop();
  forkTime.start();
  {
    ostringstream sout;
    sout << "before fork, " << __FILE__ << ": " << __LINE__;
    tree.checkBounds(2.0, sout.str().c_str());
  }
  collide_thread_fork(&tree, Fn, tau);
  {
    ostringstream sout;
    sout << "after fork, " << __FILE__ << ": " << __LINE__;
    tree.checkBounds(2.0, sout.str().c_str());
  }
  forkSoFar = forkTime.stop();
  snglTime.start();

				// Diagnostics
  unsigned error1=0, error=0;

  unsigned col1=0, col=0;	// Count number of collisions

  numcnt.clear();
  colcnt.clear();

  for (int n=0; n<nthrds; n++) {
    error1 += error1T[n];
    col1 += col1T[n];
    numcnt.insert(numcnt.end(), numcntT[n].begin(), numcntT[n].end());
    colcnt.insert(colcnt.end(), colcntT[n].begin(), colcntT[n].end());
    for (unsigned k=0; k<numdiag; k++) tdiag1[k] += tdiagT[n][k];
#ifdef DIAG
    KEtot += KEtotT[n];
    KElost += KElostT[n];
#endif
    
  }

				// For computing MFP to cell size ratio 
				// and drift ratio (diagnostic only)
  mfprat.clear();
  tsrat.clear();
  tdens.clear();
  tvolc.clear();
  for (int n=0; n<nthrds; n++) {
    mfprat.insert(mfprat.end(), mfpratT[n].begin(), mfpratT[n].end());
    tsrat. insert(tsrat.end(),   tsratT[n].begin(),  tsratT[n].end());
    tdens. insert(tdens.end(),   tdensT[n].begin(),  tdensT[n].end());
    tvolc. insert(tvolc.end(),   tvolcT[n].begin(),  tvolcT[n].end());
  }

#ifdef DIAG
  // cout << "Process " << myid << ": lost=" << KElost/KEtot << endl;
#endif
  MPI_Reduce(&col1, &col, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&error1, &error, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&ncells, &numtot, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&tdiag1[0], &tdiag0[0], numdiag, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
#ifdef DIAG
  double KEtot0=0.0, KElost0=0.0;
  MPI_Reduce(&KEtot,  &KEtot0,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&KElost, &KElost0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (myid==0) {
    cout << "Total lost=" << KElost0/KEtot0 << endl;
  }
#endif

  coltot += col;
  errtot += error;
  for (unsigned k=0; k<numdiag; k++) tdiag[k] += tdiag0[k];

  snglSoFar = snglTime.stop();

  return( col );
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

    // Skip cells with only one particle
    //
    if( number < 2 ) {
      colcntT[id].push_back(0);
      continue;  // Skip to the next cell
    }


				// Compute 1.5 maximum velocity in each MACRO cell
    double crm = 0;
    pCell *samp = c->sample;
    if (samp->state[0]>0.0) {
      for (unsigned k=0; k<3; k++) 
	crm += (samp->state[1+k] - 
		samp->state[4+k]*samp->state[4+k]/samp->state[0])/samp->state[0];
    }
    crm = 1.5*sqrt(fabs(crm));

    // Volume in the cell
    //
    double volc = c->Volume();

    // Fiducial cross section
    //
    diam = diam0;
    double cross  = M_PI*diam*diam;

    // Determine cross section based on fixed number of collisions
    //
    if (CNUM) {
      cross = 2.0*CNUM*volc/(Fn*tau*crm*number*(number-1));
      diam = sqrt(cross/M_PI);
    }

    // Diagnostic: MFP to linear cell size ratio 
    //
    mfpratT[id].push_back(pow(volc, 0.66666667)/(cross*number));
    tsratT[id].push_back(crm/1.5*tau/pow(volc,0.33333333));
    tdensT[id].push_back(number/volc);
    tvolcT[id].push_back(volc);

    // Determine number of candidate collision pairs
    // to be selected in this cell
    //
    double coeff  = 0.5*number*(number-1)*Fn/volc*cross*tau;
    double select = coeff*crm;

				// Diagnose time step in this cell
    double vmass;
    vector<double> V1, V2;
    c->Vel(vmass, V1, V2);
    double vmean = (V2[0]+V2[1]+V2[2])/vmass/3.0;
    double length = pow(c->Volume(), 0.333333);
    double taudiag = length/sqrt(vmean)/tau;
    
    int indx = (int)floor(log(taudiag)/log(10.0) + 4);
    if (indx<0) indx = 0;
    if (indx>8) indx = 8;
    tdiagT[id][indx]++;
    
				// Number of pairs to be selected
    unsigned nsel = (int)floor(select+0.5);

    collTime[id].start();
    initialize_cell(c, crm, tau, select, id);
    collCnt[id]++;

    unsigned colc = 0;
    // Loop over total number of candidate collision pairs
    //
    for (unsigned i=0; i<nsel; i++ ) {

      // Pick two particles at random out of this cell
      //
      unsigned k1 = min<int>((int)floor((*unit)()*number), number-1);
      unsigned k2 = ((int)floor((*unit)()*(number-1)) + k1 + 1) % number;
      Particle* p1 = tree->Body(c->bods[k1]); // First particle
      Particle* p2 = tree->Body(c->bods[k2]); // Second particle

      // Calculate pair's relative speed (pre-collision)
      //
      double cr = 0.0;
      for (int k=0; k<3; k++) {
	crel[k] = p1->vel[k] - p2->vel[k];
	cr += crel[k]*crel[k];
      }
      cr = sqrt(cr);

      if( cr > crm )         // If relative speed larger than crm,
        crm = cr;            // then reset crm to larger value

      // Accept or reject candidate pair according to relative speed
      //
      if( cr/crm > (*unit)() ) {
        // If pair accepted, select post-collision velocities
	//
        colc++;			// Collision counter

				// Do inelastic stuff
	error1T[id] += inelastic(tree, p1, p2, &cr, id);
				// May update relative velocity to reflect
				// excitation of internal degrees of freedom

				// Center of mass velocity
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
        for(unsigned k=0; k<3; k++ ) {
          p1->vel[k] = vcm[k] + p2->mass/tmass*vrel[k];
          p2->vel[k] = vcm[k] - p1->mass/tmass*vrel[k];
        }

	if (CBA) {

	  // Calculate pair's relative speed (post-collision)
	  cr = 0.0;
	  for (int k=0; k<3; k++) {
	    crel[k] = p1->vel[k] - p2->vel[k] - crel[k];
	    cr += crel[k]*crel[k];
	  }
	  cr = sqrt(cr);

	  // Displacement
	  if (cr>0.0) {
	    double displ;
	    for (int k=0; k<3; k++) {
	      displ = crel[k]*diam/cr;
	      if (displ > length) {
		cout << "Huge displacement, process " << myid 
		     << " id=" << id << ": displ=" << displ
		     << " len=" << length << endl;
	      }
	      p1->pos[k] += displ;
	      p2->pos[k] -= displ;
	    }
	  }
	}

      } // Loop over pairs

    }

				// Count collisions
    colcntT[id].push_back(colc);
    col1T[id] += colc;

#ifdef DIAG
    vector<double> ret = diag(); // Diagnostic output for this cell
    KEtotT[id] += ret[0];
    KElostT[id] += ret[1];
#endif
    collSoFar[id] = collTime[id].stop();
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
    return colcnt[colcnt.size()/2]; 

  } else {
    unsigned num = colcnt.size();
    MPI_Send(&num, 1, MPI_UNSIGNED, 0, 39, MPI_COMM_WORLD);
    MPI_Send(&colcnt[0], num, MPI_UNSIGNED, 0, 40, MPI_COMM_WORLD);
    return 0;
  }

}

void Collide::mfpsizeQuantile(vector<double>& quantiles, 
			      vector<double>& mfp_, vector<double>& ts_) 
{
  MPI_Status s;

  if (myid==0) {
    unsigned num;
    for (int n=1; n<numprocs; n++) {
      MPI_Recv(&num, 1, MPI_UNSIGNED, n, 39, MPI_COMM_WORLD, &s);
      vector<double> tmp(num);
      MPI_Recv(&tmp[0], num, MPI_DOUBLE, n, 40, MPI_COMM_WORLD, &s);
      mfprat.insert(mfprat.end(), tmp.begin(), tmp.end());
      MPI_Recv(&tmp[0], num, MPI_DOUBLE, n, 41, MPI_COMM_WORLD, &s);
      tsrat.insert(tsrat.end(), tmp.begin(), tmp.end());
      MPI_Recv(&tmp[0], num, MPI_DOUBLE, n, 42, MPI_COMM_WORLD, &s);
      tdens.insert(tdens.end(), tmp.begin(), tmp.end());
      MPI_Recv(&tmp[0], num, MPI_DOUBLE, n, 43, MPI_COMM_WORLD, &s);
      tvolc.insert(tvolc.end(), tmp.begin(), tmp.end());
    }

    std::sort(mfprat.begin(), mfprat.end()); 
    std::sort(tsrat.begin(),  tsrat.end()); 
    std::sort(tdens.begin(),  tdens.end()); 
    std::sort(tvolc.begin(),  tvolc.end()); 

    mfp_ = vector<double>(quantiles.size());
    ts_  = vector<double>(quantiles.size());
    for (unsigned j=0; j<quantiles.size(); j++) {
      mfp_[j] = mfprat[(unsigned)floor(quantiles[j]*mfprat.size())];
      ts_[j]  = tsrat [(unsigned)floor(quantiles[j]*tsrat.size()) ];
    }

    ofstream out("tmp.collide");
    for (unsigned j=0; j<mfprat.size(); j++)
      out << setw(8) << j 
	  << setw(18) << mfprat[j] 
	  << setw(18) << tsrat[j] 
	  << setw(18) << tdens[j] 
	  << setw(18) << tvolc[j] 
	  << endl;
    
  } else {
    unsigned num = mfprat.size();
    MPI_Send(&num, 1, MPI_UNSIGNED, 0, 39, MPI_COMM_WORLD);
    MPI_Send(&mfprat[0], num, MPI_DOUBLE, 0, 40, MPI_COMM_WORLD);
    MPI_Send(&tsrat[0],  num, MPI_DOUBLE, 0, 41, MPI_COMM_WORLD);
    MPI_Send(&tdens[0],  num, MPI_DOUBLE, 0, 42, MPI_COMM_WORLD);
    MPI_Send(&tvolc[0],  num, MPI_DOUBLE, 0, 43, MPI_COMM_WORLD);
  }
}

