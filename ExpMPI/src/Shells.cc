#include "expand.h"
#include "interp.h"

#include <Shells.H>

Shells::Shells(string& line) : PotAccel(line)
{
  nsample  = -1;

  initialize();

  rgridT = new vector<double> [nthrds];
  mgridT = new vector<double> [nthrds];
}

Shells::~Shells()
{
  delete [] rgridT;
  delete [] mgridT;
}

void Shells::initialize(void)
{
  string val;

  if (get_value("nsample", val))	nsample = atoi(val.c_str());
}

void Shells::create_grid(void)
{
  pnumber = vector<int> (numprocs);
}

void Shells::get_acceleration_and_potential(Component* C)
{
  cC = C;
  nbodies = cC->Number();

  /*======================================*/
  /* Determine potential and acceleration */
  /*======================================*/

  determine_acceleration_and_potential();
}

void Shells::determine_acceleration_and_potential(void)
{
  exp_thread_fork(false);
				// Clear external potential flag
  use_external = false;
}

void * Shells::determine_acceleration_and_potential_thread(void * arg)
{
  double rr, rfac;
  double mass, potl;
  int indx;

  unsigned nbodies = cC->levlist[mlevel].size();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  unsigned long j;		// Index of the current local particle

  for (int i=nbeg; i<nend; i++) {
    
    j = cC->levlist[mlevel][i];

				// Don't need acceleration for frozen particles
    if (cC->freeze(j)) continue;
    
    
				// Compute radius
    rr = 0.0;
    for (int k=0; k<3; k++) rr += cC->Pos(j, k) * cC->Pos(j, k);
      
    rr = sqrt(rr);

				// Compute index
    if (rr<=0.0) {
      mass = 0;
      potl = pgrid.front();
    }
    else if (rr>=rgrid.back()) {
      mass = mgrid.back();
      potl = 0.0;
    }
    else {
      indx = Vlocate(rr, rgrid);
      mass = (
	      mgrid[indx  ]*(rgrid[indx+1] - rr) +
	      mgrid[indx+1]*(rr - rgrid[indx  ])
	      ) / (rgrid[indx+1] - rgrid[indx]);
      potl = (
	      pgrid[indx  ]*(rgrid[indx+1] - rr) +
	      pgrid[indx+1]*(rr - rgrid[indx  ])
	      ) / (rgrid[indx+1] - rgrid[indx]);
    }
    
    
    // Acceleration
    rfac = mass/(rr*rr*rr);
	
    for (int k=0; k<3; k++)
      cC->AddAcc(j, k, -cC->Pos(j, k) * rfac );
      
				// Potential
    if (use_external)
      cC->AddPotExt(j, potl-mass/rr);

    else if (rr > 1.0e-16)
      cC->AddPot(j, potl-mass/rr );
  }
}

void Shells::determine_coefficients(void) 
{
  create_grid();

  for (int i=0; i<nthrds; i++) {
    rgridT[i].erase(rgridT[i].begin(), rgridT[i].end());
    mgridT[i].erase(mgridT[i].begin(), mgridT[i].end());
  }

  exp_thread_fork(true);

  if (nsample<=1) {		// Use all the particles
    rgrid1 = rgridT[0];
    mgrid1 = mgridT[0];
    for (int i=1; i<nthrds; i++) {
      rgrid1.insert(rgrid1.end(), rgridT[i].begin(), rgridT[i].end());
      mgrid1.insert(mgrid1.end(), mgridT[i].begin(), mgridT[i].end());
    }
  } else {			// Subsample the particles
    rgrid1.erase(rgrid1.begin(), rgrid1.end());
    mgrid1.erase(mgrid1.begin(), mgrid1.end());

    for (int i=0; i<nthrds; i++) {
      for (unsigned j=0; j<rgridT[i].size(); j+=nsample) {
	rgrid1.push_back(rgridT[i][j]);
	mgrid1.push_back(mgridT[i][j]);
      }
    }
  }

  int pn = rgrid1.size();

  MPI_Gather(&pn,                1, MPI_INT, 
	     &pnumber[0], numprocs, MPI_INT, 
	     0, MPI_COMM_WORLD);
  
  int nsum = 0;
  vector<int> displ(numprocs);
  if (myid==0) {
    for (int n=0; n<numprocs; n++) {
      displ[n] = nsum;
      nsum += pnumber[n];
    }
  }
  MPI_Bcast(&nsum, 1, MPI_INT, 0, MPI_COMM_WORLD);


  rgrid = vector<double>(nsum);
  mgrid = vector<double>(nsum);

  MPI_Gatherv(&rgrid1[0], pn, MPI_DOUBLE, &rgrid[0], &pnumber[0], &displ[0],
	      MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Gatherv(&mgrid1[0], pn, MPI_DOUBLE, &mgrid[0], &pnumber[0], &displ[0],
	      MPI_DOUBLE, 0, MPI_COMM_WORLD);

  grid.erase(grid.begin(), grid.end());
  for (int i=0; i<nsum; i++) 
    grid.push_back(pair<double, double>(rgrid[i], mgrid[i]));
  
  sort(grid.begin(), grid.end());

  model = vector<double>(3*(nsum+1));

  if (myid==0) {
    double rL=0.0, rC, mL=0.0, mC;
    vector<double> vval(3, 0.0);

    for (int i=0; i<nsum; i++) {
      rC = grid[i].first;
      mC = grid[i].second;

      vval[0] = rC;
      vval[1] += 0.5*(mL + mC)*(rC - rL);
      if (rL > 0.0)
	vval[2] += 0.5*(mL/rL + mC/rC)*(rC - rL);
      else
	vval[2] += 0.5*(0.0   + mC/rC)*(rC - rL);
      
      for (int k=0; k<3; k++) model[3*i+k] = vval[k];
      rL = rC;
      mL = mC;
    }

    double potlF = model[3*(nsum+1)-1];
    for (int i=0; i<=nsum; i++) model[3*i+2] -= potlF;
  }

  MPI_Bcast(&model[0], 3*(nsum+1), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  rgrid = vector<double>(nsum+1);
  mgrid = vector<double>(nsum+1);
  pgrid = vector<double>(nsum+1);

  for (int i=0; i<=nsum; i++) {
    rgrid[i] = model[3*i+0];
    mgrid[i] = model[3*i+1];
    pgrid[i] = model[3*i+2];
  }

}

void * Shells::determine_coefficients_thread(void *arg) 
{
  double rr;

  unsigned nbodies = cC->levlist[mlevel].size();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  unsigned long j;		// Index of the current local particle

  for (int i=nbeg; i<nend; i++) {
    
    j = cC->levlist[mlevel][i];

				// Don't need acceleration for frozen particles
    if (cC->freeze(j)) continue;
    
				// Compute radius
    rr = 0.0;
    for (int k=0; k<3; k++) rr += cC->Pos(j, k) * cC->Pos(j, k);
      
				// Load vectors
    rgridT[id].push_back(sqrt(rr));
    mgridT[id].push_back(cC->Mass(j));
  }
}
