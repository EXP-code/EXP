#include <Shells.H>

Shells::Shells(string& line) : PotAccel(line)
{
  nsample         = -1;
  nselect         = -1;
  self_consistent = true;
  firstime_accel  = true;

  initialize();

				// For accumulation
  rgridT = new vector<double> [nthrds];
  mgridT = new vector<double> [nthrds];
  igridT = new vector<int   > [nthrds];
  usedT  = vector<int> (nthrds);

				// For Alltoallv calls
  sdispl  = vector<int> (numprocs, 0);
  sndcts  = vector<int> (numprocs, 1);
  rcvcts  = vector<int> (numprocs);
  for (int n=0; n<numprocs; n++) rcvcts[n] = n;
  rdispl  = vector<int> (numprocs);
  snumbr  = vector<int> (numprocs);
  rnumbr  = vector<int> (numprocs);

				// For storage of samples at each level
  rgrid = new map<int, double> [multistep+1];
  mgrid = new map<int, double> [multistep+1];

  update_fr = new vector<int> [nthrds];
  update_to = new vector<int> [nthrds];
  update_ii = new vector<int> [nthrds];
}

Shells::~Shells()
{
  delete [] rgridT;
  delete [] mgridT;
  delete [] igridT;

  delete [] rgrid;
  delete [] mgrid;
}


void Shells::initialize(void)
{
  string val;

  if (get_value("nsample", val)) 
    nsample = atoi(val.c_str());

  if (get_value("nselect", val)) 
    nselect = atoi(val.c_str());

  if (get_value("self_consistent", val)) {
    if (atoi(val.c_str())) self_consistent = true; 
    else self_consistent = false;
  } else
    self_consistent = true;

}

void Shells::get_acceleration_and_potential(Component* C)
{
  cC = C;
  nbodies = cC->Number();

  if (!use_external) {

    if (firstime_accel || self_consistent || initializing) 
      determine_coefficients();
    
    firstime_accel = false;
  }

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
  int id = *((int*)arg);
  double rr, rfac;
  double mass, potl, dr;

  int indx;			// Interpolation index
  unsigned long j;		// Index of the current local particle

  unsigned nbodies;
  int nbeg, nend;


  for (int lev=mlevel; lev<=multistep; lev++) {

    nbodies = cC->levlist[lev].size();

    if (nbodies==0) continue;

    nbeg = nbodies*id/nthrds;
    nend = nbodies*(id+1)/nthrds;

    for (int i=nbeg; i<nend; i++) {
    
      j = cC->levlist[lev][i];

				// Don't need acceleration for frozen particles
      if (cC->freeze(j)) continue;
    
    
				// Compute radius
      rr = sqrt(
		cC->Pos(j, 0)*cC->Pos(j, 0) +
		cC->Pos(j, 1)*cC->Pos(j, 1) +
		cC->Pos(j, 2)*cC->Pos(j, 2) 
		);

				// Compute index
      if (rr<=0.0) {
	mass = 0;
	potl = pgrid0.front();
      }
      else if (rr>=rgrid0.back()) {
	mass = mgrid0.back();
	potl = 0.0;
    }
      else {
	indx = Vlocate(rr, rgrid0);
	dr   = rgrid0[indx+1] - rgrid0[indx];
	
	mass = (mgrid0[indx  ]*(rgrid0[indx+1] - rr) +
		mgrid0[indx+1]*(rr - rgrid0[indx  ]) ) / dr;
	
	potl = (pgrid0[indx  ]*(rgrid0[indx+1] - rr) +
		pgrid0[indx+1]*(rr - rgrid0[indx  ]) ) / dr;
      }

				// Accelerationn
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

  return (NULL);
}

void Shells::determine_coefficients(void) 
{
				// For error info
  char msg[MPI_MAX_ERROR_STRING];
  int msglen;
				// Clear the data arrays
  for (int i=0; i<nthrds; i++) {
    rgridT[i].erase(rgridT[i].begin(), rgridT[i].end());
    mgridT[i].erase(mgridT[i].begin(), mgridT[i].end());
    igridT[i].erase(igridT[i].begin(), igridT[i].end());
    usedT[i] = 0;
  }
				// Make the radius--mass lists
  exp_thread_fork(true);
  
  used1  = usedT[0];
  rgrid1 = rgridT[0];
  mgrid1 = mgridT[0];
  igrid1 = igridT[0];
  for (int i=1; i<nthrds; i++) {
    rgrid1.insert(rgrid1.end(), rgridT[i].begin(), rgridT[i].end());
    mgrid1.insert(mgrid1.end(), mgridT[i].begin(), mgridT[i].end());
    igrid1.insert(igrid1.end(), igridT[i].begin(), igridT[i].end());
    used1 += usedT[i];
  }


  MPI_Allreduce(&used1, &used, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if (used) {

    int pn = rgrid1.size();

    if (int ret=
	MPI_Alltoallv(&pn, &sndcts[0], &sdispl[0], MPI_INT, 
		      &rnumbr[0], &sndcts[0], &rcvcts[0], MPI_INT, 
		      MPI_COMM_WORLD)
	) 
      {
	MPI_Error_string(ret, msg, &msglen);
	cout << "Shells::determine_coefficients: " << msg << endl;
      }
  
    int nsum = 0;
    for (int n=0; n<numprocs; n++) {
      snumbr[n] = pn;
      rdispl[n] = nsum;
      nsum += rnumbr[n];
    }
				// Share the lists will all processes

    vector<double> rgridX(nsum);
    vector<double> mgridX(nsum);
    vector<int   > igridX(nsum);

    rgrid[mlevel].erase(rgrid[mlevel].begin(), rgrid[mlevel].end());
    mgrid[mlevel].erase(mgrid[mlevel].begin(), mgrid[mlevel].end());

    if (int ret=
	MPI_Alltoallv(&rgrid1[0], &snumbr[0], &sdispl[0], MPI_DOUBLE, 
		      &rgridX[0], &rnumbr[0], &rdispl[0], MPI_DOUBLE,
		      MPI_COMM_WORLD)
	)
      {
	MPI_Error_string(ret, msg, &msglen);
	cout << "Shells::determine_coefficients: " << msg << endl;
      }

    
    if (int ret=
	MPI_Alltoallv(&mgrid1[0], &snumbr[0], &sdispl[0], MPI_DOUBLE, 
		      &mgridX[0], &rnumbr[0], &rdispl[0], MPI_DOUBLE,
		      MPI_COMM_WORLD)
	)
      {
	MPI_Error_string(ret, msg, &msglen);
	cout << "Shells::determine_coefficients: " << msg << endl;
      }

    if (int ret=
	MPI_Alltoallv(&igrid1[0], &snumbr[0], &sdispl[0], MPI_INT, 
		      &igridX[0], &rnumbr[0], &rdispl[0], MPI_INT,
		      MPI_COMM_WORLD)
	)
      {
	MPI_Error_string(ret, msg, &msglen);
	cout << "Shells::determine_coefficients: " << msg << endl;
      }

    for (int i=0; i<nsum; i++) {
      rgrid[mlevel][igridX[i]] = rgridX[i];
      mgrid[mlevel][igridX[i]] = mgridX[i];
    }
  }

				// Make and sort the list

  double mfac = 1.0;
  if (nsample>1) mfac *= nsample;
  grid.erase(grid.begin(), grid.end());
  map<int, double>::iterator ri, mi;
  for (int m=0; m<=multistep; m++) {
    ri = rgrid[m].begin();
    mi = mgrid[m].begin();
    while (ri != rgrid[m].end() && mi != mgrid[m].end())
      grid.push_back(Dpair((ri++)->second, (mi++)->second*mfac));
  }
  
  if (grid.size()) {

    sort(grid.begin(), grid.end());

    int ntot = grid.size();
  
    rgrid0.erase(rgrid0.begin(), rgrid0.end());
    mgrid0.erase(mgrid0.begin(), mgrid0.end());
    pgrid0.erase(pgrid0.begin(), pgrid0.end());

    double rL=0.0, rC, mL=0.0, mC;
    vector<double> vval(3, 0.0);

    for (int i=0; i<ntot; i++) {
      rC = grid[i].first;
      mC = grid[i].second;

      vval[0] = rC;
      vval[1] += 0.5*(mL + mC);
      if (rL > 0.0)
	vval[2] += 0.5*(mL/rL + mC/rC);
      else
	vval[2] += 0.5*(0.0   + mC/rC);
      
      if ( (i % nselect) == 0) {
	rgrid0.push_back(vval[0]);
	mgrid0.push_back(vval[1]);
	pgrid0.push_back(vval[2]);
      }

      rL = rC;
      mL = mC;
    }

    double potlF = pgrid0.back();
    for (unsigned i=0; i<pgrid0.size(); i++) pgrid0[i] -= potlF;

  }  else {

    rgrid0.push_back(0.0);
    mgrid0.push_back(0.0);
    pgrid0.push_back(0.0);

  }
}

void * Shells::determine_coefficients_thread(void *arg) 
{
  double rr;

  unsigned nbodies = cC->levlist[mlevel].size();
  int id   = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  unsigned long j;		// Index of the current local particle

  for (int i=nbeg; i<nend; i++) {
    
    j = cC->levlist[mlevel][i];

				// Don't need acceleration for frozen particles
    if (cC->freeze(j)) continue;
    
    if (nsample>1 && (i % nsample)) continue;

				// Compute radius
    rr = sqrt(
	      cC->Pos(j, 0)*cC->Pos(j, 0) +
	      cC->Pos(j, 1)*cC->Pos(j, 1) +
	      cC->Pos(j, 2)*cC->Pos(j, 2) 
	      );
				// Load vectors
    rgridT[id].push_back(rr);
    mgridT[id].push_back(cC->Mass(j));
    igridT[id].push_back(j);
    usedT[id]++;
  }
}

//
// Clear the storage lists for the update loop to start
//
void Shells::multistep_update_begin()
{
  for (int n=0; n<nthrds; n++) {
    update_fr[n].erase(update_fr[n].begin(), update_fr[n].end());
    update_to[n].erase(update_to[n].begin(), update_to[n].end());
    update_ii[n].erase(update_ii[n].begin(), update_ii[n].end());
  }
}

void Shells::multistep_update_finish()
{
  //
  // Accumulate the change requests from each thread
  //
  for (int n=1; n<nthrds; n++) {
    update_fr[0].insert(update_fr[0].end(),
			update_fr[n].begin(), update_fr[n].end());

    update_to[0].insert(update_to[0].end(),
			update_to[n].begin(), update_to[n].end());

    update_ii[0].insert(update_ii[0].end(),
			update_ii[n].begin(), update_ii[n].end());
  }
  
  //
  // Number of changes for this process
  //
  int pn = update_fr[0].size();

  //
  // Communicate the changes per process to all processes
  //
  MPI_Alltoallv(&pn, &sndcts[0], &sdispl[0], MPI_INT, 
		&rnumbr[0], &sndcts[0], &rcvcts[0], MPI_INT, 
		MPI_COMM_WORLD);

  int nsum = 0;
  for (int n=0; n<numprocs; n++) {
    snumbr[n] = pn;
    rdispl[n] = nsum;
    nsum += rnumbr[n];
  }


  vector<int> fr(nsum);
  vector<int> to(nsum);
  vector<int> ii(nsum);

  //
  // Share exchange lists will all processes
  //
  MPI_Alltoallv(&update_fr[0][0], &snumbr[0], &sdispl[0], MPI_INT, 
		&fr[0],           &rnumbr[0], &rdispl[0], MPI_INT,
		MPI_COMM_WORLD);

  MPI_Alltoallv(&update_to[0][0], &snumbr[0], &sdispl[0], MPI_INT, 
		&to[0],           &rnumbr[0], &rdispl[0], MPI_INT,
		MPI_COMM_WORLD);

  MPI_Alltoallv(&update_ii[0][0], &snumbr[0], &sdispl[0], MPI_INT, 
		&ii[0],           &rnumbr[0], &rdispl[0], MPI_INT,
		MPI_COMM_WORLD);

  //
  // Each process performs the update
  //
  map<int, double>::iterator rr, mm;

  for (int i=0; i<nsum; i++) {

    rr = rgrid[fr[i]].find(ii[i]);
    mm = mgrid[fr[i]].find(ii[i]);

    rgrid[to[i]][ii[i]] = rr->second;
    mgrid[to[i]][ii[i]] = mm->second;

    rgrid[fr[i]].erase(rr);
    mgrid[fr[i]].erase(mm);
  }

}


//
// Save the update list for each thread to be processed at the end of
// the update loop
//
void Shells::multistep_update(int fr, int to, Component *c, int ii, int id)
{
  update_fr[id].push_back(fr);
  update_to[id].push_back(to);
  update_ii[id].push_back(ii);
}

