#include <math.h>
#include <sstream>

#include "expand.H"

#define TIMER

#include <UserDSMC.H>


UserDSMC::UserDSMC(string &line) : ExternalForce(line)
{
  id = "DSMCParticleAlgorithm";

  rcloud    = 0.001;	    // Cloud radius
  rmin      = 0.0;	    // Inner radius of cylindrical grid
  rmax      = 1.0;	    // Outer radius of cylindrical grid
  zmax      = 1.0;	    // Half-height of cylindrical grid
  vfrac     = 0.2;	    // Fraction of mean velocity for rejection
  efftol    = 0.5;	    // Minimum permitted inefficiency
  NR        = 64;	    // Number of radial bins
  NZ        = 10;	    // Number of vertical bins
  NP        = 16;	    // Number of azimuthal bins
  seed      = 11;	    // Seed for random number generator
  comp_name = "";	    // Component to apply the DSMC algorithm
  debug     = false;	    // No debugging by default
  
  initialize();

				// Look for the fiducial component for
				// stickiness
  bool found = false;
  for (auto c : comp->components) {
    if ( !comp_name.compare(c->name) ) {
      c0 = c;
      found = true;
      break;
    }
  }

  if (!found) {
    std::ostringstream sout;
    sout << "Can't find desired component <" << comp_name << ">";
    throw GenericError(sout.str(), __FILE__, __LINE__, 35, false);
  }
  
  userinfo();

				// Assign radial bins processors
  dN = (int)floor(NR/numprocs);
  for (int n=0; n<NR; n++) owners[n] = n/dN;

  rbound.resize(NR);
  owners.resize(NR);
  binlist.resize(dN*NZ*NP);

  counts.resize(numprocs);
  countsT.resize(numprocs);
  nbegpts.resize(numprocs);

  vmax.resize(dN);
  nbin.resize(dN);

  eff.resize(numprocs);
  effT.resize(numprocs);
  regrid_flag = true;

  vcom.resize(nthrds);
  vrel.resize(nthrds);
  delv.resize(nthrds);
  for (int n=0; n<nthrds; n++) {
    vcom[n].resize(3);
    vrel[n].resize(3);
    delv[n].resize(3);
  }

  gen = new ACG(seed+myid);
  unif = new Uniform(0.0, 1.0, gen);
#ifdef TIMER
  timer = new Timer;
#endif

}

UserDSMC::~UserDSMC()
{
  delete unif;
  delete gen;
#ifdef TIMER
  delete timer;
#endif
}

void UserDSMC::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine: DSMC with Rcloud=" << rcloud
       << ", Rmin=" << rmin << ", Rmax=" << rmax << ", Zmax=" << zmax
       << ", NR, NZ, NP=(" << NR << ", " << NZ << ", " << NP << ") "
       << ", seed=" << seed
       << ", efftol=" << efftol
       << ", applying to particles from component <" << comp_name << ">";

  if (debug)
    cout << ", debugging is *ON*";

  cout << endl;

  print_divider();
}

void UserDSMC::initialize()
{
  string val;

  if (get_value("compname", val))	comp_name = val;
  if (get_value("Rcloud", val))		rcloud = atof(val.c_str());
  if (get_value("Rmin", val))		rmin = atof(val.c_str());
  if (get_value("Rmax", val))		rmax = atof(val.c_str());
  if (get_value("Zmax", val))		zmax = atof(val.c_str());
  if (get_value("Vfrac", val))		vfrac = atof(val.c_str());
  if (get_value("efftol", val))		efftol = atof(val.c_str());
  if (get_value("NR", val))		NR = atoi(val.c_str());
  if (get_value("NZ", val))		NZ = atoi(val.c_str());
  if (get_value("NP", val))		NP = atoi(val.c_str());
  if (get_value("seed", val))		seed = atoi(val.c_str());
  if (get_value("debug", val))		debug = atol(val);
}


void UserDSMC::makeSort()
{
  // Exchange particles

  // Serialize redistribution list into an integer array
  // of numprocs stanzas, each with the format
  // n -- current node
  // M -- number to redistribute (may be zero)
  // index_1
  // tonode_1
  // index_2
  // tonode_2
  // .
  // .
  // index_M
  // tonode_M
  //
  // so each stanza has 2(M+1) integers

  for (int n=0; n<numprocs; n++) counts[n] = countsT[n] = 0;
  redistT.clear();
  redistT.push_back(myid);	// Index 0
  redistT.push_back(0);		// Index 1

				// Boundaries
				// 
  Particle *p;
  double R, x, y, z, phi;

  for (unsigned int k=0; k<c0->Number(); k++) {
    
    p = c0->Part(k);

    R = sqrt(p->pos[0]*p->pos[0] + p->pos[1]*p->pos[1]);
      
    iI = upper_bound(rboundI.begin(), rboundI.end(), R, BoundSort());
    if (iI == rboundI.end()) continue;

    if (owners[iI->second] != myid) {
      redistT.push_back(k);
      redistT.push_back(owners[iI->second]);
      redistT[1]++;
    }
  }

  countsT[myid] = redistT.size();
  MPI_Allreduce(&countsT[0], &counts[0], numprocs, MPI_UNSIGNED, MPI_SUM, 
		MPI_COMM_WORLD);
  int ntot = 0;
  for (int n=0; n<numprocs; n++) ntot += counts[n];
  redistT2.clear();
  redistT2.resize(ntot, 0);
  
  int ibeg;
  if (myid==0) ibeg = 0;
  else ibeg = counts[myid-1];
  for (unsigned n=0; n<redistT.size(); n++) redistT2[ibeg+n] = redistT[n];
  MPI_Allreduce(&redistT2[0], &redist[0], ntot, MPI_INT, MPI_SUM, 
		MPI_COMM_WORLD);

				// Tell the component to redistribute
				// 
  c0->redistributeByList(redist);

				// Determine efficiency
				// 
  for (int n=0; n<numprocs; n++) effT[n] = 0;
  effT[myid] = c0->Number();
  MPI_Allreduce(&effT[0], &eff[0], numprocs, MPI_UNSIGNED, MPI_SUM, 
		MPI_COMM_WORLD);

  unsigned int maxT=0, minT=std::numeric_limits<int>::max();
  for (int n=0; n<numprocs; n++) {
    minT = min<unsigned int>(minT, eff[n]);
    maxT = max<unsigned int>(maxT, eff[n]);
  }
  if ((double)minT/maxT < efftol) regrid_flag = true;


				// Clean the bin list
				// 
  for (auto i : binlist) i.clear();

  for (int i=0; i<dN; i++) {
    vmax[i] = 0.0;
    nbin[i] = 0;
  }
				// Make a new bin list
				// 
  int ir, iz, ip;
  unsigned int indx;
  double u, v, w, V;
  for (unsigned int k=0; k<c0->Number(); k++) {
    x = c0->Pos(k, 0);
    y = c0->Pos(k, 1);
    z = c0->Pos(k, 2);

    u = c0->Vel(k, 0);
    v = c0->Vel(k, 1);
    w = c0->Vel(k, 2);

    R = sqrt(x*x + y*y);

    if (R< rmin || R>rmax) continue;
    if (z<-zmax || z>zmax) continue;
    
    phi = atan2(y, x);

    V = sqrt( (u*u+v*v*w*w)/3.0 );

    iI = upper_bound(rboundI.begin(), rboundI.end(), R, BoundSort());
    if (iI == rboundI.end()) continue;
    ir = iI->second;

    if (owners[ir] != myid) continue;

    iI = upper_bound(zboundI.begin(), zboundI.end(), R, BoundSort());
    if (iI == zboundI.end()) continue;
    iz = iI->second;

    ip = (int)floor(phi/(2.0*M_PI)*NP);

    indx = ( (ir-myid*dN)*NZ + iz )*NP + ip;

    binlist[indx].push_back(k);

    vmax[ir - myid*dN] += V;
    nbin[ir - myid*dN] += 1;
  }

				// Compute mean velocity per bin
				// 
  for (int i=0; i<dN; i++) {
    if (nbin[i]>0) vmax[i] = vmax[i]*vfrac/nbin[i];
  }

}

void UserDSMC::makeGrid()
{
  if (!regrid_flag) return;

  MPI_Status  status;

  // Load and send particle numbers to master
  //
  vector<int> countT(numprocs, 0);
  countT[myid] = c0->Number();
  MPI_Reduce(&countT[0], &counts[0], numprocs, MPI_INT, MPI_SUM, 0,
	     MPI_COMM_WORLD);


  // Resize the particle buffers
  // 
  if (myid==0) {
    if (radii.size() != c0->nbodies_tot) radii.resize(c0->nbodies_tot);
    int nb = 0;
    for (int n=1; n<numprocs; n++) nb = max<int>(nb, counts[n]);
    if ((int)radiiT.size() == nb) radiiT.resize(nb);
  } else {
    if (radiiT.size() != c0->Number()) radiiT.resize(c0->Number());
  }



  // Put my particles into the buffer and the master list
  //
  nbegpts[0] = 0;
  for (int n=1; n<numprocs; n++) nbegpts[n] = counts[n-1];


  Particle *p;

  if (myid==0) {
    for (unsigned int k=0; k<c0->Number(); k++) {
      p = c0->Part(k);
      radii[k].curnode = 0;
      radii[k].curindx = k;
      radii[k].radius  = sqrt(p->pos[0]*p->pos[0] + p->pos[1]*p->pos[1]);
    }	
  } else {
    for (unsigned int k=0; k<c0->Number(); k++) {
      p = c0->Part(k);
      radiiT[k] = sqrt(p->pos[0]*p->pos[0] + p->pos[1]*p->pos[1]);
    }	
  }

  
  // Root gets the radial list
  //
  if (myid==0) {
    for (int n=1; n<numprocs; n++) {
      MPI_Recv(&radiiT[0], counts[n], MPI_FLOAT, n, 21+n, 
	       MPI_COMM_WORLD, &status);
      
      for (unsigned int k=0; k<counts[n]; k++) {
	radii[nbegpts[n]+k].curnode = n;
	radii[nbegpts[n]+k].curindx = k;
	radii[nbegpts[n]+k].radius = radiiT[k];
      }
    }

				// Check for maximum radius of cylindrical grid
				// 
    sort(radii.begin(), radii.end(), RadiusSort());
    if (rmax < radii.back().radius) {
      vector<RadiusSort>::iterator ir = 
	upper_bound(radii.begin(), radii.end(), rmax, RadiusSort());
      radii.erase(ir, radii.end());
    }

				// Check for minimum radius of cylindrical grid
				// 
    if (rmin > radii.front().radius) {
      vector<RadiusSort>::iterator ir = 
	lower_bound(radii.begin(), radii.end(), rmin, RadiusSort());
      radii.erase(radii.begin(), ir);
    }

    // Serialize redistribution list into an integer array
    // of numprocs stanzas, each with the format
    // n -- current node
    // M -- number to redistribute (may be zero)
    // index_1
    // tonode_1
    // index_2
    // tonode_2
    // .
    // .
    // index_M
    // tonode_M
    //
    // so each stanza has 2(M+1) integers

    redist.clear();
				// Boundaries
				// 
    int ibeg, iend, cntr;
    for (int n=0; n<NR; n++) {

      ibeg = (int)floor((double)(n+0)/NR*radii.size());
      iend = (int)floor((double)(n+1)/NR*radii.size());

      if (n) {
	if (ibeg+1<(int)radii.size())
	  rbound[n] = 0.5*(radii[ibeg].radius + radii[ibeg+1].radius);
	else 
	  rbound[n] = radii[ibeg].radius;
      }
      else
	rbound[n] = rmin;
				// Iterate through particles making a 
				// redistribution list
      redist.push_back(owners[n]);
      redist.push_back(0);
      cntr = redist.size()-1;

      for (int k=ibeg; k<iend; k++) {
	if (radii[k].curnode != owners[n]) {
	  redist.push_back(radii[k].curindx);
	  redist.push_back(owners[n]);
	  redist[cntr]++;
	}
      }
    }

  } else {
    int nbodies = c0->Number();
    MPI_Send(&radiiT[0], nbodies, MPI_FLOAT, 0, 21+myid, MPI_COMM_WORLD);
  }
  

				// Send boundaries to all processors
				// 
  MPI_Bcast(&rbound[0], NR, MPI_FLOAT, 0, MPI_COMM_WORLD);
  rboundI.clear();
  for (int n=0; n<NR; n++) rboundI.push_back( pair<double, int>(rbound[n], n) );

				// Send redistribution array count 
				// to all processors
  unsigned int nredist;
  if (myid==0) nredist= redist.size();
  MPI_Bcast(&nredist, numprocs, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  if (myid) redist.resize(nredist);
  MPI_Bcast(&redist, nredist, MPI_INT, 0, MPI_COMM_WORLD);

  
				// Tell the component to redistribute
				// 
  c0->redistributeByList(redist);

				// Now, finally, make the grid for each
				// node
  dPhi = 2.0*M_PI/NP;
  vector<double> vertical;
  for (unsigned int k=0; k<c0->Number(); k++) vertical.push_back(c0->Pos(k, 2));
  sort(vertical.begin(), vertical.end());
  
				// Check for min & max height
				// 
  double zbot;
  if (zmax < vertical.back()) {
    vector<double>::iterator iz = 
      upper_bound(vertical.begin(), vertical.end(), zmax);
    vertical.erase(iz, vertical.end());
  }
  if (-zmax > vertical.front()) {
    vector<double>::iterator iz = 
      lower_bound(vertical.begin(), vertical.end(), -zmax);
    vertical.erase(vertical.begin(), iz);
    zbot = -zmax;
  } else {
    zbot = vertical.front();
  }

  zboundI.clear();
  zboundI.push_back( pair<double, int>(zbot, 0) );
  for (int n=1; n<NZ; n++) {
    unsigned int iz = (int)floor((double)n/NZ*vertical.size());
    if (iz==vertical.size()-1)
      zbot = vertical[iz];
    else
      zbot = 0.5*(vertical[iz] + vertical[iz+1]);
    
    zboundI.push_back( pair<double, int>(zbot, n) );
  }

  regrid_flag = false;
}

void UserDSMC::determine_acceleration_and_potential(void)
{
#if HAVE_LIBCUDA==1		// Cuda compatibility
  getParticlesCuda(c0);
#endif

#ifdef TIMER
  if (myid==0 && debug) {
    timer->reset();
    timer->start();
  }
#endif
//----------------------------------------------------------------------
  makeGrid();
//----------------------------------------------------------------------
#ifdef TIMER
  if (debug) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) {
      std::cout << "Grid constructed in " << timer->stop() << " seconds" << std::endl;

      timer->reset();
      timer->start();
    }
  }
#endif
//----------------------------------------------------------------------
  makeSort();
//----------------------------------------------------------------------
#ifdef TIMER
  if (debug) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) {
      std::cout << "Sort completed in " << timer->stop() << " seconds" << std::endl;
      timer->reset();
      timer->start();
    }
  }
#endif
//----------------------------------------------------------------------
  exp_thread_fork(false);
//----------------------------------------------------------------------
#ifdef TIMER
  if (debug) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) {
      std::cout << "Collisions completed in " << timer->stop() << " seconds" << std::endl;
    }
  }
#endif
}


void * UserDSMC::determine_acceleration_and_potential_thread(void * arg) 
{
  unsigned int number = binlist.size();
  int id = *((int*)arg);
  int nbeg = number*id/nthrds;
  int nend = number*(id+1)/nthrds;

  double Vc, rho, Rmin, Rmax, Zmin, Zmax;
  int ir, iz, iw, ncount, Nc;

  PartMapItr it = cC->Particles().begin();
  unsigned long i;

  for (int q=0   ; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {
    i = (it++)->first;
				// Number in this cell
				// 
    Nc = binlist[i].size();
    if (Nc<2) continue;		// Can't have any collisions, next cell

				// Determine radial boundaries
				// 
    iw = i/(NZ*NP);
    ir = iw + myid*dN;
    Rmin = rboundI[ir].first;
    if (ir+1<NR)
      Rmax = rboundI[ir+1].first;
    else
      Rmax = rmax;
				// Determine vertical boundaries
				// 
    iz = (i - iw)/NP;
    Zmin = zboundI[iz].first;
    if (iz+1<NZ)
      Zmax = zboundI[iz+1].first;
    else
      Zmax = zmax;

				// Volume and density in the bin
				// 
    Vc = dPhi*(Zmin-Zmax)*0.5*(Rmax*Rmax - Rmin*Rmin);
    rho = (double)Nc/Vc;

				// Number of trial points
				// 
    ncount = (int)floor(0.5*rho*rho*M_PI*rcloud*rcloud*vmax[iw]*Vc*dtime + 0.5);

    double vmaxT = vmax[iw];

    for (int nc=0; nc<ncount; nc++) {
				// Pick two particles at random out of this cell
				// 
      int n1 = (int)floor(unif(random_gen)* Nc);
      int n2 = (1 + (int)floor(unif(random_gen)*(Nc-1)+n1) ) % Nc;

      Particle *p1 = c0->Part(binlist[i][n1]);
      Particle *p2 = c0->Part(binlist[i][n2]);

				// Calculate pair's relative speed
				// 
      double rvel=0.0;
      for (int k=0; k<3; k++) {
	vrel[id][k] = p1->vel[k] - p2->vel[k];
	rvel += vrel[id][k]*vrel[id][k];
      }
      rvel = sqrt(rvel);

				// If relative speed larger than vmax
				// then reset vmax to larger value
      if (rvel > vmaxT) vmaxT = rvel;		

				// Accept or reject candidate pair 
				// according to relative speed
      if (rvel/vmax[iw] > unif(random_gen) ) {

				// Compute post-collision velocities
	
				// Center of mass velocity
        for(int k=0; k<3; k++ )
          vcom[id][k] = 0.5*(p1->vel[k] + p2->vel[k]);

				// Random orientation

				// Cosine and sine of collision angle theta
				// 
        double cos_th = 1.0 - 2.0*unif(random_gen); 
        double sin_th = sqrt(1.0 - cos_th*cos_th);

				// Collision angle phi
				// 
        double phi = 2.0*M_PI*unif(random_gen);

				// Compute hard-sphere displacement following
				// Alexander, Garcia & Alder (1995, PRL, 74, 5212)
				// 
	delv[id][0] = rvel*cos_th          - vrel[id][0];
	delv[id][1] = rvel*sin_th*cos(phi) - vrel[id][1];
	delv[id][2] = rvel*sin_th*sin(phi) - vrel[id][2];

	double dvel = 0.0;
	for (int k=0; k<3; k++) dvel += delv[id][k]*delv[id][k];
	dvel = sqrt(dvel);
	if (dvel>0.0) {
	  for (int k=0; k<3; k++) {
	    p1->pos[k] += delv[id][k]*rcloud/dvel;
	    p2->pos[k] -= delv[id][k]*rcloud/dvel;
	  }
	}
				// Compute post-collision relative velocity
				// 
        vrel[id][0] = rvel*cos_th;
        vrel[id][1] = rvel*sin_th*cos(phi);
        vrel[id][2] = rvel*sin_th*sin(phi);
        for (int  k=0; k<3; k++ ) {
          p1->vel[k] = vcom[id][k] + 0.5*vrel[id][k];  // Update post-collision
          p2->vel[k] = vcom[id][k] - 0.5*vrel[id][k];  // velocities
        }

      }

      vmax[iw] = vmaxT;     // Update max relative speed

    }
    
  }

  return (NULL);
}


extern "C" {
  ExternalForce *makerDSMC(string& line)
  {
    return new UserDSMC(line);
  }
}

class proxysticky { 
public:
  proxysticky()
  {
    factory["usersticky"] = makerDSMC;
  }
};

static proxysticky p;
