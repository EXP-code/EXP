#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <random>
#include <algorithm>

#include "expand.H"


#include <ScatterMFP.H>

bool less_rpair(const rpair& one, const rpair& two)
{
  return (one.first < two.first);
}


ScatterMFP::ScatterMFP(const YAML::Node& conf) : ExternalForce(conf)
{
  tautab = 100;
  tauscat = 1.0;
  rmax = 100.0;
  nscat = 20;
  mfp_index = 0;

  initialize();
  
				// Look for requested id in component list
  c = NULL;
  for (auto it : comp->components) {
    if (it->id.compare(comp_id) == 0) {
      c = it;
      break;
    }
  }

				// Make sure component has been found

  if (c==NULL) {
    std::ostringstream sout;
    sout << "ScatterMFP: can not find target component <" << comp_id << ">\n";
    throw GenericError(sout.str(), __FILE__, __LINE__);
  }

				// Check for mfp in particle attribute list

  if (c->ndattrib < mfp_index+1) {
    c->ndattrib = mfp_index+1;
    PartMapItr it;
    for (it=c->particles.begin(); it!=c->particles.end(); it++)
      it->second->dattrib.resize(c->ndattrib);
    
  }
  

  if (myid==0) cerr << "Isotropic scattering on with Tau=" << tauscat << "\n";
  
  dtau  = new double [tautab];
  dtau1 = new double [tautab];
  dr = rmax/tautab;

  cntr = vector<int>(nthrds);

  cntacc = 0;
}

ScatterMFP::~ScatterMFP()
{
  delete [] dtau;
  delete [] dtau1;
}

void ScatterMFP::initialize()
{
  try {
    if (conf["tautab"])         tautab        = conf["tautab"].as<int>();
    if (conf["tauscat"])        tauscat       = conf["tauscat"].as<double>();
    if (conf["rmax"])           rmax          = conf["rmax"].as<double>();
    if (conf["nscat"])          nscat         = conf["nscat"].as<int>();
    if (conf["mfp_index"])      mfp_index     = conf["mfp_index"].as<int>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in ScatterMFP: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }
}


void ScatterMFP::get_acceleration_and_potential(Component* C)
{
  cC = C;			// "Register" component
  nbodies = cC->Number();	// And compute number of bodies

  
  /*======================================*/
  /* Determine potential and acceleration */
  /*======================================*/

  MPL_start_timer();

				// Clean table
  for (int j=0; j<tautab; j++) dtau[j] = dtau1[j] = 0.0;

				// Accumulate table
  int ind;
  double rrr, vvv, RR;
  rr2.clear();

  for (int i=1; i<=nbodies; i++) {
    RR = 0.0;
    for (int j=0; j<3; j++) 
      RR += cC->Pos(i, j) * cC->Pos(i, j);
    RR = sqrt(RR);

    ind = (int)(RR/dr);
    if (ind>=tautab) continue;
    dtau1[ind] += cC->Mass(i);
  }

				// Distribute table
  MPI_Allreduce(dtau1, dtau, tautab, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

				// Compute density from mass per bin
  for (int j=0; j<tautab; j++)
    dtau[j] /= 4.0*M_PI/3.0*(pow(dr*(j+1), 3.0) - pow(dr*j, 3.0));


				// Compute index
  sort(rr2.begin(), rr2.end(), less_rpair);

				// Reset counter
  for (int j=0; j<nthrds; j++) cntr[j] = 0;

  exp_thread_fork(false);

  int cnt0 = 0;
  for (int j=0; j<nthrds; j++) cnt0 += cntr[j];

				// Keep count of mean # of interactions
  int sum=0;
  MPI_Reduce(&cnt0, &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (myid==1) {
    
    cntacc += sum;
    if (!(this_step%nscat)) {
      string file = (string)homedir + "scatter.log";
      ofstream out(file.c_str(), ios::out | ios::app);

      if (this_step == 0)
	out << "# "
	    << setw(6) << this_step
	    << setw(15) << (double)cntacc
	    << endl;
      else
	out << "# "
	    << setw(6) << this_step
	    << setw(15) << (double)cntacc/nscat
	    << endl;

      for (int j=0; j<tautab; j++)
	out << "     " 
	    << setw(15) << dr*(j+1)
	    << setw(15) << dtau[j]
	    << endl;

      cntacc = 0;
    }
  }

#ifdef DEBUG
  // cerr.form("Process %d: rmax=%f\n", myid, rtst);
#endif

  MPL_stop_timer();
}


void * ScatterMFP::determine_acceleration_and_potential_thread(void * arg)
{

  Eigen::Vector3d vcom, vrel, vfnl;

  int i, k, ind;
  double v2;
#ifdef DEBUG  
  double rm, rp, rtst1, rtst=0.0;
#endif

  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  for (int j=nbeg; j<nend; j++) {
    
    i = rr2[j].second;
    
    if (c->freeze(i)) continue;

    ind = (int)(rr2[i].first/dr);
    if (ind>=tautab) continue;

    v2 = 0.0;
    for (int l=0; l<3; l++) 
      v2 += cC->Vel(i, l) * cC->Vel(i, l);

    Particle *p = cC->Part(i);

    p->dattrib[mfp_index] += dtau[ind] * sqrt(v2) * dtime;

    if (1.0 - exp(-p->dattrib[mfp_index]/tauscat)>unif(random_gen)) {

				// Initialize optical depth
      p->dattrib[mfp_index] = 0.0;

				// Choose a buddy
      if (j==1)
	k = rr2[2].second;
#ifdef DEBUG      
      else if ((rp=fabs(rr2[j].first-rr2[rr2[j+1].second].first)) > 
	       (rm=fabs(rr2[j].first-rr2[rr2[j-1].second].first)))
#else

      else if (fabs(rr2[j].first-rr2[rr2[j+1].second].first) > 
	       fabs(rr2[j].first-rr2[rr2[j-1].second].first) )
#endif
	k = rr2[j-1].second;
      else
	k = rr2[j+1].second;

#ifdef DEBUG
      if (j>1 && rtst<(rtst1=min<double>(rm, rp))) rtst = rtst1;
#endif

      Particle *q = cC->Part(k);

      vcom[0] = 0.5*(cC->Vel(i, 0) + cC->Vel(k, 0));
      vcom[1] = 0.5*(cC->Vel(i, 1) + cC->Vel(k, 1));
      vcom[2] = 0.5*(cC->Vel(i, 2) + cC->Vel(k, 2));

      vrel[0] = cC->Vel(k, 0) - cC->Vel(i, 0);
      vrel[1] = cC->Vel(k, 1) - cC->Vel(i, 1);
      vrel[2] = cC->Vel(k, 2) - cC->Vel(i, 2);

				// Choose a random direction for velocity
      vfnl[0] = norm(random_gen);
      vfnl[1] = norm(random_gen);
      vfnl[2] = norm(random_gen);

      vfnl *= sqrt(vrel.dot(vrel))/sqrt(vfnl.dot(vfnl));

				// To lab frame
      p->vel[0] = vcom[0] + 0.5*vfnl[0];
      p->vel[1] = vcom[1] + 0.5*vfnl[1];
      p->vel[2] = vcom[2] + 0.5*vfnl[2];

      q->vel[0] = vcom[0] - 0.5*vfnl[0];
      q->vel[1] = vcom[1] - 0.5*vfnl[1];
      q->vel[2] = vcom[2] - 0.5*vfnl[2];

      cntr[id]++;
    }
  }

  return (NULL);
}
