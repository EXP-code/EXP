#include <math.h>
#include <strstream>

#include "expand.h"

#include <ACG.h>
#include <Normal.h>

#include <UserDiffRot.H>

UserDiffRot::UserDiffRot(string &line) : ExternalForce(line)
{
  id = "Rotational randomization";

  seed = 11;			// For random number generator
  rate = 0.5;			// Rate relative to dyn time
  name = "";			// Default component name
  avoid = "";
  width = 10.0;
  maxpm = 2;
  first = true;

  initialize();

  if (name.size()>0) {
				// Look for the fiducial component
    bool found = false;
    list<Component*>::iterator cc;
    Component *c;
    for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
      c = *cc;
      if ( !name.compare(c->name) ) {
	c0 = c;
	found = true;
      break;
      }
    }

    if (!found) {
      cerr << "Process " << myid << ": can't find desired component <"
	   << name << ">" << endl;
      MPI_Abort(MPI_COMM_WORLD, 35);
    }

  }
  else
    c0 = NULL;

  if (avoid.size()>0) {
				// Look for component for particle
				// avoidance
    bool found = false;
    list<Component*>::iterator cc;
    Component *c;
    for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
      c = *cc;
      if ( !avoid.compare(c->name) ) {
	c1 = c;
	found = true;
      break;
      }
    }

    if (!found) {
      cerr << "Process " << myid << ": can't find desired component <"
	   << avoid << ">" << endl;
      c1 = NULL;
    }

  }
  else
    c1 = NULL;

  pos = vector<double>(4*maxpm);
  gen = new ACG(seed);
  normal = new Normal(0.0, 1.0, gen);
  width = width*M_PI/180.0;

  userinfo();

  // Debug

  ndyn = 25;
  dynmin = 0.001;
  dynmax = 10.0;
  ddyn = (log(dynmax) - log(dynmin))/(ndyn-1);

  vector<int> *bins = new vector<int> [nthrds];
  for (int n=0; n<nthrds; n++) bins[n] = vector<int>(ndyn, 0);
}


UserDiffRot::~UserDiffRot()
{
  delete gen;
  delete normal;
  delete [] bins;
}


void UserDiffRot::userinfo()
{
  if (myid) return;		// Return if node master node
  if (c0)
    cout << "** User routine ROTATION RANDOMIZATION initialized for component: <" << name << ">";
  else
    cout << "** User routine ROTATION RANDOMIZATION disabled: no component specified";
  
  if (avoid.size() > 0) cout << ", avoid = " << avoid;
  cout << ", maxpm = " << maxpm;
  cout << ", rate = " << rate;
  cout << ", width = " << width;
  cout << ", seed = " << seed;
  cout << endl;
  cout << "****************************************************************"
       << endl;

}

void UserDiffRot::initialize()
{
  string val;

  if (get_value("name", val))		name = val;
  if (get_value("avoid", val))		avoid = val;
  if (get_value("maxpm", val))		maxpm = atoi(val.c_str());
  if (get_value("seed", val))		seed = atoi(val.c_str());
  if (get_value("rate", val))		rate = atof(val.c_str());
  if (get_value("width", val))		width = atof(val.c_str());
  if (get_value("seed", val))		seed = atoi(val.c_str());
}


void UserDiffRot::determine_acceleration_and_potential(void)
{
  if (!c0) return;

				// Get particles to avoid
  if (c1) {

    ipm = 0;

    if (myid==0) {

      int number = -1;
      Component::Partstruct *p = c1->get_particles(&number);

      while (p) {
	
	for (int k=0; k<number && ipm<maxpm; k++, ipm++) {
	  pos[ipm*4] = sqrt(p[k].pos[0]*p[k].pos[0] +
			    p[k].pos[1]*p[k].pos[1] +
			    p[k].pos[2]*p[k].pos[2] );
	  pos[ipm*4+1] = p[k].pos[0];
	  pos[ipm*4+2] = p[k].pos[1];
	  pos[ipm*4+3] = p[k].pos[2];
	}

	if (ipm == maxpm) break;
	p = c1->get_particles(&number);
      }
      
      cout << "****************************************************************";
      cout << "UserDiffRot: avoiding " << ipm << " particles\n";
      cout << "****************************************************************";
    }

    MPI_Bcast(&ipm, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&pos[0], ipm*4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }

  exp_thread_fork(false);

  MPI_Barrier(MPI_COMM_WORLD);

  if (first) {
    
    cout << "Myid=" << myid << " in distribution diag\n";

    for (int n=1; n<nthrds; n++) {
      for (int j=0; j<ndyn; j++) bins[0][j] += bins[n][j];
    }

    vector<int> bin0(ndyn, 0);

    cout << "Myid=" << myid << " about to call Reduce\n";

    MPI_Reduce(&(bins[0][0]), &(bin0[0]), ndyn, MPI_INT, MPI_SUM, 0, 
	       MPI_COMM_WORLD);

    cout << "Myid=" << myid << " returned from Reduce\n";

    if (myid == 0) {

      cout << "****************************************************************";
      cout << "UserDiffRot: DTime distribution\n";
      for (int i=0; i<ndyn; i++)
	cout << setw(15) << exp(dynmin + ddyn*(0.5+i)) 
	     << setw(10) << bin0[i] << endl;
      cout << "****************************************************************";
    }

    first = false;
  }

}


double UserDiffRot::get_dtime(Particle& p) 
{
  double vv, E, Lx, Ly, Lz, dt;

				// Compute energy
  vv = 0.0;
  for (int k=0; k<3; k++)
    vv += p.vel[k] * p.vel[k];
      
  E = 0.5*vv + p.pot + p.potext;

  if (E>=0.0) {
    E = -1.0e-08;
    /*
    cout << "Energy error:  T=" << 0.5*vv 
	 << "  W=" << p.pot << ", " <<  p.potext << endl;
    */
  }

				// Compute angular momentum
  Lx = 
    p.pos[1]*p.vel[2]-
    p.pos[2]*p.vel[1];
  
  Ly = 
    p.pos[2]*p.vel[0]-
    p.pos[0]*p.vel[2];
  
  Lz =
    p.pos[0]*p.vel[1]-
    p.pos[1]*p.vel[0];
  
  dt = sqrt(Lx*Lx + Ly*Ly + Lz*Lz)/(-2.0*E);
  
  return dt;
}


void * UserDiffRot::determine_acceleration_and_potential_thread(void * arg) 
{
  int nbodies = particles->size();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  double phi, cosp, sinp, diffr;
  double xx, yy, uu, vv;
  
  if (first) {
				// Debug
    double dtmin = 1.0e20;
    double dtmax = 0.0;
    double dt;
    int jindx;

    indx = (*particles)[nbeg].dattrib.size();

    for (int i=nbeg; i<nend; i++) {
      dt = get_dtime((*particles)[i]);
      (*particles)[i].dattrib.push_back(tnow + dt/rate);
      dtmin = min<double>(dt, dtmin);
      dtmax = max<double>(dt, dtmax);
      if (dt<=dynmin) jindx = 0;
      else if (dt>=dynmax) jindx = ndyn-1;
      else jindx = (int)(log(dt - dynmin)/ddyn);
      bins[id][jindx]++;
    }
      
				// Debug
    
    cout << "***Size check***: Myid=" << myid << "  id=" << id
	 << "  nbeg=" << nbeg 
	 << "  nend=" << nend 
	 << "  size=" << (*particles)[nbeg].dattrib.size()
	 << "\n";


				// Debug
    /*
    cout << "****************************************************************"
	 << endl
	 << " UserDiffRot:"
	 << "  myid=" << myid 
	 << "  id=" << id 
	 << "  indx=" << indx 
	 << "  dtmin=" << dtmin 
	 << "  dtmax=" << dtmax << endl
	 << "****************************************************************"
	 << endl;
    */
    // first = false;
  }


  cout << "Myid=" << myid << "  id=" << id << " about to enter MAIN LOOP\n";

  for (int i=nbeg; i<nend; i++) {

    if ((*particles)[i].freeze()) continue;

				// Avoid?
    bool tooclose = false;

    if (c1) {
      cout << "Myid=" << myid << "  id=" << id << " in AVOID\n";
      for (int n=0; n<ipm; n++) {
	diffr = 0.0;
	for (int k=0; k<3; k++)
	  diffr += 
	    ((*particles)[i].pos[k] - pos[n*4+1+k]) *
	    ((*particles)[i].pos[k] - pos[n*4+1+k]) ;
	
	if (sqrt(diffr)<pos[n*4]) tooclose = true;
      }
    }
    
    if (tooclose) continue;

    if ((*particles)[i].dattrib.size() < indx+1) {
      cout << "***Size error***: Myid=" << myid << "  id=" << id
	   << "  check i=" << i 
	   << "  nbeg=" << nbeg 
	   << "  nend=" << nend 
	   << "  size=" << (*particles)[i].dattrib.size()
	   << "\n";
    }

    // cout << "Myid=" << myid << "  id=" << id << " n=" << i << "\n";

    if (tnow>(*particles)[i].dattrib[indx]) {

      cout << "Myid=" << myid << "  id=" << id
	   << "  time=" << (*particles)[i].dattrib[indx] << "\n";
	   
      (*particles)[i].dattrib[indx] = tnow + get_dtime((*particles)[i])/rate;

				// Do rotation
      phi =  width * (*normal)();
      cosp = cos(phi);
      sinp = sin(phi);
      
      xx = (*particles)[i].pos[0]*cosp - (*particles)[i].pos[1]*sinp;
      yy = (*particles)[i].pos[0]*sinp + (*particles)[i].pos[1]*cosp;
      uu = (*particles)[i].vel[0]*cosp - (*particles)[i].vel[1]*sinp;
      vv = (*particles)[i].vel[0]*sinp + (*particles)[i].vel[1]*cosp;

      (*particles)[i].pos[0] = xx;
      (*particles)[i].pos[1] = yy;
      (*particles)[i].vel[0] = uu;
      (*particles)[i].vel[1] = vv;
    }

  }

  cout << "Myid=" << myid << "  id=" << id << " exiting MAIN LOOP\n";

  return (NULL);
}


extern "C" {
  ExternalForce *makerDiffRot(string& line)
  {
    return new UserDiffRot(line);
  }
}

class proxy { 
public:
  proxy()
  {
    factory["userdiffrot"] = makerDiffRot;
  }
};

proxy p;
