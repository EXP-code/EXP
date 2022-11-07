#include <math.h>

#include "expand.H"

#include <UserDiffRot.H>

UserDiffRot::UserDiffRot(const YAML::Node& conf) : ExternalForce(conf)
{
  id = "Rotational randomization";

  seed = 11;			// For random number generator
  rate = 0.5;			// Rate relative to dyn time
  name = "";			// Default component name
  avoid = "";
  width = 10.0;
  maxpm = 2;
  ndyn = 25;
  dynmin = 0.001;
  dynmax = 10.0;
  first = true;

  initialize();

  if (name.size()>0) {
				// Look for the fiducial component
    bool found = false;
    for (auto c : comp->components) {
      if ( !name.compare(c->name) ) {
	c0 = c;
	found = true;
      break;
      }
    }

    if (!found) {
      std::ostringstream sout;
      sout << "Can't find desired component <" << name << ">";
      throw GenericError(sout.str(), __FILE__, __LINE__, 35, false);
    }

  }
  else
    c0 = NULL;

  if (avoid.size()>0) {
				// Look for component for particle
				// avoidance
    bool found = false;
    for (auto c : comp->components) {
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
  width = width*M_PI/180.0;

				// Dynamical time distribution initialization
  ddyn = (log(dynmax) - log(dynmin))/(ndyn-1);
  bins = new vector<int> [nthrds];
  for (int n=0; n<nthrds; n++) bins[n] = vector<int>(ndyn, 0);

  userinfo();
}


UserDiffRot::~UserDiffRot()
{
  // NADA
}


void UserDiffRot::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

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

  print_divider();
}

void UserDiffRot::initialize()
{
  try {
    if (conf["name"])    name    = conf["name"].as<string>();
    if (conf["avoid"])   avoid   = conf["avoid"].as<string>();
    if (conf["maxpm"])   maxpm   = conf["maxpm"].as<int>();
    if (conf["seed"])    seed    = conf["seed"].as<int>();
    if (conf["rate"])    rate    = conf["rate"].as<double>();
    if (conf["width"])   width   = conf["width"].as<double>();
    if (conf["seed"])    seed    = conf["seed"].as<int>();
    if (conf["ndyn"])    ndyn    = conf["ndyn"].as<int>();
    if (conf["dynmin"])  dynmin  = conf["dynmin"].as<double>();
    if (conf["dynmax"])  dynmax  = conf["dynmax"].as<double>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in UserDiffRot: "
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


void UserDiffRot::determine_acceleration_and_potential(void)
{
  if (!c0) return;

#if HAVE_LIBCUDA==1		// Cuda compatibility
  getParticlesCuda(c0);
#endif
				// Get particles to avoid
  if (c1) {

    ipm = 0;

    if (myid==0) {

      int number = -1;
      PartPtr *p = c1->get_particles(&number);

      while (p) {
	
	for (int k=0; k<number && ipm<maxpm; k++, ipm++) {
	  pos[ipm*4] = sqrt(p[k]->pos[0]*p[k]->pos[0] +
			    p[k]->pos[1]*p[k]->pos[1] +
			    p[k]->pos[2]*p[k]->pos[2] );
	  pos[ipm*4+1] = p[k]->pos[0];
	  pos[ipm*4+2] = p[k]->pos[1];
	  pos[ipm*4+3] = p[k]->pos[2];
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

  if (first) {
				// Determine next index position for
				// particle double array
    indx = cC->Part(0)->dattrib.size();
  }

  exp_thread_fork(false);

  MPI_Barrier(MPI_COMM_WORLD);

  if (first) {
    
    vector<int> bin0(ndyn, 0);
    for (int n=1; n<nthrds; n++)
      for (int j=0; j<ndyn; j++) bins[0][j] += bins[n][j];
    
    MPI_Reduce(&(bins[0][0]), &(bin0[0]), ndyn, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (myid == 0) {

      cout << "****************************************************************" << endl;
      cout << "UserDiffRot: DTime distribution\n";
      for (int i=0; i<ndyn; i++)
	cout << setw(15) << dynmin*exp(ddyn*(0.5+i)) 
	     << setw(10) << bin0[i] << endl;
      cout << "****************************************************************" << endl;
    }

				// Delete the bin storage
    delete [] bins;

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
  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  double phi, cosp, sinp, diffr;
  double xx, yy, uu, vv;
  
  if (first) {
				// Compute dynamical distribution
    double dtmin = 1.0e20;
    double dtmax = 0.0;
    double dt;
    int jindx;

    PartMapItr it = cC->Particles().begin();
    unsigned long i;

    for (int q=0   ; q<nbeg; q++) it++;
    for (int q=nbeg; q<nend; q++) {
      i = (it++)->first;

      dt = get_dtime(*(cC->Part(i)));

      cC->Part(i)->dattrib.push_back(tnow + dt/rate);

      dtmin = min<double>(dt, dtmin);
      dtmax = max<double>(dt, dtmax);

      jindx = max<int>(0, min<int>(ndyn-1, (int)((log(dt)-log(dynmin))/ddyn)));

      bins[id][jindx]++;
    }
      
  }

  PartMapItr it = cC->Particles().begin();
  unsigned long i;

  for (int q=0   ; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {
    i = (it++)->first;

    // If we are multistepping, compute accel only at or below this level
    //
    if (multistep && (cC->Part(i)->level < mlevel)) continue;
    
    // Avoid?
    bool tooclose = false;

    if (c1) {
      for (int n=0; n<ipm; n++) {
	diffr = 0.0;
	for (int k=0; k<3; k++)
	  diffr += 
	    (cC->Pos(i, k) - pos[n*4+1+k]) *
	    (cC->Pos(i, k) - pos[n*4+1+k]) ;
	
	if (sqrt(diffr)<pos[n*4]) tooclose = true;
      }
    }
    
    if (tooclose) continue;

				// Sanity check
    if ((int)cC->Part(i)->dattrib.size() < indx+1) {
      cout << "***Size error***: Myid=" << myid << "  id=" << id
	   << "  check i=" << i 
	   << "  nbeg=" << nbeg 
	   << "  nend=" << nend 
	   << "  size=" << cC->Part(i)->dattrib.size()
	   << "\n";
    }

    if (tnow>cC->Part(i)->dattrib[indx]) {

      cC->Part(i)->dattrib[indx] = tnow + get_dtime(*(cC->Part(i)))/rate;

				// Do rotation
      phi =  width * normal(random_gen);
      cosp = cos(phi);
      sinp = sin(phi);
      
      xx = cC->Pos(i, 0)*cosp - cC->Pos(i, 1)*sinp;
      yy = cC->Pos(i, 0)*sinp + cC->Pos(i, 1)*cosp;
      uu = cC->Vel(i, 0)*cosp - cC->Vel(i, 1)*sinp;
      vv = cC->Vel(i, 0)*sinp + cC->Vel(i, 1)*cosp;

      cC->Part(i)->pos[0] = xx;
      cC->Part(i)->pos[1] = yy;
      cC->Part(i)->vel[0] = uu;
      cC->Part(i)->vel[1] = vv;
    }

  }

  return (NULL);
}


extern "C" {
  ExternalForce *makerDiffRot(const YAML::Node& conf)
  {
    return new UserDiffRot(conf);
  }
}

class proxydiffrot { 
public:
  proxydiffrot()
  {
    factory["userdiffrot"] = makerDiffRot;
  }
};

static proxydiffrot p;
