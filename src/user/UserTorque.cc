#include <math.h>

#include "expand.h"

#include <localmpi.h>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>

#include <orbit.h>

class UserTorque : public ExternalForce
{
  typedef vector<float> fvector;
  SphericalModelTable *halo;
  SphericalOrbit *orb;

private:
  
  string com_name, model_name, file_name;

  Component *c0;

  int numx, numy, sgn;
  float xmin, xmax, ymin, ymax;
  double dX, dY, ton, toff, delta, boost;
  vector<fvector> array;

  int diverge;
  double diverge_rfac;

  pthread_mutex_t orb_lock;

  double interpolate(double E, double K);

  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  void userinfo();

public:

  UserTorque(string &line);
  ~UserTorque();

};


UserTorque::UserTorque(string &line) : ExternalForce(line)
{
  file_name = "dLz.dat";	// Data file
  com_name = "sphereSL";	// Default component for com
				// Model file
  model_name = "SLGridSph.model";
				// Sign for torque
  sgn = 1;
				// Use analytic cusp at small radius
  diverge = 0;
				// Cusp exponent
  diverge_rfac = 1.0;
    
  ton =-1.0e20;
  toff= 1.0e20;
  delta = 0.25;
  boost = 1.0;

  initialize();

				// Look for the fiducial component
  bool found = false;
  list<Component*>::iterator cc;
  Component *c;
  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    c = *cc;
    if ( !com_name.compare(c->id) ) {
      c0 = c;
      found = true;
      break;
    }
  }

  if (!found) {
    cerr << "Process " << myid << ": can't find desired component <"
	 << com_name << ">" << endl;
    MPI_Abort(MPI_COMM_WORLD, 35);
  }

  userinfo();

  pthread_mutex_init(&orb_lock, NULL);


}

UserTorque::~UserTorque()
{
  delete orb;
  delete halo;
  pthread_mutex_destroy(&orb_lock);
}

void UserTorque::userinfo()
{
  if (myid) return;		// Return if node master node
  print_divider();
  cout << "** User routine TORQUE GRID initialized: file=<" << file_name 
       << ">, Ton=" << ton 
       << ", Toff=" << toff
       << ", Delta=" << delta 
       << ", boost=" << boost 
       << ", sgn=" << sgn
       << endl;
  print_divider();
}

void UserTorque::initialize()
{
  string val;

  if (get_value("comname", val))	com_name = val;
  if (get_value("file_name", val))	file_name = val;
  if (get_value("model_name", val))	model_name = val;
  if (get_value("sgn", val))	        sgn = atoi(val)>0 ? 1 : -1;
  if (get_value("diverge", val))	diverge = atoi(val);
  if (get_value("diverge_rfac", val))	diverge_rfac = atof(val);
  if (get_value("ton", val))		ton = atof(val);
  if (get_value("toff", val))		toff = atof(val);
  if (get_value("delta", val))		delta = atof(val);
  if (get_value("boost", val))		boost = atof(val);
  
  ifstream in(file_name.c_str());
  if (in) {
    in.read((char *)&numx, sizeof(int));
    in.read((char *)&numy, sizeof(int));
    in.read((char *)&xmin, sizeof(float));
    in.read((char *)&xmax, sizeof(float));
    in.read((char *)&ymin, sizeof(float));
    in.read((char *)&ymax, sizeof(float));
    
    dX = (xmax - xmin)/numx;
    dY = (ymax - ymin)/numy;

    float z;

    for (int j=0; j<numy; j++) {
      fvector tmp;
      for (int i=0; i<numx; i++) {
	in.read((char *)&z, sizeof(float));
	tmp.push_back(z);
      }
      array.push_back(tmp);
    }

  } else {
    throw "UserTorque: could not open input file";
  }

  halo = new SphericalModelTable(model_name, diverge, diverge_rfac);
  orb  = new SphericalOrbit(halo);
}

double UserTorque::interpolate(double E, double K)
{
  if (E<=xmin || E>=xmax) return 0.0;
  if (K<=ymin || K>=ymax) return 0.0;

  int indx, indy;
  double ax, bx, ay, by;

  indx = (int)( (E - xmin)/dX );
  indx = max<int>(0, indx);
  indx = min<int>(numx-2, indx);
  ax = (E - dX*indx - xmin)/dX;
  bx = 1.0 - ax;

  indy = (int)( (K - ymin)/dY );
  indy = max<int>(0, indy);
  indy = min<int>(numy-2, indy);
  ay = (K - dY*indy - ymin)/dY;
  by = 1.0 - ay;

  double ret =
    ax*ay*array[indy  ][indx  ] +
    bx*ay*array[indy  ][indx+1] +
    ax*by*array[indy+1][indx  ] +
    bx*by*array[indy+1][indx+1] ;

  return ret;
}

void * UserTorque::determine_acceleration_and_potential_thread(void * arg) 
{
  double pos[3], vel[3], ll[3], rr, vv;
  double E, J, K, torque, phi;
  
  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = 1+nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  double amp = 0.25*sgn*boost*
    (1.0 + erf((tnow-ton )/delta)) *
    (1.0 + erf((toff-tnow)/delta)) ;

  PartMapItr it = cC->Particles().begin();
  unsigned long i;

  for (int q=0   ; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {
    i = (it++)->first;
				// If we are multistepping, compute accel 
				// only at or below this level

    if (multistep && (cC->Part(i)->level < mlevel)) continue;

    rr = vv = 0.0;
    for (int k=0; k<3; k++) {
      pos[k] = cC->Pos(i, k) - c0->com[k];
      vel[k] = cC->Vel(i, k);
      rr += pos[k]*pos[k];
      vv += vel[k]*vel[k];
    }

    ll[0] = pos[1]*vel[2] - pos[2]*vel[1];
    ll[1] = pos[2]*vel[0] - pos[0]*vel[2];
    ll[2] = pos[0]*vel[1] - pos[1]*vel[0];

    rr = sqrt(rr);
    E = 0.5*vv + halo->get_pot(rr);
    J = sqrt(ll[0]*ll[0] + ll[1]*ll[1] + ll[2]*ll[2]);

    if(rr > halo->get_max_radius()) continue;
	
    if(E < halo->get_pot(halo->get_min_radius())){
      cout << "illegal value of Energy: " << E
           << "  Emin[r=" << halo->get_min_radius() << "] = "
           << halo->get_pot(halo->get_min_radius());

      E = halo->get_pot(halo->get_min_radius());

      cout << ":: Now modify the Energy to E =" << E << endl;
    }
    // orbit_trans routines are not reentrant . . .
    pthread_mutex_lock(&orb_lock);
    orb->new_orbit(E, 0.5);
    K = J/orb->Jmax();
    pthread_mutex_unlock(&orb_lock);

    torque = amp*interpolate(E, K);

				// Increase/decrease tangential velocity
    phi = atan2(pos[1], pos[0]);
    
    cC->AddAcc(i, 0, -torque*sin(phi)/rr );
    cC->AddAcc(i, 1,  torque*cos(phi)/rr );
  }

  return (NULL);
}


extern "C" {
  ExternalForce *makerTorque(string& line)
  {
    return new UserTorque(line);
  }
}

class proxytorque { 
public:
  proxytorque()
  {
    factory["usertorque"] = makerTorque;
  }
};

proxytorque p;
