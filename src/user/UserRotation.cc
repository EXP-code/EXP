#include <math.h>
#include <sstream>

#include "expand.H"

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>

//
// Rotate the halo about given center
//

class UserRotation : public ExternalForce
{
private:
  
  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  double RROT, Omega;
  bool EJ;

  void userinfo();

public:

  //! Constructor
  UserRotation(string &line);

  //! Destructor
  ~UserRotation();

};

UserRotation::UserRotation(string &line) : ExternalForce(line)
{

  //  id = "Halo Rotation";	       

  Omega = 0;			// Rotation frequency
  EJ = false;                   // Default center is COM
  RROT = 0.0;			// Maximum radius in which the particle is affected by rotation 
                                // If RROT = 0.-, the influence radius is infinity 
  initialize();

  userinfo();
}

UserRotation::~UserRotation()
{
}

void UserRotation::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  std::cout << "** User routine Rotation initialized, " ;
  std::cout << "Omega=" << Omega;
  std::cout << " , the rotation center is ";
  if (EJ) 
    std::cout << " the EJ center";
  else 
    std::cout << " the COM";
  std::cout << "and R_rot is";
  if (RROT > 0.0)  
    std::cout << RROT;
  else 
    std::cout << " Infinity";
  std::cout << std::endl;
  
  print_divider();
}

void UserRotation::initialize()
{
  std::string val;

  if (get_value("Omega", val))   Omega = atof(val.c_str());
  if (get_value("RROT", val))    RROT = atof(val.c_str());
  if (get_value("EJ",val))       EJ = atol(val);

  if (RROT < 0.0) RROT =0.0;
}


void UserRotation::determine_acceleration_and_potential(void)
{
#if HAVE_LIBCUDA==1		// Cuda compatibility
  getParticlesCuda(cC);
#endif

  exp_thread_fork(false);
}


void * UserRotation::determine_acceleration_and_potential_thread(void * arg) 
{

  if (this_step) {

    unsigned nbodies = cC->Number();
    int id = *((int*)arg);
    int nbeg = nbodies*id/nthrds;
    int nend = nbodies*(id+1)/nthrds;

    double rr, angle;
    double pos0[3], vel0[3], pos1[3], vel1[3], dpos[3], dvel[3];

    for (int i=nbeg; i<nend; i++) {
      for(int k=0; k<3; k++){
	if(!EJ) {
	  pos0[k] = cC->Pos(i, k) - cC->com0[k];
	  vel0[k] = cC->Vel(i, k) - cC->cov0[k];
	} else {
	  pos0[k] = cC->Pos(i, k) - (cC->com0[k] + cC->center[k]);
	  vel0[k] = cC->Vel(i, k) - (cC->cov0[k] + cC->center[k]);
	}
      }

      rr=0.0;
      for(int k=0; k<3; k++) rr += pos0[k]*pos0[k];
      rr = sqrt(rr);

      if((rr < RROT) && (RROT > 0.0)){
	angle = Omega*dtime;
	
	pos1[0] = pos0[0]*cos(angle) - pos0[1]*sin(angle);
	pos1[1] = pos0[0]*sin(angle) + pos0[1]*cos(angle);
	pos1[2] = pos0[2];
	
	vel1[0] = vel0[0]*cos(angle) - vel0[1]*sin(angle);
	vel1[1] = vel0[0]*sin(angle) + vel0[1]*cos(angle);
	vel1[2] = vel0[2];
	
	for(int k=0; k<3; k++){
	  dpos[k] = pos1[k]-pos0[k]; 
	  dvel[k] = vel1[k]-vel0[k]; 
	}
	
	// Subtract additional artificial Coriolis force

	dvel[0] += (-2.0*vel1[1]*angle);
	dvel[1] += (2.0*vel1[0]*angle);
	dvel[2] += 0.0;

	// Rotate position and velocity
	cC->AddPos(i, dpos);
	cC->AddVel(i, dvel);
      }
    }
  }
  return (NULL);
}


extern "C" {
  ExternalForce *makerRotation(string& line)
  {
    return new UserRotation(line);
  }
}

class proxyrotate { 
public:
  proxyrotate()
  {
    factory["userrotation"] = makerRotation;
  }
};

static proxyrotate p;
