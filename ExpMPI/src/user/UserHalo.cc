#include <math.h>
#include <sstream>

#include "expand.h"

#include <UserHalo.H>

UserHalo::UserHalo(string &line) : ExternalForce(line)
{
  id = "SphericalHalo";
  model_file = "SLGridSph.model";
  diverge = 0;
  diverge_rfac = 1.0;

  initialize();

  model = new SphericalModelTable(model_file, diverge, diverge_rfac);

  userinfo();
}

UserHalo::~UserHalo()
{
  delete model;
}

void UserHalo::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine SPHERICAL HALO initialized, " ;
  cout << "Filename=" << model_file << "  diverge=" << diverge
       << "  diverge_rfac=" << diverge_rfac << endl;
  
  print_divider();
}

void UserHalo::initialize()
{
  string val;

  if (get_value("model_file", val))	model_file = val;
  if (get_value("diverge", val))	diverge = atoi(val.c_str());
  if (get_value("diverge_rfac", val))	diverge_rfac = atof(val.c_str());
}


void UserHalo::determine_acceleration_and_potential(void)
{
  exp_thread_fork(false);
}


void * UserHalo::determine_acceleration_and_potential_thread(void * arg) 
{
  int nbodies = particles->size();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  double xx, yy, zz, rr, r, pot, dpot;

  for (int i=nbeg; i<nend; i++) {

    xx = (*particles)[i].pos[0];
    yy = (*particles)[i].pos[1];
    zz = (*particles)[i].pos[2];
    rr = xx*xx + yy*yy + zz*zz;
    r = sqrt(rr);

    model->get_pot_dpot(r, pot, dpot);


    (*particles)[i].acc[0] += -dpot*xx/r;
    
    (*particles)[i].acc[1] += -dpot*yy/r;

    (*particles)[i].acc[2] += -dpot*zz/r;
    
    (*particles)[i].potext += pot;
  }

  return (NULL);
}


extern "C" {
  ExternalForce *makerHalo(string& line)
  {
    return new UserHalo(line);
  }
}

class proxyhalo { 
public:
  proxyhalo()
  {
    factory["userhalo"] = makerHalo;
  }
};

proxyhalo p;
