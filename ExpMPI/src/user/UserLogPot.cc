#include <math.h>
#include <sstream>

#include "expand.h"

#include <UserLogPot.H>

UserLogPot::UserLogPot(string &line) : ExternalForce(line)
{
  id = "LogarithmicPotential";

  initialize();

  userinfo();
}

UserLogPot::~UserLogPot()
{
}

void UserLogPot::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine LOGARITHMIC POTENTIAL initialized, " ;
  cout << "Phi = v2/2*log(R^2 + x^2 + y^2/b^2 + z^2/c^2) with R=" 
       << R << ", b=" << b << ", c=" << c << ", v2=" << v2 << endl; 
  
  print_divider();
}

void UserLogPot::initialize()
{
  string val;

  if (get_value("R", val))		R = atof(val.c_str());
  if (get_value("b", val))		b = atof(val.c_str());
  if (get_value("c", val))		c = atof(val.c_str());
  if (get_value("v2", val))		v2 = atof(val.c_str());
}


void UserLogPot::determine_acceleration_and_potential(void)
{
  exp_thread_fork(false);
}


void * UserLogPot::determine_acceleration_and_potential_thread(void * arg) 
{
  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  double xx, yy, zz, rr;

  map<unsigned long, Particle>::iterator it = cC->Particles().begin();
  unsigned long i;

  for (int q=0   ; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {
    i = (it++)->first;
				// If we are multistepping, compute accel 
				// only at or below this level

    if (multistep && (cC->Part(i)->level < mlevel)) continue;

    xx = cC->Pos(i, 0);
    yy = cC->Pos(i, 1);
    zz = cC->Pos(i, 2);
    rr = R*R + xx*xx + yy*yy/(b*b) + zz*zz/(c*c);

    cC->AddAcc(i, 0,-v2*xx/rr );
    
    cC->AddAcc(i, 1, -v2*yy/(rr*b*b) );

    cC->AddAcc(i, 2, -v2*zz/(rr*c*c) );
    
    cC->AddPotExt(i, 0.5*v2*log(rr) );
  }

  return (NULL);
}


extern "C" {
  ExternalForce *makerLogPot(string& line)
  {
    return new UserLogPot(line);
  }
}

class proxylogpot { 
public:
  proxylogpot()
  {
    factory["userlogp"] = makerLogPot;
  }
};

proxylogpot p;
