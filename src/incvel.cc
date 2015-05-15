/*
  Second part of leap-frog: step in velocity
*/

#include "expand.h"

#ifdef USE_GPTL
#include <gptl.h>
#endif

void * incr_velocity_thread(void *ptr)
{
  // Current time step
  //
  double dt = static_cast<thrd_pass_posvel*>(ptr)->dt;

  // Current level
  //
  int mlevel = static_cast<thrd_pass_posvel*>(ptr)->mlevel;

  // Thread ID
  //
  int id = static_cast<thrd_pass_posvel*>(ptr)->id;


  int nbeg, nend, indx;
  unsigned ntot;
  
  //
  // Component loop
  //
  for (auto c : comp.components) {
    
    if (mlevel>=0)		// Use a particular level
      ntot = c->levlist[mlevel].size();
    else			// Use ALL levels
      ntot = c->Number();

    if (ntot==0) continue;

    //
    // Compute the beginning and end points in the particle vector 
    // for each thread
    //
    nbeg = ntot*(id  )/nthrds;
    nend = ntot*(id+1)/nthrds;
    
    PartMapItr it = c->Particles().begin();
    
    for (int q=0   ; q<nbeg; q++) it++;
    for (int q=nbeg; q<nend; q++) {

      if (mlevel>=0)
	indx = c->levlist[mlevel][q];
      else 
	indx = (it++)->first;
      
      if (c->com_system) {
	for (int k=0; k<c->dim; k++) 
	  c->Part(indx)->vel[k] += (c->Part(indx)->acc[k] - c->acc0[k])*dt;
      } else {
	for (int k=0; k<c->dim; k++) 
	  c->Part(indx)->vel[k] += c->Part(indx)->acc[k]*dt;
      }

    }
      
  }

}

void incr_velocity(double dt, int mlevel)
{
  if (!eqmotion) return;

#ifdef USE_GPTL
    GPTLstart("incr_velocity");
#endif


  if (nthrds==1) {

    posvel_data[0].dt = dt;
    posvel_data[0].mlevel = mlevel;
    posvel_data[0].id = 0;
    
    incr_velocity_thread(&posvel_data[0]);

  } else {

    //
    // Make the <nthrds> threads
    //
    int errcode;
    void *retval;
    
    for (int i=0; i<nthrds; i++) {

      posvel_data[i].dt = dt;
      posvel_data[i].mlevel = mlevel;
      posvel_data[i].id = i;
      
      pthread_t *p = &posvel_thrd[i];
      errcode =  pthread_create(p, 0, incr_velocity_thread, &posvel_data[i]);

      if (errcode) {
	cerr << "Process " << myid
	     << " incr_velocity: cannot make thread " << i
	     << ", errcode=" << errcode << endl;
	exit(19);
      }
#ifdef DEBUG
      else {
	cout << "Process " << myid << ": thread <" << i << "> created\n";
      }
#endif
    }
    
    //
    // Collapse the threads
    //
    for (int i=0; i<nthrds; i++) {
      pthread_t p = posvel_thrd[i];
      if ((errcode=pthread_join(p, &retval))) {
	cerr << "Process " << myid
	     << " incr_velocity: thread join " << i
	     << " failed, errcode=" << errcode << endl;
	exit(20);
      }
#ifdef DEBUG    
      cout << "Process " << myid << ": incr_velocity thread <" 
	   << i << "> thread exited\n";
#endif
    }
  }

#ifdef USE_GPTL
  GPTLstop("incr_velocity");
#endif

}

void incr_com_velocity(double dt)
{
  for (auto c : comp.components) {

    if (c->com_system) {
      for (int k=0; k<c->dim; k++) c->cov0[k] += c->acc0[k]*dt;
    }
    
  }

}

