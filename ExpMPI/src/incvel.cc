/*
  Second part of leap-frog: step in velocity
*/

#include "expand.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

/// Helper class to pass info to threaded drift routine
struct thrd_pass_incvel {
  //! Time step
  double dt;

  //! Levels flag
  int mlevel;

  //! Thread counter id
  int id;
};



void * incr_velocity_thread(void *ptr)
{
  // Current time step
  double dt = static_cast<thrd_pass_incvel*>(ptr)->dt;

  // Current level
  int mlevel = static_cast<thrd_pass_incvel*>(ptr)->mlevel;

  // Thread ID
  int id = static_cast<thrd_pass_incvel*>(ptr)->id;


  list<Component*>::iterator cc;
  vector<Particle>::iterator p, pend;
  int nbeg, nend;
  Component *c;
  
  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    c = *cc;
    
    unsigned ntot = c->Number();

    //
    // Compute the beginning and end points for threads
    //
    nbeg = ntot*id    /nthrds;
    nend = ntot*(id+1)/nthrds;

    for (unsigned i=nbeg; i<nend; i++) {
      // If we are multistepping, only 
				// advance for this level ONLY
      if (multistep && (c->Part(i)->level != mlevel)) continue;

      for (int k=0; k<c->dim; k++) 
	c->Part(i)->vel[k] += (c->Part(i)->acc[k] - c->acc0[k])*dt;
    }

  }

}


void incr_velocity(double dt, int mlevel)
{
  if (!eqmotion) return;

  
  thrd_pass_incvel* td = new thrd_pass_incvel [nthrds];

  if (!td) {
    cerr << "Process " << myid
	 << ": incr_velocity: error allocating thread structures\n";
    exit(18);
  }

  pthread_t* t  = new pthread_t [nthrds];

  if (!t) {
    cerr << "Process " << myid
	 << ": incr_velocity: error allocating memory for thread\n";
    exit(18);
  }

  //
  // Make the <nthrds> threads
  //
  int errcode;
  void *retval;
  
  for (int i=0; i<nthrds; i++) {

    td[i].dt = dt;
    td[i].mlevel = mlevel;
    td[i].id = i;
    
    errcode =  pthread_create(&t[i], 0, incr_velocity_thread, &td[i]);

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
    if ((errcode=pthread_join(t[i], &retval))) {
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
  
  delete [] td;
  delete [] t;


  list<Component*>::iterator cc;
  Component *c;
  
  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    c = *cc;
    
    if (c->com_system) {	// Only do this once per multistep
      if (multistep==0 || (mstep==Mstep && mlevel==multistep))
	for (int k=0; k<c->dim; k++) c->cov0[k] += c->acc0[k]*dt;
    }
    
  }

  
}

