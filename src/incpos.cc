/*
  First part of leap-frog: step in position
*/

#include "expand.H"

#ifdef USE_GPTL
#include <gptl.h>
#endif

#ifdef HAVE_LIBCUDA
void incr_position_cuda(cuFP_t dt, int mlevel);
#endif

void * incr_position_thread(void *ptr)
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
  Component *c;
  
  //
  // Component loop
  //
  for (auto c : comp->components) {

    if (mlevel>=0)		// Use a particular level
      ntot = c->levlist[mlevel].size();
    else			// Use ALL levels
      ntot = c->Number();
      
    if (ntot==0) continue;

    //
    // Compute the beginning and end points in particle list
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

      Particle *p = c->Part(indx);
      for (int k=0; k<c->dim; k++) 
	p->pos[k] += p->vel[k]*dt;
    }
  }

  return (NULL);
}


void incr_position(double dt, int mlevel)
{
  if (!eqmotion) return;

#ifdef USE_GPTL
  GPTLstart("incr_position");
#endif

#ifdef HAVE_LIBCUDA
  if (use_cuda) {
    incr_position_cuda(static_cast<cuFP_t>(dt), mlevel);
    return;
  }
#endif

  if (nthrds==1) {

    posvel_data[0].dt = dt;
    posvel_data[0].mlevel = mlevel;
    posvel_data[0].id = 0;

    incr_position_thread(&posvel_data[0]);

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
      errcode =  pthread_create(p, 0, incr_position_thread, &posvel_data[i]);

      if (errcode) {
	std::ostringstream sout;
	sout << "Process " << myid
	     << " incr_position: cannot make thread " << i
	     << ", errcode=" << errcode;;
	throw GenericError(sout.str(), __FILE__, __LINE__, 1024, true);
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
	std::ostringstream sout;
	sout << "Process " << myid
	     << " incr_position: thread join " << i
	     << " failed, errcode=" << errcode;
	throw GenericError(sout.str(), __FILE__, __LINE__, 1024, true);
      }
#ifdef DEBUG    
      cout << "Process " << myid << ": incr_position thread <" 
	   << i << "> thread exited\n";
#endif
    }
  }
  
#ifdef USE_GPTL
  GPTLstop("incr_position");
#endif

}

void incr_com_position(double dt)
{
  for (auto c : comp->components) {
    if (c->com_system) {
      for (int k=0; k<c->dim; k++) c->com0[k] += c->cov0[k]*dt;
    }
  }
}

