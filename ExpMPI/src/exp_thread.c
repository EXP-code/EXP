
#include "expand.h"
#include <exp_thread.h>

char threading_on = 0;
pthread_mutex_t mem_lock;

void exp_thread_fork(void * (*func)(void *), char * caller)
{
  int i;
  pthread_t *thread;
  int *thr_id, errcode;
  void *retval;
  
  if (!(thread = (pthread_t *) malloc(nthrds*sizeof(pthread_t)))) {
    fprintf(stderr, "Process %d: exp_thread_fork: error allocating memory for thread\n", myid);
    exit(18);
  }
  if (!(thr_id = (int *) malloc(nthrds*sizeof(int)))) {
    fprintf(stderr, "Process %d: exp_thread_fork: error allocating memory for thread counters\n", myid);
    exit(18);
  }
				/* Make the memory mutex */
  make_mutex(&mem_lock, caller, "memory_lock");
  threading_on = 1;


				/* Make the <nthrds> threads */
  for (i=0; i<nthrds; i++) {
    thr_id[i] = i;
    if ((errcode=pthread_create(&thread[i], NULL, func, &thr_id[i]))) {
      fprintf(stderr, "Process %d, %s thread: cannot make thread %d, errcode=%d\n", 
	      myid, caller, thr_id[i], errcode);
      exit(19);
    }
  }
    
				/* Collapse the threads */
  for (i=0; i<nthrds; i++) {
    if ((errcode=pthread_join(thread[i], &retval))) {
      fprintf(stderr, "Process %d, %s thread: thread join %d failed, errcode=%d\n", 
	      myid, caller, i, errcode);
      exit(20);
    }
  }
  
				/* Remove the memory mutex */
  kill_mutex(&mem_lock, caller, "memory_lock");
  threading_on = 0;


  free(thread);
  free(thr_id);
}


void make_mutex(pthread_mutex_t *m, char *caller, char *name)
{
  int errcode;

  if ((errcode=pthread_mutex_init(m, NULL))) {
    fprintf(stderr, "Process %d, %s: mutex init %s failed, errcode=%d\n", 
	    myid, caller, name, errcode);
    exit(21);
  }
}

void kill_mutex(pthread_mutex_t *m, char * caller, char *name)
{
  int errcode;
  
  if ((errcode=pthread_mutex_destroy(m))) {
    fprintf(stderr, "%s: mutex destroy %s failed, errcode=%d\n", 
	    caller, name, errcode);
    exit(21);
  }
}
