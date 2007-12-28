
#include "expand.h"
#include <PotAccel.H>

extern "C"
void *
call_any_threads_thread_call(void *atp)
{
  thrd_pass_PotAccel *tp = (thrd_pass_PotAccel *)atp;
  PotAccel *p = (PotAccel *)tp->t;
  if (tp->coef)
    p -> determine_coefficients_thread((void*)&tp->id);
  else
    p -> determine_acceleration_and_potential_thread((void*)&tp->id);
  return NULL;
}
                   
pthread_mutex_t PotAccel::cc_lock;
int PotAccel::compute;

void PotAccel::exp_thread_fork(bool coef)
{
  int errcode;
  void *retval;
  
  td = new thrd_pass_PotAccel [nthrds];
  t = new pthread_t [nthrds];

  if (!td) {
    cerr << "Process " << myid 
	 << ": exp_thread_fork: error allocating memory for thread counters\n";
    exit(18);
  }
  if (!t) {
    cerr << "Process " << myid
	 << ": exp_thread_fork: error allocating memory for thread\n";
    exit(18);
  }
				// Enable the memory mutex
  // threading_on = 1;


				// Make the <nthrds> threads
  for (int i=0; i<nthrds; i++) {
    td[i].t = this;
    td[i].coef = coef;
    td[i].id = i;

    errcode =  pthread_create(&t[i], 0, call_any_threads_thread_call, &td[i]);
    if (errcode) {
      cerr << "Process " << myid;
      if (coef)	cerr << ", make_coefficients";
      else cerr << ", determine_acceleration";
      cerr << " thread: cannot make thread " << i
	  << ", errcode=" << errcode << endl;
      exit(19);
    }
#ifdef DEBUG
    else {
      cout << "Process " << myid << ": thread <" << i << "> created\n";
    }
#endif
  }
    
				// Collapse the threads
  for (int i=0; i<nthrds; i++) {
    if ((errcode=pthread_join(t[i], &retval))) {
      cerr << "Process " << myid;
      if (coef)	cerr << ", make_coefficients";
      else cerr << ", determine_acceleration";
      cerr << " thread: thread join " << i
	   << " failed, errcode=" << errcode << endl;
      exit(20);
    }
  }
  
				// Disable the memory mutex
  // threading_on = 0;


  delete [] td;
  delete [] t;
}


void PotAccel::make_mutex(pthread_mutex_t *m, const char *caller, 
			  const char *name)
{
  int errcode;

  if ((errcode=pthread_mutex_init(m, NULL))) {
    cerr << "Process " << myid << ", "
	 << caller << ": mutex init " 
	 << name << " failed, errcode= " << errcode << endl;
    exit(21);
  }
}

void PotAccel::kill_mutex(pthread_mutex_t *m, const char * caller, 
			  const char *name)
{
  int errcode;
  
  if ((errcode=pthread_mutex_destroy(m))) {
    cerr << "Process " << myid << ", "
	 << caller << ": mutex destroy " 
	 << name << " failed, errcode= " << errcode << endl;
    exit(21);
  }
}

void PotAccel::bomb(const string& msg)
{
  cerr << "Component [" << id << ": " << msg << endl;
  exit(-1);
}


PotAccel::PotAccel(string& line)
{
  used = 0;
  component = NULL;
  geometry = other;
  use_external = false;
  coef_dump = false;
  mlevel = 0;

  // Per thread counter
  use = new int [nthrds];
  if (!use) bomb("problem allocating <use>");

  string Line = trimLeft(trimRight(line));
  StringTok<string> tokens(Line);
  pair<string, string> datum;

				// Comma separated tokens
  string token = tokens(",");
  
  while (token.size()) {

    StringTok<string> parse(token);
    datum.first  = trimLeft(trimRight(parse("=")));
    datum.second = trimLeft(trimRight(parse("=")));
    namevalue.push_back(datum);

				// Next parameter
    token = tokens(",");
  }

}

PotAccel::~PotAccel(void)
{
  delete [] use;
}

int PotAccel::get_value(const string& name, string& value)
{
  list< pair<string, string> >::iterator it;
  for (it=namevalue.begin(); it!=namevalue.end(); it++) {
    if (it->first.compare(name) == 0) {
      value = it->second;
      return 1;
    }
  }
  return 0;
}
