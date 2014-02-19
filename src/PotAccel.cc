#include <pthread.h>
#include <sys/time.h>
#include <time.h>

#include <boost/lexical_cast.hpp>

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
  //
  // If only one thread, skip pthread call
  //
  if (nthrds==1) {

    thrd_pass_PotAccel td;

    td.t = this;
    td.coef = coef;
    td.id = 0;

    call_any_threads_thread_call(&td);

    return;
  }

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

  //
  // For determining time in threaded routines
  //
  if (comp.timing) {

    switch (comp.state) {
    case ComponentContainer::SELF:
      comp.timer_thr_acc.start();
      break;
    case ComponentContainer::INTERACTION:
      comp.timer_thr_int.start();
      break;
    case ComponentContainer::EXTERNAL:
      comp.timer_thr_ext.start();
      break;
    }

  }

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
  
  //
  // For determining time in threaded routines
  //
  if (comp.timing) {

    switch (comp.state) {
    case ComponentContainer::SELF:
      comp.timer_thr_acc.stop();
      break;
    case ComponentContainer::INTERACTION:
      comp.timer_thr_int.stop();
      break;
    case ComponentContainer::EXTERNAL:
      comp.timer_thr_ext.stop();
      break;
    }

  }

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
  used         = 0;
  component    = NULL;
  cC           = NULL;
  geometry     = other;
  use_external = false;
  coef_dump    = false;
  dof          = 3;
  mlevel       = 0;

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

  if (VERBOSE>5) {
    tv_list = vector<struct timeval>(nthrds);
    timer_list = vector<double>(2*nthrds);
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

std::map<int, std::string> PotAccel::get_value_array(const string& name)
{
  std::map<int, string> values;
  int indx;

  list< pair<string, string> >::iterator it;
  for (it=namevalue.begin(); it!=namevalue.end(); it++) {
    string key = name + "(";
    if (it->first.compare(0, key.size(), key) == 0) {
      string sindx = it->first.substr(key.size(), it->first.find(")"));
      try {
	indx = boost::lexical_cast<int>(sindx);
      } 
      catch( boost::bad_lexical_cast const& ) {
	std::cout << "PotAccel::get_value_array: input string <" 
		  << it->first << "> is not valid" << std::endl;
      }
      values[indx] = it->second;
    }
  }
  return values;
}


std::map<std::pair<int, int>, string> 
PotAccel::get_value_matrix(const string& name)
{
  std::map<std::pair<int, int>, string> values;
  std::pair<int, int> indx;

  list< pair<string, string> >::iterator it;
  for (it=namevalue.begin(); it!=namevalue.end(); it++) {
    string key = name + "(";
    if (it->first.compare(0, key.size(), key) == 0) {
      string sindx1 = it->first.substr(key.size(), it->first.find(","));
      string sindx2 = it->first.substr(key.size() + sindx1.size(), 
				       it->first.find(")"));
      try {
	indx.first  = boost::lexical_cast<int>(sindx1);
	indx.second = boost::lexical_cast<int>(sindx2);
      } 
      catch( boost::bad_lexical_cast const& ) {
	std::cout << "PotAccel::get_value_matrix: input string <" 
		  << it->first << "> is not valid" << std::endl;
      }
      values[indx] = it->second;
    }
  }
  return values;
}


void PotAccel::print_timings(const string& label) 
{
  if (VERBOSE>5) print_timings(label, timer_list);
}

void PotAccel::thread_timing_beg(int id)
{
  if (VERBOSE>5) {
    gettimeofday(&tv_list[id], 0);
    timer_list[2*id] = 
      tv_list[id].tv_usec*1.0e-6 +
      (tv_list[id].tv_sec % 1000);
  }
}

void PotAccel::thread_timing_end(int id)
{
  if (VERBOSE>5) {
    gettimeofday(&tv_list[id], 0);
    timer_list[2*id+1] = 
      tv_list[id].tv_usec*1.0e-6 +
      (tv_list[id].tv_sec % 1000);
  }
}

void PotAccel::print_timings(const string& label, vector<double>& tlist)
{
  if (VERBOSE<=5) return;

  if (myid==0) {

    NData nbeg, nend;

    nbeg.name = "beg";
    nend.name = "end";

    vector< pair<double, NData> > total_list;
    nbeg.node = nend.node = 0;
    for (int n=0; n<nthrds; n++) {
      nbeg.tid = nend.tid = n;
      total_list.push_back(pair<double, NData>(tlist[2*n+0], nbeg));
      total_list.push_back(pair<double, NData>(tlist[2*n+1], nend));
    }
    vector<double> from(2*nthrds);
    for (int np=1; np<numprocs; np++) {
      MPI_Recv(&from[0], 2*nthrds, MPI_DOUBLE, np, 37, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
      nbeg.node = nend.node = np;
      for (int n=0; n<nthrds; n++) {
	nbeg.tid = nend.tid = n;
	total_list.push_back(pair<double, NData>(tlist[2*n+0], nbeg));
	total_list.push_back(pair<double, NData>(tlist[2*n+1], nend));
      }
    }

    sort(total_list.begin(), total_list.end(), ltdub);

    string labl(label);
    if (component && cC) {	// Standard case
      if (component->name.size()>0 && cC->name.size()>0) 
	labl = labl + " [" + component->name + "==>" + cC->name + "]";
      else if (component->name.size()>0)
	labl = labl + " [" + component->name + "]";
      else if (cC->name.size()>0) 
	labl = labl + " [" + cC->name + "]";
    } else if (component) {	// I don't think this ever happens . . . 
      if (component->name.size()>0)
	labl = labl + " [" + component->name + "]";
    } else if (cC) {		// Case for external force
      if (cC->name.size()>0)
	labl = labl + " [" + cC->name + "]";
    }
      
    
    cout << labl << endl
	 << setfill('-') << setw(labl.size()) << "-" << setfill(' ') << endl;

    if (1) {
      double bmean=0.0, bdisp=0.0, emean=0.0, edisp=0.0;
      bool first_end=false, ok=true;
      int cnt=0;
      for (vector<pair<double, NData> >::iterator 
	     it=total_list.begin(); it!=total_list.end(); it++) {

	cnt++;
	if (it->second.name.compare("beg")==0) {
	  bmean += it->first;
	  bdisp += it->first * it->first;
	  if (first_end) ok = false;
	} else {
	  emean += it->first;
	  edisp += it->first * it->first;
	  first_end = true;
	}
      }

      cnt   /= 2;
      bmean /= cnt;
      bdisp = fabs(bdisp - bmean*bmean*cnt);
      if (cnt>1) bdisp /= cnt-1;

      emean /= cnt;
      edisp = fabs(edisp - emean*emean*cnt);
      if (cnt>1) edisp /= cnt-1;

      cout << "Mean(beg) = " << setw(12) << setprecision(6) << bmean
	   << "     " << setw(12) 
	   << "Stdv(beg) = " << setw(12) << setprecision(6) << sqrt(bdisp)
	   << endl
	   << "Mean(end) = " << setw(12) << setprecision(6) << emean
	   << "     " << setw(12)
	   << "Stdv(end) = " << setw(12) << setprecision(6) << sqrt(edisp)
	   << endl
	   << "DelT mean = " << setw(12) << setprecision(6) << emean - bmean
	   << "     " << setw(12) 
	   << "DelT/Stdv = " << setw(12) << (emean - bmean)/sqrt(bdisp + edisp + 1.0e-14)
	   << "  ";

      if (!ok) cout << "[OVERLAP!]" << endl;
      else     cout << endl;

    } else {

      cout << setw(5) << "Node" << setw(5) << "Tid" 
	   << setw(12) << "Time" << endl;
      for (vector<pair<double, NData> >::iterator 
	     it=total_list.begin(); it!=total_list.end(); it++) {
	cout << setw(5)  << it->second.node
	     << setw(5)  << it->second.tid
	     << setw(12) << setprecision(6) << fixed << it->first
	     << setw(5)  << it->second.name
	     << endl;
      }
    }
    cout << setw(70) << setfill('=') << "=" << endl << setfill(' ');
    
  } else {
    MPI_Send(&tlist[0], 2*nthrds, MPI_DOUBLE, 0, 37, MPI_COMM_WORLD);
  }
}

bool ltdub(const pair<double, PotAccel::NData>& A, 
	   const pair<double, PotAccel::NData>& B)
{
  return A.first < B.first;
} 

