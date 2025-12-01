#include <pthread.h>
#include <sys/time.h>
#include <time.h>

#include <Eigen/Eigen>

#include <highfive/highfive.hpp>
#include <highfive/eigen.hpp>

#include "expand.H"
#include "PotAccel.H"

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

// Map from enum to string
std::map<PotAccel::Geometry, std::string> PotAccel::geoname =
  { {PotAccel::sphere,   "sphere"  },
    {PotAccel::cylinder, "cylinder"},
    {PotAccel::cube,     "cube"    },
    {PotAccel::slab,     "slab"    },
    {PotAccel::table,    "table"   },
    {PotAccel::other,    "other"   }
  };

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
    std::ostringstream sout;
    sout << "Process " << myid 
	 << ": exp_thread_fork: error allocating memory for thread counters";
    throw GenericError(sout.str(), __FILE__, __LINE__, 1027, true);
  }
  if (!t) {
    std::ostringstream sout;
    sout << "Process " << myid
	 << ": exp_thread_fork: error allocating memory for thread\n";
    throw GenericError(sout.str(), __FILE__, __LINE__, 1027, true);
  }

  //
  // For determining time in threaded routines
  //
  if (comp->timing) {

    switch (comp->state) {
    case ComponentContainer::SELF:
      comp->timer_thr_acc.start();
      break;
    case ComponentContainer::INTERACTION:
      comp->timer_thr_int.start();
      break;
    case ComponentContainer::EXTERNAL:
      comp->timer_thr_ext.start();
      break;
    case ComponentContainer::NONE:
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
      std::ostringstream sout;
      sout << "Process " << myid;
      if (coef)	sout << ", make_coefficients";
      else sout << ", determine_acceleration";
      sout << " thread: cannot make thread " << i
	   << ", errcode=" << errcode;
      throw GenericError(sout.str(), __FILE__, __LINE__, 1027, true);
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
      std::ostringstream sout;
      sout << "Process " << myid;
      if (coef)	sout << ", make_coefficients";
      else sout << ", determine_acceleration";
      sout << " thread: thread join " << i
	   << " failed, errcode=" << errcode;
      throw GenericError(sout.str(), __FILE__, __LINE__, 1027, true);
    }
  }
  
  //
  // For determining time in threaded routines
  //
  if (comp->timing) {

    switch (comp->state) {
    case ComponentContainer::SELF:
      comp->timer_thr_acc.stop();
      break;
    case ComponentContainer::INTERACTION:
      comp->timer_thr_int.stop();
      break;
    case ComponentContainer::EXTERNAL:
      comp->timer_thr_ext.stop();
      break;
    case ComponentContainer::NONE:
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
    std::ostringstream sout;
    sout << "Process " << myid << ", "
	 << caller << ": mutex init " 
	 << name << " failed, errcode= " << errcode;
    throw GenericError(sout.str(), __FILE__, __LINE__, 1028, true);
  }
}

void PotAccel::kill_mutex(pthread_mutex_t *m, const char * caller, 
			  const char *name)
{
  int errcode;
  
  if ((errcode=pthread_mutex_destroy(m))) {
    std::ostringstream sout;
    sout << "Process " << myid << ", "
	 << caller << ": mutex destroy " 
	 << name << " failed, errcode= " << errcode;
    throw GenericError(sout.str(), __FILE__, __LINE__, 1028, true);
  }
}

PotAccel::PotAccel(Component* c0, const YAML::Node& CONF) : conf(CONF)
{
  used         = 0;
  component    = c0;
  cC           = NULL;
  geometry     = other;
  use_external = false;
  coef_dump    = false;
  play_back    = false;
  play_cnew    = false;
  compute      = false;
  dof          = 3;
  mlevel       = 0;
  scale        = 1.0;
#if HAVE_LIBCUDA==1
  cuda_aware   = false;
#endif
  
  // Per thread counter
  try {
    use.resize(nthrds);
  }
  catch (...) {
    throw GenericError("problem allocating <use>", __FILE__, __LINE__, 1027, true);
  }

  if (VERBOSE>5) {
    timer_list = vector<std::time_t>(2*nthrds);
  }

  // Add keys
  for (YAML::const_iterator it=conf.begin(); it!=conf.end(); ++it) {
    current_keys.insert(it->first.as<std::string>());
  }

}

PotAccel::~PotAccel(void)
{
}


void PotAccel::print_timings(const string& label) 
{
  if (VERBOSE>5) print_timings(label, timer_list);
}

void PotAccel::thread_timing_beg(int id)
{
  if (VERBOSE>5) {
    auto const now = std::chrono::system_clock::now();
    std::time_t newt = std::chrono::system_clock::to_time_t(now);
    timer_list[2*id] = newt;
  }
}

void PotAccel::thread_timing_end(int id)
{
  if (VERBOSE>5) {
    auto const now = std::chrono::system_clock::now();
    std::time_t newt = std::chrono::system_clock::to_time_t(now);
    timer_list[2*id+1] = newt;
  }
}

void PotAccel::print_timings(const string& label, TList& tlist)
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
      for (auto it : total_list) {
	cnt++;
	if (it.second.name.compare("beg")==0) {
	  bmean += it.first;
	  bdisp += it.first * it.first;
	  if (first_end) ok = false;
	} else {
	  emean += it.first;
	  edisp += it.first * it.first;
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
      for (auto it : total_list) {
	cout << setw(5)  << it.second.node
	     << setw(5)  << it.second.tid
	     << setw(12) << setprecision(6) << fixed << it.first
	     << setw(5)  << it.second.name
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

