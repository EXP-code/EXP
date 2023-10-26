#include <sys/timeb.h>
#include <math.h>
#include <sstream>

#include "expand.H"

#include <PeriodicBC.H>

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

const std::set<std::string>
PeriodicBC::valid_keys = {
  "compname",
  "sx",
  "sy",
  "sz",
  "cx",
  "cy",
  "cz",
  "dT",
  "nbin",
  "tcol",
  "vunit",
  "btype"
};

PeriodicBC::PeriodicBC(const YAML::Node& conf) : ExternalForce(conf)
{
  (*barrier)("PeriodicBC: BEGIN construction", __FILE__, __LINE__);

  id = "PeriodicBC";		// Periodic boundary condition ID

				// Sizes in each dimension
  L = vector<double>(3, 2.0);	
				// Center offset in each dimension
  offset = vector<double>(3, 1.0);

  bc = "ppp";			// Periodic BC in all dimensions

  comp_name = "";		// Default component (must be specified)

  nbin = 0;			// Number of bins in trace (0 means no trace)
  dT = 1.0;			// Interval for trace
  tcol = -1;			// Column for temperture info (ignored if <0)
  trace = false;		// Tracing off until signalled

  initialize();

				// Look for the fiducial component
  bool found = false;
  for (auto c : comp->components) {
    if ( !comp_name.compare(c->name) ) {
      c0 = c;
      found = true;
      break;
    }
  }

  if (!found) {
    std::ostringstream sout;
    sout << "PeriodicBC: "
	 << "can't find fiducial component <" << comp_name << ">";
    throw GenericError(sout.str(), __FILE__, __LINE__, 35, false);
  }
  
  // 
  // Initialize structures for trace data
  //
  if (nbin) {
    cbinT = vector< vector<unsigned> >(nthrds);
    mbinT = vector< vector<double>   >(nthrds);
    tbinT = vector< vector<double>   >(nthrds);
    for (int n=0; n<nthrds; n++) {
      cbinT[n] = vector<unsigned>(nbin, 0);
      mbinT[n] = vector<double>  (nbin, 0);
      tbinT[n] = vector<double>  (nbin, 0);
    }
    cbin = vector<unsigned>(nbin);
    mbin = vector<double>  (nbin);
    tbin = vector<double>  (nbin);
    Tnext = tnow;
    dX = L[0]/nbin;
  }

  userinfo();

  atomic_weights.resize(13, -1.0);

  atomic_weights[0]  = 0.000548579909; // Mass of electron
  atomic_weights[1]  = 1.0079;
  atomic_weights[2]  = 4.0026;
  atomic_weights[3]  = 6.941;
  atomic_weights[4]  = 9.0122;
  atomic_weights[5]  = 10.811;
  atomic_weights[6]  = 12.011;
  atomic_weights[7]  = 14.007;
  atomic_weights[8]  = 15.999;
  atomic_weights[9]  = 18.998;
  atomic_weights[10] = 20.180;
  atomic_weights[11] = 22.990;
  atomic_weights[12] = 24.305;

  (*barrier)("Periodic: END construction", __FILE__, __LINE__);
}

PeriodicBC::~PeriodicBC()
{
}

void PeriodicBC::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine PERIODIC BOUNDARY CONDITION initialized"
       << " using component <" << comp_name << ">";
  if (nbin) cout << " with gas trace, dT=" << dT 
		 << ", nbin=" << nbin << ", tcol=" << tcol;
  cout << endl;

  cout << "   Cube sides (x , y , z) = (" 
       << L[0] << " , " 
       << L[1] << " , " 
       << L[2] << " ) " << endl; 

  cout << "Center offset (x , y , z) = (" 
       << offset[0] << " , " 
       << offset[1] << " , " 
       << offset[2] << " ) " << endl; 

  cout << "Boundary type (x , y , z) = (" 
       << bc[0] << " , " 
       << bc[1] << " , " 
       << bc[2] << " ) " << endl;

  print_divider();
}

void PeriodicBC::initialize()
{
  // Check for unmatched keys
  //
  auto unmatched = YamlCheck(conf, valid_keys);
  if (unmatched.size())
    throw YamlConfigError("PeriodicBC", "parameter", unmatched,
			  __FILE__, __LINE__);

  // Assign values from YAML
  //
  try {
    if (conf["compname"])       comp_name          = conf["compname"].as<string>();
    
    if (conf["sx"])             L[0]               = conf["sx"].as<double>();
    if (conf["sy"])             L[1]               = conf["sy"].as<double>();
    if (conf["sz"])             L[2]               = conf["sz"].as<double>();
    
    if (conf["cx"])             offset[0]          = conf["cx"].as<double>();
    if (conf["cy"])             offset[1]          = conf["cy"].as<double>();
    if (conf["cz"])             offset[2]          = conf["cz"].as<double>();
    
    if (conf["dT"])             dT                 = conf["dT"].as<double>();
    if (conf["nbin"])           nbin               = conf["nbin"].as<int>();
    if (conf["tcol"])           tcol               = conf["tcol"].as<int>();
    
    if (conf["vunit"])          vunit              = conf["vunit"].as<double>();
    if (conf["temp"])           temp               = conf["temp"].as<double>();
    
    
    if (conf["btype"]) {
      std::string val = conf["btype"].as<std::string>();
      if (strlen(val.c_str()) >= 3) {
	for (int k=0; k<3; k++) {
	  switch (val.c_str()[k]) {
	  case 'p':
	    bc[k] = 'p';		// Periodic
	    break;
	  case 'r':
	    bc[k] = 'r';		// Reflection
	    break;
	  default:
	    bc[k] = 'v';		// Vacuum
	    break;
	  }
	}
      }
    }
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in PeriodicBC: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }
    
}


void PeriodicBC::determine_acceleration_and_potential(void)
{
  if (cC != c0) return;
#if HAVE_LIBCUDA==1		// Cuda compatibility
  getParticlesCuda(cC);
#endif
  if (nbin && tnow>=Tnext) trace = true;
  exp_thread_fork(false);
  if (trace) write_trace();
  print_timings("PeriodicBC: thread timings");
}

void PeriodicBC::write_trace()
{
  for (int n=1; n<nthrds; n++) {
    for (int k=0; k<nbin; k++) {
      cbinT[0][k] += cbinT[n][k];
      mbinT[0][k] += mbinT[n][k];
      tbinT[0][k] += tbinT[n][k];
    }
  }

  MPI_Reduce(&cbinT[0][0], &cbin[0], nbin, MPI_UNSIGNED, MPI_SUM, 0, 
	     MPI_COMM_WORLD);

  MPI_Reduce(&mbinT[0][0], &mbin[0], nbin, MPI_DOUBLE,   MPI_SUM, 0, 
	     MPI_COMM_WORLD);

  MPI_Reduce(&tbinT[0][0], &tbin[0], nbin, MPI_DOUBLE,   MPI_SUM, 0, 
	     MPI_COMM_WORLD);

  if (myid==0) {
    ostringstream fout;
    fout << outdir << runtag << ".shocktube_trace";
    ofstream out(fout.str().c_str(), ios::app);
    for (int k=0; k<nbin; k++)
      out << setw(18) << tnow
	  << setw(18) << dX*(0.5+k) - offset[0]
	  << setw(18) << mbin[k]/(dX*L[1]*L[2])
	  << setw(18) << tbin[k]/(mbin[k]+1.0e-10)
	  << setw(10) << cbin[k]
	  << endl;
    out << endl;
  }

  //
  // Clean data structures for next call
  //
  for (int n=0; n<nthrds; n++) {
    for (int k=0; k<nbin; k++) {
      cbinT[n][k] = 0;
      mbinT[n][k] = tbinT[n][k] = 0.0;
    }
  }

  trace = false;		// Tracing off until
  Tnext += dT;			// tnow>=Tnext
}


void * PeriodicBC::determine_acceleration_and_potential_thread(void * arg) 
{
  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  
  thread_timing_beg(id);

  double pos, delta;
  PartMapItr it = cC->Particles().begin();
  
  std::advance(it, nbeg);

  for (int q=nbeg; q<nend; q++) {
    
				// Index for the current particle
    unsigned long i = (it++)->first;
    
    Particle *p = cC->Part(i);
    double   mi = 0.0;

    for (int k=0; k<3; k++) {

      // Ignore vacuum boundary dimensions
      //
      if (bc[k] == 'v') continue;

      // Increment so that the positions range
      // between 0 and L[k]
      //
      pos = p->pos[k] + offset[k];

      //
      // Reflection BC
      //
      if (bc[k] == 'r') {
	if (pos < 0.0) {
	  delta = -pos - L[k]*floor(-pos/L[k]);
	  p->pos[k] = delta - offset[k];
	  p->vel[k] *= -1.0;
	} 
	if (pos >= L[k]) {
	  delta = pos - L[k]*floor(pos/L[k]);
	  p->pos[k] =  L[k] - delta - offset[k];
	  p->vel[k] *= -1.0;
	}
      }

      //
      // Periodic BC
      //
      if (bc[k] == 'p') {
	if (pos < 0.0) {
	  p->pos[k] += L[k]*floor(1.0+fabs(pos/L[k]));
	  
	}
	if (pos >= L[k]) {
	  p->pos[k] += - L[k]*floor(fabs(pos/L[k]));
	}
      }

      //
      // Sanity check
      //
      if (p->pos[k] < -offset[k] || p->pos[k] >= L[k]-offset[k]) {
	cout << "Process " << myid << " id=" << id 
	     << ": Error in pos[" << k << "]=" << p->pos[k] << endl;
      }
    }

    //
    // Acccumlate data for shocktube trace
    //
    if (trace) {
      int indx = static_cast<int>(floor((p->pos[0]+offset[0])/dX));
      if (indx>=0 && indx<nbin) {
	cbinT[id][indx]++;
	mbinT[id][indx] += p->mass;
	if (tcol>=0 && tcol<static_cast<int>(p->dattrib.size()))
	  tbinT[id][indx] += p->mass*p->dattrib[tcol];
      }
    }

  }
  
  thread_timing_end(id);

  return (NULL);
}

