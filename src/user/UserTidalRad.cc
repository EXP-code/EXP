#include <sstream>
#include <numeric>
#include <cmath>

#include "expand.H"
#include <localmpi.H>
#include <gaussQ.H>

#include <UserTidalRad.H>


UserTidalRad::UserTidalRad(const YAML::Node &conf) : ExternalForce(conf)
{
  id = "TidalRadius";
				// Output file name
  filename = outdir + runtag + ".TidalRad";

  comp_name = "";		// Default component for com
  rtrunc    = 1.0;		// Default tidal truncation
  rfactor   = 1.0;		// Fraction of rtruc for setting the scale
  rtorig    = 0.5;		// Original percentile radius target
  dtTrunc   = 0.0;		// Truncation setting interval (zero: skip)
  dtScale   = 0.0;		// Force scale setting interval (zero: skip)
  debug     = false;		// Print mean percentile radii (testing)
  pctile    = 0.9;              // Percentile for radial average

  pcnbin    = 0;		// Number of points on either side of
				// percential target for averaging. Set
				// from simulation, if zero.

  boxcar    = 1;		// Number of steps in time series average
  diag      = 0;		// Diagnostic output stride (0 means never)
  
  initialize();

  boxcar = std::max<unsigned>(1, boxcar);

  if (comp_name.size()>0) {
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
      cerr << "Process " << myid << ": can't find desired component <"
	   << comp_name << ">" << endl;
      MPI_Abort(MPI_COMM_WORLD, 35);
    }

  } else {
    if (myid==0) {
      std:: cerr << "UserTidalRad: desired component name must be specified"
	<< std::endl;
	   MPI_Abort(MPI_COMM_WORLD, 36);
    }
  }

  if (myid==0) {

    rt_cur = rtrunc;

    if (restart) {

      // Backup up old file
      std::string backupfile = filename + ".bak";
      std::string command("cp ");
      command += filename + " " + backupfile;
      if (system(command.c_str()) == -1) {
	std::cerr << "UserTidalRad: error in executing <"
		  << command << ">" << endl;
      }
	
      // Open new output stream for writing
      std::ofstream out(filename);
      if (!out) {
	throw FileCreateError(filename,
			      "UserTidalRad: error opening new log file",
			      __FILE__, __LINE__);
      }
	
      // Open old file for reading
      ifstream in(backupfile.c_str());
      if (!in) {
	throw FileOpenError(backupfile, "UserTidalRad: error opening original log file",
			    __FILE__, __LINE__);
      }
      
      const int linesize = 1024;
      char line[linesize];
      
      in.getline(line, linesize); // Discard header
      in.getline(line, linesize); // Next line
      
      double lasttime, radius;
      
      while (in) {
	istringstream ins(line);
	
	ins >> lasttime;
	ins >> radius;

	radsol.push_back(radius);

	if (lasttime >= tnow) break;

	out << line;

	in.getline(line, linesize); // Next line
      }

      cout << "UserTidalRad: restart at T=" << lasttime 
	   << " with rtrunc=" << radius
	   << std::endl;
      
      rt_cur = radius;

    } else {

      ofstream out(filename.c_str(), ios::out);

      out << std::left
	  << setw(16) << "# Time"
	  << setw(16) << "Radius"
	  << setw(16) << "Energy"
	  << setw(16) << "u"
	  << setw(16) << "v"
	  << setw(16) << "w"
	  << setw(16) << "E_min"
	  << setw(16) << "E_max"
	  << setw(10) << "Number"
	  << std::endl;
    }

    // Radial averages for debugging
    //
    if (debug) {
      // Open output stream for writing
      //
      std::ofstream dbg("TidalRadDebug." + runtag, ios::out);
      if (dbg.good()) {
	const std::vector<double> pct = {0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95};
	// Labels
	//
	dbg << "# " << std::setw(14) << std::left << "Time |";
	for (auto p : pct) {
	  std::ostringstream sout1, sout2;
	  sout1 << "E(" << p << ") |";
	  sout2 << "R(" << p << ") |";
	  dbg << std::setw(16) << std::left << sout1.str()
	      << std::setw(16) << std::left << sout2.str();
	}
	{
	  std::ostringstream sout1, sout2;
	  sout1 << "E(" << 1.0 << ") |";
	  sout2 << "R(" << 1.0 << ") |";
	  dbg << std::setw(16) << std::left << sout1.str()
	      << std::setw(16) << std::left << sout2.str()
	      << std::endl;
	}
	dbg << std::endl;

	// Column numbers
	//
	int coln = 1;
	std::ostringstream sout;
	sout << "[" << coln++ << "] |";
	dbg << "# " << std::setw(14) << std::left << sout.str();
	for (size_t n=0; n<pct.size(); n++) {
	  std::ostringstream sout1, sout2;
	  sout1 << "[" << coln++ << "] |";
	  sout2 << "[" << coln++ << "] |";
	  dbg << std::setw(16) << std::left << sout1.str()
	      << std::setw(16) << std::left << sout2.str();
	}
	{
	  std::ostringstream sout1, sout2;
	  sout1 << "[" << coln++ << "] |";
	  sout2 << "[" << coln++ << "] |";
	  dbg << std::setw(16) << std::left << sout1.str()
	      << std::setw(16) << std::left << sout2.str()
	      << std::endl;
	}
	dbg << std::endl;

	// Seperators
	//
	dbg << "# " << std::setw(14) << std::setfill('-') << '+';
	for (size_t n=0; n<pct.size(); n++) {
	  dbg << std::setw(16) << std::setfill('-') << '+'
	      << std::setw(16) << std::setfill('-') << '+';
	}
	{
	  dbg << std::setw(16) << std::setfill('-') << '+'
	      << std::setw(16) << std::setfill('-') << '+';
	}
	dbg << std::endl;
      }
    }
    // Debug radial percentile file stanza
  }
  // File initialization loop

  MPI_Bcast(&rt_cur, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Set next truncation time
  //
  if (dtTrunc>0.0) tnextT = tnow + dtTrunc;

  // Set next scale change time
  //
  if (dtScale>0.0) tnextS = tnow + dtScale;

  userinfo();

}

UserTidalRad::~UserTidalRad()
{
}

void UserTidalRad::userinfo()
{
  if (myid) return;		// Return if node master node


  print_divider();

  cout << "** User routine TIDAL RADIUS initialized "
       << "using component <" << comp_name << "> "
       << "with boxcar=" << boxcar << ", dtTrunc=" << dtTrunc
       << ", dtScale=" << dtScale << ", diag=" << diag
       << std::endl;

  print_divider();
}

void UserTidalRad::initialize()
{
  try {
    if (conf["compname"]) comp_name = conf["compname"].as<std::string>();
    if (conf["rtrunc"])   rtrunc    = conf["rtrunc"].  as<double>();
    if (conf["rfactor"])  rfactor   = conf["rfactor"]. as<double>();
    if (conf["rtorig"])   rtorig    = conf["rtorig"].  as<double>();
    if (conf["pctile"])   pctile    = conf["pctile"].  as<double>();
    if (conf["pcnbin"])   pcnbin    = conf["pcnbin"].  as<int>();
    if (conf["boxcar"])   boxcar    = conf["boxcar"].  as<unsigned>();
    if (conf["dtTrunc"])  dtTrunc   = conf["dtTrunc"]. as<double>();
    if (conf["dtScale"])  dtScale   = conf["dtScale"]. as<double>();
    if (conf["diag"])     diag      = conf["diag"].    as<int>();
    if (conf["debug"])    debug     = conf["debug"].   as<bool>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in UserTidalRad: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }

  cov.resize(3*nthrds);
  mas.resize(nthrds);
}


double UserTidalRad::radavg()
{
  if (radsol.size()>boxcar) {
    size_t old = radsol.size() - boxcar;
    radsol.erase(radsol.begin(), radsol.begin() + old);
  }

  double ret = std::accumulate(radsol.begin(), radsol.end(), 0.0);
  if (radsol.size()) ret /= radsol.size();
  return ret;
}


void UserTidalRad::determine_acceleration_and_potential(void)
{
  if (cC != c0) return;		// Check that this component is the target

#if HAVE_LIBCUDA==1		// Cuda compatibility
  getParticlesCuda(cC);
#endif
				// Only compute for top level
  if (multistep && mlevel>0) return;

				// Resize storage vectors
  erg.resize(c0->Number());
  rad.resize(c0->Number());
  std::fill(erg.begin(), erg.end(), 0.0);
  std::fill(rad.begin(), rad.end(), 0.0);

				// Initialize storage with zeros
  std::fill(cov.begin(), cov.end(), 0.0);
  std::fill(mas.begin(), mas.end(), 0.0);

  // Center of velocity pass
  //
  pass = 0;

  exp_thread_fork(false);


  // Compute center of velocity
  //
  for (int n=1; n<nthrds; n++) {
    for (int k=0; k<3; k++) cov[k] += cov[3*n+k];
    mas[0] += mas[n];
  }

  MPI_Allreduce(MPI_IN_PLACE, &cov[0], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &mas[0], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (mas[0]>0.0) {
    for (int k=0; k<3; k++) cov[k] /= mas[0];
  }

  // Binding energy pass
  //
  pass = 1;

  exp_thread_fork(false);

  // Compbine at root node
  //
  std::vector<int> Ns, Nd;
  std::vector<double> totErg, totRad;
  int Ntot = 0.0, N = c0->Number();
  double min_erg, max_erg, rt_erg;

  if (myid) {
				// Gather number of particles on each
				// node to the root
    MPI_Gather(&N, 1, MPI_INT, NULL, 1, MPI_INT, 0, MPI_COMM_WORLD);
				// Gather energy and radius vectors
				// from each node to the root
    MPI_Gatherv(&erg[0], N, MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gatherv(&rad[0], N, MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  } else {
				// Gather number of particles on each
				// node
    Ns.resize(numprocs);
    Nd.resize(numprocs);
    MPI_Gather(&N, 1, MPI_INT, &Ns[0], 1, MPI_INT, 0, MPI_COMM_WORLD);

				// Compute totals and displacements
				// 
    for (int n=0; n<numprocs; n++) {
      Nd[n] = Ntot;
      Ntot += Ns[n];
    }
    totErg.resize(Ntot);
    totRad.resize(Ntot);

				// Gather total energy and radius vectors
				//
    MPI_Gatherv(&erg[0], N, MPI_DOUBLE, &totErg[0], &Ns[0], &Nd[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Gatherv(&rad[0], N, MPI_DOUBLE, &totRad[0], &Ns[0], &Nd[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
				// Find maximum radius with bound energy
				//

				// Energy, radius pairs for sorting
    std::vector<std::pair<double, double>> erg_rad(Ntot);

    for (int n=0; n<Ntot; n++) {
      erg_rad[n] = {totErg[n], totRad[n]};
    }
    
				// Sort the energy-radius pairs and
				// find the energy zero index
    std::sort(erg_rad.begin(), erg_rad.end());
    std::pair<double, double> zero(0.0, 0.0);
    auto it = std::lower_bound(erg_rad.begin(), erg_rad.end(), zero);
    if (it != erg_rad.begin()) it--;
    if (it == erg_rad.end())   it--;

				// Compute the average radius for the
				// target percentile

    int nbin = 100;		// [The default minimum number per bin]
    if (pcnbin) nbin = std::max<int>(nbin, pcnbin);
    else        nbin = std::max<int>(nbin, sqrt(c0->Number()));
    int edge = std::distance(erg_rad.begin(), it);
    size_t k = static_cast<size_t>(edge*pctile);
    size_t kmin = k-nbin/2, kmax = k+nbin/2;
    kmin = std::max<size_t>(kmin, 0);
    kmax = std::min<size_t>(kmax, erg_rad.size());
    rt_cur = 0.0;
    for (size_t j=kmin; j<kmax; j++) rt_cur += erg_rad[j].second;
    rt_erg = erg_rad[k].first;
    if (kmax > kmin) rt_cur /= kmax - kmin;
    else std::cout << "UserTidalRad logic error: kmax=" << kmax
		   << " !> kmin=" << kmin << " erg=" << rt_erg
		   <<  std::endl;

    radsol.push_back(rt_cur);
    rt_cur = radavg();

    // Append radial averages to file for debugging
    //
    if (debug) {
      // Open output stream for writing
      //
      std::ofstream dbg("TidalRadDebug." + runtag, ios::out | ios::app);
      if (dbg.good()) {
	const std::vector<double> pct = {0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95};
	int edge = std::distance(erg_rad.begin(), it);
	dbg << std::setw(16) << tnow;
	for (auto p : pct) {
	  size_t k = static_cast<size_t>(edge*p);
	  size_t kmin = k-nbin/2, kmax = k+nbin/2;
	  kmin = std::max<size_t>(kmin, 0);
	  kmax = std::min<size_t>(kmax, erg_rad.size());
	  double rcur=0;
	  for (size_t j=kmin; j<kmax; j++) rcur += erg_rad[j].second;
	  dbg << std::setw(16) << erg_rad[k].first
	      << std::setw(16) << rcur/(kmax - kmin);
	}
	dbg << std::setw(16) << it->first
	    << std::setw(16) << it->second
	    << std::endl;
      }
    }

    // Diagnostic values
    //
    min_erg = erg_rad.begin() ->first;
    max_erg = erg_rad.rbegin()->first;
  }

  // Finally, communicate the rt_cur value to all nodes
  //
  MPI_Bcast(&rt_cur, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Tidal radius diagnostics
  //
  if (diag>0 and myid==0 and this_step % diag==0) {
    std::ofstream out("TidalRadDiag." + runtag, ios::out | ios::app);
    if (out.good())
      out << std::setw(18) << tnow
	  << std::setw(18) << rt_cur
	  << std::setw(18) << rt_cur/rtorig
	  << std::setw(18) << rt_cur/rtorig * rfactor
	  << std::setw(18) << rt_cur/rtorig * rtrunc * rfactor
	  << std::endl;
  }

  // Update rtrunc in the component
  //
  if (dtTrunc>0.0 and tnow >= tnextT) {
    c0->rtrunc = rt_cur/rtorig * rtrunc * rfactor;
    tnextT += dtTrunc;
  }

  // Change force scale to account for change in radial scale
  //
  if (dtScale>0.0 and tnow >= tnextS) {
    c0->force->setScale(rt_cur/rtorig * rfactor);
    tnextS += dtScale;
  }

  if (myid==0) {
    // Open output stream for writing
    //
    std::ofstream out(filename.c_str(), ios::out | ios::app);
    if (out.good()) {
      out << std::left
	  << setw(16) << tnow
	  << setw(16) << rt_cur
	  << setw(16) << rt_erg
	  << setw(16) << cov[0]
	  << setw(16) << cov[1]
	  << setw(16) << cov[2]
	  << setw(16) << min_erg
	  << setw(16) << max_erg
	  << setw(10) << c0->Number()
	  << std::endl;
    } else {
      std::cerr << "UserTidalRad: error opening <" << filename
		<< "> for append" << std::endl;
    }
  }
}


void * UserTidalRad::determine_acceleration_and_potential_thread(void * arg) 
{
  int nbodies = c0->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  thread_timing_beg(id);

  PartMapItr it = cC->Particles().begin();

  for (int q=0   ; q<nbeg; q++) it++;

  for (int q=nbeg; q<nend; q++) {

    unsigned long i = (it++)->first;

    if (pass == 0) {

      double rr = 0.0;
      for (int k=0; k<3; k++) {
	double pos = c0->Pos(i, k, Component::Local | Component::Centered);
	rr += pos*pos;
      }
				// Radius for index q
      rad[q] = rr = sqrt(rr);
      
      if (rr < rt_cur) {
	mas[id] += c0->Part(i)->mass;
	for (int k=0; k<3; k++)
	  cov[3*id+k] += c0->Part(i)->mass * c0->Part(i)->vel[k];
      }

    } else {
      double v2 = 0.0;
      for (int k=0; k<3; k++) {
	double dv = c0->Part(i)->vel[k] - cov[k];
	v2 += dv*dv;
      }
				// Energy for index q
      erg[q] = 0.5*v2 + c0->Part(i)->pot;
    }
  }

  thread_timing_end(id);

  return (NULL);
}


extern "C" {
  ExternalForce *makerTidalRad(const YAML::Node& conf)
  {
    return new UserTidalRad(conf);
  }
}

class proxytidalrad { 
public:
  proxytidalrad()
  {
    factory["usertidalrad"] = makerTidalRad;
  }
};

static proxytidalrad p;
