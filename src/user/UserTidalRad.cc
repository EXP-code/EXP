#include <math.h>
#include <sstream>

#include "expand.h"
#include <localmpi.h>
#include <gaussQ.h>

#include <UserTidalRad.H>

UserTidalRad::UserTidalRad(const YAML::Node &conf) : ExternalForce(conf)
{
  id = "TidalRadius";
				// Output file name
  filename = outdir + runtag + ".TidalRad";

  comp_name = "";		// Default component for com
  rtrunc    = 1.0;		// Default tidal truncation
  rfactor   = 1.0;		// Fraction of rtruc for setting the scale
  firsttime = true;		// Used to set fiducial scale on first pass
  
  initialize();

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
      ofstream out(filename);
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
  }

  MPI_Bcast(&rt_cur, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

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
       << "using component <" << comp_name << ">" << std::endl;

  print_divider();
}

void UserTidalRad::initialize()
{
  try {
    if (conf["compname"]) comp_name = conf["compname"].as<std::string>();
    if (conf["rtrunc"])   rtrunc    = conf["rtrunc"].as<double>();
    if (conf["rfactor"])  rfactor   = conf["rfactor"].as<double>();
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


void UserTidalRad::determine_acceleration_and_potential(void)
{
  if (cC != c0) return;		// Check that this component is the target

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
  if (myid) {
    MPI_Reduce(&rad[0], NULL, c0->Number(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&erg[0], NULL, c0->Number(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  } else {
    MPI_Reduce(MPI_IN_PLACE, &rad[0], c0->Number(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &erg[0], c0->Number(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }
  

  // Find maximum radius with bound energy
  //
  double min_erg, max_erg;

  if (myid==0) {

    std::vector<std::pair<double, double>> erg_rad(c0->Number());

    for (size_t n=0; n<c0->Number(); n++) {
      erg_rad[n] = {erg[n], rad[n]};
    }

    std::sort(erg_rad.begin(), erg_rad.end());
    std::pair<double, double> zero(0.0, 0.0);
    auto it = std::lower_bound(erg_rad.begin(), erg_rad.end(), zero);
    if (it != erg_rad.begin()) it--;
    if (it == erg_rad.end())   it--;
    
    rt_erg = it->first;
    rt_cur = it->second;

    min_erg = erg_rad.begin()->first;
    max_erg = erg_rad.rbegin()->first;
  }

  MPI_Bcast(&rt_erg, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&rt_cur, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Update rtrunc in the component
  //
  c0->rtrunc = rt_cur*rfactor;

  // Update scale in force
  //
  if (firsttime) {
    rt_cur0 = rt_cur;
    firsttime = false;
  }

  // Change force scale to account for change in radial scale
  //
  c0->force->setScale(rt_cur/rt_cur0);

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

proxytidalrad p;
