#include "expand.H"
#include <localmpi.H>

#include <SatelliteOrbit.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>
#include <ResPot.H>
#include <biorth.H>
#include <sphereSL.H>
#include <UserRPtest.H>
#include <BarForcing.H>

#include <sstream>
#include <memory>
#include <cmath>

const std::set<std::string>
UserRPtest::valid_keys = {
  "L0",
  "M0",
  "L1",
  "L2",
  "rmin",
  "rmax",
  "scale",
  "NUMX",
  "NUME",
  "RECS",
  "with_ps",
  "npart",
  "model",
  "ctrname",
  "filename"
};

UserRPtest::UserRPtest(const YAML::Node& conf) : ExternalForce(conf)
{
  LMAX = 2;
  NMAX = 20;
  NUMR = 1000;
  L0 = 2;
  M0 = 2;
  L1 = -1;
  L2 =  2;
  rmin = 1.0e-4;
  rmax = 1.95;
  scale = 0.067;

  NUMX = 400;			// Points in Ang mom grid
  NUME = 100;			// Points in Energy grid
  RECS = 100;			// Points in Angle grid

  with_ps = false;		// Don't print phase space for each particle
  npart = 5;			// Number of particles to trace

  first = true;

				// Tabled spherical model
  model_file = "SLGridSph.model";
  ctr_name = "";		// Default component for com is none

				// Log file name
  filename = outdir + "RPtest." + runtag;

  initialize();

  if (ctr_name.size()>0) {
				// Look for the fiducial component for
				// centering
    bool found = false;
    for (auto c : comp->components) {
      if ( !ctr_name.compare(c->name) ) {
	c0 = c;
	found = true;
      break;
      }
    }

    if (!found) {
      cerr << "Process " << myid << ": can't find desired component <"
	   << ctr_name << ">" << endl;
    }

  }
  else
    c0 = NULL;


				// Perturbation
  BarForcing::L0 = L0;
  BarForcing::M0 = M0;
  auto bar = std::make_shared<BarForcing>(NMAX, 0.1, 0.1, 1.0);
  bar->compute_quad_parameters(0.2, 0.2);

				// Set up for resonance potential
  auto hm = std::make_shared<SphericalModelTable>(model_file);
  halo_model = hm;

  ResPot::NUMX = NUMX;
  ResPot::NUME = NUME;
  ResPot::RECS = RECS;
  respot = std::make_shared<ResPot>(halo_model, bar, L0, M0, L1, L2);

  userinfo();
}

UserRPtest::~UserRPtest()
{
  // Nothing
}

void UserRPtest::userinfo()
{
  if (myid) return;		// Return if node master node
  print_divider();
  cout << "** User routine RESONANCE POTENTIAL TEST initialized";
  cout << ", L=" << L0
       << ", M=" << M0
       << ", l_1=" << L1
       << ", l_2=" << L2
       << "\n";
  print_divider();
}

void UserRPtest::initialize()
{
  // Check for unmatched keys
  //
  auto unmatched = YamlCheck(conf, valid_keys);
  if (unmatched.size())
    throw YamlConfigError("UserRPtest", "parameter", unmatched,
			  __FILE__, __LINE__);

  // Assign values from YAML
  //
  try {
    if (conf["L0"])             L0                 = conf["L0"].as<int>();
    if (conf["M0"])             M0                 = conf["M0"].as<int>();
    if (conf["L1"])             L1                 = conf["L1"].as<int>();
    if (conf["L2"])             L2                 = conf["L2"].as<int>();
    
    if (conf["rmin"])           rmin               = conf["rmin"].as<double>();
    if (conf["rmax"])           rmax               = conf["rmax"].as<double>();
    if (conf["scale"])          scale              = conf["scale"].as<double>();
    
    if (conf["NUMX"])           NUMX               = conf["NUMX"].as<int>();
    if (conf["NUME"])           NUME               = conf["NUME"].as<int>();
    if (conf["RECS"])           RECS               = conf["RECS"].as<int>();
    
    if (conf["with_ps"])        with_ps            = conf["with_ps"].as<bool>();
    if (conf["npart"])          npart              = conf["npart"].as<int>();
    
    if (conf["model"])          model_file         = conf["model"].as<string>();
    if (conf["ctrname"])        ctr_name           = conf["ctrname"].as<string>();
    if (conf["filename"])       filename           = conf["filename"].as<string>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in UserRPtest: "
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

void UserRPtest::determine_acceleration_and_potential(void)
{
  if (first) {
    
#if HAVE_LIBCUDA==1		// Cuda compatibility
  getParticlesCuda(cC);
#endif

    if (restart) {
      
      if (myid == 0) {
	// Backup up old file
	string curfile = outdir + filename;
	string backupfile = curfile + ".bak";
	string command("cp ");
	command += curfile + " " + backupfile;
	if (system(command.c_str()) == -1) {
	  std::cerr << "UserRPtest: error in executing <"
		    << command << ">" << endl;
	}
	
	// Open new output stream for writing
	ofstream out(curfile.c_str());
	if (!out) {
	  std::ostringstream sout;
	  sout << "UserRPtest: error opening new log file <" 
	       << curfile << "> for writing";
	  throw GenericError(sout.str(), __FILE__, __LINE__);
	}
	
	// Open old file for reading
	ifstream in(backupfile.c_str());
	if (!in) {
	  std::ostringstream sout;
	  sout << "UserRPtest: error opening original log file <" 
	       << backupfile << "> for reading";
	  throw GenericError(sout.str(), __FILE__, __LINE__);
	}
	
	const int linesize = 1024;
	char line[linesize];
	
	in.getline(line, linesize); // Discard header
	in.getline(line, linesize); // Next line
	
	double tlast1;
	bool firstline = true;

	while (in) {
	  istringstream ins(line);

	  ins >> tlast1;

	  if (tlast1 >= tnow) {
	    if (firstline) {
	      std::string msg = "UserRPtest: can't read log file, aborting";
	      throw GenericError(msg, __FILE__, __LINE__, 123, false);
	    }
	    break;
	  }

	  firstline = false;

	  out << line << "\n";

	  in.getline(line, linesize); // Next line
	}
      }

    }

    if (!restart) {

      if (myid == 0) {		// Write header
	ofstream out(string(outdir+filename).c_str(), ios::out | ios::app);
	out.setf(ios::left);
	out << setw(15) << "# Time";

	const int nlabels = 8;
	string labels[] = {"E", "K", "w1", "w2", "w3", "f", "beta", "psi"};

	const int nlabels2 = 6;
	string labels2[] = {"X", "Y", "Z", "U", "V", "W"};

	for (unsigned i=1; i<=npart; i++) {
	  for (int j=0; j<nlabels; j++) {
	    ostringstream sout;
	    sout << labels[j] << "(" << i << ")";
	    out << "| " << setw(13) << sout.str();
	  }
	  if (with_ps) {
	    for (int j=0; j<nlabels2; j++) {
	      ostringstream sout;
	      sout << labels2[j] << "(" << i << ")";
	      out << "| " << setw(13) << sout.str();
	    }
	  }
	}
	out << endl;


	
	out << setw(15) << "# 1";
	int cntr = 2;

	for (unsigned i=1; i<=npart; i++) {
	  for (int j=0; j<nlabels; j++) {
	    out << "| " << setw(13) << cntr++;
	  }
	  if (with_ps) {
	    for (int j=0; j<nlabels2; j++) {
	    out << "| " << setw(13) << cntr++;
	    }
	  }
	}
	cout << endl;
	
      }
      
    }

    first = false;

  }

  exp_thread_fork(false);
}

void * UserRPtest::determine_acceleration_and_potential_thread(void * arg) 
{
  double pos[3], vel[3];
  double R2, R, pot, dpot;
  double E, K, I1, J, O1, O2, w1, w2, w3, f, beta, psi;
  
  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  ofstream out;

  if (myid==0 && id==0) {
    out.open(string(outdir+filename).c_str(), ios::out | ios::app);
    out.setf(ios::left);
  }


  PartMapItr it = cC->Particles().begin();
  unsigned long i;

  for (int q=0   ; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {
    i = (it++)->first;
				// If we are multistepping, compute accel 
				// only at or below this level

    if (multistep && (cC->Part(i)->level < mlevel)) continue;

    R2 = 0.0;
    for (int k=0; k<3; k++) {
      pos[k] = cC->Pos(i, k);
      if (c0) pos[k] -= c0->com[k];
      vel[k] = cC->Vel(i, k);
      R2 += pos[k]*pos[k];
    }
    R = sqrt(R2);

    halo_model->get_pot_dpot(R, pot, dpot);

    if (myid==0 && id==0 && mlevel==0 && i<npart) {
      if (i==0) out << setw(15) << tnow;
      respot->coord(pos, vel, E, K, I1, J, O1, O2, w1, w2, w3, f, beta, psi);
      out << setw(15) << E
	  << setw(15) << K
	  << setw(15) << w1
	  << setw(15) << w2
	  << setw(15) << w3
	  << setw(15) << f
	  << setw(15) << beta
	  << setw(15) << psi;

      if (with_ps) {
	for (int k=0; k<3; k++) out << setw(15) << pos[k];
	for (int k=0; k<3; k++) out << setw(15) << vel[k];
      }
      if (i==npart-1) out << endl;
    }

    for (int k=0; k<3; k++) cC->AddAcc(i, k, -dpot*pos[k]/R );
    
    cC->AddPotExt(i, cC->Mass(i) * pot );

  }

  if (myid==0 && id==0) {
    out.close();
  }

  return (NULL);
}


extern "C" {
  ExternalForce *makerRPtest(const YAML::Node& conf)
  {
    return new UserRPtest(conf);
  }
}

class proxyP { 
public:
  proxyP()
  {
    factory["userrptest"] = makerRPtest;
   }
};

static proxyP p;
