#include <math.h>
#include "expand.H"
#include <localmpi.H>

#include <SatelliteOrbit.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>
#include <ResPot.H>
#include <biorth.H>
#include <sphereSL.H>
#include <UserResPot.H>
#include <BarForcing.H>
#include <CircularOrbit.H>

#include <sstream>

#include <pthread.h>  
#ifdef DEBUG
static pthread_mutex_t iolock = PTHREAD_MUTEX_INITIALIZER;
#endif

string respot_mpi_id()
{
  static bool first = true;
  static string id;

  if (first) {
    ostringstream sout;
    sout << outdir << runtag << ".respot_dbg." << myid;
    id = sout.str();
    first = false;
  }

  return id;
}

UserResPot::UserResPot(const YAML::Node& conf) : ExternalForce(conf)
{
  LMAX = 2;
  NMAX = 20;
  NUMR = 800;
  L0 = 2;
  M0 = 2;
  L1 = -1;
  L2 =  2;
  rmin = 1.0e-3;
  rmax = 1.98;
  Klim = 1.0;
  scale = 0.067;
  drfac = 0.05;

  ton = -20.0;			// Turn on time
  toff = 1.0e20;		// Turn off time
  delta = 1.0;			// Turn on duration
  toffset = 0.0;		// Time offset for orbit
  omega = -1.0;			// Patern speed
  phase0 = 0.0;			// Initial phase

  NUMX = 400;			// Points in Ang mom grid
  NUME = 200;			// Points in Energy
  RECS = 100;			// Points in Angle grid
  ITMAX = 50;			// Number of iterations for mapping solution

  MASS = -1.0;			// Bar mass
  MFRAC = 0.05;			// Fraction of enclosed mass
  LENGTH = 0.067;		// Bar length
  AMP = 1.0;			// Mass prefactor
  COROT = 10;			// Corotation factor
  A21 = 0.2;			// Major to semi-minor ratio
  A32 = 0.05;			// Semi-minor to minor ratio

  domega = 0.0;			// Rate of forced slow down
  tom0 = 1.0;			// Midpoint of forced bar slow down
  dtom = -1.0;			// Width of forced bar slow down
  fileomega = "";		// File containing Omega vs T

  usebar = true;		// Use BarForcing and 
  useorb = false;		// Not CircularOrbit

  first = true;

  usetag = -1;			// Flag not used unless explicitly defined

				// Tabled spherical model
  model_file = "SLGridSph.model";
  ctr_name = "";		// Default component for com is none

				// Log file name
  filename = outdir + "ResPot." + runtag;

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

				// Set up for resonance potential
  auto hm = std::make_shared<SphericalModelTable>(model_file);
  halo_model = hm;

				// Perturbation
  if (MASS < 0.0) MASS = hm->get_mass(LENGTH);

  if (usebar) {
    BarForcing::L0 = L0;
    BarForcing::M0 = M0;
    auto bar = std::make_shared<BarForcing>(NMAX, MFRAC*MASS, LENGTH, COROT);

    bar->set_model(halo_model);
    bar->compute_quad_parameters(A21, A32);
    omega = omega0 = bar->Omega();
    Iz = bar->get_Iz();

    pert = bar;
  } else {
    auto orb = std::make_shared<CircularOrbit>(NMAX, L0, M0, MASS, LENGTH);

    if(omega <= 0.0)
      omega = omega0 = sqrt(MASS/(LENGTH*LENGTH*LENGTH));
    else
      omega0 = omega;

    Iz = MASS*LENGTH*LENGTH;

    pert = orb;
  }

  ResPot::NUMX = NUMX;
  ResPot::NUME = NUME;
  ResPot::RECS = RECS;
  ResPot::ITMAX = ITMAX;
  respot = std::make_shared<ResPot>(halo_model, pert, L0, M0, L1, L2);

				// Construct debug file name
  ostringstream sout;
  sout << outdir << runtag << ".respot_dbg." << myid;
  respot->set_debug_file(sout.str());

  btotn = vector<int>(ResPot::NumDesc-1);
  bcount = vector< vector<int> >(nthrds);
  difLz = vector<double>(nthrds);
  for (int i=0; i<nthrds; i++)
    bcount[i] = vector<int>(ResPot::NumDesc-1);
  
				// Read omega file
  if (fileomega.size()) {
    ifstream in(fileomega.c_str());
    const int sizebuf = 1024;
    char linebuf[sizebuf];

    double t, om;
    if (in) {

      while (in) {
	in.getline(linebuf, sizebuf);
	if (!in) break;
	if (linebuf[0]=='#') continue;

	istringstream sin(linebuf);
	sin >> t;
	sin >> om;
	if (sin) {
	  Time.push_back(t);
	  Omega.push_back(om);
	}
      }

    } else {
      std::ostringstream sout;
      sout << "UserResPot could not open <" << fileomega << ">";
      throw GenericError(sout.str(), __FILE__, __LINE__, 102, false);
    }
    
  }

  userinfo();
}

UserResPot::~UserResPot()
{
  // Nothing
}

void UserResPot::userinfo()
{
  if (myid) return;		// Return if node master node
  print_divider();
  cout << "** User routine RESONANCE POTENTIAL initialized";
  if (usebar) cout << " using bar";
  else cout << " using satellite";
  cout << " with Length=" << LENGTH 
       << ", Mass=" << MASS 
       << ", Mfrac=" << MFRAC 
       << ", Amp=" << AMP
       << ", Iz=" << Iz
       << ", Omega=" << omega 
       << ", Domega=" << domega 
       << ", tom0=" << tom0
       << ", dtom=" << dtom
       << ", b/a=" << A21
       << ", c/b=" << A32
       << ", Ton=" << ton
       << ", Toff=" << toff
       << ", Delta=" << delta
       << ", L=" << L0
       << ", M=" << M0
       << ", L1=" << L1
       << ", L2=" << L2
       << ", Klim=" << Klim;
  if (fileomega.size())
    cout << ", using table <" << fileomega << "> for Omega(t)";
  if (dtom>0) cout << ", T_om=" << tom0 << ", dT_om=" << dtom;
  cout << ", Domega=" << domega;
  if (usetag>=0)
    cout << ", with bad value tagging";
  cout << ", ITMAX=" << ITMAX << "\n";
  print_divider();
}

void UserResPot::initialize()
{
  try {
    if (conf["LMAX"])           LMAX               = conf["LMAX"].as<int>();
    if (conf["NMAX"])           NMAX               = conf["NMAX"].as<int>();
    if (conf["NUMR"])           NUMR               = conf["NUMR"].as<int>();
    
    if (conf["L0"])             L0                 = conf["L0"].as<int>();
    if (conf["M0"])             M0                 = conf["M0"].as<int>();
    if (conf["L1"])             L1                 = conf["L1"].as<int>();
    if (conf["L2"])             L2                 = conf["L2"].as<int>();
    
    if (conf["rmin"])           rmin               = conf["rmin"].as<double>();
    if (conf["rmax"])           rmax               = conf["rmax"].as<double>();
    if (conf["Klim"])           Klim               = conf["Klim"].as<double>();
    if (conf["scale"])          scale              = conf["scale"].as<double>();
    if (conf["drfac"])          drfac              = conf["drfac"].as<double>();
    
    if (conf["ton"])            ton                = conf["ton"].as<double>();
    if (conf["toff"])           toff               = conf["toff"].as<double>();
    if (conf["delta"])          delta              = conf["delta"].as<double>();
    if (conf["toffset"])        toffset            = conf["toffset"].as<double>();
    if (conf["phase0"])         phase0             = conf["phase0"].as<double>();
    
    if (conf["MASS"])           MASS               = conf["MASS"].as<double>();
    if (conf["MFRAC"])          MFRAC              = conf["MFRAC"].as<double>();
    if (conf["LENGTH"])         LENGTH             = conf["LENGTH"].as<double>();
    if (conf["AMP"])            AMP                = conf["AMP"].as<double>();
    if (conf["COROT"])          COROT              = conf["COROT"].as<double>();
    if (conf["A21"])            A21                = conf["A21"].as<double>();
    if (conf["A32"])            A32                = conf["A32"].as<double>();
    
    if (conf["NUMX"])           NUMX               = conf["NUMX"].as<int>();
    if (conf["NUME"])           NUME               = conf["NUME"].as<int>();
    if (conf["RECS"])           RECS               = conf["RECS"].as<int>();
    if (conf["ITMAX"])          ITMAX              = conf["ITMAX"].as<int>();
    
    if (conf["omega"])          omega              = conf["omega"].as<double>();
    if (conf["domega"])         domega             = conf["domega"].as<double>();
    if (conf["tom0"])           tom0               = conf["tom0"].as<double>();
    if (conf["dtom"])           dtom               = conf["dtom"].as<double>();
    
    
    if (conf["model"])          model_file         = conf["model"].as<string>();
    if (conf["ctrname"])        ctr_name           = conf["ctrname"].as<string>();
    if (conf["filename"])       filename           = conf["filename"].as<string>();
    if (conf["fileomega"])      fileomega          = conf["fileomega"].as<string>();
    if (conf["usetag"])         usetag             = conf["usetag"].as<int>();
    
    if (conf["usebar"])
      {
	usebar = conf["usebar"].as<bool>();
	useorb = !usebar;
      }
    if (conf["useorb"])   
      {
	useorb = conf["useorb"].as<bool>();
	usebar = !useorb;
      }
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in UserResPot: "
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
  
  
double UserResPot::get_omega(double t)
{
  if (t<Time.front()) return Omega.front();
  if (t>Time.back())  return Omega.back();

  return odd2(t, Time, Omega, 0);
}


void UserResPot::determine_acceleration_and_potential(void)
{
  if (multistep and mlevel>0) return;

#if HAVE_LIBCUDA==1		// Cuda compatibility
  getParticlesCuda(cC);
#endif

  if (first) {

    if (restart) {

      if (myid == 0) {
				// Backup up old file
	string backupfile = outdir + filename + ".bak";
	string command("cp ");
	command += outdir + filename + " " + backupfile;
	if (system(command.c_str()) == -1) {
	  std::cerr << "UserResPot: error in executing <"
		    << command << ">" << endl;
	}
	
				// Open new output stream for writing
	ofstream out(string(outdir+filename).c_str());
	if (!out) {
	  throw FileCreateError(outdir+filename, "UserResPotOrb: error opening new log file",
				__FILE__, __LINE__);
	}
	
				// Open old file for reading
	ifstream in(backupfile.c_str());
	if (!in) {
	  throw FileCreateError(backupfile, "UserResPotOrb: error opening original log file",
				__FILE__, __LINE__);
	}

	const int linesize = 1024;
	char line[linesize];
	
	in.getline(line, linesize); // Discard header
	in.getline(line, linesize); // Next line

	double phase1, omlast1, tlast1;
	bool firstline = true;

	while (in) {

	  if (line[0] != '#' && line[0] != '!') {

	    istringstream ins(line);

	    ins >> tlast1;
	    ins >> phase1;
	    ins >> omlast1;

	    if (tlast1 >= tnow) {
	      if (firstline) {
		std::ostringstream sout;
		sout << "UserResPot: can't read log file, aborting" << endl;
		sout << "UserResPot: line=" << line;
		throw GenericError(sout.str(), __FILE__, __LINE__, 123, false);
	      }
	      break;
	    }

	    firstline = false;
	    tlast = tlast1;
	    phase = phase1;
	    omlast = omlast1;

	  }

	  out << line << "\n";
	  
	  in.getline(line, linesize); // Next line
	}
				// Trapezoidal rule step
	phase += (tnow - tlast)*0.5*(omega + omlast);
      }

      MPI_Bcast(&tlast, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&phase, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&omega, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    if (!restart) {

      if (myid==0) {		// Write header
	ofstream out(string(outdir+filename).c_str(), ios::out | ios::app);
	out.setf(ios::left);
	out << setw(15) << "# Time"
	    << setw(15) << "Phase"
	    << setw(15) << "Omega"
	    << setw(15) << "dOmega(tot)"
	    << setw(15) << "Bounds"
	    << endl;
	for (int j=1; j<ResPot::NumDesc; j++)
	  out << setw(15) << ResPot::ReturnDesc[j];
	out << endl;
	
	char c = out.fill('-');
	int ncnt=1;
	out << "# " << setw(13) << ncnt++;
	out << "| " << setw(13) << ncnt++;
	out << "| " << setw(13) << ncnt++;
	out << "| " << setw(13) << ncnt++;
	out << "| " << setw(13) << ncnt++;
      	for (int j=1; j<ResPot::NumDesc; j++)
	  out << "| " << setw(13) << ncnt++;
	out << endl;
	out.fill(c);
      }
      
      phase = phase0;	// Initial phase 
    }

  } else {			
				// Trapezoidal rule integration
    phase += (tnow - tlast)*0.5*(omega + omlast);
  }

				// Store current state
  tlast = tnow;
  omlast = omega;

				// Clear bounds counter
  for (int n=0; n<nthrds; n++) {
    for (int j=0; j<ResPot::NumDesc-1; j++) bcount[n][j] = 0;
  }

				// Clear difLz array
  for (int n=0; n<nthrds; n++) difLz[n] = 0.0;
  difLz0 = 0.0;

  // -----------------------------------------------------------
  // Compute the mapping
  // -----------------------------------------------------------

  exp_thread_fork(false);

  // -----------------------------------------------------------

				// Get total number out of bounds
  for (int n=1; n<nthrds; n++) {
    for (int j=0; j<ResPot::NumDesc-1; j++) {
      bcount[0][j] += bcount[n][j];
      btotn[j] = 0;
    }
  }
  MPI_Reduce(&(bcount[0][0]), &btotn[0], ResPot::NumDesc-1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

				// Get total change in angular momentum
  difLz0 = 0.0;
  for (int n=1; n<nthrds; n++) difLz[0] += difLz[n];
  MPI_Allreduce(&difLz[0], &difLz0, 1,
		MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (fileomega.size())
    omega = get_omega(tnow);
  else if (dtom>0.0)
    omega = omega0*(1.0 + domega*0.5*(1.0 + erf( (tnow - tom0)/dtom )));
  else
    omega = omega0*(1.0 + domega*(tnow - tom0*0.5));

				// Write diagnostic log
  if (myid==0) {
    int btot=0;
    for (int j=0; j<ResPot::NumDesc-1; j++) btot += btotn[j];
    ofstream out(string(outdir+filename).c_str(), ios::out | ios::app);
    out.setf(ios::left);
    out << setw(15) << tnow
	<< setw(15) << phase
	<< setw(15) << omega
	<< setw(15) << -difLz0/Iz
	<< setw(15) << btot;
    for (int j=0; j<ResPot::NumDesc-1; j++)
      out << setw(15) << btotn[j];
    out << endl;
  }

  first = false;

  print_timings("UserResPot: acceleration timings");
}
void * UserResPot::determine_acceleration_and_potential_thread(void * arg) 
{
  double amp, R2, R;
  double posI[3], posO[3], velI[3], velO[3], Lz0, Lz1;
  
  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  thread_timing_beg(id);

  amp = AMP *
    0.5*(1.0 + erf( (tnow - ton) /delta )) *
    0.5*(1.0 + erf( (toff - tnow)/delta )) ;
    
  
  vector<double> Phase(3);
  Phase[0] = phase;
  Phase[1] = phase + omega*0.5*dtime;
  Phase[2] = phase + omega*dtime;

				// Check for nan (can get rid of this
				// eventually)
  bool updated;
  bool found_nan = false;
  ResPot::ReturnCode ret;
  double dpot;

  for (int i=nbeg; i<nend; i++) {

				// If we are multistepping, compute accel 
				// only at or below this level

    if (multistep && (cC->Part(i)->level < mlevel)) continue;

    ret = ResPot::OK;		// Reset error flags
    updated = false;

    if (usetag>=0 && cC->Part(i)->iattrib[usetag]) continue;

				// Initial conditions
    for (int k=0; k<3; k++) {
      posI[k] = cC->Pos(i, k);
      if (c0) posI[k] -= c0->com[k];
      velI[k] = cC->Vel(i, k);
    }
    
				// Initial Lz
    Lz0 = posI[0]*velI[1] - posI[1]*velI[0];
    
    
    R2 = 0.0;
    for (int k=0; k<3; k++) {
      R2 += posI[k]*posI[k];
    }
    R = sqrt(R2);

    if (R>rmin && R<rmax) {

      if ((ret=respot-> 
	   Update(dtime, Phase, amp, posI, velI, posO, velO)) == ResPot::OK) {
	
				// Current ang mom
	Lz1 = posO[0]*velO[1] - posO[1]*velO[0];

				// Accumulate change in Lz for each resonance
	if (respot->K()<Klim)
	  difLz[id] += cC->Mass(i)*(Lz1 - Lz0);
	
	updated = true;

      } else {

	if (usetag>=0) cC->Part(i)->iattrib[usetag] = 1;
#ifdef DEBUG
	pthread_mutex_lock(&iolock);
	cout << "Process " << myid << " id=" << id << ":"
	     << " i=" << myid << " Error=" << ResPot::ReturnDesc[ret] << endl;
	pthread_mutex_unlock(&iolock);
#endif
      }

    }
     
    if (!updated) {
      dpot = halo_model->get_dpot(R);
      for (int k=0; k<3; k++) {
	posO[k] = velI[k] * dtime;
	if (R>0.01*rmin) velO[k] = -dpot*posI[k]/R * dtime;
      }
    }

    for (int k=0; k<3; k++) {
      cC->Part(i)->pos[k] = posO[k];
      cC->Part(i)->vel[k] = velO[k];
      cC->Part(i)->acc[k] = (velO[k] - velI[k])/dtime;
      if (!found_nan) {
	if ( std::isnan(cC->Pos(i, k)) ||
	     std::isnan(cC->Vel(i, k)) ||
	     std::isnan(cC->Acc(i, k)) ) found_nan = true; 
      }
    }
    cC->Part(i)->potext = halo_model->get_pot(R);
    
    if (found_nan) {
      cout << "Process " << myid << ": found nan\n";
      for (int k=0; k<3; k++) cout << setw(15) << cC->Pos(i, k);
      for (int k=0; k<3; k++) cout << setw(15) << cC->Vel(i, k);
      for (int k=0; k<3; k++) cout << setw(15) << cC->Acc(i, k);
      cout << endl << flush;
      found_nan = false;
    }
    
  }

  thread_timing_end(id);

  return (NULL);
}


extern "C" {
  ExternalForce *makerResPot(const YAML::Node& conf)
  {
    return new UserResPot(conf);
  }
}

class proxyrespot { 
public:
  proxyrespot()
  {
    factory["userrespot"] = makerResPot;
  }
};

static proxyrespot p;
