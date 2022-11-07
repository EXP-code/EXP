#include <math.h>
#include "expand.H"
#include <localmpi.H>

#include <SatelliteOrbit.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>
#include <ResPotOrb.H>
#include <biorth.H>
#include <sphereSL.H>
#include <UserResPotOrb.H>
#include <BarForcing.H>

#include <sstream>

#include <pthread.h>  

#ifdef DEBUG
static pthread_mutex_t iolock = PTHREAD_MUTEX_INITIALIZER;
#endif

string respotorb_mpi_id()
{
  static bool first = true;
  static string id;

  if (first) {
    ostringstream sout;
    sout << outdir << runtag << ".respotorb_dbg." << myid;
    id = sout.str();
    first = false;
  }

  return id;
}

UserResPotOrb::UserResPotOrb(const YAML::Node& conf) : ExternalForce(conf)
{
  LMAX = 2;
  NUMR = 800;
  L0 = 2;
  M0 = 2;
  Klim = 1.0;

  ton = -20.0;			// Turn on time
  toff = 1.0e20;		// Turn off time
  delta = 1.0;			// Turn on duration

  NUMX = 400;			// Points in Ang mom grid
  NUME = 200;			// Points in Energy
  RECS = 100;			// Points in Angle grid
  ITMAX = 50;			// Number of iterations for mapping solution
  DELE = 0.001;			// Fractional offset in E grid
  DELK = 0.001;			// Offset in Kappa
  DELB = 0.001;			// Offset in Beta
  ALPHA = 0.25;			// Power law index for J distribution

  MASS = 0.05;			// Satellite mass
  AMP = 1.0;			// Mass prefactor

  first = true;
  debug = false;		// Diagnostic output

  usetag = -1;			// Flag not used unless explicitly defined

  pmass = -1.0;			// Mass for two-body diffusion
  diffuse = 0;

				// Tabled spherical model
  model_file = "SLGridSph.model";
  ctr_name = "";		// Default component for com is none

				// Orbit data file
  data_file = "orbit.data";

				// Log file name
  filename = outdir + "ResPot." + runtag;

  initialize();

  if (numRes==0)  {
    throw GenericError("You must specify at least one resonance!", __FILE__, __LINE__);
  }

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

  bool ok = true;
  if (data_file.size()>0) {
				// Open satellite orbit file
    ifstream in(data_file.c_str());
    if (!in) ok = false;
    
    double t, r, x[3], v[3];
    while (in) {
      in >> t;
      r = 0.0;
      for (int k=0; k<3; k++) {
	in >> x[k];
	r += x[k]*x[k];
      }
      r = sqrt(r);
      
      for (int k=0; k<3; k++) in >> v[k];

      if (in) {
	Time.push_back(t);
	Phase.push_back(atan2(x[1], x[0]));
	Radius.push_back(r);

	if (r>0.0)
	  Omega.push_back( (x[0]*v[1] - x[1]*v[0])/(r*r) );
	else
	  Omega.push_back( 0.0 );
      }

    }

  }
  else
    ok = false;

  if (!ok) {
    std::ostringstream sout;
    sout << "Process " << myid << ": can't open or read from <"
	 << data_file << ">";
    throw GenericError(sout.str(), __FILE__, __LINE__, 124, false);
  }

				// Set up for resonance potential
  auto hm = std::make_shared<SphericalModelTable>(model_file);
  halo_model = hm;

  double r = get_radius(tnow);
  omega = get_omega(tnow);
  Iz = MASS*r*r;

  ResPotOrb::NUMX = NUMX;
  ResPotOrb::NUME = NUME;
  ResPotOrb::RECS = RECS;
  ResPotOrb::ALPHA = ALPHA;
  ResPotOrb::ITMAX = ITMAX;
  ResPotOrb::DELTA_E = DELE;
  ResPotOrb::DELTA_K = DELK;
  ResPotOrb::DELTA_B = DELB;
				// Instantiate one for each resonance
  for (int i=0; i<numRes; i++) {
    respot.push_back(std::make_shared<ResPotOrb>(halo_model, MASS, L0, M0, L1[i], L2[i]));
  }
				// Construct debug file names
  for (int i=0; i<numRes; i++) {
    ostringstream sout;
    sout << outdir << runtag 
	 << ".respot_dbg." << L1[i] << "_" << L2[i] << "." << myid;
    respot[i]->set_debug_file(sout.str());
  }

  // Initialize two-body diffusion
  if (pmass>0.0) {
    diffuse = std::make_shared<TwoBodyDiffuse>(pmass);
    if (debug) {
      ostringstream file;
      file << outdir << "diffusion_grid." << runtag << "." << myid;
      ofstream out(file.str().c_str());
      if (out) diffuse->dump_grid(&out);
    }
  }

  btotn = vector<int>(ResPotOrb::NumDesc-1);
  difLz0 = vector<double>(numRes);
  bcount = vector< vector<int> >(nthrds);
  difLz = vector< vector<double> >(nthrds);
  for (int i=0; i<nthrds; i++) {
    difLz[i] = vector<double>(numRes);
    bcount[i] = vector<int>(ResPotOrb::NumDesc-1);
  }

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
      sout << "UserResPotOrbOrb could not open <" << fileomega << ">";
      throw GenericError(sout.str(), __FILE__, __LINE__, 103, false);
    }
    
  }

  userinfo();
}

UserResPotOrb::~UserResPotOrb()
{
  // Nothing
}

void UserResPotOrb::userinfo()
{
  if (myid) return;		// Return if node master node
  print_divider();
  cout << "** User routine RESONANCE POTENTIAL initialized";

  cout << " using satellite";
  cout << " with Mass=" << MASS 
       << ", Amp=" << AMP
       << ", Iz=" << Iz
       << ", Omega=" << omega 
       << ", Ton=" << ton
       << ", Toff=" << toff
       << ", Delta=" << delta
       << ", L=" << L0
       << ", M=" << M0
       << ", Klim=" << Klim
       << ", model=" << model_file
       << ", orbit=" << data_file;

  if (usetag>=0)
    cout << ", with bad value tagging";
  for (int ir=0; ir<numRes; ir++)
    cout << ", (l_1,l_2)_" << ir << "=(" << L1[ir] << "," << L2[ir] << ")";
  cout << ", ITMAX=" << ITMAX;

  if (pmass>0.0) cout << ", using two-body diffusion with logL=5.7 and mass=" 
		      << pmass;
  cout << endl;
  print_divider();
}

void UserResPotOrb::initialize()
{
  try {
    if (conf["LMAX"])           LMAX               = conf["LMAX"].as<int>();
    if (conf["NMAX"])           NMAX               = conf["NMAX"].as<int>();
    if (conf["NUMR"])           NUMR               = conf["NUMR"].as<int>();
    
    if (conf["L0"])             L0                 = conf["L0"].as<int>();
    if (conf["M0"])             M0                 = conf["M0"].as<int>();
    
    for (numRes=0; numRes<1000; numRes++) {
      ostringstream countl1, countl2;
      countl1 << "L1(" << numRes+1 << ")";
      countl2 << "L2(" << numRes+1 << ")";
      if (conf[countl1.str()]) {
	L1.push_back(conf[countl1.str()].as<int>());
	if (conf[countl2.str()])
	  L2.push_back(conf[countl2.str()].as<int>());
	else break;
      } else break;
    }
    
    if (L1.size() != L2.size() || numRes != (int)L1.size()) {
      std::ostingstream sout;
      sout << "UserResPotOrb: error parsing resonances, "
	   << "  Size(L1)=" << L1.size() << "  Size(L2)=" << L2.size() 
	   << "  numRes=" << numRes;
      throw GenericError(sout.str(), __FILE__, __LINE__, 119, false);
    }
    
    if (conf["Klim"])           Klim               = conf["Klim"].as<double>();
    
    if (conf["ton"])            ton                = conf["ton"].as<double>();
    if (conf["toff"])           toff               = conf["toff"].as<double>();
    if (conf["delta"])          delta              = conf["delta"].as<double>();
    
    if (conf["MASS"])           MASS               = conf["MASS"].as<double>();
    if (conf["AMP"])            AMP                = conf["AMP"].as<double>();
    if (conf["COROT"])          COROT              = conf["COROT"].as<double>();
    if (conf["A21"])            A21                = conf["A21"].as<double>();
    if (conf["A32"])            A32                = conf["A32"].as<double>();
    
    if (conf["NUMX"])           NUMX               = conf["NUMX"].as<int>();
    if (conf["NUME"])           NUME               = conf["NUME"].as<int>();
    if (conf["RECS"])           RECS               = conf["RECS"].as<int>();
    if (conf["ITMAX"])          ITMAX              = conf["ITMAX"].as<int>();
    if (conf["DELE"])           DELE               = conf["DELE"].as<double>();
    if (conf["DELK"])           DELK               = conf["DELK"].as<double>();
    if (conf["DELB"])           DELB               = conf["DELB"].as<double>();
    if (conf["ALPHA"])          ALPHA              = conf["ALPHA"].as<double>();
    
    if (conf["pmass"])          pmass              = conf["pmass"].as<double>();
    
    if (conf["model"])          model_file         = conf["model"].as<string>();
    if (conf["data"])           data_file          = conf["data"].as<string>();
    if (conf["ctrname"])        ctr_name           = conf["ctrname"].as<string>();
    if (conf["filename"])       filename           = conf["filename"].as<string>();
    if (conf["debug"])          debug              = conf["debug"].as<bool>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in UserResPotOrb: "
			   << error.what() << std::endl;
    MPI_Finalize();
    exit(-1);
  }
}

double UserResPotOrb::get_radius(double t)
{
  if (t<Time.front()) return Radius.front();
  if (t>Time.back())  return Radius.back();
  
  return odd2(t, Time, Radius, 0);
}

double UserResPotOrb::get_phase(double t)
{
  if (t<Time.front()) return Phase.front();
  if (t>Time.back())  return Phase.back();
  
  return odd2(t, Time, Phase, 0);
}

double UserResPotOrb::get_omega(double t)
{
  if (t<Time.front()) return Omega.front();
  if (t>Time.back())  return Omega.back();
  
  return odd2(t, Time, Omega, 0);
}


void UserResPotOrb::determine_acceleration_and_potential(void)
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
	  std::cerr << "UserResPotOrb: error in executing <"
		    << command << ">" << endl;
	}
	
				// Open new output stream for writing
	ofstream out(string(outdir + filename).c_str());
	if (!out) {
	  throw FileCreateError(outdir+filename, "UserResPotOrb: error opening new log file",
				__FILE__, __LINE__);
	}
	
				// Open old file for reading
	ifstream in(backupfile.c_str());
	if (!in) {
	  throw FileOpenError(backupfile, "UserResPotOrb: error opening original log file",
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
		std::ostingstream sout;
		sout << "UserResPotOrb: can't read log file, aborting" << endl;
		sout << "UserResPotOrb: line=" << line;
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
	ofstream out(string(outdir + filename).c_str(), ios::out | ios::app);
	out.setf(ios::left);
	out << setw(15) << "# Time"
	    << setw(15) << "Phase"
	    << setw(15) << "Omega"
	    << setw(15) << "dOmega(tot)"
	    << setw(15) << "Bounds(tot)";
	for (int ir=0; ir<numRes; ir++) {
	  ostringstream olab;
	  olab << "dOmega(" << L1[ir] << "," << L2[ir] << ")";
	  out << setw(15) << olab.str().c_str();
	}
	for (int j=1; j<ResPotOrb::NumDesc; j++)
	  out << setw(15) << ResPotOrb::ReturnDesc[j];
	out << endl;

	char c = out.fill('-');
	int ncnt=1;
	out << "# " << setw(13) << ncnt++;
	out << "| " << setw(13) << ncnt++;
	out << "| " << setw(13) << ncnt++;
	out << "| " << setw(13) << ncnt++;
	out << "| " << setw(13) << ncnt++;
	for (int ir=0; ir<numRes; ir++)
	  out << "| " << setw(13) << ncnt++;
      	for (int j=1; j<ResPotOrb::NumDesc; j++)
	  out << "| " << setw(13) << ncnt++;
	out << endl;
	out.fill(c);
      }
      
    }

  }

				// Store current state
  tlast = tnow;
  omlast = omega;

  phase = get_phase(tnow);
  omega = get_omega(tnow);

				// Clear bounds counter
  for (int n=0; n<nthrds; n++) {
    for (int j=0; j<ResPotOrb::NumDesc-1; j++) bcount[n][j] = 0;
  }

				// Clear difLz array
  for (int n=0; n<nthrds; n++) {
    for (int j=0; j<numRes; j++) difLz[n][j] = 0.0;
  }
  for (int j=0; j<numRes; j++) difLz0[j] = 0.0;

  // -----------------------------------------------------------
  // Compute the mapping
  // -----------------------------------------------------------

  exp_thread_fork(false);

  // -----------------------------------------------------------

				// Get total number out of bounds
  for (int j=0; j<ResPotOrb::NumDesc-1; j++) {
    for (int n=1; n<nthrds; n++)  bcount[0][j] += bcount[n][j];
    btotn[j] = 0;
  }
  MPI_Reduce(&(bcount[0][0]), &btotn[0], ResPotOrb::NumDesc-1, 
	     MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

				// Get total change in angular momentum
  for (int ir=0; ir<numRes; ir++) {
    for (int n=1; n<nthrds; n++) difLz[0][ir] += difLz[n][ir];
  }
  MPI_Allreduce(&difLz[0][0], &difLz0[0], numRes, 
		MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


  double difLzT = 0.0;
  for (int ir=0; ir<numRes; ir++) difLzT += difLz0[ir];

  double r = get_radius(tnow);
  omega = get_omega(tnow);
  Iz = MASS*r*r;
				// Write diagnostic log
  if (myid==0) {
    int btot=0;
    for (int j=0; j<ResPotOrb::NumDesc-1; j++) btot += btotn[j];
    ofstream out(string(outdir+filename).c_str(), ios::out | ios::app);
    out.setf(ios::left);
    out << setw(15) << tnow
	<< setw(15) << phase
	<< setw(15) << omega
	<< setw(15) << -difLzT/Iz
	<< setw(15) << btot;
    for (int ir=0; ir<numRes; ir++) out << setw(15) << -difLz0[ir]/Iz;
    for (int j=0; j<ResPotOrb::NumDesc-1; j++)
      out << setw(15) << btotn[j];
    out << endl;
  }

  first = false;

  print_timings("UserResPotOrb: acceleration timings");
}


void * UserResPotOrb::determine_acceleration_and_potential_thread(void * arg) 
{
  double amp, R0, R1;
  double posI[3], posO[3], velI[3], velO[3], vdif[3], Lz0, Lz1;
  
  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  thread_timing_beg(id);

  amp = AMP *
    0.5*(1.0 + erf( (tnow - ton) /delta )) *
    0.5*(1.0 + erf( (toff - tnow)/delta )) ;
    
  
  double rsat = get_radius(tnow);
  double phase = get_phase(tnow);
  double omega = get_omega(tnow);

  vector<double> Phase(3);
  Phase[0] = phase;
  Phase[1] = phase + omega*0.5*dtime;
  Phase[2] = phase + omega*dtime;

				// Check for nan (can get rid of this
				// eventually)
  bool updated;
  bool found_nan = false;
  ResPotOrb::ReturnCode ret;
  double dpot;
  int ir;

  PartMapItr it = cC->Particles().begin();
  unsigned long i;

  for (int q=0   ; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {
    i = (it++)->first;
				// If we are multistepping, compute accel 
				// only at or below this level

    if (multistep && (cC->Part(i)->level < mlevel)) continue;

    ret = ResPotOrb::OK;	// Reset error flags
    updated = false;

    if (usetag>=0 && cC->Part(i)->iattrib[usetag]) continue;

				// Initial conditions
    R0 = 0.0;
    for (int k=0; k<3; k++) {
      posI[k] = cC->Pos(i, k);
      if (c0) posI[k] -= c0->com[k];
      velI[k] = cC->Vel(i, k);

      R0 += posI[k]*posI[k];
    }
    R0 = sqrt(R0);
    
				// Initial Lz
    Lz0 = posI[0]*velI[1] - posI[1]*velI[0];
    
    ir = i % numRes;
      
    if ((ret=respot[ir]-> 
	 Update(dtime, Phase, rsat, amp, posI, velI, posO, velO)) == ResPotOrb::OK) {
	
				// Apply two-body diffusion
      if (pmass>0.0) {
	diffuse->get_diffusion(dtime, posO, velO, vdif);
	for (int k=0; k<3; k++) velO[k] += vdif[k];
      }
				// Current ang mom
      Lz1 = posO[0]*velO[1] - posO[1]*velO[0];

				// Accumulate change in Lz for each resonance
      if (respot[ir]->K()<Klim)
	difLz[id][ir] += cC->Mass(i)*(Lz1 - Lz0);
	
      updated = true;

    } else {

      bcount[id][ret-1]++;

      if (usetag>=0) cC->Part(i)->iattrib[usetag] = 1;
#ifdef DEBUG
      pthread_mutex_lock(&iolock);
      cout << "Process " << myid << " id=" << id << ":"
	   << " i=" << myid << " Error=" << ResPotOrb::ReturnDesc[ret] << endl;
      pthread_mutex_unlock(&iolock);
#endif
    }



    if (!updated) {
				// Try zero amplitude update
      if ((ret=respot[ir]-> 
	   Update(dtime, Phase, rsat, 0.0, posI, velI, posO, velO)) != ResPotOrb::OK) {
	
	dpot = halo_model->get_dpot(R0);
	for (int k=0; k<3; k++) {
	  posO[k] = posI[k] + velI[k] * dtime;
	  velO[k] = velI[k] - dpot*posI[k]/R0 * dtime;
	}
      }
    }

    R1 = 0.0;
    for (int k=0; k<3; k++) {
      cC->Part(i)->pos[k] = posO[k];
      cC->Part(i)->vel[k] = velO[k];
      cC->Part(i)->acc[k] = (velO[k] - velI[k])/dtime;
      if (!found_nan) {
	if ( std::isnan(cC->Pos(i, k)) ||
	     std::isnan(cC->Vel(i, k)) ||
	     std::isnan(cC->Acc(i, k)) ) found_nan = true; 
      }
      R1 += posO[k]*posO[k];
    }
    R1 = sqrt(R1);

    cC->Part(i)->potext = halo_model->get_pot(R1);
    
    if (found_nan) {
      cout << "Process " << myid << ": found nan\n";
      for (int k=0; k<3; k++) cout << setw(15) << cC->Pos(i, k);
      for (int k=0; k<3; k++) cout << setw(15) << cC->Vel(i, k);
      for (int k=0; k<3; k++) cout << setw(15) << cC->Acc(i, k);
      cout << endl << flush;
      found_nan = false;
    }
    
#if 1
    if (i==10) {
      ostringstream sout;
      sout << outdir << "test_orbit.respot." << myid;
      ofstream out(sout.str().c_str(), ios::app | ios::out);
      if (out) {
	out << setw(15) << tnow;
	for (int k=0; k<3; k++) out << setw(15) << posI[k];
	double v2 = 0.0;
	for (int k=0; k<3; k++) {
	  out << setw(15) << velI[k];
	  v2 += velI[k]*velI[k];
	}
	out << setw(15) << 0.5*v2 + halo_model->get_pot(R0);
	
	for (int k=0; k<3; k++) out << setw(15) << posO[k];
	v2 = 0.0;
	for (int k=0; k<3; k++) {
	  out << setw(15) << velO[k];
	  v2 += velO[k]*velO[k];
	}
	out << setw(15) << 0.5*v2 + halo_model->get_pot(R1) << endl;
      }
    }
#endif

  }

  thread_timing_end(id);

  return (NULL);
}


extern "C" {
  ExternalForce *makerResPotOrb(const YAML::Node& conf)
  {
    return new UserResPotOrb(conf);
  }
}

class proxyrporb { 
public:
  proxyrporb()
  {
    factory["userrespotorb"] = makerResPotOrb;
  }
};

static proxyrporb p;
