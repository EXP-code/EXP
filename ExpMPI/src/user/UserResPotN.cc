#include <math.h>
#include "expand.h"
#include <localmpi.h>

#include <SatelliteOrbit.h>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>
#include <ResPot.H>
#include <biorth.h>
#include <sphereSL.h>
#include <UserResPotN.H>
#include <BarForcing.H>

#include <sstream>

#include <pthread.h>  
static pthread_mutex_t iolock = PTHREAD_MUTEX_INITIALIZER;


static int respot_mpi_id_var;
int respot_mpi_id()
{
  return respot_mpi_id_var;
}

UserResPotN::UserResPotN(string &line) : ExternalForce(line)
{
  LMAX = 2;
  NUMR = 800;
  L0 = 2;
  M0 = 2;
  rmin = 1.0e-3;
  rmax = 1.98;
  Klim = 1.0;
  scale = 0.067;
  drfac = 0.05;

  ton = -20.0;			// Turn on time
  toff = 1.0e20;		// Turn off time
  delta = 1.0;			// Turn on duration
  toffset = 0.0;		// Time offset for orbit
  omega = 18.9;			// Patern speed
  phase0 = 0.0;			// Initial phase

  NUMX = 400;			// Points in Ang mom grid
  NUME = 200;			// Points in Energy
  RECS = 100;			// Points in Angle grid
  ITMAX = 50;			// Number of iterations for mapping solution

  MASS = 0.05;			// Bar mass
  LENGTH = 0.067;		// Bar length
  AMP = 1.0;			// Mass prefactor
  COROT = 10;			// Corotation factor
  A21 = 0.2;			// Major to semi-minor ratio
  A32 = 0.05;			// Semi-minor to minor ratio

  self = true;			// Self consistent slow down
  domega = 0.0;			// Rate of forced slow down
  tom0 = 1.0;			// Midpoint of forced bar slow down
  dtom = -1.0;			// Width of forced bar slow down


  first = true;
  debug = false;

  usetag = -1;			// Flag not used unless explicitly defined

				// Tabled spherical model
  model_file = "SLGridSph.model";
  ctr_name = "";		// Default component for com is none

				// Log file name
  filename = outdir + "ResPot." + runtag;

  initialize();

  if (numRes==0)  {
    if (myid==0) cerr << "You must specify at least one resonance!\n";
    MPI_Abort(MPI_COMM_WORLD, 120);
    exit(0);
  }

  if (ctr_name.size()>0) {
				// Look for the fiducial component for
				// centering
    bool found = false;
    list<Component*>::iterator cc;
    Component *c;
    for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
      c = *cc;
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
  SphericalModelTable *hm = new SphericalModelTable(model_file);
  halo_model = hm;

				// Perturbation
  BarForcing::L0 = L0;
  BarForcing::M0 = M0;
  BarForcing *bar = new BarForcing(NMAX, MASS, LENGTH, COROT);
  bar->set_model(halo_model);
  bar->compute_quad_parameters(A21, A32);
  omega = omega0 = bar->Omega();
  Iz = bar->get_Iz();

  pert = bar;


  ResPot::NUMX = NUMX;
  ResPot::NUME = NUME;
  ResPot::RECS = RECS;
  ResPot::ITMAX = ITMAX;

  for (int i=0; i<numRes; i++) {
    respot.push_back(new ResPot(halo_model, pert, L0, M0, L1[i], L2[i]));
  }

  bcount = vector<int>(nthrds);
  difLz = vector< vector<double> >(nthrds);
  for (int i=0; i<nthrds; i++) difLz[i] = vector<double>(numRes);
  difLz0 = vector<double>(numRes);
  difLNP = vector<double>(nthrds);

  userinfo();

  respot_mpi_id_var = myid;
}

UserResPotN::~UserResPotN()
{
  for (int i=0; i<numRes; i++) delete respot[i];
  delete halo_model;
  delete pert;
}

void UserResPotN::userinfo()
{
  if (myid) return;		// Return if node master node
  print_divider();
  cout << "** User routine RESONANCE POTENTIAL initialized";
  cout << " with Length=" << LENGTH 
       << ", Mass=" << MASS 
       << ", Amp=" << AMP
       << ", Iz=" << Iz
       << ", Omega=" << omega 
       << ", b/a=" << A21
       << ", c/b=" << A32
       << ", Ton=" << ton
       << ", Toff=" << toff
       << ", Delta=" << delta
       << ", L=" << L0
       << ", M=" << M0
       << ", Klim=" << Klim;
  if (self)  cout << ", with self-consistent slow down";
  else {
    if (dtom>0) cout << ", T_om=" << tom0 << ", dT_om=" << dtom;
    cout << ", Domega=" << domega;
  }
  if (debug) cout << ", with debug ang mom";
  if (usetag>=0)
    cout << ", with bad value tagging";
  for (int ir=0; ir<numRes; ir++)
    cout << ", (l_1,l_2)_" << ir << "=(" << L1[ir] << "," << L2[ir] << ")";
  cout << ", ITMAX=" << ITMAX << "\n";
  print_divider();
}

void UserResPotN::initialize()
{
  string val;

  if (get_value("LMAX", val))     LMAX = atoi(val.c_str());
  if (get_value("NMAX", val))     NMAX = atoi(val.c_str());
  if (get_value("NUMR", val))     NUMR = atoi(val.c_str());

  if (get_value("L0", val))       L0 = atoi(val.c_str());
  if (get_value("M0", val))       M0 = atoi(val.c_str());

  for (numRes=0; numRes<1000; numRes++) {
    ostringstream countl1, countl2;
    countl1 << "L1(" << numRes+1 << ")";
    countl2 << "L2(" << numRes+1 << ")";
    if (get_value(countl1.str(), val)) {
      L1.push_back(atoi(val.c_str()));
      if (get_value(countl2.str(), val))
	L2.push_back(atoi(val.c_str()));
      else break;
    } else break;
  }

  if (L1.size() != L2.size() || numRes != (int)L1.size()) {
    cerr << "UserResPotN: error parsing resonances, "
	 << "  Size(L1)=" << L1.size() << "  Size(L2)=" << L2.size() 
	 << "  numRes=" << numRes << endl;
    MPI_Abort(MPI_COMM_WORLD, 119);
  }

  if (get_value("rmin", val))     rmin = atof(val.c_str());
  if (get_value("rmax", val))     rmax = atof(val.c_str());
  if (get_value("Klim", val))     Klim = atof(val.c_str());
  if (get_value("scale", val))    scale = atof(val.c_str());
  if (get_value("drfac", val))    drfac = atof(val.c_str());

  if (get_value("ton", val))      ton = atof(val.c_str());
  if (get_value("toff", val))     toff = atof(val.c_str());
  if (get_value("delta", val))    delta = atof(val.c_str());
  if (get_value("toffset", val))  toffset = atof(val.c_str());
  if (get_value("phase0", val))   phase0 = atof(val.c_str());

  if (get_value("MASS", val))     MASS = atof(val.c_str());
  if (get_value("LENGTH", val))   LENGTH = atof(val.c_str());
  if (get_value("AMP", val))      AMP = atof(val.c_str());
  if (get_value("COROT", val))    COROT = atof(val.c_str());
  if (get_value("A21", val))      A21 = atof(val.c_str());
  if (get_value("A32", val))      A32 = atof(val.c_str());

  if (get_value("NUMX", val))     NUMX = atoi(val.c_str());
  if (get_value("NUME", val))     NUME = atoi(val.c_str());
  if (get_value("RECS", val))     RECS = atoi(val.c_str());
  if (get_value("ITMAX", val))    ITMAX = atoi(val.c_str());
  
  if (get_value("self", val))     self = atoi(val.c_str());
  if (get_value("domega", val))   domega = atof(val.c_str());
  if (get_value("tom0", val))     tom0 = atof(val.c_str());
  if (get_value("dtom", val))     dtom = atof(val.c_str());


  if (get_value("model", val))    model_file = val;
  if (get_value("ctrname", val))  ctr_name = val;
  if (get_value("filename", val)) filename = val;
  if (get_value("usetag", val))   usetag = atoi(val.c_str());
  if (get_value("debug", val))    debug = atoi(val.c_str());
}

void UserResPotN::determine_acceleration_and_potential(void)
{

  if (first) {

    if (restart) {

      if (myid == 0) {
				// Backup up old file
	string backupfile = filename + ".bak";
	string command("cp ");
	command += filename + " " + backupfile;
	system(command.c_str());
	
				// Open new output stream for writing
	ofstream out(filename.c_str());
	if (!out) {
	  cout << "UserResPotN: error opening new log file <" 
	       << filename << "> for writing\n";
	  MPI_Abort(MPI_COMM_WORLD, 121);
	  exit(0);
	}
	
				// Open old file for reading
	ifstream in(backupfile.c_str());
	if (!in) {
	  cout << "UserResPotN: error opening original log file <" 
	       << backupfile << "> for reading\n";
	  MPI_Abort(MPI_COMM_WORLD, 122);
	  exit(0);
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

	    if (tlast1 >= tpos) {
	      if (firstline) {
		cerr << "UserResPotN: can't read log file, aborting" << endl;
		cerr << "UserResPotN: line=" << line << endl;
		MPI_Abort(MPI_COMM_WORLD, 123);
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
	ofstream out(filename.c_str(), ios::out | ios::app);
	out.setf(ios::left);
	out << setw(15) << "# Time"
	    << setw(15) << "Phase"
	    << setw(15) << "Omega"
	    << setw(15) << "dOmega(tot)"
	    << setw(15) << "Bounds";
	for (int ir=0; ir<numRes; ir++) {
	  ostringstream olab;
	  olab << "dOmega(" << L1[ir] << "," << L2[ir] << ")";
	  out << setw(15) << olab.str().c_str();
	}
	if (debug) out << setw(15) << "dLz(NP)";
	out << endl;

	char c = out.fill('-');
	int ncnt=1;
	out << "# " << setw(13) << ncnt++ 
	    << "| " << setw(13) << ncnt++
	    << "| " << setw(13) << ncnt++
	    << "| " << setw(13) << ncnt++
	    << "| " << setw(13) << ncnt++;
	for (int ir=0; ir<numRes; ir++) {
	  out << "| " << setw(13) << ncnt++;
	}
	if (debug) {
	  out << "| " << setw(13) << ncnt++;
	}
	out << endl;
	out.fill(c);
      }
      
      phase = phase0;	// Initial phase 
    }

  } else {			// Trapezoidal rule integration
    phase += (tnow - tlast)*0.5*(omega + omlast);
  }

				// Store current state
  tlast = tnow;
  omlast = omega;

				// Clear bounds counter
  for (int n=0; n<nthrds; n++) bcount[n] = 0;

				// Clear difLz array
  for (int n=0; n<nthrds; n++) {
    for (int j=0; j<numRes; j++) difLz[n][j] = 0.0;
    difLNP[n] = 0.0;
  }
  for (int j=0; j<numRes; j++) difLz0[j] = 0.0;

  // -----------------------------------------------------------
  // Compute the mapping
  // -----------------------------------------------------------

  exp_thread_fork(false);

  // -----------------------------------------------------------

				// Get total number out of bounds
  int btot=0;
  for (int n=1; n<nthrds; n++) bcount[0] += bcount[n];
  MPI_Reduce(&bcount[0], &btot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

				// Get total change in angular momentum
  for (int ir=0; ir<numRes; ir++) {
    for (int n=1; n<nthrds; n++) difLz[0][ir] += difLz[n][ir];
  }
  MPI_Allreduce(&difLz[0][0], &difLz0[0], numRes, 
		MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


				// Total change in ang mom *without* pert
				// for debugging
  double difLNP0 = 0.0;
  if (debug) {
    for (int n=1; n<nthrds; n++) difLNP[0] += difLNP[n];
    MPI_Reduce(&difLNP[0], &difLNP0, 1, 
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }

  double difLzT = 0.0;
  for (int ir=0; ir<numRes; ir++) difLzT += difLz0[ir];

  if (self)
    omega -= difLzT/Iz;
  else {
    if (dtom>0.0)
      omega = omega0*(1.0 + domega*0.5*(1.0 + erf( (tnow - tom0)/dtom )));
    else
      omega += omega0*domega*dtime*
	0.5*(1.0 + erf( (tnow - ton) /delta )) *
	0.5*(1.0 + erf( (toff - tnow)/delta )) ;
  }

				// Write diagnostic log
  if (myid==0) {
    ofstream out(filename.c_str(), ios::out | ios::app);
    out.setf(ios::left);
    out << setw(15) << tnow
	<< setw(15) << phase
	<< setw(15) << omega
	<< setw(15) << -difLzT/Iz
	<< setw(15) << btot;
    for (int ir=0; ir<numRes; ir++) out << setw(15) << -difLz0[ir]/Iz;
    if (debug)
      out << setw(15) << difLNP0;
    out << endl;
  }

  first = false;
}


void * UserResPotN::determine_acceleration_and_potential_thread(void * arg) 
{
  double amp, R2, R, res;
  double posI[3], posO[3], velI[3], velO[3], Lz0, Lz1, Lz2;
  double pos1[3], vel1[3], dpos[3], dvel[3];
  
  unsigned nbodies = cC->Number();
  int id = *((int*)arg);
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  amp = AMP *
    0.5*(1.0 + erf( (tnow - ton) /delta )) *
    0.5*(1.0 + erf( (toff - tnow)/delta )) ;
    
  

				// Check for nan (can get rid of this
				// eventually)
  bool found_nan = false;
  ResPot::ReturnCode ret;
  int ir;

  for (int i=nbeg; i<nend; i++) {

    if (usetag>=0 && cC->Part(i)->iattrib[usetag]) continue;

				// Initial conditions
    for (int k=0; k<3; k++) {
      posI[k] = cC->Pos(i, k);
      if (c0) posI[k] -= c0->com[k];
      velI[k] = cC->Vel(i, k);
    }
    
				// Initial Lz
    Lz0 = posI[0]*velI[1] - posI[1]*velI[0];
    
    
				// Update without perturbation
    ret = respot[0]->
      Update(dtime, phase, omega, 0.0, posI, velI, pos1, vel1, &res);
    
    Lz2 = pos1[0]*vel1[1] - pos1[1]*vel1[0];

    R2 = 0.0;
    for (int k=0; k<3; k++) {
				// Accumulated difference
      dpos[k] = 0.0;
      dvel[k] = 0.0;

      R2 += posI[k]*posI[k];
    }
    R = sqrt(R2);

    if (R>rmin && R<rmax && ret == ResPot::OK) {

      ir = i % numRes;
      
      if ((ret=respot[ir]-> 
	   Update(dtime, phase, omega, amp,
		  posI, velI, posO, velO, &res)) == ResPot::OK) {
	
				// Current ang mom
	Lz1 = posO[0]*velO[1] - posO[1]*velO[0];

	  
				// Accumulate change in Lz for each resonance
	if (respot[ir]->K()<Klim)
	  difLz[id][ir] += cC->Mass(i)*(Lz1 - Lz0);
	
	if (debug && ir==0 && respot[ir]->K()<Klim) 
	  difLNP[id] += cC->Mass(i)*(Lz2 - Lz0);

	/*
	if (fabs(Lz1-Lz0) > 1.0e-12) {
	  pthread_mutex_lock(&iolock);
	  cout << setw(15) << tnow
	       << setw(5) << myid
	       << setw(5) << id
	       << setw(10) << i 
	       << setw(15) << R
	       << setw(15) << Lz0
	       << setw(15) << Lz1
	       << setw(15) << Lz2
	       << setw(15) << Lz1 - Lz0
	       << setw(15) << Lz2 - Lz0
	       << endl;
	  pthread_mutex_unlock(&iolock);
	}
	*/

				// Accumulate changes in PS due to resonances
	for (int k=0; k<3; k++) {
	  dpos[k] += posO[k] - pos1[k];
	  dvel[k] += velO[k] - vel1[k];
	}

      }
	
    }

    if (ret != ResPot::OK) {
      bcount[id]++;
      if (usetag>=0) cC->Part(i)->iattrib[usetag] = 1;
#ifdef DEBUG
      pthread_mutex_lock(&iolock);
      cout << "Process " << myid << " id=" << id << ":"
	   << " i=" << myid << " Error=" << ResPot::ReturnDesc[ret] << endl;
      pthread_mutex_unlock(&iolock);
#endif
    }

    for (int k=0; k<3; k++) {
      cC->Part(i)->pos[k] = pos1[k] + dpos[k];
      cC->Part(i)->vel[k] = vel1[k] + dvel[k];
      cC->Part(i)->acc[k] = 0.0;
      if (!found_nan) {
	if ( isnan(cC->Pos(i, k)) ||
	     isnan(cC->Vel(i, k)) ||
	     isnan(cC->Acc(i, k)) ) found_nan = true; 
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
	

    /*
    Lz1 = 
      cC->Part(i)->pos[0] * cC->Part(i)->vel[1] - 
      cC->Part(i)->pos[1] * cC->Part(i)->vel[0] ;

    if (respot[0]->K()<Klim) 
      difLz[id][0] += cC->Mass(i)*(Lz1 - Lz0);

    if (debug && respot[0]->K()<Klim) 
      difLNP[id] += cC->Mass(i)*(Lz2 - Lz0);
    */
    
  }

  /*
  pthread_mutex_lock(&iolock);
  cout << setw(15) << tnow
       << setw(15) << amp
       << setw(15) << ResPot::ReturnDesc[ret]
       << setw(5) << myid
       << setw(5) << id
       << setw(8) << nbeg
       << setw(8) << nend;
  for (int ir=0; ir<numRes; ir++)
    cout << setw(15) << difLz[id][ir];
  cout << endl;
  pthread_mutex_unlock(&iolock);
  */

  return (NULL);
}


extern "C" {
  ExternalForce *makerResPotN(string& line)
  {
    return new UserResPotN(line);
  }
}

class proxyN { 
public:
  proxyN()
  {
    factory["userrespot2"] = makerResPotN;
  }
};

static proxyN p;
