#include <sys/timeb.h>
#include <stdlib.h>
#include <sys/types.h>
#include <getopt.h>
#include <time.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <boost/lexical_cast.hpp>

#include <expand.h>
#include <ExternalCollection.H>
#include <UserTreeDSMC.H>
#include <CollideLTE.H>
#include <CollideIon.H>

// #define DEBUG

#ifdef USE_GPTL
#include <gptl.h>
#endif

using namespace std;

//
// Debugging check
//
static bool sampcel_debug = false;

//
// Physical units
//

static double pc = 3.08567758e18;       // cm
static double a0 = 0.052917721092e-7;	// cm (Bohr radius)
static double boltz = 1.381e-16;	// cgs
// static double year = 365.242*24*3600;	// seconds
static double mp = 1.67262178e-24;      // g
static double msun = 1.9891e33;		// g

double UserTreeDSMC::Lunit = 3.0e5*pc;
double UserTreeDSMC::Munit = 1.0e12*msun;
double UserTreeDSMC::Tunit = sqrt(Lunit*Lunit*Lunit/(Munit*6.67384e-08));
double UserTreeDSMC::Vunit = Lunit/Tunit;
double UserTreeDSMC::Eunit = Munit*Vunit*Vunit;
bool   UserTreeDSMC::use_effort = true;

//std::map<int, double> UserTreeDSMC::atomic_weights;
std::set<std::string> UserTreeDSMC:: colltypes;

UserTreeDSMC::UserTreeDSMC(string& line) : ExternalForce(line)
{
  (*barrier)("TreeDSMC: BEGIN construction", __FILE__, __LINE__);

  id = "TreeDSMC";		// ID string

				// Default parameter values
  ncell      = 7;		// 
  Ncell      = 64;
  cnum       = 0;
  madj       = 512;		// No tree pruning by default
  epsm       = -1.0;
  diamfac    = 1.0;
  boxsize    = 1.0;
  boxratio   = 1.0;
  comp_name  = "gas disk";
  ctype      = "Ion";
  nsteps     = -1;
  msteps     = -1;
  use_temp   = -1;
  use_dens   = -1;
  use_delt   = -1;
  use_Kn     = -1;
  use_St     = -1;
  use_vol    = -1;
  use_exes   = -1;
  coolfrac   = 0.1;
  enhance    = 1.0;
  frontier   = false;
  tsdiag     = false;
  voldiag    = false;
  tspow      = 4;
  mfpstat    = false;
  cbadiag    = false;
  dryrun     = false;
  nocool     = false;
  use_multi  = false;
  use_pullin = false;
  esol       = false;
  ntc        = true;
  cba        = true;
  tube       = false;
  slab       = false;
  sub_sample = true;
  treechk    = false;
  mpichk     = false;

				// static initialization
  initialize_colltypes();
				// Initialize using input parameters
  initialize();

  // initialize the atomic_weights map hardcode the atomic weight map
  // for use in collFrac

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

  // add in atomic weights for any other higher, more tracer, species

				// Look for the fiducial component
  bool found = false;
  list<Component*>::iterator cc;
  Component *c;
  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    c = *cc;
    if ( !comp_name.compare(c->name) ) {
      c0 = c;
      found = true;
      break;
    }
  }

  if (!found) {
    cerr << "UserTreeDSMC: process " << myid 
	 << ": can't find fiducial component <" << comp_name << ">" << endl;
    MPI_Abort(MPI_COMM_WORLD, 35);
  }
  
  (*barrier)("TreeDSMC: BEFORE use checks", __FILE__, __LINE__);

  //
  // Sanity check on excess attribute if excess calculation is
  // desired
  //
  if (use_exes>=0) {

    int ok1 = 1, ok;

    PartMapItr p    = c0->Particles().begin();
    PartMapItr pend = c0->Particles().end();
    for (; p!=pend; p++) {
      if (use_exes >= static_cast<int>(p->second.dattrib.size())) {
	ok1 = 0;
	break;
      }
    }

    MPI_Allreduce(&ok1, &ok, 1, MPI_INT, MPI_PROD, MPI_COMM_WORLD);

    // Turn off excess calculation
    // if particles have incompatible attributes
    //
    if (ok==0) {
      if (myid==0) {
	cout << "UserTreeDSMC: excess calculation requested but some" << endl
	     << "particles have incompatible float attribute counts." << endl
	     << "Attribute #" << use_exes << ". Continuing without excess."
	     << endl;
      }
      use_exes = -1;
    }
  }

  //
  // Sanity check on Knudsen number calculation
  //
  if (use_Kn>=0) {

    int ok1 = 1, ok;

    PartMapItr p    = c0->Particles().begin();
    PartMapItr pend = c0->Particles().end();
    for (; p!=pend; p++) {
      if (use_Kn >= static_cast<int>(p->second.dattrib.size())) {
	ok1 = 0;
	break;
      }
    }

    MPI_Allreduce(&ok1, &ok, 1, MPI_INT, MPI_PROD, MPI_COMM_WORLD);

    // Turn off excess calculation
    // if particles have incompatible attributes
    //
    if (ok==0) {
      if (myid==0) {
	cout << "UserTreeDSMC: Knudsen number calculation requested but some" << endl
	     << "particles have incompatible float attribute counts." << endl
	     << "Attribute #" << use_Kn << ". Continuing without Knudsen number computation."
	     << endl;
      }
      use_Kn = -1;
    }
  }

  //
  // Sanity check on Strouhal number calculation
  //
  if (use_St>=0) {

    int ok1 = 1, ok;

    PartMapItr p    = c0->Particles().begin();
    PartMapItr pend = c0->Particles().end();
    for (; p!=pend; p++) {
      if (use_St >= static_cast<int>(p->second.dattrib.size())) {
	ok1 = 0;
	break;
      }
    }

    MPI_Allreduce(&ok1, &ok, 1, MPI_INT, MPI_PROD, MPI_COMM_WORLD);

    // Turn off Strouhal number calculation
    // if particles have incompatible attributes
    //
    if (ok==0) {
      if (myid==0) {
	cout << "UserTreeDSMC: Strouhal number calculation requested but some" << endl
	     << "particles have incompatible float attribute counts." << endl
	     << "Attribute #" << use_St << ". Continuing without Strouhal number computation."
	     << endl;
      }
      use_St = -1;
    }
  }

  (*barrier)("TreeDSMC: BEFORE species map", __FILE__, __LINE__);

  //
  // Make the initial species map
  //
  {
    int ok1 = 1, ok;

    PartMapItr p    = c0->Particles().begin();
    PartMapItr pend = c0->Particles().end();

    for (; p!=pend; p++) {
      int Zi = static_cast<int>(p->second.Z);
      for (int i=1; i<=Zi+1; i++) {
	speciesKey indxi(p->second.Z, i);
	if (spec1.find(indxi) == spec1.end()) spec1[indxi] = 0;
      }
      speciesKey indx(Zi, p->second.C);
      if (spec1.find(indx) == spec1.end()) spec1[indx] = 1;
      else                                 spec1[indx]++;
    }


    MPI_Allreduce(&ok1, &ok, 1, MPI_INT, MPI_PROD, MPI_COMM_WORLD);

    // cout << ok << "\t" << ok1 << endl;
    
    if (ok) {

      int sizm;
      speciesKey indx;
      unsigned long cnts;
      std::map<speciesKey, unsigned long>::iterator it, it2;

      spec = spec1;
      for (int i=0; i<numprocs; i++) {
	if (i == myid) {
	  sizm = spec1.size();
	  MPI_Bcast(&sizm, 1, MPI_INT, i, MPI_COMM_WORLD);

	  for (it=spec1.begin(); it != spec1.end(); it++) {
	    indx = it->first;
	    cnts = it->second;
	    MPI_Bcast(&indx.first, 1, MPI_UNSIGNED, i, MPI_COMM_WORLD);
	    MPI_Bcast(&indx.second, 1, MPI_UNSIGNED, i, MPI_COMM_WORLD);
	    MPI_Bcast(&cnts, 1, MPI_UNSIGNED_LONG, i, MPI_COMM_WORLD);
	  }
	} else {
	  MPI_Bcast(&sizm, 1, MPI_INT, i, MPI_COMM_WORLD);
	  for (int j=0; j<sizm; j++) {
	    MPI_Bcast(&indx.first, 1, MPI_UNSIGNED, i, MPI_COMM_WORLD);
	    MPI_Bcast(&indx.second, 1, MPI_UNSIGNED, i, MPI_COMM_WORLD);
	    MPI_Bcast(&cnts, 1, MPI_UNSIGNED_LONG, i, MPI_COMM_WORLD);
	    if (spec.find(indx) == spec.end()) spec[indx]  = cnts;
	    else                               spec[indx] += cnts;
	  }
	}
      }
      
      for (it=spec.begin(); it != spec.end(); it++)  {
	indx = it->first;
	spec_list.insert(indx);
        if (collFrac.find(indx) == collFrac.end()) 
	  collFrac[indx]  = 1.0/(atomic_weights[indx.first]);
	else {			// Ask Brandt about this . . . 
	  collFrac[indx] *= 1.0/(atomic_weights[indx.first]);
	  std::cout << "Weird stuff!" << std::endl;
	}
      }

      std::map<speciesKey, unsigned long> check = spec;
      for (it=check.begin(); it!=check.end(); it++) it->second = 0;

      if (myid==0) {
	cout << endl
	     << "--------------" << endl
	     << "Species counts" << endl
	     << "--------------" << endl
	     << endl;

	cout << setw(4) << right << "#";
	for (it=spec.begin(); it != spec.end(); it++)
	  cout << setw(8) << right 
	       << "(" << it->first.first << "," << it->first.second << ")";
	cout << endl;

	cout << setw(4) << right << "---";
	for (it=spec.begin(); it != spec.end(); it++)
	  cout << setw(12) << right << "--------";
	cout << endl;

	it2 = spec1.begin();
	cout << setw(4) << right << myid;
	for (it=spec.begin(); it != spec.end(); it++) {
	  if (it->first == it2->first) {
	    cout << setw(12) << right << it2->second;
	    check[it->first] += it2->second;
	    it2++;
	  } else {
	    cout << setw(12) << right << 0;
	  }
	}
	cout << endl;
      }
      
      unsigned val;
      for (int i=1; i<numprocs; i++) {
	if (i == myid) {
	  it2 = spec1.begin();
	  for (it=spec.begin(); it != spec.end(); it++) {
	    if (it->first == it2->first) val = (it2++)->second;
	    else                         val = 0;
	    MPI_Send(&val, 1, MPI_UNSIGNED, 0, 142, MPI_COMM_WORLD);
	  }
	}
	if (myid == 0) {
	  cout << setw(4) << right << i;
	  for (it=spec.begin(); it != spec.end(); it++) {
	    MPI_Recv(&val, 1, MPI_UNSIGNED, i, 142, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    cout << setw(12) << right << val;
	    check[it->first] += val;
	  }
	  cout << endl;
	}
      }
      if (myid==0) {
	cout << setw(4) << right << "---";
	for (it=spec.begin(); it != spec.end(); it++)
	  cout << setw(12) << right << "--------";
	cout << endl;
	cout << setw(4) << right << "TOT";
	for (it=spec.begin(); it != spec.end(); it++)
	  cout << setw(12) << right << it->second;
	cout << endl;
	cout << setw(4) << right << "CHK";
	for (it=check.begin(); it != check.end(); it++)
	  cout << setw(12) << right << it->second;
	cout << endl << endl;
      }
    }
  }

  (*barrier)("TreeDSMC: AFTER species map", __FILE__, __LINE__);


  Vunit = Lunit/Tunit;
  Eunit = Munit*Vunit*Vunit;

				// Number of protons per mass unit
  for (std::map<speciesKey, double>::iterator 
	 it=collFrac.begin(); it!=collFrac.end(); it++) it->second *= Munit/mp;

  pHOT::sub_sample = sub_sample;

  c0->HOTcreate(spec_list);

  if (tube) {
    c0->Tree()->setSides (boxsize*boxratio, boxsize, boxsize);
    c0->Tree()->setOffset(0.0,              0.0,     0.0);
  } 
  else if (slab) {
    c0->Tree()->setSides (boxsize, boxsize, boxsize*boxratio);
    c0->Tree()->setOffset(0.0,              0.0,     0.0);
  } 
  else {
    c0->Tree()->setSides (2.0*boxsize, 2.0*boxsize, 2.0*boxsize*boxratio);
    c0->Tree()->setOffset(    boxsize,     boxsize,     boxsize*boxratio);
  }

  pCell::bucket = ncell;
  pCell::Bucket = Ncell;

  volume = pHOT::sides[0] * pHOT::sides[1] * pHOT::sides[2];

  //
  // Set collision parameters
  //
  Collide::NTC     = ntc;
  Collide::CBA     = cba;
  Collide::CBADIAG = cbadiag;
  Collide::PULLIN  = use_pullin;
  Collide::ESOL    = esol;
  Collide::EPSMratio = epsm;
  Collide::DRYRUN  = dryrun;
  Collide::NOCOOL  = nocool;
  Collide::TSDIAG  = tsdiag;
  Collide::VOLDIAG = voldiag;
  Collide::TSPOW   = tspow;
  Collide::MFPDIAG = mfpstat;
  Collide::EFFORT  = use_effort;
  Collide::ENHANCE = enhance;

  //
  // Create the collision instance from the allowed list
  //
  if (ctype.compare("LTE") == 0)
    collide = new CollideLTE(this, diamfac, nthrds);
  if (ctype.compare("Ion") == 0)
    collide = new CollideIon(this, diamfac, nthrds);
  else {
    std::cout << "No such Collide type: " << ctype << std::endl;
    exit(-1);
  }

  collide->set_temp_dens(use_temp, use_dens);
  if (esol) collide->set_timestep(-1);
  else      collide->set_timestep(use_delt);
  collide->set_Kn(use_Kn);
  collide->set_St(use_St);
  collide->set_excess(use_exes);
  ElostTotCollide = ElostTotEPSM = 0.0;

  //
  // Timers: set precision to microseconds
  //
  
  partnTime.Microseconds();
  tree1Time.Microseconds();
  tradjTime.Microseconds();
  tcellTime.Microseconds();
  tstepTime.Microseconds();
  llistTime.Microseconds();
  clldeTime.Microseconds();
  clldeWait.Microseconds();
  partnWait.Microseconds();
  tree1Wait.Microseconds();
  tree2Wait.Microseconds();
  timerDiag.Microseconds();

  //
  // Quantiles for distribution diagnstic
  //
  quant.push_back(0.0);		// 0
  quant.push_back(0.01);	// 1
  quant.push_back(0.05);	// 2
  quant.push_back(0.1);		// 3
  quant.push_back(0.2);		// 4
  quant.push_back(0.5);		// 5
  quant.push_back(0.8);		// 6
  quant.push_back(0.9);		// 7
  quant.push_back(0.95);	// 8
  quant.push_back(0.99);	// 9
  quant.push_back(1.0);		// 10

  userinfo();

  (*barrier)("TreeDSMC: END construction", __FILE__, __LINE__);
}

UserTreeDSMC::~UserTreeDSMC()
{
  delete collide;
}

void UserTreeDSMC::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine TreeDSMC initialized, "
       << "Lunit=" << Lunit << ", Tunit=" << Tunit << ", Munit=" << Munit
       << ", Vunit=" << Vunit << ", Eunit=" << Eunit
       << ", cnum=" << cnum << ", diamfac=" << diamfac
       << ", madj=" << madj << ", epsm=" << epsm << ", boxsize=" << boxsize 
       << ", ncell=" << ncell << ", Ncell=" << Ncell 
       << ", boxratio=" << boxratio << ", compname=" << comp_name;
  if (msteps>=0) 
    cout << ", with diagnostic output at levels <= " << msteps;
  else if (nsteps>0) 
    cout << ", with diagnostic output every " << nsteps << " steps";
  if (use_temp>=0) cout << ", temp at pos="   << use_temp;
  if (use_dens>=0) cout << ", dens at pos="   << use_dens;
  if (use_Kn>=0)   cout << ", Kn at pos="     << use_Kn;
  if (use_St>=0)   cout << ", St at pos="     << use_St;
  if (use_vol>=0)  cout << ", cell volume at pos=" << use_vol;
  if (use_exes>=0) cout << ", excess at pos=" << use_exes;
  if (use_pullin)  cout << ", Pullin algorithm enabled";
  if (dryrun)      cout << ", collisions disabled";
  if (nocool)      cout << ", cooling disabled";
  if (epsm>0.0)    cout << ", using EPSM";
  else             cout << ", EPSM disabled";
  if (ntc)         cout << ", using NTC";
  else             cout << ", NTC disabled";
  if (cba)         cout << ", using CBA";
  else             cout << ", CBA disabled";
  if (cba && cbadiag)     
                   cout << " with diagnostics";
  if (tube)        cout << ", using TUBE mode";
  else if (slab)   cout << ", using THIN SLAB mode";
  if (use_effort)  cout << ", with effort-based load";
  else             cout << ", with uniform load";
  if (fabs(enhance-1.0)>1.0e-6)
                   cout << ", with enhanced cooling of " << enhance;
  if (use_multi) {
    cout << ", multistep enabled";
    if (use_delt>=0) 
      cout << ", time step at pos=" << use_delt << ", coolfrac=" << coolfrac;
  }
  cout << endl;

  print_divider();
}

void UserTreeDSMC::initialize()
{
  string val;

  if (get_value("Lunit", val))		Lunit      = atof(val.c_str());
  if (get_value("Tunit", val))		Tunit      = atof(val.c_str());
  if (get_value("Munit", val))		Munit      = atof(val.c_str());
  if (get_value("cnum", val))		cnum       = atoi(val.c_str());
  if (get_value("madj", val))		madj       = atoi(val.c_str());
  if (get_value("epsm", val))		epsm       = atof(val.c_str());
  if (get_value("diamfac", val))	diamfac    = atof(val.c_str());
  if (get_value("boxsize", val))	boxsize    = atof(val.c_str());
  if (get_value("boxratio", val))	boxratio   = atof(val.c_str());
  if (get_value("coolfrac", val))	coolfrac   = atof(val.c_str());
  if (get_value("enhance", val))	enhance    = atof(val.c_str());
  if (get_value("nsteps", val))		nsteps     = atoi(val.c_str());
  if (get_value("msteps", val))		msteps     = atoi(val.c_str());
  if (get_value("ncell", val))		ncell      = atoi(val.c_str());
  if (get_value("Ncell", val))		Ncell      = atoi(val.c_str());
  if (get_value("compname", val))	comp_name  = val;
  if (get_value("use_temp", val))	use_temp   = atoi(val.c_str());
  if (get_value("use_dens", val))	use_dens   = atoi(val.c_str());
  if (get_value("use_delt", val))	use_delt   = atoi(val.c_str());
  if (get_value("use_Kn", val))		use_Kn     = atoi(val.c_str());
  if (get_value("use_St", val))		use_St     = atoi(val.c_str());
  if (get_value("use_vol", val))	use_vol    = atoi(val.c_str());
  if (get_value("use_exes", val))	use_exes   = atoi(val.c_str());
  if (get_value("frontier", val))	frontier   = atol(val);
  if (get_value("tspow", val))		tspow      = atoi(val.c_str());
  if (get_value("tsdiag", val))		tsdiag     = atol(val);
  if (get_value("voldiag", val))	voldiag    = atol(val);
  if (get_value("mfpstat", val))	mfpstat    = atol(val);
  if (get_value("cbadiag", val))	cbadiag    = atol(val);
  if (get_value("dryrun", val))		dryrun     = atol(val);
  if (get_value("nocool", val))		nocool     = atol(val);
  if (get_value("use_multi", val))	use_multi  = atol(val);
  if (get_value("use_pullin", val))	use_pullin = atol(val);
  if (get_value("use_effort", val))	use_effort = atol(val);
  if (get_value("esol", val))		esol       = atol(val);
  if (get_value("cba", val))		cba        = atol(val);
  if (get_value("ntc", val))		ntc        = atol(val);
  if (get_value("tube", val))		tube       = atol(val);
  if (get_value("slab", val))		slab       = atol(val);
  if (get_value("sub_sample", val))	sub_sample = atol(val);
  if (get_value("treechk", val))	treechk    = atol(val);
  if (get_value("mpichk", val))		mpichk     = atol(val);

  if (get_value("ctype", val)) {
    if (check_ctype(val)) ctype = val;
    else {
      if (myid==0) {
	std::cerr << "UserTreeDSMC: invalid ctype <" << ctype << ">" 
		  << std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      exit(-1);
    }
  }

  /**
     Look for array values in the parameter string of the form
     spc(1,2)=3.1, spc(3,4)=5.6, etc.
  */
  
  if (ctype.compare("LTE")==0) {
    std::map<std::pair<int, int>, string> vals;
    if ((vals = get_value_matrix("spc")).size()) {
      std::map<std::pair<int, int>, string>::iterator it=vals.begin();
      while (it != vals.end()) {
	try {
	  speciesKey p(it->first.first, it->first.second);
	  collFrac[p] = boost::lexical_cast<double>(it->second);
	} 
	catch( boost::bad_lexical_cast const& ) {
	  std::cout << "UserTreeDSMC::initialize: bad double value, "
		    << "input string was: " << it->second << std::endl;
	}
      }
      it++;
    } else {
      // The default key is defined in pCell.H
      collFrac[defaultKey] = 1.0;
    }
  }
}


//  __  __    _    ___ _   _   ____   ___  _   _ _____ ___ _   _ _____ 
// |  \/  |  / \  |_ _| \ | | |  _ \ / _ \| | | |_   _|_ _| \ | | ____|
// | |\/| | / _ \  | ||  \| | | |_) | | | | | | | | |  | ||  \| |  _|  
// | |  | |/ ___ \ | || |\  | |  _ <| |_| | |_| | | |  | || |\  | |___ 
// |_|  |_/_/   \_\___|_| \_| |_| \_\\___/ \___/  |_| |___|_| \_|_____|
//                                                                    

void UserTreeDSMC::determine_acceleration_and_potential(void)
{
  static unsigned cat = 0;	// For debugging sync
  {
    std::ostringstream sout;
    sout << "TreeDSMC: call=" << ++cat;
    (*barrier)(sout.str(), __FILE__, __LINE__);
  }

#ifdef USE_GPTL
  GPTLstart("UserTreeDSMC::determine_acceleration_and_potential");
#endif

  static bool firstime = true;
  static unsigned nrep = 0;

  //
  // Only compute DSMC when passed the fiducial component
  //

  if (cC != c0) {
#ifdef USE_GPTL
    GPTLstop("UserTreeDSMC::determine_acceleration_and_potential");
#endif
    return;
  }

  //
  // Get timing for entire step so far to load balancing the partition
  //
  double pot_time = c0->get_time_sofar();


  (*barrier)("TreeDSMC: after initialization", __FILE__, __LINE__);


  //
  // Make the cells
  //

  if (firstime) {
    //
    // This is a full repartition tree build
    //
    c0->Tree()->setWeights(use_effort);
    c0->Tree()->Repartition(0); nrep++;
    c0->Tree()->makeTree();
    c0->Tree()->makeCellLevelList();
#ifdef DEBUG
    if (!c0->Tree()->checkBodycell()) {
      cout << "Process " << myid << ": "
	   << "makeTree completed: body cell check FAILED!" << endl;
    }    
    if (!c0->Tree()->checkParticles(cout)) {
      cout << "Process " << myid << ": "
	   << "makeTree completed: initial particle check FAILED!" << endl;
    }    
    if (!c0->Tree()->checkFrontier(cout)) {
      cout << "Process " << myid << ": "
	   << "makeTree completed: frontier check FAILED!" << endl;
    }

    if (!c0->Tree()->checkKeybods()) {
      cout << "Process " << myid 
	   << ": makeTree: ERROR particle key not in keybods AFTER makeTree(), T=" 
	   << tnow << endl;
    }
#endif
    if (sampcel_debug) {
      ostringstream sout;
      sout << "after MAKE tree, first time, "
	   << "[" << 0 << ", " << tnow  << "]";
      c0->Tree()->checkSampleCells(sout.str().c_str());
      c0->Tree()->logFrontierStats();
    }

    if (use_temp || use_dens || use_vol) assignTempDensVol();

    stepnum = 0;
    curtime = tnow;

#ifdef DEBUG
    cout << "Computed partition and tree [firstime on #" 
	 << setw(4) << left << myid << "]" << endl;
#endif
  } else {
    
    if (tnow-curtime < 1.0e-14) {
      if (myid==0) {
	cout << "UserTreeDSMC: attempt to redo step at T=" << tnow << endl;
      }
#ifdef USE_GPTL
      GPTLstop("UserTreeDSMC::determine_acceleration_and_potential");
#endif
      return; 			// Don't do this time step again!
    }

    stepnum++;
    curtime = tnow;
  }

#ifdef DEBUG
  c0->Tree()->densCheck();
#endif
  
#ifdef DEBUG
  if (!c0->Tree()->checkParticles(cout)) {
    cout << "After init only: Particle check FAILED [" << right
	 << setw(3) << mlevel << ", " << setw(3) << myid << "]" << endl;
  }
#endif

  (*barrier)("TreeDSMC: after cell computation", __FILE__, __LINE__);

  //
  // Only run diagnostics every nsteps
  //
  bool diagstep = (nsteps>0 && mstep%nsteps == 0);

  //
  // Diagnostics run at levels <= msteps (takes prececdence over nsteps)
  //
  if (msteps>=0) 
    diagstep = (mlevel <= static_cast<unsigned>(msteps)) ? true : false;

  //
  // Compute time step
  //
  // Now, DSMC is computed on the smallest step, every step
  double tau = dtime*mintvl[multistep]/Mstep;


  MPI_Bcast(&tau, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  TimeElapsed partnSoFar, tree1SoFar, tradjSoFar, tcellSoFar, tstepSoFar;
  TimeElapsed waitcSoFar, waitpSoFar, wait1SoFar, wait2SoFar, timerSoFar;
  TimeElapsed collideSoFar;

  overhead.Start();

  //
  // Sort the particles into cells
  //
  if (mlevel<=madj) {

#ifdef USE_GPTL
    GPTLstart("UserTreeDSMC::pHOT");
    GPTLstart("UserTreeDSMC::waiting");
    (*barrier)("TreeDSMC: pHOT waiting");
    GPTLstop ("UserTreeDSMC::waiting");
    GPTLstart("UserTreeDSMC::repart");
#endif

    (*barrier)("TreeDSMC: after pHOT wait", __FILE__, __LINE__);

    partnTime.start();
    c0->Tree()->Repartition(mlevel); nrep++;
    partnSoFar = partnTime.stop();

    partnWait.start();
    (*barrier)("TreeDSMC: after repartition", __FILE__, __LINE__);
    waitpSoFar = partnWait.stop();

#ifdef USE_GPTL
    GPTLstop ("UserTreeDSMC::repart");
    GPTLstart("UserTreeDSMC::makeTree");
#endif
    tree1Time.start();
    c0->Tree()->makeTree();
    c0->Tree()->makeCellLevelList();
    tree1Time.stop();
    tree1Wait.start();
    (*barrier)("TreeDSMC: after makeTree", __FILE__, __LINE__);
    wait1SoFar = tree1Wait.stop();
#ifdef USE_GPTL
    GPTLstop ("UserTreeDSMC::makeTree");
    GPTLstart("UserTreeDSMC::pcheck");
#endif
    tree1Time.start();
#ifdef DEBUG
    cout << "Made partition, tree and level list [" << mlevel << "]" << endl;
    if (!c0->Tree()->checkParticles(cout)) {
      cout << "Particle check on new tree FAILED [" << mlevel << "]" << endl;
    }
#endif
    if (sampcel_debug) {
      ostringstream sout;
      sout << "after MAKE tree at level "
	   << "[" << mlevel << ", " << tnow  << "]";
      c0->Tree()->checkSampleCells(sout.str().c_str());
      c0->Tree()->logFrontierStats();
    }
    tree1SoFar = tree1Time.stop();
    
#ifdef USE_GPTL
    GPTLstop("UserTreeDSMC::pcheck");
    GPTLstop("UserTreeDSMC::pHOT");
#endif

  } else {

#ifdef USE_GPTL
    GPTLstart("UserTreeDSMC::pHOT_2");
    GPTLstart("UserTreeDSMC::adjustTree");
#endif

#ifdef DEBUG
    if (myid==0)
      cout << "About to adjust tree [" << mlevel << "]" << endl;
#endif
    tradjTime.start();
    c0->Tree()->adjustTree(mlevel);
    tradjSoFar = tradjTime.stop();
    tcellTime.start();
    c0->Tree()->adjustCellLevelList(mlevel);
    tcellSoFar = tcellTime.stop();
    tree2Wait.start();
    (*barrier)("TreeDSMC: after adjustTree", __FILE__, __LINE__);
    wait2SoFar = tree2Wait.stop();

    if (sampcel_debug) {
      ostringstream sout;
      sout << "after ADJUST tree at level "
	   << "[" << mlevel << ", " << tnow  << "]";
      c0->Tree()->checkSampleCells(sout.str().c_str());
    }

#ifdef USE_GPTL
    GPTLstop("UserTreeDSMC::adjustTree");
    GPTLstop("UserTreeDSMC::pHOT_2");
#endif

  }

  overhead.Stop();
  pot_time += overhead.getTime();
  pot_time /= max<unsigned>(1, c0->Number());
  
  if (use_effort) {
    PartMapItr pitr = c0->Particles().begin(), pend = c0->Particles().end();
    for (; pitr!= pend; pitr++) pitr->second.effort = pot_time;
  }

  //
  // Evaluate collisions among the particles
  //

  clldeTime.start();
  
  if (0) {
    ostringstream sout;
    sout << "before Collide [" << nrep << "], " 
	 << __FILE__ << ": " << __LINE__;
    c0->Tree()->checkBounds(2.0, sout.str().c_str());
  }

#ifdef USE_GPTL
  GPTLstart("UserTreeDSMC::collide");
#endif

  //
  // So far, all computations have been about repartition and
  // tessellation.  All of the collision stuff is done by the current
  // Collide class instance.
  //

  collide->collide(*c0->Tree(), collFrac, tau, mlevel, diagstep);
    
  collideSoFar = clldeTime.stop();

  clldeWait.start();
  (*barrier)("TreeDSMC: after collide", __FILE__, __LINE__);

#ifdef USE_GPTL
  GPTLstop("UserTreeDSMC::collide");
#endif

  waitcSoFar = clldeWait.stop();

  // -----------------
  // Time step request
  // -----------------
  //
  // New timesteps are selected for the cells based on the collision
  // diagnostics from the current step.
  //

#ifdef USE_GPTL
  GPTLstart("UserTreeDSMC::collide_timestep");
#endif


  (*barrier)("TreeDSMC: before collide timestep", __FILE__, __LINE__);

  tstepTime.start();
  if (use_multi) collide->compute_timestep(c0->Tree(), coolfrac);
  tstepSoFar = tstepTime.stop();
  
#ifdef USE_GPTL
  GPTLstop("UserTreeDSMC::collide_timestep");
#endif

  //
  // Repartition and remake the tree after first step to adjust load
  // balancing for work queue effort method
  //
  if (firstime && use_effort && mlevel==0) {
    c0->Tree()->Repartition(0); nrep++;
    c0->Tree()->makeTree();
    c0->Tree()->makeCellLevelList();
  }

  firstime = false;

  //
  // Periodically display the current progress
  //
  // Lots of diagnostics are computed and emitted here . . .
  //
  if (diagstep) {
#ifdef USE_GPTL
    GPTLstart("UserTreeDSMC::collide_diag");
#endif

				// Uncomment for debug
    // collide->Debug(tnow);

    unsigned medianNumb = collide->medianNumber();
    unsigned collnum=0, coolnum=0;
    if (mfpstat) {
      collide->collQuantile(quant, coll_);
      collide->mfpsizeQuantile(quant, mfp_, ts_, nsel_, cool_, rate_,
			       collnum, coolnum);
    }

    double ExesCOLL, ExesEPSM;
    if (use_exes>=0) collide->energyExcess(ExesCOLL, ExesEPSM);
      
    if (frontier) {
      ostringstream sout;
      sout << outdir << runtag << ".DSMC_frontier";
      string filen = sout.str();
      c0->Tree()->testFrontier(filen);
    }

    vector<unsigned> ncells, bodies;
    c0->Tree()->countFrontier(ncells, bodies);

    if (mfpstat && myid==0) {

				// Generate the file name
      ostringstream sout;
      sout << outdir << runtag << ".DSMC_mfpstat";
      string filen = sout.str();

				// Check for existence
      ifstream in(filen.c_str());
      if (in.fail()) {
				// Write a new file
	ofstream out(filen.c_str());
	if (out) {
	  out << "# " << right
	      << setw(12) << "Time"
	      << setw(14) << "Quantiles" 
	      << setw(14) << "Cells"
	      << setw(14) << "MFP/size"
	      << setw(14) << "Flight/size"
	      << setw(14) << "Collions/cell"
	      << setw(14) << "Nsel/Number"
	      << setw(14) << "Energy ratio"
	      << setw(14) << "Excess ratio"
	      << endl;
	}
      }
      in.close();

				// Open old file to write a stanza
				// 
      ofstream out(filen.c_str(), ios::app);
      if (out) {
	for (unsigned nq=0; nq<quant.size(); nq++) {
	  out << setw(14) << tnow
	      << setw(14) << quant[nq]
	      << setw(14) << collnum
	      << setw(14) << mfp_[nq] 
	      << setw(14) << ts_[nq] 
	      << setw(14) << coll_[nq] 
	      << setw(14) << nsel_[nq] 
	      << setw(14) << cool_[nq] 
	      << setw(14) << rate_[nq] 
	      << endl;
	}
	out << endl;
      }
    }

    (*barrier)("TreeDSMC: after mfp stats", __FILE__, __LINE__);

				// Overall statistics
				// 
    double KEtotl=collide->Etotal(), KEtot=0.0;
    double Mtotal=collide->Mtotal(), Mtotl=0.0;
    double Elost1, Elost2, ElostC=0.0, ElostE=0.0;

    collide->Elost(&Elost1, &Elost2);

    MPI_Reduce(&KEtotl, &KEtot,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Mtotal, &Mtotl,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Elost1, &ElostC, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Elost2, &ElostE, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

				// Computing mass-weighted temperature
    const double f_H = 0.76;
    double mm = f_H*mp + (1.0-f_H)*4.0*mp;
    double meanT = 0.0;
    if (Mtotl>0.0) meanT = 2.0*KEtot/Mtotl*Eunit/3.0 * mm/Munit/boltz;

    unsigned cellBods = c0->Tree()->checkNumber();
    unsigned oobBods  = c0->Tree()->oobNumber();

    double Mass = 0.0;
    unsigned Counts = 0;
    c0->Tree()->totalMass(Counts, Mass);

				// Check frontier for mass at or below 
				// current level
    double cmass1=0.0, cmass=0.0;
    pHOT_iterator pit(*(c0->Tree()));

    (*barrier)("TreeDSMC: checkAdjust", __FILE__, __LINE__);

    while (pit.nextCell()) {
      pCell *cc = pit.Cell();
      if (cc->maxplev >= mlevel && cc->ctotal > 1) cmass1 += cc->stotal[0];
    }

				// Collect up info from all processes
    timerDiag.start();
    MPI_Reduce(&cmass1, &cmass, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    c0->Tree()->CollectTiming();
    collide->CollectTiming();
    timerSoFar = timerDiag.stop();

    const unsigned nf = 12;
    const unsigned nt = pHOT::ntile+2;
    if (tt.size() != nf) tt = vector<TimeElapsed>(nf);
    vector<double> in(nf), IN(nf);
    vector< vector<double> > out(nt), OUT(nt);
    for (unsigned i=0; i<nt; i++) {
      out[i] = vector<double>(nf);
      if (mlevel==0) OUT[i] = vector<double>(nf);
    }

    in[ 0] = partnSoFar();
    in[ 1] = tree1SoFar();
    in[ 2] = tradjSoFar();
    in[ 3] = tcellSoFar();
    in[ 4] = tstepSoFar();
    in[ 5] = llistTime.getTime()();
    in[ 6] = collideSoFar();
    in[ 7] = timerSoFar();
    in[ 8] = waitpSoFar();
    in[ 9] = waitcSoFar();
    in[10] = wait1SoFar();
    in[11] = wait2SoFar();

    tt[ 0] += partnSoFar;
    tt[ 1] += tree1SoFar;
    tt[ 2] += tradjSoFar;
    tt[ 3] += tcellSoFar;
    tt[ 4] += tstepSoFar;
    tt[ 5] += llistTime.getTime();
    tt[ 6] += collideSoFar;
    tt[ 7] += timerSoFar;
    tt[ 8] += waitpSoFar;
    tt[ 9] += waitcSoFar;
    tt[10] += wait1SoFar;
    tt[11] += wait2SoFar;

    if (mlevel==0) {
      for (unsigned k=0; k<nf; k++) {
	IN[k] = tt[k]();
	tt[k].zero();
      }
    }

    // Get the timing info from each process
    vector<double> valu(numprocs);
    for (unsigned j=0; j<nf; j++) {
      MPI_Gather(&in[j], 1, MPI_DOUBLE, &valu[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      // Sort the array
      sort(valu.begin(), valu.end());

      // Select the quantiles
      out[0][j]     = valu.front();
      for (unsigned k=0; k<nt-2; k++)
	out[k+1][j] = valu[static_cast<int>(floor(valu.size()*0.01*pHOT::qtile[k]))];
      out[nt-1][j]  = valu.back();

      if (mlevel==0) {

	MPI_Gather(&IN[j], 1, MPI_DOUBLE, &valu[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	// Sort the array
	sort(valu.begin(), valu.end());

	// Select the quantiles
	OUT[0][j]     = valu.front();
	for (unsigned k=0; k<nt-2; k++)
	  OUT[k+1][j] = valu[static_cast<int>(floor(valu.size()*0.01*pHOT::qtile[k]))];
	OUT[nt-1][j]  = valu.back();
      }
    }

    vector<double> tot(nt, 0.0), TOT(nt, 0.0);
    for (unsigned k=0; k<nt; k++) {
      for (unsigned i=0; i<nf; i++) {
	tot[k] += out[k][i];
	if (mlevel==0) TOT[k] += OUT[k][i];
      }
    }

    int pCellTot;
    MPI_Reduce(&pCell::live, &pCellTot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    collide->EPSMtimingGather();
    collide->CPUHogGather();

    const int ebins = 1000;
    vector<unsigned> efrt(ebins, 0);
    double minEff, maxEff;

    if (use_effort) {
      PartMapItr pitr, pend = c0->Particles().end();
      double minEff1 = 1.0e20, maxEff1 = 0.0;
      for (pitr=c0->Particles().begin(); pitr!= pend; pitr++) {
	minEff1 = min<double>(pitr->second.effort, minEff1);
	maxEff1 = max<double>(pitr->second.effort, maxEff1);
      }

      MPI_Allreduce(&minEff1, &minEff, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&maxEff1, &maxEff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

      if (minEff>0.0) {
	minEff = log(minEff);
	maxEff = log(maxEff);

	int indx;
	double Lvalue;
	vector<unsigned> efrt1(ebins, 0);
	for (pitr=c0->Particles().begin(); pitr!= pend; pitr++) {
	  Lvalue = log(pitr->second.effort);
	  indx   = floor((Lvalue - minEff) / (maxEff - minEff) * ebins);
	  if (indx<0)      indx = 0;
	  if (indx>=ebins) indx = ebins-1;
	  efrt1[indx]++;
	}

	MPI_Reduce(&efrt1[0], &efrt[0], ebins, MPI_UNSIGNED, MPI_SUM, 0, 
		   MPI_COMM_WORLD);
      }
    }

    collide->printCollGather();

    if (myid==0) {

      unsigned sell_total = collide->select();
      unsigned coll_total = collide->total();
      unsigned coll_error = collide->errors();
      unsigned epsm_total = collide->EPSMtotal();
      unsigned epsm_cells = collide->EPSMcells();
      collide->printCollSummary();

      vector<double> disp;
      collide->dispersion(disp);
      double dmean = (disp[0]+disp[1]+disp[2])/3.0;

      ostringstream sout;
      sout << outdir << runtag << ".DSMC_log";
      ofstream mout(sout.str().c_str(), ios::app);

      mout << "Summary:" << endl << left << "--------" << endl << scientific
	   << setw(6) << " " << setw(20) << tnow       << "current time" << endl
	   << setw(6) << " " << setw(20) << mlevel     << "current level" << endl
	   << setw(6) << " " << setw(20) << Counts     << "total counts" << endl
	   << setw(6) << " " << setw(20) << Mass       << "total mass" << endl
	   << setw(6) << " " << setw(20) << meanT      << "mass-weighted temperature" << endl
	   << setw(6) << " " << setw(20) << Mtotl      << "accumulated mass" << endl
	   << setw(6) << " " << setw(20) << cmass      << "mass at this level" << endl
	   << setw(6) << " " << setw(20) << mstep      << "step number" << endl
	   << setw(6) << " " << setw(20) << stepnum    << "step count" << endl
	   << setw(6) << " " << setw(20) << sell_total << "targets" << endl
	   << setw(6) << " " << setw(20) << coll_total << "collisions" << endl
	   << setw(6) << " " << setw(20) << coll_error << "collision errors (" 
	   << setprecision(2) << fixed 
	   << 100.0*coll_error/(1.0e-08+coll_total) << "%)" << endl
	   << setw(6) << " " << setw(20) << oobBods << "out-of-bounds" << endl
	   << endl;

      collide->colldeTime(mout);

      if (epsm>0) mout << setw(6) << " " << setw(20) << epsm_total 
		       << "EPSM particles ("
		       << 100.0*epsm_total/c0->nbodies_tot << "%)" 
		       << scientific << endl;
      mout << setw(6) << " " << setw(20) << medianNumb << "number/cell" << endl
	   << setw(6) << " " << setw(20) << c0->Tree()->TotalNumber() 
	   << "total # cells" << endl;

      if (epsm>0) mout << setw(6) << " " << setw(20) << epsm_cells 
		       << "EPSM cells (" << setprecision(2) << fixed 
		       << 100.0*epsm_cells/c0->Tree()->TotalNumber() 
		       << "%)" << scientific << endl;

      if (mfpstat) {
	mout << setw(6) << " " << setw(20) << "--------" << "--------------------" << endl
	     << setw(6) << " " << setw(20) << nsel_[0] << "collision/body @  0%" << endl
	     << setw(6) << " " << setw(20) << nsel_[2] << "collision/body @  5%" << endl
	     << setw(6) << " " << setw(20) << nsel_[5] << "collision/body @ 50%" << endl
	     << setw(6) << " " << setw(20) << nsel_[8] << "collision/body @ 95%" << endl
	     << setw(6) << " " << setw(20) << nsel_[10] << "collision/body @100%" << endl;

	mout << fixed << setprecision(0)
	     << setw(6) << " " << setw(20) << "--------" << "--------------------" << endl

	     << setw(6) << " " << setw(20) << coll_[0] << "collision/cell @  0%" << endl
	     << setw(6) << " " << setw(20) << coll_[2] << "collision/cell @  5%" << endl
	     << setw(6) << " " << setw(20) << coll_[5] << "collision/cell @ 50%" << endl
	     << setw(6) << " " << setw(20) << coll_[8] << "collision/cell @ 95%" << endl
	     << setw(6) << " " << setw(20) << coll_[10] << "collision/cell @100%" << endl
	     << setw(6) << " " << setw(20) << "--------" << "--------------------" << endl

	     << setw(6) << " " << setw(20) << c0->Tree()->CellCount(0.0) 
	     << "occupation @  0%" << endl
	     << setw(6) << " " << setw(20) << c0->Tree()->CellCount(0.05) 
	     << "occupation @  5%" << endl
	     << setw(6) << " " << setw(20) << c0->Tree()->CellCount(0.50) 
	     << "occupation @ 50%" << endl
	     << setw(6) << " " << setw(20) << c0->Tree()->CellCount(0.95) 
	     << "occupation @ 95%" << endl
	     << setw(6) << " " << setw(20) << c0->Tree()->CellCount(1.0) 
	     << "occupation @100%" << endl
	     << setw(6) << " " << setw(20) << "--------" << "--------------------" << endl
	     << setw(6) << " " << setw(20) << cellBods
	     << "total number in cells" << endl
	     << endl << setprecision(4);
      }
	
      ElostTotCollide += ElostC;
      ElostTotEPSM    += ElostE;

      mout << scientific << "Energy (system):" << endl 
	   << "----------------" << endl 
	   << " Lost collide = " << ElostC << endl;
      if (epsm>0) mout << "    Lost EPSM = " << ElostE << endl;
      mout << "   Total loss = " << ElostTotCollide+ElostTotEPSM << endl;
      if (epsm>0) mout << "   Total EPSM = " << ElostTotEPSM << endl;
      mout << "     Total KE = " << KEtot << endl;
      if (use_exes>=0) {
	mout << "  COLL excess =" << ExesCOLL << endl;
	if (epsm>0) mout << "  EPSM excess = " << ExesEPSM << endl;
      }
      if (KEtot<=0.0) mout << "         Ratio= XXXX" << endl;
      else mout << "   Ratio lost = " << (ElostC+ElostE)/KEtot << endl;
      mout << "     3-D disp = " << disp[0] << ", " << disp[1] 
	   << ", " << disp[2] << endl;
      if (dmean>0.0) {
	mout << "   Disp ratio = " << disp[0]/dmean << ", " 
	     << disp[1]/dmean << ", " << disp[2]/dmean << endl << endl;
      }
	
      unsigned sumcells=0, sumbodies=0;
      mout << endl;
      mout << "-----------------------------------------------------" << endl;
      mout << "-----Cell/body diagnostics---------------------------" << endl;
      mout << "-----------------------------------------------------" << endl;
      mout << right << setw(8) << "Level" << setw(15) << "Scale(x)"
	   << setw(10) << "Cells" << setw(10) << "Bodies" << endl;
      mout << "-----------------------------------------------------" << endl;
      for (unsigned n=0; n<ncells.size(); n++) {
	mout << setw(8) << n << setw(15) << pHOT::sides[0]/(1<<n)
	     << setw(10) << ncells[n] << setw(10) << bodies[n]
	     << endl;
	sumcells  += ncells[n];
	sumbodies += bodies[n];
      }
      mout << "-----------------------------------------------------" << endl;
      mout << setw(8) << "TOTALS" 
	   << setw(15) << "**********"
	   << setw(10) << sumcells << setw(10) << sumbodies << endl;
      mout << "-----------------------------------------------------" << endl;
      mout << left << endl;
      
      vector<float>    keymake, xchange, convert, overlap, prepare;
      vector<float>    cupdate, scatter, repartn, tadjust, keycall;
      vector<float>    keycomp, keybods, waiton0, waiton1, waiton2;
      vector<float>    keynewc, keyoldc, treebar, diagdbg;
      vector<unsigned> numbods;

      c0->Tree()->Timing(keymake, xchange, convert, overlap, prepare,
			 cupdate, scatter, repartn, tadjust, keycall,
			 keycomp, keybods, waiton0, waiton1, waiton2,
			 keynewc, keyoldc, treebar, diagdbg, numbods);

      mout << "-----------------------------" << endl
	   << "Timing (secs) at mlevel="      << mlevel << endl
	   << "-----------------------------" << endl;

      outHeader0(mout);

      outHelper0(mout, "partition",    0, out, tot);
      outHelper0(mout, "partn wait",   8, out, tot);
      outHelper0(mout, "make tree",    1, out, tot);
      outHelper0(mout, "make wait",   10, out, tot);
      outHelper0(mout, "adjust tree",  2, out, tot);
      outHelper0(mout, "adjust cell",  3, out, tot);
      outHelper0(mout, "adjust wait", 11, out, tot);
      mout << endl;

      outHeader1(mout);

      outHelper1<float>(mout, "keymake", keymake);
      outHelper1<float>(mout, "xchange", xchange);
      outHelper1<float>(mout, "convert", convert);
      outHelper1<float>(mout, "overlap", overlap);
      outHelper1<float>(mout, "prepare", prepare);
      outHelper1<float>(mout, "cupdate", cupdate);
      outHelper1<float>(mout, "scatter", scatter);
      outHelper1<float>(mout, "repartn", repartn);
      outHelper1<float>(mout, "tadjust", tadjust);
      outHelper1<float>(mout, "keycall", keycall);
      outHelper1<float>(mout, "keycomp", keycomp);
      outHelper1<float>(mout, "keybods", keybods);
      outHelper1<float>(mout, "new key", keynewc);
      outHelper1<float>(mout, "old key", keyoldc);
      outHelper1<float>(mout, "diagnos", diagdbg);

      if (mpichk) {
	outHelper1<float>(mout, "wait #0", waiton0);
	outHelper1<float>(mout, "wait #1", waiton1);
	outHelper1<float>(mout, "wait #2", waiton2);
	outHelper1<float>(mout, "barrier", treebar);
      }
      outHelper1<unsigned int>(mout, "numbods", numbods);
      mout << endl;

      outHeader0(mout);

      outHelper0(mout, "timesteps", 4, out, tot);
      outHelper0(mout, "step list", 5, out, tot);
      outHelper0(mout, "collide  ", 6, out, tot);
      outHelper0(mout, "coll wait", 9, out, tot);
      outHelper0(mout, "overhead ", 7, out, tot);

      collide->tsdiag(mout);
      collide->voldiag(mout);

      //
      // Debugging usage
      //
      collide->EPSMtiming(mout);
      collide->CPUHog (mout);

      //
      // Collision effort
      //
      if (use_effort) {
	const int nlst = 15;
	const int lst[] = {0  ,   4,   9,  19,  49,  99, 199, 
			   499, 799, 899, 949, 979, 989, 994, 998};

				// Cumulate
	for (int i=1; i<ebins; i++) efrt[i] += efrt[i-1];
	if (efrt[ebins-1]>0) {
	  double norm = 1.0/efrt[ebins-1];
	  mout << "-----------------------------" << endl
	       << "Collision effort, mlevel="     << mlevel << endl
	       << "-----------------------------" << endl << left
	       << setw(12) << "Effort"
	       << setw(10) << "Count"
	       << setw(10) << "Sum"
	       << setw(12) << "Fraction"
	       << setw(12) << "Cumulate" << endl
	       << setw(12) << "------"
	       << setw(10) << "-----"
	       << setw(10) << "---"
	       << setw(12) << "--------"
	       << setw(12) << "--------" << endl;

	  
	  mout << setprecision(3) << scientific;

	  double ptile = exp(minEff + (maxEff - minEff)*(0.5+lst[0])/ebins);

	  mout << setw(12) << ptile
	       << setw(10) << efrt[lst[0]]
	       << setw(10) << efrt[lst[0]]
	       << setw(12) << norm*efrt[lst[0]]
	       << setw(12) << norm*efrt[lst[0]]
	       << endl;
	  
	  for (int i=1; i<nlst; i++) {
	    ptile = exp(minEff + (maxEff - minEff)*(0.5+lst[i])/ebins);
	    mout << setw(12) << ptile
		 << setw(10) << efrt[lst[i]] - efrt[lst[i-1]]
		 << setw(10) << efrt[lst[i]]
		 << setw(12) << norm*(efrt[lst[i]] - efrt[lst[i-1]])
		 << setw(12) << norm*efrt[lst[i]]
		 << endl;
	  }

	} else {
	  mout << "-----------------------------" << endl
	       << "Empty effort list, mlevel="    << mlevel << endl
	       << "-----------------------------" << endl;
	}
      }

      
      if (mlevel==0) {

	mout << "-----------------------------" << endl
	     << "Totals for the entire step   " << endl
	     << "-----------------------------" << endl;

	outHeader0(mout);

	outHelper0(mout, "partition",    0, OUT, TOT);
	outHelper0(mout, "partn wait",   8, OUT, TOT);
	outHelper0(mout, "make tree",    1, OUT, TOT);
	outHelper0(mout, "make wait",   10, OUT, TOT);
	outHelper0(mout, "adjust tree",  2, OUT, TOT);
	outHelper0(mout, "adjust cell",  3, OUT, TOT);
	outHelper0(mout, "adjust wait", 11, OUT, TOT);
	outHelper0(mout, "timesteps",    4, OUT, TOT);
	outHelper0(mout, "step list",    5, OUT, TOT);
	outHelper0(mout, "collide  ",    6, OUT, TOT);
	outHelper0(mout, "coll wait",    9, OUT, TOT);
	outHelper0(mout, "overhead ",    7, OUT, TOT);
	mout << endl;
      }

      //
      // Cell instance diagnostics
      //
      mout << "-----------------------------" << endl
	   << "Object counts at mlevel="      << mlevel << endl
	   << "-----------------------------" << endl
	   << " pCell # = " << pCellTot       << endl
	   << "-----------------------------" << endl;
      
    }

    //
    // Reset the collide counters (CollideIon)
    //
    collide->resetColls();


    //
    // Reset the timers
    //
    partnTime.reset();
    tree1Time.reset();
    tradjTime.reset();
    tcellTime.reset();
    tstepTime.reset();
    llistTime.reset();
    clldeTime.reset();
    timerDiag.reset();
    partnWait.reset();
    tree1Wait.reset();
    tree2Wait.reset();
    clldeWait.reset();

#ifdef USE_GPTL
    GPTLstop("UserTreeDSMC::collide_diag");
#endif

  }

#ifdef DEBUG
  if (!c0->Tree()->checkParticles(cout)) {
    cout << "Before level list: Particle check FAILED [" << right
	 << setw(3) << mlevel << ", " << setw(3) << myid << "]" << endl;
  }
#endif
  if (sampcel_debug) {
    ostringstream sout;
    sout << "before LEVEL LIST at level "
	 << "[" << mlevel << ", " << tnow  << "]";
    c0->Tree()->checkSampleCells(sout.str().c_str());
  }

  (*barrier)("TreeDSMC: after collision diags", __FILE__, __LINE__);

				// Remake level lists because particles
				// will (usually) have been exchanged 
				// between nodes
  llistTime.start();
  c0->reset_level_lists();
  llistTime.stop();

#ifdef DEBUG
  if (!c0->Tree()->checkParticles(cout)) {
    cout << "After level list: Particle check FAILED [" << right
	 << setw(3) << mlevel << ", " << setw(3) << myid << "]" << endl;

  }
#endif

#ifdef USE_GPTL
  GPTLstop("UserTreeDSMC::determine_acceleration_and_potential");
#endif

				// Debugging disk
				//
  // triggered_cell_body_dump(0.01, 0.002);
  // TempHisto();

  (*barrier)("TreeDSMC: end of accel routine", __FILE__, __LINE__);

}

void UserTreeDSMC::triggered_cell_body_dump(double time, double radius)
{
  static bool done = false;
  if (tnow<time) return;
  if (done) return;

  unsigned cnt=0;
  vector<double> p(3);
  pHOT_iterator c(*c0->Tree());

  for (unsigned n=0; n<c0->Tree()->Number(); n++) {
    c.nextCell();
    double r2 = 0.0;

    c.Cell()->MeanPos(p);
    for (unsigned k=0; k<3; k++) r2 += p[k]*p[k];

    if (r2 < radius*radius) {
      ostringstream ostr;
      ostr << outdir << runtag << ".testcell." << myid << "." << cnt++;
      ofstream out(ostr.str().c_str());

      // for (set<unsigned long>::iterator j=c.Cell()->bods.begin();
      for (vector<unsigned long>::iterator j=c.Cell()->bods.begin();
	   j!=c.Cell()->bods.end(); j++) {
	for (unsigned k=0; k<3; k++) 
	  out << setw(18) << c0->Tree()->Body(*j)->pos[k];
	for (unsigned k=0; k<3; k++) 
	  out << setw(18) << c0->Tree()->Body(*j)->vel[k];
	out << endl;
      }
    }
  }

  done = true;
}


void UserTreeDSMC::assignTempDensVol()
{
  const double f_H = 0.76;
  double mm = f_H*mp + (1.0-f_H)*4.0*mp;
  double KEtot, KEdsp, T;
  double Tfac;
  if (ctype == "LTE") Tfac = 2.0*UserTreeDSMC::Eunit/3.0 * mm/UserTreeDSMC::Munit/boltz;
  pCell *cell;
#ifdef DEBUG
  unsigned nbod=0, ntot, zbod=0, ztot, pcel=0, ptot;
  unsigned sing=0, stot, ceqs=0, qtot, zero=0, none;
  double n2=0.0, n1=0.0, N2, N1;
  double MinT, MaxT, MeanT, VarT;
  double minT=1e20, maxT=0.0, meanT=0.0, varT=0.0;
#endif

  pHOT_iterator it(*c0->Tree());

  while (it.nextCell()) {
    cell = it.Cell();
    //cell->sample->KE(KEtot, KEdsp);
    
    //T = KEdsp * Tfac;

    // Assign temp and/or density to particles
    //
    if (use_temp>=0 || use_dens>=0 || use_vol>=0) {
#ifdef DEBUG
      unsigned ssz = cell->sample->ctotal;
#endif
      unsigned csz = cell->ctotal;
      double  volm = cell->Volume();
      double  dens = cell->Mass()/volm;

      vector<unsigned long>::iterator j = cell->bods.begin();
      while (j != cell->bods.end()) {
	if (*j == 0) {
	  cout << "proc=" << myid << " id=" << id 
	       << " ptr=" << hex << cell << dec
	       << " indx=" << *j << "/" << csz << endl;
	  j++;
	  continue;
	}
        std::pair<int, int> sKey(cell->Body(j)->Z, cell->Body(j)->C);
	cell->sample->KE(sKey, KEtot, KEdsp);
	double mi = mp*atomic_weights[cell->Body(j)->Z];
	Tfac = 2.0*UserTreeDSMC::Eunit/3.0 * mi/UserTreeDSMC::Munit/boltz;
	T = KEdsp* Tfac;
        //if( T ==0) { cout << "UserTree T = 0 at j = " << *j << " KE = " << KEdsp << " Tfrac = " << Tfac << endl; }
	
	int sz = cell->Body(j)->dattrib.size();
	if (use_temp>=0 && use_temp<sz) 
	  cell->Body(j)->dattrib[use_temp] = T;
	if (use_dens>=0 && use_dens<sz) 
	  cell->Body(j)->dattrib[use_dens] = dens;
	if ((use_vol>=0) && use_vol<sz)
	  cell->Body(j)->dattrib[use_vol]  = volm;
	j++;
      }
#ifdef DEBUG
      if (T>0.0) {
	nbod += csz;
	minT = min<double>(T, minT);
	maxT = max<double>(T, maxT);
	meanT += csz*T;
	varT  += csz*T*T;
      } else {
	zbod += csz;
	if (ssz>1) {	// Sample cell has more than 1?
	  n1 += ssz;
	  n2 += ssz * ssz;
	  pcel++;
	  // Singletons sample cells
	} else if (ssz==1) {
	  sing++;
	  if (cell->sample == cell) ceqs++;
	} else {		// Empty cells
	  zero++;
	}
      }
#endif
    }
  }

#ifdef DEBUG
  MPI_Reduce(&nbod,  &ntot,  1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&zbod,  &ztot,  1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&pcel,  &ptot,  1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sing,  &stot,  1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&ceqs,  &qtot,  1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&zero,  &none,  1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&n1,    &N1,    1, MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&n2,    &N2,    1, MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&meanT, &MeanT, 1, MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&varT,  &VarT,  1, MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&minT,  &MinT,  1, MPI_DOUBLE,   MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&maxT,  &MaxT,  1, MPI_DOUBLE,   MPI_MAX, 0, MPI_COMM_WORLD);

  if (myid==0) {
    cout << setfill('-') << setw(70) << '-' << setfill(' ') << endl;
    cout << "Non-zero temperature assigned for " << ntot << " bodies"  << endl
	 << stot << " cells are singletons and " << qtot << " are roots" << endl
	 << none << " cells are empty" << endl
	 << "Zero temperature assigned for " << ztot << " bodies" << endl;
    if (ptot>1)
      cout << ", mean(N) = " << N1/ptot << endl
	   << " stdev(N) = " << sqrt( (N2 - N1*N1/ptot)/(ptot-1) ) << endl;
    cout << "MinT = " << MinT << endl << "MaxT = " << MaxT << endl;
    if (ntot>0)
      cout << " mean(T) = " << MeanT/ntot << endl;
    if (ntot>1) 
      cout << "stdev(T) = " 
	   << sqrt( (VarT - MeanT*MeanT/ntot)/(ntot-1) ) << endl;
    cout << setfill('-') << setw(70) << '-' << setfill(' ') << endl;
  }

  TempHisto();
#endif
}

void UserTreeDSMC::TempHisto()
{
  if (use_temp<0) return;

  pCell *cell;
  double Tlog, T, M, V, totalM1 = 0.0, totalM0;
  const int numT = 40;
  const double TlogMin = 3.0, TlogMax = 8.0;
  vector<double> td1(numT+2, 0), td0(numT+2);
  vector<double> vd1(numT+2, 0), vd0(numT+2);

  pHOT_iterator pit(*(c0->Tree()));

  while (pit.nextCell()) {
    cell = pit.Cell();
    // set<unsigned long>::iterator j = cell->bods.begin();
    vector<unsigned long>::iterator j = cell->bods.begin();
    V = cell->Volume();
    while (j != cell->bods.end()) {
      T = cell->Body(j)->dattrib[use_temp];
      if (T>0.0) {
	M = cell->Body(j)->mass;
	totalM1 += M;
	Tlog = log(T)/log(10.0);
	if (Tlog<TlogMin) {
	  td1[0] += M;
	  vd1[0] += V*M;
	} else if (Tlog>=TlogMax) {
	  td1[numT+1] += M;
	  vd1[numT+1] += V*M;
	} else {
	  int indx = floor((log(T)/log(10.0) - TlogMin) /
			   (TlogMax - TlogMin)*numT) + 1;
	  td1[indx] += M;
	  vd1[indx] += V*M;
	}
      }
      j++;
    }
  }

  
  MPI_Reduce(&totalM1, &totalM0, 1,      MPI_DOUBLE, MPI_SUM, 0, 
	     MPI_COMM_WORLD);

  MPI_Reduce(&td1[0],  &td0[0],  numT+2, MPI_DOUBLE, MPI_SUM, 0, 
	     MPI_COMM_WORLD);

  MPI_Reduce(&vd1[0],  &vd0[0],  numT+2, MPI_DOUBLE, MPI_SUM, 0, 
	     MPI_COMM_WORLD);

  if (myid==0) {

    for (int i=0; i<numT+2; i++) {
      if (td0[i]>0.0) vd0[i] /= td0[i];
    }

    cout << "----------------" << endl
	 << "Temperature dist" << endl
	 << "Time=" << tnow    << endl
	 << "----------------" << endl
	 << setw(10) << "<1000" 
	 << setw(10) << setprecision(2) << td0[0]/totalM0
	 << setw(10) << setprecision(2) << vd0[0] << endl;
    for (int i=0; i<numT; i++) 
      cout << setw(10) << setprecision(2) 
	   << pow(10.0, TlogMin + (TlogMax-TlogMin)/numT*(0.5+i))
	   << setw(10) << setprecision(2) << td0[i+1]/totalM0
	   << setw(10) << setprecision(2) << vd0[i+1] << endl;
    cout << setw(10) << ">1e8" << setw(10) 
	 << setprecision(2) << td0[numT+1]/totalM0
	 << setw(10) << setprecision(2) << vd0[numT+1] << endl;
  }
}

extern "C" {
  ExternalForce *makerTreeDSMC(string& line)
  {
    return new UserTreeDSMC(line);
  }
}

class proxytreedsmc { 
public:
  proxytreedsmc()
  {
    factory["usertreedsmc"] = makerTreeDSMC;
  }
};

proxytreedsmc p;
