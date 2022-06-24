#include <sys/timeb.h>
#include <stdlib.h>
#include <sys/types.h>
#include <getopt.h>
#include <time.h>

#include <iostream>
#include <iomanip>
#include <memory>
#include <random>
#include <fstream>
#include <sstream>

#include <expand.H>
#include <ExternalCollection.H>
#include <TreeDSMC.H>
#include <CollideLTE.H>
#include <CollideIon.H>

// #define DEBUG

#ifdef USE_GPTL
#include <gptl.h>
#endif

#include <NVTX.H>

//
// Debugging check
//
static bool sampcel_debug = false;
static bool levelst_debug = false;

//
// Version string for log file stamping
//
const std::string TreeDSMC::version = "0.49 [06/03/22 cuda trace RC]";

//
// Simulation units
//
double TreeDSMC::Lunit = 3.0e5*pc;
double TreeDSMC::Munit = 1.0e12*msun;
double TreeDSMC::Tunit = sqrt(Lunit*Lunit*Lunit/(Munit*6.67384e-08));
double TreeDSMC::Vunit = Lunit/Tunit;
double TreeDSMC::Eunit = Munit*Vunit*Vunit;
bool   TreeDSMC::use_effort = true;

std::vector<double> TreeDSMC::atomic_weights;

std::set<std::string> TreeDSMC:: colltypes;

TreeDSMC::TreeDSMC(const YAML::Node& conf) : ExternalForce(conf)
{
  (*barrier)("TreeDSMC: BEGIN construction", __FILE__, __LINE__);
  
  id = "TreeDSMC";		// ID string
  
				// Default parameter values
  ncell      = 7;		// 
  Ncell      = 64;
  cnum       = 0;
  madj       = 512;		// No tree pruning by default
  seed       = 11;		// Default seed for Collide
  epsm       = -1.0;
  hsdiam     = 1.0;
  crossfac   = 1.0;
  boxsize    = 1.0;
  boffset    = 0.0;
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
  use_Eint   = -1;
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
  sub_sample = true;
  treechk    = false;
  mpichk     = false;
  mfpts      = false;
  hybrid     = false;
  dumpHOT    = false;
  ageout     = 0;
  cntHOT     = 0;
  
  bSiz.resize(3, boxsize);
  bOff.resize(3, boffset);

  // Static initialization
  initialize_colltypes();

  // Initialize using input parameters
  initialize();
  
  // Set Collide random seed
  if (seed==0) {
    std::ifstream rin("/dev/urandom");
    rin.read((char *)&seed, sizeof(seed));
  }
  Collide::seed = seed;

  // Set cuda state
  //
  Collide::use_cuda = use_cuda;

  // Update derived units from Munit, Lunit, Tunit
  Vunit = Lunit/Tunit;
  Eunit = Munit*Vunit*Vunit;
  
  // Report TreeDSMC version to stdout log
  //
  if (myid==0) {
    print_divider();
    cout << "** User routine TreeDSMC, version " << version << std::endl;
    print_divider();
  }
  
  // Initialize the atomic_weights map hardcode the atomic weight map
  // for use in collFrac.  Weights in atomic mass units (amu)
  
  atomic_weights.resize(15, -1.0);

  atomic_weights[0]  = 0.000548579909; // Electron
  atomic_weights[1]  = 1.0079;	       // Hydrogen
  atomic_weights[2]  = 4.0026;	       // Helium
  atomic_weights[3]  = 6.941;	       // Lithium
  atomic_weights[4]  = 9.0122;	       // Beryllium
  atomic_weights[5]  = 10.811;	       // Boron
  atomic_weights[6]  = 12.011;	       // Carbon
  atomic_weights[7]  = 14.007;	       // Nitrogen
  atomic_weights[8]  = 15.999;	       // Oxygen
  atomic_weights[9]  = 18.998;	       // Florine
  atomic_weights[10] = 20.180;	       // Neon
  atomic_weights[11] = 22.990;	       // Sodium
  atomic_weights[12] = 24.305;	       // Magnesium
  atomic_weights[13] = 26.982;	       // Aluminium
  atomic_weights[14] = 28.085;	       // Silicon
  
  // add in atomic weights for any other higher, more tracer, species
  
  // Look for the fiducial component
  bool found = false;
  for (auto cc : comp->components) {
    if ( !comp_name.compare(cc->name) ) {
      c0 = cc;
      found = true;
      break;
    }
  }
  
  if (!found) {
    if (myid==0) {
      std::cerr << "TreeDSMC: process " << myid 
		<< ": can't find fiducial component <" << comp_name << ">"
		<< endl;
    }
    MPI_Finalize();
    exit(35);
  }
  
  (*barrier)("TreeDSMC: BEFORE use checks", __FILE__, __LINE__);
  
  //
  // Get use_key postion index from the fiducial component; this will
  // be checked for consistency in the Collide constructor
  //
  use_key = c0->keyPos;
  
  //
  // Sanity check on excess attribute if excess calculation is
  // desired
  //
  if (use_exes>=0) {
    
    int ok1 = 1, ok;
    
    PartMapItr p    = c0->Particles().begin();
    PartMapItr pend = c0->Particles().end();
    for (; p!=pend; p++) {
      if (use_exes >= static_cast<int>(p->second->dattrib.size())) {
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
	cout << "TreeDSMC: excess calculation requested but some" << endl
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
      if (use_Kn >= static_cast<int>(p->second->dattrib.size())) {
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
	std::cout 
	  << "TreeDSMC: Knudsen number calculation requested but some" 
	  << std::endl
	  << "particles have incompatible float attribute counts." 
	  << std::endl
	  << "Attribute #" << use_Kn 
	  << ". Continuing without Knudsen number computation."
	  << std::endl;
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
      if (use_St >= static_cast<int>(p->second->dattrib.size())) {
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
	std::cout 
	  << "TreeDSMC: Strouhal number calculation requested but some" 
	  << std::endl
	  << "particles have incompatible float attribute counts." 
	  << std::endl
	  << "Attribute #" << use_St 
	  << ". Continuing without Strouhal number computation."
	  << std::endl;
      }
      use_St = -1;
    }
  }
  
  //
  // Make the initial species map
  //
  if (use_key>=0) {
    
    (*barrier)("TreeDSMC: BEFORE species map construction", __FILE__, __LINE__);
    makeSpeciesMap();		// Compute fractions in each species
    (*barrier)("TreeDSMC: AFTER species map construction",  __FILE__, __LINE__);
    
    typedef std::map<speciesKey, unsigned long> spCountMap;
    typedef spCountMap::iterator  spCountMapItr;
    
    //
    // Get element list
    //
    typedef std::set<unsigned short> elemL;
    elemL elems;
    for (spCountMapItr it=spec.begin(); it != spec.end(); it++)
      elems.insert(it->first.first);
    
    //
    // Create all possible species
    //
    for (elemL::iterator it=elems.begin(); it != elems.end(); it++)  {
      if (hybrid) {
	speciesKey indx(*it, 0);
	spec_list.insert(indx);
	collFrac[indx] = 1.0;
      } else {
	for (unsigned short C=0; C<*it+2; C++) {
	  speciesKey indx(*it, C);
	  spec_list.insert(indx);
	  collFrac[indx] = 1.0;
	}
      }
    }
    
    if (myid==0) {
      //
      // Make total count array
      //
      std::vector< std::vector<unsigned long> > specG(numprocs);
      std::vector<unsigned long> specT(spec.size(), 0);
      for (int n=0; n<numprocs; n++) specG[n].resize(spec.size(), 0);
      
      // Add root's counts
      //
      size_t j = 0;
      spCountMapItr it2 = spec1.begin();
      for (spCountMapItr it=spec.begin(); it!=spec.end(); it++) {
	if (it->first == it2->first) specG[0][j] = it2++->second;
	j++;
      }
      
      // Get data from nodes
      //
      for (int n=1; n<numprocs; n++) {
	MPI_Recv(&specG[n][0], spec.size(), MPI_UNSIGNED_LONG, n, 8,
		 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
      
      // Make output table
      //
      cout << endl
	   << "--------------" << endl
	   << "Species counts" << endl
	   << "--------------" << endl
	   << endl;
      
      cout << setw(4) << right << "#";
      for (spCountMapItr it=spec.begin(); it != spec.end(); it++)
	cout << setw(8) << right 
	     << "(" << it->first.first << "," << it->first.second << ")";
      cout << endl;
      
      cout << setw(4) << right << "---";
      for (spCountMapItr it=spec.begin(); it != spec.end(); it++)
	cout << setw(12) << right << "--------";
      cout << endl;
      
      for (int n=0; n<numprocs; n++) {
	cout << setw(4) << right << n;
	for (size_t j=0; j<spec.size(); j++) {
	  cout << setw(12) << right << specG[n][j];
	  specT[j] += specG[n][j];
	}
	cout << endl;
      }
      
      cout << setw(4) << right << "---";
      for (spCountMapItr it=spec.begin(); it != spec.end(); it++)
	cout << setw(12) << right << "--------";
      cout << endl;
      
      cout << setw(4) << right << "TOT";
      for (size_t j=0; j<spec.size(); j++) 
	cout << setw(12) << right << specT[j];
      cout << endl << endl;
      
    } else {
      
      std::vector<unsigned long> specT(spec.size(), 0);
      
      // Pack node's data spec counts to send to root
      //
      size_t j = 0;
      spCountMapItr it2 = spec1.begin();
      for (spCountMapItr it=spec.begin(); it!=spec.end(); it++) {
	if (it->first == it2->first) specT[j] = it2++->second;
	j++;
      }
      
      // Send . . . 
      //
      MPI_Send(&specT[0], spec.size(), MPI_UNSIGNED_LONG, 0, 8,
	       MPI_COMM_WORLD);
    }
    
  } else {		   
    spec_list.insert(Particle::defaultKey);
    collFrac[Particle::defaultKey] = 1.0;
  }
  
  // Number of protons per mass unit
  for (auto it=collFrac.begin(); it!=collFrac.end(); it++) 
    it->second *= Munit/amu;
  
  pHOT::sub_sample = sub_sample;
  
  c0->HOTcreate(spec_list);
  
  c0->Tree()->setSides (bSiz[0], bSiz[1], bSiz[2]);
  c0->Tree()->setOffset(bOff[0], bOff[1], bOff[2]);
  
  pCell::bucket = ncell;
  pCell::Bucket = Ncell;
  
  volume = pHOT::sides[0] * pHOT::sides[1] * pHOT::sides[2];
  
  //
  // Set collision parameters
  //
  Collide::NTC            = ntc;	   
  Collide::PULLIN    	  = use_pullin; 
  Collide::ESOL      	  = esol;	   
  Collide::EPSMratio 	  = epsm;	   
  Collide::DRYRUN    	  = dryrun;	   
  Collide::NOCOOL    	  = nocool;	   
  Collide::TSDIAG    	  = tsdiag;	   
  Collide::VOLDIAG   	  = voldiag;	   
  Collide::TSPOW     	  = tspow;	   
  Collide::MFPDIAG   	  = mfpstat;	   
  Collide::EFFORT    	  = use_effort; 
  Collide::ENHANCE   	  = enhance;
  
  //
  // Create the collision instance from the allowed list
  //
  if (ctype.compare("LTE") == 0)
    collide = new CollideLTE(this, c0, hsdiam, crossfac, nthrds);
  if (ctype.compare("Ion") == 0)
    collide = new CollideIon(this, c0, hsdiam, crossfac, nthrds);
  else {
    std::ostringstream sout;
    sout << "No such Collide type: " << ctype;
    throw GenericError(sout.str(), __FILE__, __LINE__);
  }
  
  //
  // Set additional parameters
  //
  
  collide->set_temp_dens(use_temp, use_dens);
  
  if (esol) collide->set_timestep(-1);
  else      collide->set_timestep(use_delt);
  
  collide->set_key   (use_key);
  collide->set_Eint  (use_Eint);
  collide->set_Kn    (use_Kn);
  collide->set_St    (use_St);
  collide->set_excess(use_exes);
  collide->set_MFPTS (mfpts);
  
  ElostTotCollide = ElostTotEPSM = 0.0;
  
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

TreeDSMC::~TreeDSMC()
{
  delete collide;
}

void TreeDSMC::userinfo()
{
  if (myid) return;		// Return if node master node
  
  print_divider();
  
  // Report version and key parameters
  //
  cout << "** User routine TreeDSMC initialized, " 
       << "Lunit=" << Lunit << ", Tunit=" << Tunit << ", Munit=" << Munit
       << ", Vunit=" << Vunit << ", Eunit=" << Eunit
       << ", cnum=" << cnum << ", hsdiam=" << hsdiam 
       << ", crossfac=" << crossfac << ", madj=" << madj 
       << ", epsm=" << epsm << ", ncell=" << ncell << ", Ncell=" << Ncell 
       << ", boxsize=[" << bSiz[0] << ", " << bSiz[1] << ", " << bSiz[2] << "]"
       << ", boffset=[" << bOff[0] << ", " << bOff[1] << ", " << bOff[2] << "]"
       << ", compname=" << comp_name
       << ", keyPos=" << use_key << ", vol(tree)=" << volume;
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
  if (use_Eint>=0) cout << ", internal energy at pos=" << use_Eint;
  if (use_pullin)  cout << ", Pullin algorithm enabled";
  if (dryrun)      cout << ", collisions disabled";
  if (nocool)      cout << ", cooling disabled";
  if (epsm>0.0)    cout << ", using EPSM";
  else             cout << ", EPSM disabled";
  if (ntc)         cout << ", using NTC";
  else             cout << ", NTC disabled";
  if (cba && cbadiag)     
    cout << " with diagnostics";
  if (use_effort)  cout << ", with effort-based load";
  else             cout << ", with uniform load";
  if (fabs(enhance-1.0)>1.0e-6)
    cout << ", with enhanced cooling of " << enhance;
  if (use_multi) {
    cout << ", multistep enabled";
    if (use_delt>=0) 
      cout << ", time step at pos=" << use_delt << ", coolfrac=" << coolfrac;
  }
  
  cout << ", collsion type=" << ctype;
  
  if (ctype.compare("Ion") == 0) {
    cout << ", rr type=" << Ion::getRRtype();
    if (CollideIon::equiptn)
      cout << ", e-i equipartition is ON";
    else
      cout << ", e-i equipartition is OFF";
    if (CollideIon::scatter)
      cout << ", scattering only is ON";
    else
      cout << ", scattering only is OFF";
  }
  
  cout << endl;
  
  print_divider();
}

template<typename Out>
void split(const std::string &s, char delim, Out result) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    *(result++) = item;
  }
}

std::vector<std::string> split(const std::string &s, char delim)
{
  std::vector<std::string> elems;
  split(s, delim, std::back_inserter(elems));
  return elems;
}

std::pair<int, int> splitPair(const std::string &s, char delim)
{
  std::vector<std::string> i = split(s, delim);
  if (i.size() != 2) {
    std::cerr << "Error spliting indices: found string <" << s << "> and delimiter=" << delim << std::endl;
    exit(-1);
  }
  return std::pair<int, int>(std::stoi(i[0]), std::stoi(i[1]));
}

void TreeDSMC::initialize()
{
  try {
    if (conf["Lunit"])          Lunit              = conf["Lunit"].as<double>();
    if (conf["Tunit"])          Tunit              = conf["Tunit"].as<double>();
    if (conf["Munit"])          Munit              = conf["Munit"].as<double>();
    if (conf["cnum"])           cnum               = conf["cnum"].as<int>();
    if (conf["madj"])           madj               = conf["madj"].as<int>();
    if (conf["epsm"])           epsm               = conf["epsm"].as<double>();
    if (conf["hsdiam"])         hsdiam             = conf["hsdiam"].as<double>();
    if (conf["crossfac"])       crossfac           = conf["crossfac"].as<double>();
    if (conf["coolfrac"])       coolfrac           = conf["coolfrac"].as<double>();
    if (conf["enhance"])        enhance            = conf["enhance"].as<double>();
    if (conf["nsteps"])         nsteps             = conf["nsteps"].as<int>();
    if (conf["msteps"])         msteps             = conf["msteps"].as<int>();
    if (conf["ncell"])          ncell              = conf["ncell"].as<int>();
    if (conf["Ncell"])          Ncell              = conf["Ncell"].as<int>();
    if (conf["compname"])       comp_name          = conf["compname"].as<string>();
    if (conf["use_temp"])       use_temp           = conf["use_temp"].as<int>();
    if (conf["use_dens"])       use_dens           = conf["use_dens"].as<int>();
    if (conf["use_delt"])       use_delt           = conf["use_delt"].as<int>();
    if (conf["use_Kn"])         use_Kn             = conf["use_Kn"].as<int>();
    if (conf["use_St"])         use_St             = conf["use_St"].as<int>();
    if (conf["use_vol"])        use_vol            = conf["use_vol"].as<int>();
    if (conf["use_exes"])       use_exes           = conf["use_exes"].as<int>();
    if (conf["use_Eint"])       use_Eint           = conf["use_Eint"].as<int>();
    if (conf["frontier"])       frontier           = conf["frontier"].as<bool>();
    if (conf["tspow"])          tspow              = conf["tspow"].as<int>();
    if (conf["tsdiag"])         tsdiag             = conf["tsdiag"].as<bool>();
    if (conf["voldiag"])        voldiag            = conf["voldiag"].as<bool>();
    if (conf["mfpstat"])        mfpstat            = conf["mfpstat"].as<bool>();
    if (conf["cbadiag"])        cbadiag            = conf["cbadiag"].as<bool>();
    if (conf["dryrun"])         dryrun             = conf["dryrun"].as<bool>();
    if (conf["nocool"])         nocool             = conf["nocool"].as<bool>();
    if (conf["use_multi"])      use_multi          = conf["use_multi"].as<bool>();
    if (conf["use_pullin"])     use_pullin         = conf["use_pullin"].as<bool>();
    if (conf["use_effort"])     use_effort         = conf["use_effort"].as<bool>();
    if (conf["esol"])           esol               = conf["esol"].as<bool>();
    if (conf["cba"])            cba                = conf["cba"].as<bool>();
    if (conf["ntc"])            ntc                = conf["ntc"].as<bool>();
    if (conf["sub_sample"])     sub_sample         = conf["sub_sample"].as<bool>();
    if (conf["treechk"])        treechk            = conf["treechk"].as<bool>();
    if (conf["mpichk"])         mpichk             = conf["mpichk"].as<bool>();
    if (conf["mfpts"])          mfpts              = conf["mfpts"].as<bool>();
    if (conf["hybrid"])         hybrid             = conf["hybrid"].as<bool>();
    if (conf["dumpHOT"])        dumpHOT            = conf["dumpHOT"].as<bool>();
    if (conf["ageout"])         ageout             = conf["ageout"].as<int>();

    NTC::NTCdb::maxAge = ageout;

    if (conf["ntc_chkpt"])      NTC::NTCdb::intvl  = conf["ntc_chkpt"].as<int>();
    if (conf["ntc_verbose"])    NTC::NTCdb::chatty = conf["ntc_verbose"].as<bool>();
  
    if (conf["rrtype"]) {
      std::string val = conf["rrtype"].as<std::string>();
      if (Ion::setRRtype(val)) {
	std::ostringstream sout;
	sout << "TreeDSMC: invalid rrtype <" << val << ">";
	throw GenericError(sout.str(), __FILE__, __LINE__);
      }
    }
    
    if (conf["ctype"]) {
      std::string val = conf["ctype"].as<std::string>();
      if (check_ctype(val)) ctype = val;
      else {
	std::ostringstream sout;
	sout << "TreeDSMC: invalid ctype <" << ctype << ">";
	throw GenericError(sout.str(), __FILE__, __LINE__);
      }
    }
  
    if (conf["boxsize"]) {      // Check for scalar
      if (conf["boxsize"].IsScalar()) {
	double v = conf["boxsize"].as<double>();
	for (size_t k=0; k<3; k++) bSiz[k] = v;
      }	                        // Check for vector
      else if (conf["boxsize"].IsSequence()) {
	auto vz = conf["boxsize"].as<std::vector<double>>();
	if (vz.size()) {
	  for (size_t k=0; k<3; k++) {
	    size_t indx = std::min<size_t>(k, vz.size()-1);
	    bSiz[k] = vz[indx];
	  }
	}
      }	                        // Unexpected type
      else {
	std::ostringstream sout;
	sout << "TreeDSMC: invalid boxsize";
	throw GenericError(sout.str(), __FILE__, __LINE__);
      }
    }


    if (conf["boffset"]) {	// Check for scalar
      if (conf["boffset"].IsScalar()) {
	double v = conf["boffset"].as<double>();
	for (size_t k=0; k<3; k++) bOff[k] = v;
      }	                        // Check for vector
      else if (conf["boffset"].IsSequence()) {
	auto vz = conf["boffset"].as<std::vector<double>>();
	if (vz.size()) {
	  for (size_t k=0; k<3; k++) {
	    size_t indx = std::min<size_t>(k, vz.size()-1);
	    bOff[k] = vz[indx];
	  }
	}
      } else {                  // Unexpected type
	std::ostringstream sout;
	sout << "TreeDSMC: invalid boffset";
	throw GenericError(sout.str(), __FILE__, __LINE__);
      }
    }

    // Look for array values in the parameter string of the form
    // spc(1,2)=3.1, spc(3,4)=5.6, etc.
    //
    if (ctype.compare("LTE")==0) {
    
      for (YAML::const_iterator it=conf.begin(); it!=conf.end(); it++) {
	std::string key = it->first.as<std::string>();
	
	size_t pos = key.find("spc(");
	
	if (pos == 0) {
	  size_t bpos = key.find("(");
	  size_t epos = key.find(")");
	  if (bpos != std::string::npos and epos != std::string::npos) {
	    std::pair<int, int> pr = splitPair(key.substr(bpos+1, epos-bpos-1), ',');
	    speciesKey p(pr.first, pr.second);
	    collFrac[p] = it->second.as<double>();
	  }
	} else {
	  collFrac[Particle::defaultKey] = 1.0;
	}
      }
    }
    
    // Ion method specific parameters
    //
    if (ctype.compare("Ion")==0) {
      if (conf["equiptn"]) CollideIon::equiptn  = conf["equiptn"].as<bool>();
      if (conf["scatter"]) CollideIon::scatter  = conf["scatter"].as<bool>();
      if (conf["nodelC"])  CollideIon::NoDelC   = conf["nodelC"].as<int>();
      if (conf["logL"])    CollideIon::logL     = conf["logL"].as<double>();
      if (conf["config"])  CollideIon::config   = conf["config"];
    }
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in TreeDSMC: "
			   << error.what()         << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                << std::endl
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }
}


//  __  __    _    ___ _   _   ____   ___  _   _ _____ ___ _   _ _____ 
// |  \/  |  / \  |_ _| \ | | |  _ \ / _ \| | | |_   _|_ _| \ | | ____|
// | |\/| | / _ \  | ||  \| | | |_) | | | | | | | | |  | ||  \| |  _|  
// | |  | |/ ___ \ | || |\  | |  _ <| |_| | |_| | | |  | || |\  | |___ 
// |_|  |_/_/   \_\___|_| \_| |_| \_\\___/ \___/  |_| |___|_| \_|_____|
//                                                                    

void TreeDSMC::determine_acceleration_and_potential(void)
{
  static unsigned cat = 0;	// For debugging sync
  {
    std::ostringstream sout;
    sout << "TreeDSMC: call=" << ++cat;
    (*barrier)(sout.str(), __FILE__, __LINE__);
  }
  
#ifdef USE_GPTL
  GPTLstart("TreeDSMC::determine_acceleration_and_potential");
#endif
  
  static bool firstime = true;
  static unsigned nrep = 0;
  
  nvTracerPtr tPtr, tPtr2;
  if (cuda_prof) tPtr = nvTracerPtr(new nvTracer("TreeDSMC: init"));

  //
  // Only compute DSMC when passed the fiducial component
  //
  
  if (cC != c0) {
#ifdef USE_GPTL
    GPTLstop("TreeDSMC::determine_acceleration_and_potential");
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
    
    //
    // Extreme debugging
    //
    if (levelst_debug) {
      
      std::ostringstream sout;
      sout << __FILE__ << ":" << __LINE__ << ", Node " << myid 
	   << ": ERROR bodycell after tree build , T=" << tnow 
	   << " at level [" << mlevel << "]";
      
      if (!c0->Tree()->checkBodycell(sout.str())) {
	cout << "Process " << myid << ": "
	     << "makeTree completed: body cell check FAILED!" << endl;
      }    
      if (!c0->Tree()->checkParticles(cout, sout.str())) {
	cout << "Process " << myid << ": "
	     << "makeTree completed: initial particle check FAILED!" << endl;
      }    
      if (!c0->Tree()->checkFrontier(cout, sout.str())) {
	cout << "Process " << myid << ": "
	     << "makeTree completed: frontier check FAILED!" << endl;
      }
      
      if (!c0->Tree()->checkKeybods(sout.str())) {
	cout << "Process " << myid 
	     << ": makeTree: ERROR particle key not in keybods AFTER makeTree(), T=" 
	     << tnow << endl;
      }
    }
    
    if (sampcel_debug) {
      std::ostringstream sout;
      sout << __FILE__ << ":" << __LINE__ << ", Node " << myid 
	   << ": AFTER makeTree(), first time, "
	   << "[" << 0 << ", " << tnow  << "]";
      
      c0->Tree()->checkSampleCells(sout.str().c_str());
      c0->Tree()->logFrontierStats();
    }
    
    if (use_temp || use_dens || use_vol) assignTempDensVol();
    
    stepnum = 0;
    tM0 = curtime = tnow;
    
#ifdef DEBUG
    cout << "Computed partition and tree [firstime on #" 
	 << setw(4) << left << myid << "]" << endl;
#endif
  } else {
    
    if (tnow-curtime < 1.0e-14) {
      if (myid==0) {
	cout << "TreeDSMC: attempt to redo step at T=" << tnow << endl;
      }
#ifdef USE_GPTL
      GPTLstop("TreeDSMC::determine_acceleration_and_potential");
#endif
      return; 			// Don't do this time step again!
    }
    
    stepnum++;
    curtime = tnow;
  }
  
				// Time the entire step
  fullstep.Start();

#ifdef DEBUG
  c0->Tree()->densCheck();
#endif
  
#ifdef DEBUG
  if (!c0->Tree()->checkParticles(cout, "TreeDSMC: after init")) {
    cout << "After init only: Particle check FAILED [" << right
	 << setw(3) << mlevel << ", " << setw(3) << myid << "]" << endl;
  }
#endif
  
  (*barrier)("TreeDSMC: after cell computation", __FILE__, __LINE__);
  
  //
  // Only run diagnostics every nsteps
  //
  bool diagstep = nsteps>0 && (stepnum % nsteps == 0);
  
  //
  // Diagnostics run at levels <= msteps (takes prececdence over nsteps)
  //
  if (msteps>=0 && diagstep) 
    diagstep = (mlevel <= static_cast<unsigned>(msteps)) ? true : false;
  
  overhead.Start();
  
  if (cuda_prof) {
    tPtr.reset();
    tPtr = nvTracerPtr(new nvTracer("TreeDSMC: (re)partition"));
  }

  //
  // Sort the particles into cells
  //
  if (mlevel<=madj) {
    
#ifdef USE_GPTL
    GPTLstart("TreeDSMC::pHOT");
    GPTLstart("TreeDSMC::waiting");
    (*barrier)("TreeDSMC: pHOT waiting");
    GPTLstop ("TreeDSMC::waiting");
    GPTLstart("TreeDSMC::repart");
#endif
    
    (*barrier)("TreeDSMC: after pHOT wait", __FILE__, __LINE__);
    
    partnTime.start();
    c0->Tree()->Repartition(mlevel); nrep++;
    partnSoFar += partnTime.stop();
    
    partnWait.start();
    (*barrier)("TreeDSMC: after repartition", __FILE__, __LINE__);
    waitpSoFar += partnWait.stop();
    
#ifdef USE_GPTL
    GPTLstop ("TreeDSMC::repart");
    GPTLstart("TreeDSMC::makeTree");
#endif
    tree1Time.start();
    if (cuda_prof) {
      tPtr2.reset();
      tPtr2 = nvTracerPtr(new nvTracer("TreeDSMC: makeTree"));
    }
    c0->Tree()->makeTree();
    if (cuda_prof) {
      tPtr2.reset();
      tPtr2 = nvTracerPtr(new nvTracer("TreeDSMC: makeCellLevelList"));
    }
    c0->Tree()->makeCellLevelList();
    tree1Time.stop();
    tree1Wait.start();
    (*barrier)("TreeDSMC: after makeTree", __FILE__, __LINE__);
    wait1SoFar += tree1Wait.stop();
#ifdef USE_GPTL
    GPTLstop ("TreeDSMC::makeTree");
    GPTLstart("TreeDSMC::pcheck");
#endif
    tree1Time.start();
#ifdef DEBUG
    cout << "Made partition, tree and level list [" << mlevel << "]" << endl;
    if (!c0->Tree()->checkParticles(cout, "TreeDSMC: after partition")) {
      cout << "Particle check on new tree FAILED [" << mlevel << "]" << endl;
    }
#endif
    if (sampcel_debug) {
      std::ostringstream sout;
      sout << __FILE__ << ":" << __LINE__ << ", Node " << myid 
	   << ": after makeTree(), T=" << tnow 
	   << " at level [" << mlevel << "]";
      
      c0->Tree()->checkSampleCells(sout.str().c_str());
      c0->Tree()->logFrontierStats();
    }
    tree1SoFar += tree1Time.stop();
    
#ifdef USE_GPTL
    GPTLstop("TreeDSMC::pcheck");
    GPTLstop("TreeDSMC::pHOT");
#endif
    
  } else {
    
#ifdef USE_GPTL
    GPTLstart("TreeDSMC::pHOT_2");
    GPTLstart("TreeDSMC::adjustTree");
#endif
    
#ifdef DEBUG
    if (myid==0)
      cout << "About to adjust tree [" << mlevel << "]" << endl;
#endif
    tradjTime.start();
    if (cuda_prof) {
      tPtr2.reset();
      tPtr2 = nvTracerPtr(new nvTracer("TreeDSMC: adjustTree"));
    }
    c0->Tree()->adjustTree(mlevel);
    tradjSoFar += tradjTime.stop();
    tcellTime.start();
    if (cuda_prof) {
      tPtr2.reset();
      tPtr2 = nvTracerPtr(new nvTracer("TreeDSMC: adjustCellLevelList"));
    }
    c0->Tree()->adjustCellLevelList(mlevel);
    tcellSoFar += tcellTime.stop();
    tree2Wait.start();
    (*barrier)("TreeDSMC: after adjustTree", __FILE__, __LINE__);
    wait2SoFar += tree2Wait.stop();
    
    if (sampcel_debug) {
      std::ostringstream sout;
      sout << __FILE__ << ":" << __LINE__ << ", Node " << myid 
	   << ", after adjustTree() at level "
	   << "[" << mlevel << ", " << tnow  << "]";
      
      c0->Tree()->checkSampleCells(sout.str().c_str());
    }
    
#ifdef USE_GPTL
    GPTLstop("TreeDSMC::adjustTree");
    GPTLstop("TreeDSMC::pHOT_2");
#endif
    
  }
  
  if (cuda_prof) tPtr2.reset();

  overhead.Stop();
  pot_time += overhead.getTime();
  pot_time /= max<unsigned>(1, c0->Number());
  
  if (use_effort) {
    PartMapItr pitr = c0->Particles().begin(), pend = c0->Particles().end();
    for (; pitr!= pend; pitr++) pitr->second->effort = pot_time;
  }
  
  //
  // Evaluate collisions among the particles
  //
  
  clldeTime.start();
  if (cuda_prof) {
    tPtr.reset();
    tPtr = nvTracerPtr(new nvTracer("TreeDSMC: collide"));
  }
  
  if (0) {
    ostringstream sout;
    sout << "before Collide [" << nrep << "], " 
	 << __FILE__ << ": " << __LINE__;
    c0->Tree()->checkBounds(2.0, sout.str().c_str());
  }
  
#ifdef USE_GPTL
  GPTLstart("TreeDSMC::collide");
#endif
  
  //
  // So far, all computations have been about repartition and
  // tessellation.  All of the collision stuff is done by the current
  // Collide class instance.
  //
  
  (*barrier)("TreeDSMC: BEFORE Collide::collide", __FILE__, __LINE__);
  
  const Collide::UU&
    CC = collide->collide(*c0->Tree(), collFrac, mlevel, diagstep);
  
  (*barrier)("TreeDSMC: AFTER Collide::collide", __FILE__, __LINE__);
  
  // Collide:UU is a std::tuple<unsigned, unsigned>
  //
  // Tuple_0 is the body count at this level
  // Tuple_1 is the collision count
  
  collideSoFar += clldeTime.stop();
  
  clldeWait.start();
  (*barrier)("TreeDSMC: after collide", __FILE__, __LINE__);
  
#ifdef USE_GPTL
  GPTLstop("TreeDSMC::collide");
#endif
  
  waitcSoFar += clldeWait.stop();
  
  // -----------------
  // Time step request
  // -----------------
  //
  // New timesteps are selected for the cells based on the collision
  // diagnostics from the current step.
  //
  
  if (cuda_prof) {
    tPtr.reset();
    tPtr = nvTracerPtr(new nvTracer("TreeDSMC: finish"));
  }
#ifdef USE_GPTL
  GPTLstart("TreeDSMC::collide_timestep");
#endif
  
  (*barrier)("TreeDSMC: before collide timestep", __FILE__, __LINE__);
  
  if (cuda_prof) tPtr2 = nvTracerPtr(new nvTracer("TreeDSMC: timesteps"));

  tstepTime.start();
  if (use_multi) collide->compute_timestep(coolfrac);
  tstepSoFar += tstepTime.stop();
  
#ifdef USE_GPTL
  GPTLstop("TreeDSMC::collide_timestep");
#endif
  
  if (cuda_prof) tPtr2.reset();

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
  
  
  // Remake level lists because particles will (usually) have been
  // exchanged between nodes
  //
  // Begin with extreme debugging, if enabled
  //
  if (levelst_debug) {
    
    std::ostringstream sout;
    sout << __FILE__ << ":" << __LINE__ << ", Node " << myid 
	 << ": ERROR bodycell after remaking level lists , T=" << tnow 
	 << " at level [" << mlevel << "]";
    
    if (!c0->Tree()->checkParticles(cout, sout.str())) {
      cout << "Before level list: Particle check FAILED [" << right
	   << setw(3) << mlevel << ", " << setw(3) << myid << "]" << endl;
    }
  }
  
  if (sampcel_debug) {
    std::ostringstream sout;
    sout << __FILE__ << ":" << __LINE__ << ", Node " << myid 
	 << ", before LEVEL LIST at level "
	 << "[" << mlevel << ", " << tnow  << "]";
    c0->Tree()->checkSampleCells(sout.str().c_str());
  }
  
  llistTime.start();
  c0->reset_level_lists();
  llistTime.stop();
  
  // Finish with extreme debugging, if enabled
  //
  if (levelst_debug) {
    std::ostringstream sout;
    sout << __FILE__ << ":" << __LINE__ << ", Node " << myid 
	 << ": ERROR bodycell after resetting level lists , T=" << tnow 
	 << " at level [" << mlevel << "]";
    
    if (!c0->Tree()->checkParticles(cout, sout.str())) {
      cout << "After level list: Particle check FAILED [" << right
	   << setw(3) << mlevel << ", " << setw(3) << myid << "]" << endl;
      
    }
  }
  
  {
    std::ostringstream sout;
    sout << "TreeDSMC: BEFORE diagstep stanza, T=" << tnow
	 << ", bods=" << std::get<0>(CC);
    (*barrier)(sout.str(), __FILE__, __LINE__);
  }
  
  if (cuda_prof) {
    tPtr2.reset();
    tPtr2 = nvTracerPtr(new nvTracer("TreeDSMC: diagnostics"));
  }

  //
  // Periodically display the current progress
  //
  // Only makes sense to perform the following diagnostics and updates
  // if there are particles in cells whose interactions will be
  // computed at this level
  //
  // Lots of diagnostics are computed and emitted here . . .
  //
  if (diagstep && std::get<0>(CC)>0) {
    
    (*barrier)("TreeDSMC: ENTERING diagstep stanza", __FILE__, __LINE__);
    
#ifdef USE_GPTL
    GPTLstart("TreeDSMC::collide_diag");
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
	      << setw(14) << "Coll/cell"
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
    
    double mm = amu * collide->molWeight();
    double meanT = 0.0;
    if (Mtotl>0.0) meanT = 2.0*KEtot/Mtotl*Eunit/3.0 * mm/Munit/boltz;
    
    unsigned cellBods = c0->Tree()->checkNumber();
    unsigned oobBods  = c0->Tree()->oobNumber();
    
    double Mass = 0.0;
    unsigned Counts = 0;
    c0->Tree()->totalMass(Counts, Mass);
    c0->Tree()->computeCellStates();
    
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
    timerSoFar += timerDiag.stop();
    
    const unsigned nf = 12;
    const unsigned nt = pHOT::ntile+2;
    if (tt.size() != nf) tt = vector<double>(nf);
    vector<double> in(nf), IN(nf);
    vector< vector<double> > out(nt), OUT(nt);
    for (unsigned i=0; i<nt; i++) {
      out[i] = vector<double>(nf);
      if (mlevel==0) OUT[i] = vector<double>(nf);
    }
    
    in[ 0] = partnSoFar;
    in[ 1] = tree1SoFar;
    in[ 2] = tradjSoFar;
    in[ 3] = tcellSoFar;
    in[ 4] = tstepSoFar;
    in[ 5] = llistTime.getTime();
    in[ 6] = collideSoFar;
    in[ 7] = timerSoFar;
    in[ 8] = waitpSoFar;
    in[ 9] = waitcSoFar;
    in[10] = wait1SoFar;
    in[11] = wait2SoFar;
    
    tt[ 0] = partnSoFar;
    tt[ 1] = tree1SoFar;
    tt[ 2] = tradjSoFar;
    tt[ 3] = tcellSoFar;
    tt[ 4] = tstepSoFar;
    tt[ 5] = llistTime.getTime();
    tt[ 6] = collideSoFar;
    tt[ 7] = timerSoFar;
    tt[ 8] = waitpSoFar;
    tt[ 9] = waitcSoFar;
    tt[10] = wait1SoFar;
    tt[11] = wait2SoFar;
    
    if (mlevel==0) {
      for (unsigned k=0; k<nf; k++) {
	IN[k] = tt[k];
	tt[k] = 0.0;
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
    collide->KElossGather();
    
    const int ebins = 1000;
    vector<unsigned> efrt(ebins, 0);
    double minEff, maxEff;
    
    if (use_effort) {
      PartMapItr pitr, pend = c0->Particles().end();
      double minEff1 = 1.0e20, maxEff1 = 0.0;
      for (pitr=c0->Particles().begin(); pitr!= pend; pitr++) {
	minEff1 = min<double>(pitr->second->effort, minEff1);
	maxEff1 = max<double>(pitr->second->effort, maxEff1);
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
	  Lvalue = log(pitr->second->effort);
	  if (maxEff>minEff) {
	    indx = floor((Lvalue - minEff) / (maxEff - minEff) * ebins);
	    if (indx<0)      indx = 0;
	    if (indx>=ebins) indx = ebins-1;
	    efrt1[indx]++;
	  } else efrt1[0]++;
	}
	
	MPI_Reduce(&efrt1[0], &efrt[0], ebins, MPI_UNSIGNED, MPI_SUM, 0, 
		   MPI_COMM_WORLD);
      }
    }
    
    (*barrier)("TreeDSMC: BEFORE Collide::printCollGather",  __FILE__, __LINE__);
    collide->printCollGather();
    (*barrier)("TreeDSMC: AFTER Collide::printCollGather",  __FILE__, __LINE__);
    
    (*barrier)("TreeDSMC: BEFORE Collide::auxGather",  __FILE__, __LINE__);
    collide->auxGather();
    (*barrier)("TreeDSMC: AFTER Collide::auxGather",  __FILE__, __LINE__);
    
    if (use_key>=0) {
      (*barrier)("TreeDSMC: BEFORE species map update", __FILE__, __LINE__);
      makeSpeciesMap();
      (*barrier)("TreeDSMC: AFTER species map update",  __FILE__, __LINE__);
    } 
    
    // Compute fractions in trace counters and temperature for direct and trace
    //
    (*barrier)("TreeDSMC: BEFORE Collide::gatherSpecies",  __FILE__, __LINE__);
    collide->gatherSpecies();
    (*barrier)("TreeDSMC: AFTER Collide::gatherSpecies",   __FILE__, __LINE__);
    
    // Get NTC statistics
    //
    (*barrier)("TreeDSMC: BEFORE Collide::NTCgather",  __FILE__, __LINE__);
    collide->NTCgather();
    (*barrier)("TreeDSMC: AFTER Collide::NTCgather",   __FILE__, __LINE__);
    
    // Get level statistics from tree
    //
    (*barrier)("TreeDSMC: BEFORE pHOT::gatherCellLevelList",  __FILE__, __LINE__);
    c0->Tree()->gatherCellLevelList();
    (*barrier)("TreeDSMC: AFTER pHOT::gatherCellLevelList",   __FILE__, __LINE__);
    
    if (myid==0) {
      
      unsigned sell_total = collide->select();
      unsigned coll_total = collide->total();
      unsigned coll_error = collide->errors();
      unsigned epsm_total = collide->EPSMtotal();
      unsigned epsm_cells = collide->EPSMcells();
      
      collide->printSpecies(spec, TempTot);
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
      
      collide->NTCstats(mout);
      
      collide->colldeTime(mout);
      
      if (epsm>0) mout << setw(6) << " " << setw(20) << epsm_total 
		       << "EPSM particles ("
		       << 100.0*epsm_total/c0->CurTotal() << "%)" 
		       << scientific << endl;
      mout << setw(6) << " " << setw(20) << mlevel     << "multi level" << endl
	   << setw(6) << " " << setw(20) << medianNumb << "number/cell" << endl
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
      
      collide->auxPrint(mout);
      
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
      vector<float>    keynewc, keyoldc, treebar, diagdbg, keysort;
      vector<float>    keygenr;
      vector<unsigned> numbods;
      
      c0->Tree()->Timing(keymake, xchange, convert, overlap, prepare,
			 cupdate, scatter, repartn, tadjust, keycall,
			 keycomp, keybods, keysort, keygenr, waiton0,
			 waiton1, waiton2, keynewc, keyoldc, treebar,
			 diagdbg, numbods);
      
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
      outHelper1<float>(mout, "keysort", keysort);
      outHelper1<float>(mout, "keygenr", keygenr);
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

      mout << endl;
      
      collide->tsdiag(mout);
      collide->voldiag(mout);
      
      //
      // Debugging usage
      //
      collide->EPSMtiming(mout);
      collide->CPUHog (mout);
      collide->KEloss (mout);
      
      //
      // Collision effort
      //
      if (use_effort) {
	const int nlst = 15;
	const int lst[] = {0  ,   4,   9,  19,  49,  99, 199, 
			   499, 799, 899, 949, 979, 989, 994, 998};
	
	// Cumulate
	//
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
      
      
      //
      // Maximum levels for each cell
      //
      c0->Tree()->printCellLevelList(mout, "Cell time step levels");
    }
    
    //
    // Reset the collide counters (CollideIon)
    //
    (*barrier)("TreeDSMC: BEFORE collide counters", __FILE__, __LINE__);
    collide->resetColls();
    (*barrier)("TreeDSMC: AFTER collide counters",  __FILE__, __LINE__);
    
    //
    // Clean NTC catch
    //
    if (ageout>0) {
      (*barrier)("TreeDSMC: BEFORE NTC expire", __FILE__, __LINE__);
      collide->ageout();
      (*barrier)("TreeDSMC: AFTER NTC expire",  __FILE__, __LINE__);
    }
    
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
    
    //
    // Zero the elapsed time counters
    //
    partnSoFar   = 0.0;
    tree1SoFar   = 0.0;
    tradjSoFar   = 0.0;
    tcellSoFar   = 0.0;
    tstepSoFar   = 0.0;
    collideSoFar = 0.0;
    timerSoFar   = 0.0;
    waitpSoFar   = 0.0;
    waitcSoFar   = 0.0;
    wait1SoFar   = 0.0;
    wait2SoFar   = 0.0;

#ifdef USE_GPTL
    GPTLstop("TreeDSMC::collide_diag");
#endif
    
    (*barrier)("TreeDSMC: after collision diags", __FILE__, __LINE__);

  } // END: diagstep
  
  if (cuda_prof) tPtr2.reset();

  
#ifdef USE_GPTL
  GPTLstop("TreeDSMC::determine_acceleration_and_potential");
#endif
  
  // pHOT frontier dump
  //
  if (dumpHOT) {
    std::ofstream out;
    if (myid==0) {
      std::ostringstream sout;
      sout << outdir << runtag << "." << cntHOT++ << ".frontier";
      out.open(sout.str());
    }
    c0->Tree()->dumpFrontier(out);
  }

  // Debugging disk
  //
  // triggered_cell_body_dump(0.01, 0.002);
  // TempHisto();
  
  (*barrier)("TreeDSMC: end of accel routine", __FILE__, __LINE__);
  
  if (myid==0 and mlevel==0 and tnow>tM0) {
    tM0 = tnow;
    fullstep.Stop();

    ostringstream sout;
    sout << outdir << runtag << ".DSMC_log";
    ofstream mout(sout.str().c_str(), ios::app);
    
    mout << endl
	 << "------------------------------------" << endl
	 << "Timing (microsec) at T="              << tnow << endl
	 << "------------------------------------" << endl
	 << "----Running average = "               << fullstep.getTavg()*1.0e+06 << endl
	 << "----Last full step  = "               << fullstep.getTime()*1.0e+06 << endl
	 << "----Last time intvl = "               << fullstep.getLast()*1.0e+06 << endl
	 << "------------------------------------" << endl;
      
  } else {
    fullstep.Stop(false);
  }
}

void TreeDSMC::triggered_cell_body_dump(double time, double radius)
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


void TreeDSMC::assignTempDensVol()
{
  const double f_H = 0.76;
  // This should be generalized and
  // computed dynamically rather than
  // hardwired
  double mm = f_H*mp + (1.0-f_H)*4.0*mp;
  double KEtot, KEdsp, T;
  double Tfac;
  if (ctype == "LTE") Tfac = 2.0*TreeDSMC::Eunit/3.0 * mm/TreeDSMC::Munit/boltz;
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
	
        speciesKey sKey = Particle::defaultKey;
	if (use_key>=0) {
	  KeyConvert kc(cell->Body(j)->iattrib[use_key]);
	  sKey = kc.getKey();
	}
	
	cell->sample->KE(sKey, KEtot, KEdsp);
	double mi = amu;
	// For non-trace or multiple species
	// computations
	if (sKey != Particle::defaultKey) mi = mp*atomic_weights[sKey.first];
	
	Tfac = 2.0*TreeDSMC::Eunit/3.0 * mi/TreeDSMC::Munit/boltz;
	T    = KEdsp* Tfac;
	
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
	nbod  += csz;
	minT   = min<double>(T, minT);
	maxT   = max<double>(T, maxT);
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

void TreeDSMC::TempHisto()
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
    
    std::cout << "----------------" << std::endl
	      << "Temperature dist" << std::endl
	      << "Time=" << tnow    << std::endl
	      << "----------------" << std::endl
	      << std::setw(10) << "<1000" 
	      << std::setw(10) << std::setprecision(2) << td0[0]/totalM0
	      << std::setw(10) << std::setprecision(2) << vd0[0] << std::endl;
    for (int i=0; i<numT; i++) 
      std::cout << std::setw(10) << std::setprecision(2) 
		<< exp10(TlogMin + (TlogMax-TlogMin)/numT*(0.5+i))
		<< std::setw(10) << std::setprecision(2) << td0[i+1]/totalM0
		<< std::setw(10) << std::setprecision(2) << vd0[i+1] << std::endl;
    std::cout << std::setw(10) << ">1e8" << std::setw(10) 
	      << std::setprecision(2) << td0[numT+1]/totalM0
	      << std::setw(10) << std::setprecision(2) << vd0[numT+1] << std::endl;
  }
}

// Make species map
//
void TreeDSMC::makeSpeciesMap()
{
  // Clean the local and global maps
  //
  spec .clear();
  spec1.clear();
  std::vector<double> T1(2, 0), T0(2, 0);
  
  // Particle loop
  //
  PartMapItr p    = c0->Particles().begin();
  PartMapItr pend = c0->Particles().end();
  for (; p!=pend; p++) {
    //
    // Species accumulation
    //
    KeyConvert kc  (p->second->iattrib[use_key]);
    speciesKey indx(kc.getKey());
    unsigned short Z = indx.first;
    for (unsigned short i=1; i<=Z+1; i++) {
      speciesKey indxi(Z, i);
      if (spec1.find(indxi) == spec1.end()) spec1[indxi] = 0;
    }
    if (spec1.find(indx) == spec1.end()) spec1[indx] = 1;
    else                                 spec1[indx]++;
    //
    // Temperature accumulation
    //
    if (use_temp>=0) {
      T1[0] += p->second->mass;
      T1[1] += p->second->mass * p->second->dattrib[use_temp];
    }
  }
  
  
  // Send the local map to other nodes
  //
  int sizm;
  speciesKey indx;
  unsigned long cnts;
  std::map<speciesKey, unsigned long>::iterator it, it2;
  
  
  spec = spec1;			// Copy local map to global
  for (int i=0; i<numprocs; i++) {
    if (i == myid) {
      sizm = spec1.size();
      // Local map size
      MPI_Bcast(&sizm, 1, MPI_INT, i, MPI_COMM_WORLD);
      // Send local map
      for (it=spec1.begin(); it != spec1.end(); it++) {
	indx = it->first;
	cnts = it->second;
	MPI_Bcast(&indx.first, 1, MPI_UNSIGNED, i, MPI_COMM_WORLD);
	MPI_Bcast(&indx.second, 1, MPI_UNSIGNED, i, MPI_COMM_WORLD);
	MPI_Bcast(&cnts, 1, MPI_UNSIGNED_LONG, i, MPI_COMM_WORLD);
      }
    } else {
      // Receive from Node i
      MPI_Bcast(&sizm, 1, MPI_INT, i, MPI_COMM_WORLD);
      for (int j=0; j<sizm; j++) {
	MPI_Bcast(&indx.first, 1, MPI_UNSIGNED, i, MPI_COMM_WORLD);
	MPI_Bcast(&indx.second, 1, MPI_UNSIGNED, i, MPI_COMM_WORLD);
	MPI_Bcast(&cnts, 1, MPI_UNSIGNED_LONG, i, MPI_COMM_WORLD);
	// Update global map
	if (spec.find(indx) == spec.end()) spec[indx]  = cnts;
	else                               spec[indx] += cnts;
      }
    }
  }
  
  MPI_Allreduce(&T1[0], &T0[0], 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  if (T0[0] > 0.0) TempTot = T0[1]/T0[0];
  else TempTot = 0.0;
}

void TreeDSMC::finish()
{ 
  collide->finish(); 
}

