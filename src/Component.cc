#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>
#include <memory>
#include <map>

#include <Component.H>
#include <Bessel.H>
#include <CBrock.H>
#include <CBrockDisk.H>
#include <Hernquist.H>
#include <Sphere.H>
#include <EJcom.H>
#include <Cylinder.H>
#include <Cube.H>
#include <Slab.H>
#include <SlabSL.H>
#include <Direct.H>
#include <Shells.H>
#include <NoForce.H>
#include <Orient.H>
#include <pHOT.H>

#include "expand.H"

// For sort algorithm below
bool less_loadb(const loadb_datum& one, const loadb_datum& two)
{
  return (one.top < two.top);
}

// Constructor
Component::Component(YAML::Node& CONF)
{
  // Make a copy
  conf = CONF;

  try {
    name = conf["name"].as<std::string>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << __FILE__ << ": " << __LINE__ << std::endl
			   << "Error parsing component 'name': "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;

    MPI_Finalize();
    exit(-1);
  }

  try {
    cconf = conf["parameters"];
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing 'parameters' for Component <"
			   << name << ">: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf
			   << std::string(60, '-') << std::endl;

    MPI_Finalize();
    exit(-2);
  }
  
  pfile = conf["bodyfile"].as<std::string>();

  YAML::Node force;
  try {
    force = conf["force"];
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing 'force' for Component <"
			   << name << ">: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;

    MPI_Finalize();
    exit(-3);
  }

  id = force["id"].as<std::string>();

  try {
    fconf = force["parameters"];
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing force 'parameters' for Component <"
			   << name << ">: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << force                << std::endl
			   << std::string(60, '-') << std::endl;

    MPI_Finalize();
    exit(-4);
  }

  EJ          = 0;
  nEJkeep     = 100;
  nEJwant     = 500;
  EJkinE      = true;
  EJext       = false;
  EJdiag      = false;
  EJdryrun    = false;
  EJx0        = 0.0;
  EJy0        = 0.0;
  EJz0        = 0.0;
  EJu0        = 0.0;
  EJv0        = 0.0;
  EJw0        = 0.0;
  EJdT        = 0.0;
  EJlinear    = false;
  EJdamp      = 1.0;

  binary      = false;

  adiabatic   = false;
  ton         = -1.0e20;
  toff        =  1.0e20;
  twid        = 0.1;

  rtrunc      = 1.0e20;
  rcom        = 1.0e20;
  consp       = false;
  tidal       = -1;

  com_system  = false;
  com_log     = false;

#if HAVE_LIBCUDA==1
  bunchSize   = 100000;
#endif
  
  timers      = false;
				// Null out pointers
  orient      = 0;

  com         = 0;
  cov         = 0;
  coa         = 0;
  center      = 0;
  angmom      = 0;
  ps          = 0;

  com0        = 0;
  cov0        = 0;
  acc0        = 0;

  indexing    = true;
  aindex      = false;
  umagic      = true;

  nlevel      = -1;
  keyPos      = -1;
  top_seq     = 0;

  pBufSiz     = 100000;		// Default number particles in MPI-IO buffer
  blocking    = false;		// Default for MPI_File_write blocking
  buffered    = true;		// Use buffered writes for POSIX binary

  set_default_values();

  mdt_ctr = vector< vector<unsigned> > (multistep+1);
  for (unsigned n=0; n<=multistep; n++) mdt_ctr[n] = vector<unsigned>(mdtDim, 0);

  angmom_lev  = vector<double>(3*(multistep+1), 0);
  com_lev     = vector<double>(3*(multistep+1), 0);
  cov_lev     = vector<double>(3*(multistep+1), 0);
  coa_lev     = vector<double>(3*(multistep+1), 0);
  com_mas     = vector<double>(multistep+1, 0);

  tree = 0;

  pbuf.resize(PFbufsz);

  // Enter unset defaults in YAML conf
  //
  if (CONF["parameters"]) CONF["parameters"] = cconf;

  configure();

  initialize_cuda();

  read_bodies_and_distribute_ascii();

  reset_level_lists();

  modified = 0;
}

void Component::set_default_values()
{
  if (!cconf["EJ"])              cconf["EJ"]          = EJ;
  if (!cconf["nEJkeep"])         cconf["nEJkeep"]     = nEJkeep;
  if (!cconf["nEJwant"])         cconf["nEJwant"]     = nEJwant;
  if (!cconf["EJkinE"])          cconf["EJkinE"]      = EJkinE;
  if (!cconf["EJext"])           cconf["EJext"]       = EJext;
  if (!cconf["EJdiag"])          cconf["EJdiag"]      = EJdiag;
  if (!cconf["EJdryrun"])        cconf["EJdryrun"]    = EJdryrun;
  if (!cconf["EJx0"])            cconf["EJx0"]        = EJx0;
  if (!cconf["EJy0"])            cconf["EJy0"]        = EJy0;
  if (!cconf["EJz0"])            cconf["EJz0"]        = EJz0;
  if (!cconf["EJu0"])            cconf["EJu0"]        = EJu0;
  if (!cconf["EJv0"])            cconf["EJv0"]        = EJv0;
  if (!cconf["EJw0"])            cconf["EJw0"]        = EJw0;
  if (!cconf["EJdT"])            cconf["EJdT"]        = EJdT;
  if (!cconf["EJlinear"])        cconf["EJlinear"]    = EJlinear;
  if (!cconf["EJdamp"])          cconf["EJdamp"]      = EJdamp;
  if (!cconf["binary"])          cconf["binary"]      = binary;
  if (!cconf["adiabatic"])       cconf["adiabatic"]   = adiabatic;
  if (!cconf["ton"])             cconf["ton"]         = ton;
  if (!cconf["toff"])            cconf["toff"]        = toff;
  if (!cconf["twid"])            cconf["twid"]        = twid;
  if (!cconf["rtrunc"])          cconf["rtrunc"]      = rtrunc;
  if (!cconf["rcom"])            cconf["rcom"]        = rcom;
  if (!cconf["consp"])           cconf["consp"]       = consp;
  if (!cconf["tidal"])           cconf["tidal"]       = tidal;
  if (!cconf["comlog"])          cconf["comlog"]      = com_log;
#if HAVE_LIBCUDA==1
  if (!cconf["bunch"])           cconf["bunch"]       = bunchSize;
#endif
  if (!cconf["timers"])          cconf["timers"]      = timers;
  if (!cconf["com_system"])      cconf["com_system"]  = com_system;
  if (!cconf["com"])             cconf["com"]         = com_system;
  if (!cconf["indexing"])        cconf["indexing"]    = indexing;
  if (!cconf["aindex"])          cconf["aindex"]      = aindex;
  if (!cconf["umagic"])          cconf["umagic"]      = umagic;
  if (!cconf["nlevel"])          cconf["nlevel"]      = nlevel;
  if (!cconf["keyPos"])          cconf["keyPos"]      = keyPos;
  if (!cconf["pBufSiz"])         cconf["pBufSiz"]     = pBufSiz;
  if (!cconf["blocking"])        cconf["blocking"]    = blocking;
  if (!cconf["buffered"])        cconf["buffered"]    = buffered;
}


void Component::HOTcreate(std::set<speciesKey> spec_list)
{
  delete tree;
  tree = new pHOT(this, spec_list);
}


void Component::HOTdelete()
{
  delete tree;
}


class thrd_pass_reset
{
public:
  int id;
  Component *c;
  vector< vector<int> > newlist;
};

static vector<thrd_pass_reset> td;
static pthread_t* t  = 0;

void * reset_level_lists_thrd(void *ptr)
{
  // Thread ID
  int id = static_cast<thrd_pass_reset*>(ptr)->id;

  // Component
  Component *c = static_cast<thrd_pass_reset*>(ptr)->c;


  int nbodies = c->Number();
  int nbeg = nbodies*(id  )/nthrds;
  int nend = nbodies*(id+1)/nthrds;
  
  PartMapItr it = c->Particles().begin();

  std::advance(it, nbeg);

  vector< vector<int> > *v = &static_cast<thrd_pass_reset*>(ptr)->newlist;

  for (int n=nbeg; n<nend; n++, it++) {
    (*v)[it->second->level].push_back(it->first);
  }
  
  return (NULL);
}

void Component::reset_level_lists()
{
  if (td.size()==0) {
    td = vector<thrd_pass_reset>(nthrds);

    t  = new pthread_t [nthrds];

    if (!t) {
      std::ostringstream sout;
      sout << "Process " << myid
	   << ": reset_level_lists: error allocating memory for thread";
      throw GenericError(sout.str(), __FILE__, __LINE__);
    }
  }
  
  if (nthrds==1) {

    td[0].id = 0;
    td[0].c  = this;
    td[0].newlist  = vector< vector<int> >(multistep+1);

    reset_level_lists_thrd(&td[0]);

  } else {

    //
    // Make the <nthrds> threads
    //
    int errcode;
    void *retval;
  
    for (int i=0; i<nthrds; i++) {
    
      td[i].id = i;
      td[i].c  = this;
      td[i].newlist  = vector< vector<int> >(multistep+1);
      
      errcode =  pthread_create(&t[i], 0, reset_level_lists_thrd, &td[i]);

      if (errcode) {
	std::ostringstream sout;
	sout << "Process " << myid
	     << " reset_level_lists: cannot make thread " << i
	     << ", errcode=" << errcode;;
	throw GenericError(sout.str(), __FILE__, __LINE__);
      }
#ifdef DEBUG
      else {
	cout << "Process " << myid << ": thread <" << i << "> created\n";
      }
#endif
    }
    
    //
    // Collapse the threads
    //
    for (int i=0; i<nthrds; i++) {
      if ((errcode=pthread_join(t[i], &retval))) {
	std::ostringstream sout;
	sout << "Process " << myid
	     << " reset_level_lists: thread join " << i
	     << " failed, errcode=" << errcode;
	throw GenericError(sout.str(), __FILE__, __LINE__);
      }
#ifdef DEBUG    
      cout << "Process " << myid << ": multistep thread <" << i << "> thread exited\n";
#endif
    }
  }
				// Particle list per level.
				// Begin with empty lists . . .
  levlist = vector< vector<int> > (multistep+1);
  for (int i=0; i<nthrds; i++) {
    for (unsigned n=0; n<=multistep; n++) {
      levlist[n].insert(levlist[n].end(),
			td[i].newlist[n].begin(), 
			td[i].newlist[n].end());
    }
  }
  
  if (VERBOSE>10) {
				// Level creation check
    for (int n=0; n<numprocs; n++) {
      if (n==myid) {
	if (myid==0) 
	  cout << endl
	       << "----------------------------------------------" << endl
	       << "Level creation in Component <" << name << ">:" << endl 
	       << "----------------------------------------------" << endl
	       << setw(4) << left << "ID" << setw(4) << "lev"
	       << setw(12) << "first" << setw(12) << "last" 
	       << setw(12) << "count" << endl;
	for (unsigned j=0; j<=multistep; j++) {
	  cout << left << setw(4) << myid << setw(4) << j;
	  if (levlist[j].size())
	    cout << left
		 << setw(12) << levlist[j].front()
		 << setw(12) << levlist[j].back() 
		 << setw(12) << levlist[j].size()
		 << endl;
	  else
	    cout << left
		 << setw(12) << (int)(-1)
		 << setw(12) << (int)(-1) 
		 << setw(12) << levlist[j].size()
		 << endl;
	}
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << endl;
  }

}

void Component::print_level_lists(double T)
{
				// Print out level info

  if (nlevel>0 && (this_step % nlevel == 0)) {

#if HAVE_LIBCUDA==1
    if (use_cuda) {		// Call the CUDA version
      print_level_lists_cuda(T);
      return;
    }
#endif

    std::vector< std::vector<unsigned> > cntr(multistep+1);
    for (unsigned n=0; n<=multistep; n++) {
      cntr[n] = std::vector<unsigned>(mdtDim, 0);
      MPI_Reduce(&mdt_ctr[n][0], &cntr[n][0], mdtDim, MPI_UNSIGNED,
		 MPI_SUM, 0, MPI_COMM_WORLD);
      
      for (int k=0; k<mdtDim; k++) mdt_ctr[n][k] = 0;
    }
    
    if (myid==0) {

      unsigned tot=0;
      for (unsigned n=0; n<=multistep; n++) tot += cntr[n][mdtDim-1];

      if (!tot && myid==0)
	std::cout << "print_level_lists [" << name 
		  << ", T=" << tnow << "]: tot=" << tot << std::endl;
      
      if (tot) {

	std::ostringstream ofil;
	ofil << outdir << runtag << ".levels";
	std::ofstream out(ofil.str().c_str(), ios::app);

	unsigned curn, dtcnt, sum=0;
	out << setw(90) << std::setfill('-') << '-' << std::endl;
	std::ostringstream sout;
	sout << "--- Component <" << name 
	     << ", " << id  << ">, T=" << T;
	out << std::setw(90) << std::left << sout.str().c_str() << std::endl;
	out << std::setw(90) << '-' << std::endl << std::setfill(' ');
	out << std::setw(3)  << "L" 
	    << std::setw(10) << "Number" 
	    << std::setw(10) << "dN/dL" 
	    << std::setw(10) << "N(<=L)";
	if (DTold) {
	  out << std::setw(10) << "f(r/v)"
	      << std::setw(10) << "f(s/v)"
	      << std::setw(10) << "f(v/a)"
	      << std::setw(10) << "f(r/a)";
	  dtcnt = 5;
	} else {
	  out << std::setw(10) << "f(q/v)"
	      << std::setw(10) << "f(v/a)"
	      << std::setw(10) << "f(s/v)"
	      << std::setw(10) << "f(r/v)" 
	      << std::setw(10) << "f(r/a)";
	  dtcnt = 6;
	}
	out << std::setw(10) << "f(int)" << std::endl;
	out << std::setw(90) << std::setfill('-') << '-' << std::endl
	    << std::setfill(' ');
	for (unsigned n=0; n<=multistep; n++) {
	  curn = cntr[n][mdtDim-1];
	  sum += curn;
	  out << std::setw(3)  << n 
	      << std::setw(10) << curn << setprecision(3) << std::fixed
	      << std::setw(10) << static_cast<double>(curn)/tot
	      << std::setw(10) << static_cast<double>(sum) /tot;
	  for (unsigned k=0; k<dtcnt; k++) {
				// If there are counts at this level:
	    if (curn) out << std::setw(10)
			  << static_cast<double>(cntr[n][k])/curn;
				// No counts at this level:
	    else      out << std::setw(10) << "*";
	  }
	  out << std::endl;
	}
	out << std::endl << std::setw(3) << "T" << std::setw(10) << tot
	    << std::endl << std::endl << std::right;
      }
    }
  }
}

Component::Component(YAML::Node& CONF, istream *in, bool SPL) : conf(CONF)
{
  // Make a copy
  conf = CONF;

  try {
    name = conf["name"].as<std::string>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << __FILE__ << ": " << __LINE__ << std::endl
			   << "Error parsing component 'name': "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;

    MPI_Finalize();
    exit(-5);
  }

  try {
    cconf = conf["parameters"];
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing 'parameters' for Component <"
			   << name << ">: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf
			   << std::string(60, '-') << std::endl;

    MPI_Finalize();
    exit(-6);
  }
  
  pfile = conf["bodyfile"].as<std::string>();

  YAML::Node cforce;
  try {
    cforce = conf["force"];
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing 'force' for Component <"
			   << name << ">: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;

    MPI_Finalize();
    exit(-7);
  }

  id = cforce["id"].as<std::string>();

  try {
    fconf = cforce["parameters"];
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing force 'parameters' for Component <"
			   << name << ">: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << cforce                << std::endl
			   << std::string(60, '-') << std::endl;

    MPI_Finalize();
    exit(-8);
  }

  // Defaults
  //
  EJ          = 0;
  nEJkeep     = 100;
  nEJwant     = 500;
  EJkinE      = true;
  EJext       = false;
  EJdiag      = false;
  EJdryrun    = false;
  EJx0        = 0.0;
  EJy0        = 0.0;
  EJz0        = 0.0;
  EJu0        = 0.0;
  EJv0        = 0.0;
  EJw0        = 0.0;
  EJdT        = 0.0;
  EJlinear    = false;
  EJdamp      = 1.0;

  binary      = true;

  adiabatic   = false;
  ton         = -1.0e20;
  toff        =  1.0e20;
  twid        = 0.1;

  rtrunc      = 1.0e20;
  rcom        = 1.0e20;
  consp       = false;
  tidal       = -1;

  com_system  = false;
  com_log     = false;
  com_restart = 0;

#if HAVE_LIBCUDA==1
  bunchSize   = 100000;
#endif

  timers      = false;

  force       = 0;		// Null out pointers
  orient      = 0;

  com         = 0;
  cov         = 0;
  coa         = 0;
  center      = 0;
  angmom      = 0;
  ps          = 0;

  com0        = 0;
  cov0        = 0;
  acc0        = 0;

  indexing    = true;
  aindex      = false;
  umagic      = true;

  keyPos      = -1;
  nlevel      = -1;
  top_seq     = 0;

  pBufSiz     = 100000;
  blocking    = false;

  configure();

  initialize_cuda();

  if (SPL) read_bodies_and_distribute_binary_spl(in);
  else     read_bodies_and_distribute_binary_out(in);

  mdt_ctr = vector< vector<unsigned> > (multistep+1);
  for (unsigned n=0; n<=multistep; n++) mdt_ctr[n] = vector<unsigned>(mdtDim, 0);

  angmom_lev  = vector<double>(3*(multistep+1), 0);
  com_lev     = vector<double>(3*(multistep+1), 0);
  cov_lev     = vector<double>(3*(multistep+1), 0);
  coa_lev     = vector<double>(3*(multistep+1), 0);
  com_mas     = vector<double>(multistep+1, 0);

  reset_level_lists();

  tree = 0;

  pbuf.resize(PFbufsz);
}

void Component::configure(void)
{
  // Load parameters from YAML configuration node
  try {
    if (cconf["com"     ])  com_system = cconf["com"     ].as<bool>();
    if (cconf["comlog"  ])     com_log = cconf["comlog"  ].as<bool>();
    if (cconf["timers"  ])      timers = cconf["timers"  ].as<bool>();
  
    if (cconf["tidal"]) {
      tidal = cconf["tidal"].as<int>();
      consp = true;
    }

    if (cconf["EJ"      ])         EJ  = cconf["EJ"].as<int>();
    if (cconf["eEJ0"    ] and myid==0)
      std::cout << "Component: eEJ0 is no longer used, Ecurr is computed from the bodies using the expansion directly" << std::endl;
    if (cconf["nEJkeep" ])    nEJkeep  = cconf["nEJkeep" ].as<int>();
    if (cconf["nEJwant" ])    nEJwant  = cconf["nEJwant" ].as<int>();
    if (cconf["EJx0"    ])       EJx0  = cconf["EJx0"    ].as<double>();
    if (cconf["EJy0"    ])       EJy0  = cconf["EJy0"    ].as<double>();
    if (cconf["EJz0"    ])       EJz0  = cconf["EJz0"    ].as<double>();
    if (cconf["EJu0"    ])       EJu0  = cconf["EJu0"    ].as<double>();
    if (cconf["EJv0"    ])       EJv0  = cconf["EJv0"    ].as<double>();
    if (cconf["EJw0"    ])       EJw0  = cconf["EJw0"    ].as<double>();
    if (cconf["EJdT"    ])       EJdT  = cconf["EJdT"    ].as<double>();
    if (cconf["EJkinE"  ])     EJkinE  = cconf["EJkinE"  ].as<bool>();
    if (cconf["EJext"   ])      EJext  = cconf["EJext"   ].as<bool>();
    if (cconf["EJdiag"  ])     EJdiag  = cconf["EJdiag"  ].as<bool>();
    if (cconf["EJdryrun"])   EJdryrun  = cconf["EJdryrun"].as<bool>();
    if (cconf["EJlinear"])   EJlinear  = cconf["EJlinear"].as<bool>();
    if (cconf["EJdamp"  ])     EJdamp  = cconf["EJdamp"  ].as<double>();
    if (cconf["rmax"    ])       rmax  = cconf["rmax"    ].as<double>();
    if (cconf["rtrunc"  ])     rtrunc  = cconf["rtrunc"  ].as<double>();
    if (cconf["rcom"    ])       rcom  = cconf["rcom"    ].as<double>();
    if (cconf["magic"   ])     umagic  = cconf["magic"   ].as<bool>();
    if (cconf["indexing"])   indexing  = cconf["indexing"].as<bool>();
    if (cconf["aindex"  ])     aindex  = cconf["aindex"  ].as<bool>();
    if (cconf["nlevel"  ])     nlevel  = cconf["nlevel"  ].as<int>();
    if (cconf["keypos"  ])     keyPos  = cconf["keypos"  ].as<int>();
    if (cconf["pbufsiz" ])    pBufSiz  = cconf["pbufsiz" ].as<int>();
    if (cconf["blocking"])   blocking  = cconf["blocking"].as<bool>();
#if HAVE_LIBCUDA==1
    if (cconf["bunch"   ])  bunchSize  = cconf["bunch"   ].as<int>();
#endif
    
    if (cconf["ton"]) {
      ton = cconf["ton"].as<double>();
      adiabatic = true;
    }

    if (cconf["toff"]) {
      toff = cconf["toff"].as<double>();
      adiabatic = true;
    }

    if (cconf["twid"]) {
      twid = cconf["twid"].as<double>();
      adiabatic = true;
    }
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters for Component <"
			   << name << ">: "
			   << error.what()         << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << cconf                << std::endl
			   << std::string(60, '-') << std::endl;

    MPI_Finalize();
    exit(-9);
  }


  // Instantiate the force ("reflection" by hand)
  //
  if ( !id.compare("bessel") ) {
    force = new Bessel(this, fconf);
  }
  else if ( !id.compare("c_brock") ) {
    force = new CBrock(this, fconf);
  }
  else if ( !id.compare("c_brock_disk") ) {
    force = new CBrockDisk(this, fconf);
  }
  else if ( !id.compare("hernq") ) {
    force = new Hernquist(this, fconf);
  }
  else if ( !id.compare("sphereSL") ) {
    force = new Sphere(this, fconf);
  }
  else if ( !id.compare("EJcom") ) {
    force = new EJcom(this, fconf);
  }
  else if ( !id.compare("cube") ) {
    force = new Cube(this, fconf);
  }
  else if ( !id.compare("slab") ) {
    force = new Slab(this, fconf);
  }
  else if ( !id.compare("slabSL") ) {
    force = new SlabSL(this, fconf);
  }
  else if ( !id.compare("cylinder") ) {
    force = new Cylinder(this, fconf);
  }
  else if ( !id.compare("direct") ) {
    force = new Direct(this, fconf);
  }
  else if ( !id.compare("shells") ) {
    force = new Shells(this, fconf);
  }
  else if ( !id.compare("noforce") ) {
    force = new NoForce(this, fconf);
  }
  else {
    string msg("I don't know about the force: ");
    msg += id;
    bomb(msg);
  }

  dim = force->dof;
}


void Component::initialize(void)
{
  com    = new double [3];
  center = new double [3];
  cov    = new double [3];
  coa    = new double [3];
  angmom = new double [3];
  ps     = new double [6];

				// For COM system
  com0   = new double[3];
  cov0   = new double[3];
  acc0   = new double[3];

  for (int k=0; k<3; k++) {
    com[k]  = center[k] = cov[k]  = coa[k]    = 0.0;
    com0[k] = cov0[k]   = acc0[k] = angmom[k] = 0.0;
  }  

  if (com_system) {

    if (consp) {
      comE_lev = vector<double>(3*(multistep+1), 0);
      covE_lev = vector<double>(3*(multistep+1), 0);
      comE_mas = vector<double>(multistep+1, 0);
    }

    initialize_com_system();

    if (myid==0) {
      cout << "---- Component <" <<  name 
	   << ">: center of mass system is *ON*, rtrunc=" << rtrunc;
      if (consp) cout << ", conserving com momentum [iattr #=" << tidal << "]";
      cout << ", computed COM system:";
      cout << endl << "\t\t(x, y, z)=("
	   << setw(15) << com0[0] << ", "
	   << setw(15) << com0[1] << ", "
	   << setw(15) << com0[2] << ") "
	   << endl << "\t\t"
	   << "(u, v, w)=("
	   << setw(15) << cov0[0] << ", "
	   << setw(15) << cov0[1] << ", "
	   << setw(15) << cov0[2] << ") "
	   << endl;
      
      if (com_log) {

	comfile = outdir + name + ".comlog." + runtag;
	bool newfile = true;

	if (restart) {

	  // Open old file for reading
	  ifstream in(comfile.c_str());

	  if (in) {
	    
	    // Backup up old file
	  
	    string backupfile = comfile + ".bak";
	    if (rename(comfile.c_str(), backupfile.c_str())) {
	      ostringstream message;
	      message << "Component: error making backup file <" 
		      << backupfile << ">\n";
	      bomb(message.str().c_str());
	    }

	    // Open new output stream for writing
	  
	    ofstream out(comfile.c_str());
	    if (!out) {
	      ostringstream message;
	      message << "Component: error opening new log file <" 
		      << comfile << "> for writing\n";
	      bomb(message.str().c_str());
	    }
	  
	    const int cbufsiz = 16384;
	    char *cbuffer = new char [cbufsiz];
	    double ttim, ttim0;
	    int tarrow = 1;
	    bool first_data = true;

	    while (in) {
	      in.getline(cbuffer, cbufsiz);
	      if (!in) break;
	      string line(cbuffer);
	      
	      if (line.find_first_of("#") != string::npos) {

		// Header/comment lines

		out << cbuffer << "\n";
		
	      } else {
		
		// Data lines
	      
		StringTok<string> toks(line);
		ttim  = atof(toks(" ").c_str());

		if (first_data) {
		  istringstream istr(line);
		  istr >> ttim0;
		  for (int k=0; k<3; k++) istr >> com0[k];
		  for (int k=0; k<3; k++) istr >> cov0[k];
		  for (int k=0; k<3; k++) istr >> acc0[k];
		  for (int k=0; k<3; k++) istr >> center[k];
		  first_data = false;
		}

		// Compute the direction of time
		
		if (ttim != ttim0) tarrow = ttim - ttim0 ? 1 : -1;

		if ( (tnow - ttim)*tarrow < 1.0e-8 ) {
		  istringstream istr(line);

		  istr >> ttim;

		  for (int k=0; k<3; k++) istr >> com0[k];
		  for (int k=0; k<3; k++) istr >> cov0[k];
		  for (int k=0; k<3; k++) istr >> acc0[k];
		  for (int k=0; k<3; k++) istr >> center[k];
	    
		  cout << "\t\tRead com log for Component <" << name 
		       << "> at T=" << ttim << ", using:";

		  cout << endl << "\t\t(x, y, z)=("
		       << setw(15) << com0[0] << ", "
		       << setw(15) << com0[1] << ", "
		       << setw(15) << com0[2] << ") "
		       << endl << "\t\t"
		       << "(u, v, w)=("
		       << setw(15) << cov0[0] << ", "
		       << setw(15) << cov0[1] << ", "
		       << setw(15) << cov0[2] << ") "
		       << endl;

		  newfile = false;
		  com_restart = 1;

		  break;
		}
		out << cbuffer << "\n";
	      }
	    }

	    delete [] cbuffer;

	    if (newfile) {
	      cout << "Component: time=" << tnow << " not found in <"
		   << comfile << ">, starting new log file\n";
	    }

	  } else {
	    cout << "Component: error opening original log file <" 
		 << comfile << "> for reading, starting new log file\n";
	  }
	}

	if (newfile) {
	  ofstream out(comfile.c_str());
	  if (!out) {
	    ostringstream message;
	    message << "Component: error opening new log file <" 
		    << comfile << "> for writing\n";
	    bomb(message.str().c_str());
	  }
	  
	  out.setf(ios::left);
	  out << setw(15) << "#\n";
	  out << setw(15) << "# Time"
	      << setw(15) << "X"
	      << setw(15) << "Y"
	      << setw(15) << "Z"
	      << setw(15) << "U"
	      << setw(15) << "V"
	      << setw(15) << "W"
	      << setw(15) << "aX"
	      << setw(15) << "aY"
	      << setw(15) << "aZ"
	      << setw(15) << "cX"
	      << setw(15) << "cY"
	      << setw(15) << "cZ"
	      << endl;
	  out << "#\n";
	}
      }
      cout << "\n";		// Close off info line
    }
				// Send com to all processes
    restart_com_system();
  }


  if (EJ) {

    if (EJdiag) cout << "Process " << myid << ": about to create Orient with"
		     << " nkeep="  << nEJkeep
		     << " nwant="  << nEJwant
		     << " EJkinE=" << EJkinE
		     << " EJext="  << EJext;
    
    if (myid==0) {
      if (EJ & Orient::CENTER) {
	cout << "---- Component <" << name << ">: EJ center finding is *ON*";
	cout << " with damping=" << EJdamp;
	if (EJkinE)   cout << ", using particle kinetic energy";
	if (EJext)    cout << ", using external potential";
	if (EJdryrun) cout << ", dryrun";
	cout << endl;
      }
      if (EJ & Orient::AXIS) {
	cout << "---- Component <" << name << ">: AXIS orientation is *ON*";
	if (EJdryrun) cout << ", dryrun";
	cout << endl;
      }
    }
      
    string EJlogfile = outdir + name + ".orient." + runtag;

    unsigned EJctl = 0;
    if (EJdiag)		EJctl |= Orient::DIAG;
    if (EJkinE)		EJctl |= Orient::KE;
    if (EJext)		EJctl |= Orient::EXTERNAL;

    orient = new Orient(nEJkeep, nEJwant, EJ, EJctl, EJlogfile, EJdT, EJdamp);
    
    if (restart && (EJ & Orient::CENTER)) {
      Eigen::VectorXd::Map(&center[0], 3) = orient->currentCenter();
    } else {
      if (EJlinear) orient -> set_linear();
      if (not com_system) {
	orient -> set_center(EJx0, EJy0, EJz0);
	orient -> set_cenvel(EJu0, EJv0, EJw0);
	center[0] = EJx0;
	center[1] = EJy0;
	center[2] = EJz0;
      }
    }

    if (EJdiag) cout << "Process " << myid << ": Orient successful\n";
  }

  if (myid == 0) {		// Center status
    cout << "---- Component <" << name;
    if (restart)
      cout << ">: current center on restart: x, y, z: " 
	   << center[0] << ", " 
	   << center[1] << ", " 
	   << center[2] << std::endl;
    else if (not com_system) {
      cout << ">: user specified initial center: x, y, z: " 
	   << EJx0 << ", " 
	   << EJy0 << ", " 
	   << EJz0 << std::endl;
    } else {
      cout << ">: default center relative to com: x, y, z: " 
	   << 0.0 << ", " 
	   << 0.0 << ", " 
	   << 0.0 << std::endl;
    }

    cout << "---- Component <" << name << ">: ";

    if (nlevel<0)
      cout << "no multistep level reporting";
    else
      cout << "multistep level reporting every " << nlevel << " steps";

    cout << endl << endl;
  }
  
}

void Component::initialize_cuda()
{
#if HAVE_LIBCUDA==1
  cudaDevice = -1;

  if (use_cuda) {

    // Get device count; exit on failure
    //
    cuda_safe_call_mpi(cudaGetDevice(&cudaDevice), __FILE__, __LINE__,
		       myid, "cudaGetDevice failure");
      
    // Check device context
    //
    if (cudaDevice>=0) {

      std::cout << "---- Component <" << name << ">: "
		<< "on CUDA device on Rank [" << myid
		<< "] on [" << processor_name << "]"
		<< std::endl;

      cuda_initialize();
    } else {
      std::cout << "---- Component <" << name << ">: "
		<< "could not find an initialized CUDA device on Rank ["
		<< myid	<< "] on [" << processor_name << "] . . . "
		  << "this will cause a failure" << std::endl;
    }
  }
#endif

}

Component::~Component(void)
{
  delete force;

  delete orient;

  delete [] com;
  delete [] center;
  delete [] cov;
  delete [] coa;
  delete [] angmom;
  delete [] ps;

  delete [] com0;
  delete [] cov0;
  delete [] acc0;

  delete tree;
}

void Component::bomb(const string& msg)
{
  std::ostringstream sout;
  sout << "Component <" << name << ", " << id << ">: " << msg;
  throw GenericError(sout.str(), __FILE__, __LINE__);
}

void Component::read_bodies_and_distribute_ascii(void)
{
				// Open file
  std::ifstream fin;
  const int nline = 2048;
  char line[nline];
  
  if (myid == 0) {
    fin.open(pfile);

    if (fin.fail()) {
      std::ostringstream sout;
      sout << "Couldn't open " << pfile << " . . . quitting";
      throw GenericError(sout.str(), __FILE__, __LINE__);
    }

    fin.getline(line, nline);
    istringstream ins(line);
    
    ins >> nbodies_tot;		
    if (!ins) {
      std::ostringstream sout;
      sout << "Error reading nbodies_tot . . . quitting";
      throw GenericError(sout.str(), __FILE__, __LINE__);
    }
    ins >> niattrib;
    if (!ins) {
      std::ostringstream sout;
      sout << "Error reading integer attribute # . . . quitting\n";
      throw GenericError(sout.str(), __FILE__, __LINE__);
    }
    ins >> ndattrib;
    if (!ins) {
      std::ostringstream sout;
      sout << "Error reading double attribute # . . . quitting";
      throw GenericError(sout.str(), __FILE__, __LINE__);
    }
  }
				// Broadcast attributes for this
				// phase-space component
  MPI_Bcast(&nbodies_tot, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&niattrib,    1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ndattrib,    1, MPI_INT, 0, MPI_COMM_WORLD);

  double rmax1=0.0, r2;

  is_init = 1;
  setup_distribution();
  is_init = 0;
				// Initialize the particle ferry
				// instance with dynamic attribute
				// sizes
  if (not pf) pf = ParticleFerryPtr(new ParticleFerry(niattrib, ndattrib));

  if (myid==0) {
				// Read in Node 0's particles
    for (unsigned i=1; i<=nbodies_table[0]; i++) {

      PartPtr part = std::make_shared<Particle>(niattrib, ndattrib);
      
      part->readAscii(aindex, i, &fin);
				// Get the radius
      double r2 = 0.0;
      for (int j=0; j<3; j++) r2 += part->pos[j]*part->pos[j];
      rmax1 = max<double>(r2, rmax1);
      
				// Load the particle
      particles[part->indx] = part;

				// Record top_seq
      top_seq = std::max<unsigned long>(part->indx, top_seq);
    }

    nbodies = nbodies_table[0];

    unsigned icount, ibufcount;
    for (int n=1; n<numprocs; n++) {

      pf->ShipParticles(n, 0, nbodies_table[n]);

      icount = 0;
      ibufcount = 0;
      while (icount < nbodies_table[n]) {

	PartPtr part = std::make_shared<Particle>(niattrib, ndattrib);

	int i = nbodies_index[n-1] + 1 + icount;
	part->readAscii(aindex, i, &fin);

	r2 = 0.0;
	for (int k=0; k<3; k++) r2 += part->pos[k]*part->pos[k];
	rmax1 = max<double>(r2, rmax1);

	pf->SendParticle(part);
	icount++;

				// Record top_seq
	top_seq = std::max<unsigned long>(part->indx, top_seq);
      }

    }

  } else {

    pf->ShipParticles(myid, 0, nbodies);
      
#ifdef DEBUG
    int icount = 0;
#endif

    while (PartPtr part=pf->RecvParticle()) {
      particles[part->indx] = part;
#ifdef DEBUG
      if (icount<5) {
	cout << "Process " << myid << ": received ";
	cout << setw(14) << part->mass;
	for (int k=0; k<3; k++) cout << setw(14) << part->pos[k];
	cout << endl;
      }
      icount++;
#endif
    }
  }
				// Default: set to max radius
				// can be overriden by parameter

  rmax = sqrt(fabs(rmax1));
  MPI_Bcast(&rmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

				// Send top_seq to all nodes
  MPI_Bcast(&top_seq, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

				// COM HERE?
  if (myid==0) {
    try {
      fin.close();
    }
    catch (const ofstream::failure& e) {
      std::cout << "Component: exception closing file <" << pfile
		<< ": " << e.what() << std::endl;
    }
  }

  initialize();

#ifdef DEBUG
  if (particles.size()) {
    unsigned long imin = std::numeric_limits<unsigned long>::max();
    unsigned long imax = 0, kmin = imax, kmax = 0;
    for (auto p : particles) {
      imin = std::min<unsigned long>(imin, p.first);
      imax = std::max<unsigned long>(imax, p.first);
      kmin = std::min<unsigned long>(kmin, p.second->indx);
      kmax = std::max<unsigned long>(kmax, p.second->indx);
    }
    cout << "read_bodies_and_distribute_ascii: process " << myid 
	 << " name=" << name << " bodies [" << kmin << ", "
	 << kmax << "], [" << imin << ", " << imax << "]"
	 << " #=" << particles.size() << endl;
  } else {
    cout << "read_bodies_and_distribute_ascii: process " << myid 
	 << " name=" << name
	 << " #=" << particles.size() << endl;
  }
#endif
}

void Component::read_bodies_and_distribute_binary_out(istream *in)
{
				// Get component header
  ComponentHeader header;
				// Node local parameter buffer
  int ninfochar;
  std::shared_ptr<char> info;
  
  if (myid == 0) {

    rsize = sizeof(double);

    if (umagic) {
      unsigned long cmagic;
      in->read((char*)&cmagic, sizeof(unsigned long));
      if ( (cmagic & nmask) != magic ) {
	std::string msg("Error identifying new PSP.  Is this an old PSP?");
	throw GenericError(msg, __FILE__, __LINE__);
      }
      rsize = cmagic & mmask;
    }

    if (!header.read(in)) {
      std::string msg("Error reading component header");
      throw GenericError(msg, __FILE__, __LINE__);
    }

    nbodies_tot = header.nbod;
    niattrib    = header.niatr;
    ndattrib    = header.ndatr;
    ninfochar   = header.ninfochar;

    // Use this as of C++17
    // info = std::make_shared<char[]>(ninfochar+1);

    // C++14 Workaround
    info = std::shared_ptr<char>(new char[ninfochar+1],
				 std::default_delete<char[]>());

				// Zero fill array
    std::fill(info.get(), info.get()+ninfochar+1, 0);
				// Copy into array
    memcpy(info.get(), header.info.get(), ninfochar);
  }

  if (umagic)
    MPI_Bcast(&rsize, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

				// Broadcast attributes for this
				// phase-space component
  MPI_Bcast(&nbodies_tot, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&niattrib,    1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ndattrib,    1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ninfochar,   1, MPI_INT, 0, MPI_COMM_WORLD);
  if (myid) {
    // Use this as of C++17
    // info = std::make_shared<char[]>(ninfochar+1);

    // C++14 workaround:
    info = std::shared_ptr<char>(new char[ninfochar+1],
				 std::default_delete<char[]>());

				// Zero fill array
    std::fill(info.get(), info.get()+ninfochar+1, 0);
  }
  MPI_Bcast(info.get(), ninfochar, MPI_CHAR, 0, MPI_COMM_WORLD);

				// Parse info field to get 
				// id and parameter strings
  YAML::Node config;

  if (ignore_info and VERBOSE>3) {		// Ignore parameter info
    if (myid==0) std::cout << std::string(60, '-') << std::endl
			   << "ignore_info debug"  << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    config = conf;
  } else {			// Use parameter info
    std::istringstream sin(info.get());
    try {
      config = YAML::Load(sin);
    }
    catch (YAML::Exception& error) {
      if (myid==0)
	std::cerr << "YAML: error parsing <" << info.get() << "> "
		  << "in " << __FILE__ << ":" << __LINE__ << std::endl
		  << "YAML error: " << error.what() << std::endl;
      MPI_Finalize();
      exit(-10);
    }

    try {
      name  = config["name"].as<std::string>();
      cconf = config["parameters"];
      pfile = config["bodyfile"].as<std::string>();
    }
    catch (YAML::Exception & error) {
      if (myid==0) std::cout << "Error parsing YAML in PSP file: "
			     << error.what() << std::endl
			     << std::string(60, '-') << std::endl
			     << "Config node"        << std::endl
			     << std::string(60, '-') << std::endl
			     << config               << std::endl
			     << std::string(60, '-') << std::endl;
      MPI_Finalize();
      exit(-11);
    }

    YAML::Node force;
    
    try {
      force = config["force"];
      id    = force["id"].as<std::string>();
      fconf = force["parameters"];
    }
    catch (YAML::Exception & error) {
      if (myid==0) std::cout << "Error parsing YAML force stanza in PSP file: "
			     << error.what() << std::endl
			     << std::string(60, '-') << std::endl
			     << "Config node"        << std::endl
			     << std::string(60, '-') << std::endl
			     << config               << std::endl
			     << std::string(60, '-') << std::endl;
      
      MPI_Finalize();
      exit(-12);
    }

    // Assign local conf
    //
    conf["name"]       = name;
    conf["parameters"] = cconf;
    conf["bodyfile"]   = pfile;
    conf["force"]      = force;
				// Informational output
    if (myid==0)  {
      cconf.SetStyle(YAML::EmitterStyle::Flow);
      fconf.SetStyle(YAML::EmitterStyle::Flow);

      cout << std::string(60, '-') << endl
	   << "--- New Component"  << endl
	   << setw(20) << " name   :: " << name        << endl
	   << setw(20) << " id     :: " << id          << endl
	   << setw(20) << " cparam :: " << cconf       << endl
	   << setw(20) << " fparam :: " << fconf       << endl
	   << std::string(60, '-') << endl;
    }
  } // END: parse and assign parameter info from PSP
  
  double rmax1=0.0, r2;

  is_init = 1;
  setup_distribution();
  is_init = 0;
				// Initialize the particle ferry
				// instance with dynamic attribute
				// sizes
  if (not pf) pf = ParticleFerryPtr(new ParticleFerry(niattrib, ndattrib));

				// Form cumulative and differential
				// bodies list
  unsigned int ipart=0;

  if (myid==0) {
				// Read root node particles
    seq_cur = 0;

    rmax1 = 0.0;
    for (unsigned i=1; i<=nbodies_table[0]; i++)
    {
      PartPtr part = std::make_shared<Particle>(niattrib, ndattrib);
      
      part->readBinary(rsize, indexing, ++seq_cur, in);

      r2 = 0.0;
      for (int j=0; j<3; j++) r2 += part->pos[j]*part->pos[j];
      rmax1 = max<double>(r2, rmax1);

				// Load the particle
      particles[part->indx] = part;

				// Record top_seq
      top_seq = std::max<unsigned long>(part->indx, top_seq);
    }

    nbodies = nbodies_table[0];


				// Now load the other nodes
    unsigned icount;
    for (int n=1; n<numprocs; n++) {

      cout << "---- Component [" << name << "]: loading node <" << n << ">\n";

      pf->ShipParticles(n, 0, nbodies_table[n]);

      icount = 0;
      while (icount < nbodies_table[n]) {
	PartPtr part = std::make_shared<Particle>(niattrib, ndattrib);

	part->readBinary(rsize, indexing, ++seq_cur, in);

	r2 = 0.0;
	for (int k=0; k<3; k++) 
	  r2 += part->pos[k]*part->pos[k];

	rmax1 = max<double>(r2, rmax1);

	icount++;
	pf->SendParticle(part);

				// Record top_seq
	top_seq = std::max<unsigned long>(part->indx, top_seq);
      }

    }

  } else {

    pf->ShipParticles(myid, 0, nbodies);
      
    int icount = 0;
    PartPtr part;
    while (part=pf->RecvParticle()) {
      particles[part->indx] = part;
      icount++;
    }
  }


				// Default: set to max radius
  rmax = sqrt(fabs(rmax1));
  MPI_Bcast(&rmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

				// Send top_seq to all nodes
  MPI_Bcast(&top_seq, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

  initialize();

#ifdef DEBUG
  unsigned long imin = std::numeric_limits<unsigned long>::max();
  unsigned long imax = 0, kmin = imax, kmax = 0;
  for (auto p : particles) {
    imin = std::min<unsigned long>(imin, p.first);
    imax = std::max<unsigned long>(imax, p.first);
    kmin = std::min<unsigned long>(kmin, p.second->indx);
    kmax = std::max<unsigned long>(kmax, p.second->indx);
  }
  cout << "read_bodies_and_distribute_binary_out: process " << myid 
       << " name=" << name << " bodies [" << kmin << ", "
       << kmax << "], [" << imin << ", " << imax << "]"
       << " #=" << particles.size() << endl;
#endif
}


void Component::openNextBlob(std::ifstream& in,
			     std::list<std::string>::iterator& fit,
			     int& N)
{
  in.close();			// Close current file

  std::string curfile;

  try {
    if (outdir.back() != '/')	// Check whether directory has trailing '/'
      curfile = outdir + '/' + *fit;
    else
      curfile = outdir + *fit;
    
    in.open(curfile);
  } catch (...) {
    std::ostringstream sout;
    sout << "Could not open SPL blob <" << curfile << ">";
    throw std::runtime_error(sout.str());
  }

				// Double check file status
  if (not in.good()) {
    std::ostringstream sout;
    sout << "Could not open SPL blob <" << curfile << ">";
    throw std::runtime_error(sout.str());
  }

				// Get particle count
  try {
    in.read((char*)&N, sizeof(unsigned int));
  } catch (...) {
    std::ostringstream sout;
    sout << "Could not get particle count from <" << curfile << ">";
    throw std::runtime_error(sout.str());
  }

  // Advance filename iterator
  //
  fit++;
}


void Component::read_bodies_and_distribute_binary_spl(istream *in)
{
  // Will contain the component header
  //
  ComponentHeader header;

  // Node local parameter buffer
  //
  int ninfochar;
  std::unique_ptr<char[]> info;
  
  // Number of split files
  int number = 0;

  if (myid == 0) {

    rsize = sizeof(double);

    unsigned long cmagic;
    try {
      in->read((char*)&cmagic, sizeof(unsigned long));
      in->read((char*)&number, sizeof(int));
    } catch (...) {
      std::ostringstream sout;
      sout << "Error reading magic info and file count from master";
      throw std::runtime_error(sout.str());
    }

    if (umagic) {
      if ( (cmagic & nmask) != magic ) {
	std::string msg("Error identifying new PSP.  Is this an old PSP?");
	throw GenericError(msg, __FILE__, __LINE__);
      }
      rsize = cmagic & mmask;
    }

    if (!header.read(in)) {
      std::string msg("Error reading component header");
      throw GenericError(msg, __FILE__, __LINE__);
    }

    nbodies_tot = header.nbod;
    niattrib    = header.niatr;
    ndattrib    = header.ndatr;
    ninfochar   = header.ninfochar;

    info = std::make_unique<char[]>(ninfochar+1);

    // Zero fill array
    //
    std::fill(info.get(), info.get() + (ninfochar+1), 0);

    // Copy into array
    //
    memcpy(info.get(), header.info.get(), ninfochar);
  }

  if (umagic)
    MPI_Bcast(&rsize, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

  // Broadcast attributes for this phase-space component
  //
  MPI_Bcast(&nbodies_tot, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&niattrib,    1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ndattrib,    1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ninfochar,   1, MPI_INT, 0, MPI_COMM_WORLD);
  if (myid) {
    info = std::make_unique<char[]>(ninfochar+1);

    // Zero fill array
    //
    std::fill(info.get(), info.get() + (ninfochar+1), 0);
  }
  MPI_Bcast(info.get(), ninfochar, MPI_CHAR, 0, MPI_COMM_WORLD);

  // Parse info field to get id and parameter strings
  //
  YAML::Node config;

  if (ignore_info and VERBOSE>3) { // Ignore parameter info
    if (myid==0) std::cout << std::string(60, '-') << std::endl
			   << "ignore_info debug"  << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    config = conf;
  } else {			// Use parameter info
    std::istringstream sin(info.get());
    try {
      config = YAML::Load(sin);
    }
    catch (YAML::Exception& error) {
      std::cerr << "YAML: error parsing <" << info.get() << "> "
		<< "in " << __FILE__ << ":" << __LINE__ << std::endl
		<< "YAML error: " << error.what() << std::endl;
      throw error;
    }

    try {
      name  = config["name"].as<std::string>();
      cconf = config["parameters"];
      pfile = config["bodyfile"].as<std::string>();
    }
    catch (YAML::Exception & error) {
      if (myid==0) std::cout << "Error parsing YAML in PSP file: "
			     << error.what() << std::endl
			     << std::string(60, '-') << std::endl
			     << "Config node"        << std::endl
			     << std::string(60, '-') << std::endl
			     << config               << std::endl
			     << std::string(60, '-') << std::endl;
      MPI_Finalize();
      exit(-13);
    }

    YAML::Node force;
    
    try {
      force = config["force"];
      id    = force["id"].as<std::string>();
      fconf = force["parameters"];
    }
    catch (YAML::Exception & error) {
      if (myid==0) std::cout << "Error parsing YAML force stanza in PSP file: "
			     << error.what() << std::endl
			     << std::string(60, '-') << std::endl
			     << "Config node"        << std::endl
			     << std::string(60, '-') << std::endl
			     << config               << std::endl
			     << std::string(60, '-') << std::endl;
      
      MPI_Finalize();
      exit(-14);
    }

    // Assign local conf
    //
    conf["name"]       = name;
    conf["parameters"] = cconf;
    conf["bodyfile"]   = pfile;
    conf["force"]      = force;

    // Informational output
    //
    if (myid==0)  {
      cconf.SetStyle(YAML::EmitterStyle::Flow);
      fconf.SetStyle(YAML::EmitterStyle::Flow);

      cout << std::string(60, '-') << endl
	   << "--- New Component"  << endl
	   << setw(20) << " name   :: " << name        << endl
	   << setw(20) << " id     :: " << id          << endl
	   << setw(20) << " cparam :: " << cconf       << endl
	   << setw(20) << " fparam :: " << fconf       << endl
	   << std::string(60, '-') << endl;
    }
  }
  // END: parse and assign parameter info from PSP
  
  // Get file names for split PSP parts
  //
  std::list<std::string> parts;

  if (myid==0) {
    const size_t PBUF_SIZ = 1024;
    char buf [PBUF_SIZ];

    for (int n=0; n<number; n++) {
      in->read((char *)buf, PBUF_SIZ);
      parts.push_back(buf);
    }
  }

  double rmax1=0.0, r2;

  is_init = 1;
  setup_distribution();
  is_init = 0;
				// Initialize the particle ferry
				// instance with dynamic attribute
				// sizes
  if (not pf) pf = ParticleFerryPtr(new ParticleFerry(niattrib, ndattrib));

				// Form cumulative and differential
				// bodies list
  unsigned int ipart=0;

  if (myid==0) {

				// Set iterator to beginning of split list
    auto fit = parts.begin();

    std::ifstream fin;		// Open file stream
    int N;			// Number of particles in this stream

				// Get next blob (which is the first blob here)
    openNextBlob(fin, fit, N);

    int fcount = 0;		// Count number of particles read from this stream

				// Read root node particles
    seq_cur = 0;

    rmax1 = 0.0;
    for (unsigned i=1; i<=nbodies_table[0]; i++)
    {
      PartPtr part = std::make_shared<Particle>(niattrib, ndattrib);
      
      part->readBinary(rsize, indexing, ++seq_cur, &fin);

				// Check particle count in current split file
      if (++fcount == N) {
	fcount = N = 0;		// Maybe last particle . . .
	if (fit != parts.end()) openNextBlob(fin, fit, N);
      }

      r2 = 0.0;
      for (int j=0; j<3; j++) r2 += part->pos[j]*part->pos[j];
      rmax1 = max<double>(r2, rmax1);

				// Load the particle
      particles[part->indx] = part;

				// Record top_seq
      top_seq = std::max<unsigned long>(part->indx, top_seq);
    }

    nbodies = nbodies_table[0];


				// Now load the other nodes
    unsigned icount;
    for (int n=1; n<numprocs; n++) {

      cout << "---- Component [" << name << "]: loading node <" << n << ">\n";

      pf->ShipParticles(n, 0, nbodies_table[n]);

      icount = 0;
      while (icount < nbodies_table[n]) {
	PartPtr part = std::make_shared<Particle>(niattrib, ndattrib);

	part->readBinary(rsize, indexing, ++seq_cur, &fin);

				// Check particle count in current split file
	if (++fcount == N) {
	  fcount = N = 0;	// May be last particle . . .
	  if (fit != parts.end()) openNextBlob(fin, fit, N);
	}

	r2 = 0.0;
	for (int k=0; k<3; k++) 
	  r2 += part->pos[k]*part->pos[k];

	rmax1 = max<double>(r2, rmax1);

	icount++;		// Send the particle
	pf->SendParticle(part);

				// Record top_seq
	top_seq = std::max<unsigned long>(part->indx, top_seq);
      }

    }

  } else {

    pf->ShipParticles(myid, 0, nbodies);
      
    int icount = 0;
    PartPtr part;
    while (part=pf->RecvParticle()) {
      particles[part->indx] = part;
      icount++;
    }
  }


				// Default: set to max radius
  rmax = sqrt(fabs(rmax1));
  MPI_Bcast(&rmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

				// Send top_seq to all nodes
  MPI_Bcast(&top_seq, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

  initialize();

#ifdef DEBUG
  unsigned long imin = std::numeric_limits<unsigned long>::max();
  unsigned long imax = 0, kmin = imax, kmax = 0;
  for (auto p : particles) {
    imin = std::min<unsigned long>(imin, p.first);
    imax = std::max<unsigned long>(imax, p.first);
    kmin = std::min<unsigned long>(kmin, p.second->indx);
    kmax = std::max<unsigned long>(kmax, p.second->indx);
  }
  cout << "read_bodies_and_distribute_binary_spl: process " << myid 
       << " name=" << name << " bodies [" << kmin << ", "
       << kmax << "], [" << imin << ", " << imax << "]"
       << " #=" << particles.size() << endl;
#endif
}


PartPtr * Component::get_particles(int* number)
{
  static std::vector<unsigned> totals;
  static unsigned counter = 0;	// Counter for all
  static int      node    = 0;	// Current node
  
				// Reset
  if (*number < 0) {
    counter = 0;
    node    = 0;
  }
				// Done?
  if (counter == nbodies_tot) {
    *number = 0;
    return 0;
  }

  int curcount = 0;		// Counter for this bunch

				// Do this on first call ONLY
  if (counter == 0) {
				// Add new particles to sequence
    seq_new_particles();

#ifdef DEBUG
    std::cout << "get_particles: Node " << myid
	      << " # part=" << particles.size() << " nbodies=" << nbodies
		  << std::endl << std::flush;
#endif    

				// Every process report particle numbers
    totals.resize(numprocs);
    MPI_Allgather(&nbodies, 1, MPI_UNSIGNED, &totals[0], 1, MPI_UNSIGNED,
		  MPI_COMM_WORLD);
				// Cumulate
    for (int n=1; n<numprocs; n++) totals[n] += totals[n-1];
  }


  std::map<unsigned long, PartPtr> tlist;

  unsigned icount;
  unsigned beg = counter;
  unsigned end = counter + PFbufsz;

  bool complete = false;
  if (end >= totals[node]) {
    complete = true;
    end = totals[node];
  }


  if (myid==0) {
				// Do root's particle first
    if (node==0) {

	auto itb = particles.begin();
	auto ite = particles.begin();
	
	std::advance(itb, beg);
	if (complete) ite = particles.end();
	else          std::advance(ite, end);

	icount = 0;
	for (auto it=itb; it!=ite; it++) pbuf[icount++] = it->second;

#ifdef DEBUG
	std::cout << "get_particles: master loaded " 
		  << icount << " of its own particles"
		  << ", beg=" << beg << ", iend=" << end
		  << ", dist=" << std::distance(itb, ite)
		  << ", expected=" << totals[0]
		  << std::endl << std::flush;
#endif    

    } else {
      
      unsigned number;
      pf->ShipParticles(0, node, number);
      
      icount = 0;
      while (PartPtr part=pf->RecvParticle()) pbuf[icount++] = part;
#ifdef DEBUG
      std::cout << "Process " << myid 
	   << ": received " << icount << " particles from Slave " << node
		<< ", expected " << number << ", total=" << totals[node]
		<< std::endl << std::flush;
#endif    
    }

    // Load the array
    for (unsigned n=0; n<icount; n++) {
      tlist[pbuf[n]->indx] = pbuf[n];
      curcount++;
      counter++;
    }
				// Nodes send particles to master
  } else if (myid == node) {
      
    auto itb = particles.begin();
    auto ite = particles.begin();

    std::advance(itb, beg - totals[node-1]);
    if (complete) ite = particles.end();
    else          std::advance(ite, end - totals[node-1]);
    

    icount = std::distance(itb, ite);
    pf->ShipParticles(0, myid, icount);

    for (auto it=itb; it!=ite; it++) pf->SendParticle(it->second);

#ifdef DEBUG
      std::cout << "Process " << myid 
		<< ": sent " << icount << " particles from Slave " << node
		<< std::endl << std::flush;
#endif    
  }

  MPI_Bcast(&counter, 1, MPI_INT, 0, MPI_COMM_WORLD);

				// Return values
  *number = curcount;

#ifdef DEBUG
  if (myid==0) {
    std::cout << "get_particles: master size of tlist=" << tlist.size() 
	      << " current count=" << curcount << std::endl;
  }
#endif

  int n = 0;
  for (auto cur : tlist) pbuf[n++] = cur.second;
  
  // Move to next node?
  //
  if (complete) node++;

#ifdef DEBUG
  std::cout << "Process " << myid 
       << ": received next counter=" << counter;
  if (counter > nbodies_tot) cout << " [this means we are done]";
  std::cout << std::endl << std::flush;
  MPI_Barrier(MPI_COMM_WORLD);
#endif    

  return &pbuf[0];
}


void Component::write_binary(ostream* out, bool real4)
{
  ComponentHeader header;

  if (myid == 0) {

    header.nbod  = nbodies_tot;
    header.niatr = niattrib;
    header.ndatr = ndattrib;
  
    std::ostringstream outs;
    if (conf.Type() != YAML::NodeType::Null) outs << conf << std::endl;

    // Resize info string, if necessary
    size_t infosz = outs.str().size() + 4;
    if (header.ninfochar < outs.str().size()) {
      header.ninfochar = outs.str().size();
      // Use this as of C++17
      // header.info = std::make_shared<char[]>(header.ninfochar+1);

      // C++14 workaround:
      header.info = std::shared_ptr<char>(new char[header.ninfochar+1],
					  std::default_delete<char[]>());
    }

    // Copy to info string
    strncpy(header.info.get(), outs.str().c_str(), header.ninfochar);

    // DEBUGGING
    if (false and myid==0) {
      std::cout << std::string(72, '-') << std::endl
		<< "Serialized YAML header looks like this:" << std::endl
		<< std::string(72, '-') << std::endl
		<< outs.str() << std::endl
		<< "Cur size=" << outs.str().size()
		<< " max size=" << header.ninfochar << std::endl
		<< std::string(72, '-') << std::endl;
    }

    if (real4) rsize = sizeof(float);
    else       rsize = sizeof(double);
    unsigned long cmagic = magic + rsize;

    out->write((const char*)&cmagic, sizeof(unsigned long));

    if (!header.write(out)) {
      std::string msg("Component::write_binary: Error writing particle header");
      throw GenericError(msg, __FILE__, __LINE__);
    }
  }

				// First bunch of particles
  int number = -1;
  PartPtr *p = get_particles(&number);

  float tf;
  double pot0, pv;
  while (p) {

    if (myid == 0) {
      for (int k=0; k<number; k++) {
	p[k]->writeBinary(rsize, indexing, out);
      }
    }
				// Next bunch of particles
    p = get_particles(&number);

  }
    
}

void Component::write_binary_header(ostream* out, bool real4, const std::string prefix, int nth)
{
  ComponentHeader header;

  if (myid == 0) {

    header.nbod  = nbodies_tot;
    header.niatr = niattrib;
    header.ndatr = ndattrib;
  
    std::ostringstream outs;
    if (conf.Type() != YAML::NodeType::Null) outs << conf << std::endl;

    // Resize info string, if necessary
    size_t infosz = outs.str().size() + 4;
    if (header.ninfochar < outs.str().size()) {
      header.ninfochar = outs.str().size();
      // Use this as of C++17
      // header.info = std::make_shared<char[]>(header.ninfochar+1);

      // C++14 workaround:
      header.info = std::shared_ptr<char>(new char[header.ninfochar+1],
					  std::default_delete<char[]>());
    }

    // Copy to info string
    strncpy(header.info.get(), outs.str().c_str(), header.ninfochar);

    // DEBUGGING
    if (false and myid==0) {
      std::cout << std::string(72, '-') << std::endl
		<< "Serialized YAML header looks like this:" << std::endl
		<< std::string(72, '-') << std::endl
		<< outs.str() << std::endl
		<< "Cur size=" << outs.str().size()
		<< " max size=" << header.ninfochar << std::endl
		<< std::string(72, '-') << std::endl;
    }

    if (real4) rsize = sizeof(float);
    else       rsize = sizeof(double);
    unsigned long cmagic = magic + rsize;

    int nfiles = numprocs*nth;

    out->write((const char*)&cmagic,  sizeof(unsigned long));
    out->write((const char*)&nfiles,  sizeof(int));

    if (!header.write(out)) {
      std::string msg("Component::write_binary: Error writing particle header");
      throw GenericError(msg, __FILE__, __LINE__);
    }

    const size_t PBUF_SIZ = 1024;
    char buf [PBUF_SIZ];

    for (int n=0; n<nfiles; n++) {
      std::ostringstream sout;
      sout << prefix << "-" << n << '\0';
      sout.str().copy(buf, sout.str().size());
      out->write((const char*)buf, PBUF_SIZ);
    }
  }

}

void Component::write_binary_particles(std::ostream* out, bool real4)
{
  unsigned int N = particles.size();
  out->write((const char*)&N, sizeof(unsigned int));

  if (real4) rsize = sizeof(float);
  else       rsize = sizeof(double);

  // Use buffered writes
  //
  if (buffered) {
    ParticleBuffer buf(rsize, indexing, particles.begin()->second.get());
    for (auto p : particles)
      p.second->writeBinaryBuffered(rsize, indexing, out, buf);
    buf.writeBuffer(out, true);	// Complete the write
  }
  // Unbuffered, direct writes
  //
  else {
    for (auto p : particles) p.second->writeBinary(rsize, indexing, out);
  }
}

void Component::write_binary_particles
(std::vector<std::shared_ptr<std::ofstream>>& out, bool real4)
{
  // Set number of OpenMP threads equal to number of streams.  It's up
  // to the user to make that cores are available.
  //
  size_t nthrds = out.size();
  omp_set_num_threads(nthrds);

  // Divide the particles between the streams
  //
  std::vector<unsigned> NN(nthrds, 0);

  // Move to write position
  //
#pragma omp parallel for
  for (size_t n=0; n<nthrds; n++) {
    out[n]->seekp(sizeof(unsigned int)); // Skip size of particle
					 // number
  }

  // Assign data size
  //
  if (real4) rsize = sizeof(float);
  else       rsize = sizeof(double);

  // Use buffered writes (default)
  //
  if (buffered) {

    // Make a buffer for each stream
    //
    std::vector<std::shared_ptr<ParticleBuffer>> buf(nthrds);
#pragma omp parallel for
    for (size_t n=0; n<nthrds; n++) {
      buf[n] = std::make_shared<ParticleBuffer>
	(rsize, indexing, particles.begin()->second.get());
    }

    // Write the particles in parallel
    //
#pragma omp parallel for
    for (size_t b=0; b<particles.bucket_count(); b++) {
      auto tid = omp_get_thread_num();
      for (auto p=particles.begin(b); p!=particles.end(b); p++) {
	p->second->writeBinaryBuffered(rsize, indexing, 
				       out[tid].get(), *buf[tid]);
	// Tally number of particles per stream
	NN[tid]++;
      }
    }

    // Flush the buffer
    //
#pragma omp parallel for
    for (size_t n=0; n<nthrds; n++) {
      buf[n]->writeBuffer(out[n].get(), true);
    }
  }
  // Unbuffered, direct writes
  //
  else {
#pragma omp parallel for
    for (size_t b=0; b<particles.bucket_count(); b++) {
      auto tid = omp_get_thread_num();
      for (auto p=particles.begin(b); p!=particles.end(b); p++) {
	p->second->writeBinary(rsize, indexing, out[tid].get());
	// Tally number of particles per stream
	NN[tid]++;
      }
    }
  }

  // Write the particle counts into each file blob
  //
#pragma omp parallel for
  for (size_t n=0; n<nthrds; n++) {
    out[n]->seekp(0);		// Position at beginning
    out[n]->write((const char *)&NN[n], sizeof(unsigned int));
  }
}


// Helper class that manages two buffers that can be swapped to
// support non-blocking MPI-IO writes
//
class DoubleBuf
{
private:
  // The storage
  std::vector<char> src1, src2;

  // Pointers to the char buffers
  char *curr, *next;

public:

  // Initialize with size chars
  DoubleBuf(int size)
  {
    src1.resize(size);
    src2.resize(size);
    curr = &src1[0];
    next = &src2[0];
  }

  // Get the current buffer
  char* operator()()
  {
    return curr;
  }

  // Swap buffers and return the new current buffer
  char* swap()
  {
    char* temp = curr;
    curr = next;
    next = temp;
    return curr;
  }

};


void Component::write_binary_mpi_b(MPI_File& out, MPI_Offset& offset, bool real4)
{
  ComponentHeader header;
  MPI_Status status;
  char err[MPI_MAX_ERROR_STRING];
  int len;

  if (real4) rsize = sizeof(float);
  else       rsize = sizeof(double);

  if (myid == 0) {

    header.nbod  = nbodies_tot;
    header.niatr = niattrib;
    header.ndatr = ndattrib;
  
    std::ostringstream outs;
    outs << conf << std::endl;
    strncpy(header.info.get(), outs.str().c_str(), header.ninfochar);

    // DEBUGGING
    if (false and myid==0) {
      std::cout << std::string(72, '-') << std::endl
		<< "Serialized YAML header looks like this:" << std::endl
		<< std::string(72, '-') << std::endl
		<< outs.str() << std::endl
		<< "Cur size=" << outs.str().size()
		<< " max size=" << header.ninfochar << std::endl
		<< std::string(72, '-') << std::endl;
    }

    unsigned long cmagic = magic + rsize;

    int ret =
      MPI_File_write_at(out, offset, &cmagic, 1, MPI_UNSIGNED_LONG, &status);

    if (ret != MPI_SUCCESS) {
      MPI_Error_string(ret, err, &len);
      std::cout << "Component::write_binary_mpi_b: " << err
		<< " at line " << __LINE__ << std::endl;
    }

    offset += sizeof(unsigned long);

    if (!header.write_mpi(out, offset)) {
      std::string msg("Component::write_binary_mpi_b: Error writing particle header");
      throw GenericError(msg, __FILE__, __LINE__);
    }

  } else {
    offset += sizeof(unsigned long) + header.getSize();
  }

  unsigned N = particles.size();
  std::vector<unsigned> numP(numprocs, 0);

  MPI_Allgather(&N, 1, MPI_UNSIGNED, &numP[0], 1, MPI_UNSIGNED,	MPI_COMM_WORLD);
  
  for (int i=1; i<numprocs; i++) numP[i] += numP[i-1];
  unsigned bSiz = particles.begin()->second->getMPIBufSize(rsize, indexing);
  if (myid) offset += numP[myid-1] * bSiz;
  
  std::vector<char> buffer(pBufSiz*bSiz);
  size_t count = 0;
  char *buf = &buffer[0], *bufl;

  for (auto & p : particles) {
    buf += p.second->writeBinaryMPI(buf, rsize, indexing);
    count++;

    if (count==pBufSiz) {
      int ret =
	MPI_File_write_at(out, offset, &buffer[0], bSiz*count, MPI_CHAR, &status);

      if (ret != MPI_SUCCESS) {
	MPI_Error_string(ret, err, &len);
	std::cout << "Component::write_binary_mpi_b: " << err
		  << " at line " << __LINE__ << std::endl;
      }

      offset += bSiz*count;
      count   = 0;
      buf     = &buffer[0];
    }
  }

  if (count) {
    int ret = MPI_File_write_at(out, offset, &buffer[0], bSiz*count, MPI_CHAR, &status);

    if (ret != MPI_SUCCESS) {
      MPI_Error_string(ret, err, &len);
      std::cout << "Component::write_binary_mpi_b: " << err
		<< " at line " << __LINE__ << std::endl;
    }

    offset += bSiz*count;
  }

  // Position file offset at end of particles
  //
  offset += (numP[numprocs-1] - numP[myid]) * bSiz;
}


void Component::write_binary_mpi_i(MPI_File& out, MPI_Offset& offset, bool real4)
{
  ComponentHeader header;
  MPI_Request request = MPI_REQUEST_NULL;
  MPI_Status status;
  char err[MPI_MAX_ERROR_STRING];
  int len;

  if (real4) rsize = sizeof(float);
  else       rsize = sizeof(double);

  if (myid == 0) {

    header.nbod  = nbodies_tot;
    header.niatr = niattrib;
    header.ndatr = ndattrib;
  
    std::ostringstream outs;
    outs << conf << std::endl;
    strncpy(header.info.get(), outs.str().c_str(), header.ninfochar);

    // DEBUGGING
    if (true and myid==0) {
      std::cout << std::string(72, '-') << std::endl
		<< "Serialized YAML header looks like this:" << std::endl
		<< std::string(72, '-') << std::endl
		<< outs.str() << std::endl
		<< "Cur size=" << outs.str().size()
		<< " max size=" << header.ninfochar << std::endl
		<< std::string(72, '-') << std::endl;
    }

    unsigned long cmagic = magic + rsize;

    int ret =
      MPI_File_write_at(out, offset, &cmagic, 1, MPI_UNSIGNED_LONG, &status);

    if (ret != MPI_SUCCESS) {
      MPI_Error_string(ret, err, &len);
      std::cout << "Component::write_binary_mpi_i: " << err
		<< " at line " << __LINE__ << std::endl;
    }

    offset += sizeof(unsigned long);

    if (!header.write_mpi(out, offset)) {
      std::string msg("Component::write_binary_mpi_i: Error writing particle header");
      throw GenericError(msg, __FILE__, __LINE__);
    }

  } else {
    offset += sizeof(unsigned long) + header.getSize();
  }

  unsigned N = particles.size();
  std::vector<unsigned> numP(numprocs, 0);

  MPI_Allgather(&N, 1, MPI_UNSIGNED, &numP[0], 1, MPI_UNSIGNED,	MPI_COMM_WORLD);
  
  for (int i=1; i<numprocs; i++) numP[i] += numP[i-1];
  unsigned bSiz = particles.begin()->second->getMPIBufSize(rsize, indexing);
  if (myid) offset += numP[myid-1] * bSiz;
  
  DoubleBuf buffer(pBufSiz*bSiz);
  char *buf = buffer();
  size_t count = 0;

  for (auto & p : particles) {
    buf += p.second->writeBinaryMPI(buf, rsize, indexing);
    count++;

    if (count==pBufSiz) {
      // Check for completion of last write
      //
      int ret = MPI_Wait(&request, &status);
      
      if (ret != MPI_SUCCESS) {
	MPI_Error_string(ret, err, &len);
	std::cout << "Component::write_binary_mpi_i: " << err
		  << " at line " << __LINE__ << std::endl;
      }
      
      // Non-blocking write allows next buffer to be filled
      //
      ret =
	MPI_File_iwrite_at(out, offset, buffer(), bSiz*count, MPI_CHAR, &request);

      if (ret != MPI_SUCCESS) {
	MPI_Error_string(ret, err, &len);
	std::cout << "Component::write_binary_mpi_i: " << err
		  << " at line " << __LINE__ << std::endl;
      }

      offset += bSiz*count;
      count   = 0;
      buf     = buffer.swap();
    }
  }

  // Check for completion of last write
  //
  int ret = MPI_Wait(&request, &status);
      
  if (ret != MPI_SUCCESS) {
    MPI_Error_string(ret, err, &len);
    std::cout << "Component::write_binary_mpi_i: " << err
	      << " at line " << __LINE__ << std::endl;
  }

  if (count) {

    // Block on final write
    //
    ret = MPI_File_write_at(out, offset, buffer(), bSiz*count, MPI_CHAR, &status);

    if (ret != MPI_SUCCESS) {
      MPI_Error_string(ret, err, &len);
      std::cout << "Component::write_binary_mpi_i: " << err
		<< " at line " << __LINE__ << std::endl;
    }

    offset += bSiz*count;
  }

  // Position file offset at end of particles
  //
  offset += (numP[numprocs-1] - numP[myid]) * bSiz;
}


void Component::write_ascii(ostream* out, bool accel)
{
  int number = -1;
  PartPtr *p = get_particles(&number);

  while (p) {
    if (myid == 0) {
      for (int k=0; k<number; k++) {
	p[k]->writeAscii(indexing, accel, out);
      }
    }

    p = get_particles(&number);
  }
    
}


void Component::initialize_com_system()
{
  double mtot1;
  double *com1 = new double [3];
  double *cov1 = new double [3];
  
  
  PartMapItr p, pend;

				// Zero stuff out
  mtot0 = mtot1 = 0.0;
  for (int k=0; k<dim; k++) com1[k] = cov1[k] = 0.0;

				// Particle loop
  pend = particles.end();
  for (p=particles.begin(); p != pend; p++) {
    
    mtot1 += p->second->mass;

    for (int k=0; k<dim; k++) com1[k] += p->second->mass*p->second->pos[k];
    for (int k=0; k<dim; k++) cov1[k] += p->second->mass*p->second->vel[k];
    
  }
  
  MPI_Allreduce(&mtot1, &mtot0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(com1, com0, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(cov1, cov0, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  if (mtot0 > 0.0) {
    for (int k=0; k<dim; k++) com0[k] /= mtot0;
    for (int k=0; k<dim; k++) cov0[k] /= mtot0;
  }

  for (int k=0; k<dim; k++) {
    center[k] = 0.0;
  }

  delete [] com1;
  delete [] cov1;
}

void Component::restart_com_system()
{
  MPI_Bcast(&com_restart, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (com_restart) {
    MPI_Bcast(&com0[0],   3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cov0[0],   3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&acc0[0],   3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&center[0], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }

}


struct thrd_pass_posn
{
  int id;
  Component *c;
  int  tidal;
  bool consp;
  bool com_system;
  unsigned mlevel;
  vector<double> com,  cov,  coa,  mtot;
  vector<double> comE, covE, mtotE;
};



void * fix_positions_thread(void *ptr)
{
  int id          =   static_cast<thrd_pass_posn*>(ptr)->id;
  Component *c    =   static_cast<thrd_pass_posn*>(ptr)->c;

  int  tidal      =   static_cast<thrd_pass_posn*>(ptr)->tidal;
  bool consp      =   static_cast<thrd_pass_posn*>(ptr)->consp;
  bool com_system =   static_cast<thrd_pass_posn*>(ptr)->com_system;

  unsigned mlevel =   static_cast<thrd_pass_posn*>(ptr)->mlevel;

  double *com     = &(static_cast<thrd_pass_posn*>(ptr)->com[0]);
  double *cov     = &(static_cast<thrd_pass_posn*>(ptr)->cov[0]);
  double *coa     = &(static_cast<thrd_pass_posn*>(ptr)->coa[0]);
  double *mtot    = &(static_cast<thrd_pass_posn*>(ptr)->mtot[0]);


  double *comE, *covE, *mtotE;

  if (consp && com_system) {
    comE          = &(static_cast<thrd_pass_posn*>(ptr)->comE[0]);
    covE          = &(static_cast<thrd_pass_posn*>(ptr)->covE[0]);
    mtotE         = &(static_cast<thrd_pass_posn*>(ptr)->mtotE[0]);
  }

  for (unsigned mm=mlevel; mm<=multistep; mm++) {

    int nbodies = c->levlist[mm].size();
    int nbeg    = nbodies*(id  )/nthrds;
    int nend    = nbodies*(id+1)/nthrds;

				// Particle loop
    for (int q=nbeg; q<nend; q++) {
    
      unsigned long n = c->levlist[mm][q];
      Particle     *p = c->Part(n);

      if (consp and tidal>=0) {
	if (c->escape_com(*p) && p->iattrib[tidal]==0) {
				// Set flag indicating escaped particle
	  p->iattrib[tidal] = 1;

	  if (com_system) {	// Conserve momentum of center of mass
				// and compute center of acceleration
	    mtotE[mm] += p->mass;
	    for (unsigned k=0; k<3; k++) {
	      comE[3*mm+k] += p->mass*p->pos[k]; 
	      covE[3*mm+k] += p->mass*p->vel[k]; 
	    }
	  }
	  continue;
	}
	
	if (p->iattrib[tidal]==1) continue;
      }

      if (c->freeze(n)) continue;

      mtot[mm] += p->mass;

      // Compute new center of mass quantities
      //
      for (int k=0; k<c->dim; k++) {
	com[3*mm+k] += p->mass*p->pos[k];
	cov[3*mm+k] += p->mass*p->vel[k];
	coa[3*mm+k] += p->mass*p->acc[k];
      }
    }
  }

  return (NULL);
}
  

void Component::fix_positions_cpu(unsigned mlevel)
{
				// Zero center
  for (int i=0; i<3; i++) center[i] = 0.0;

  				// Zero variables
  mtot = 0.0;
  for (int k=0; k<dim; k++) com[k] = cov[k] = coa[k] = 0.0;

				// Zero multistep counters at and
				// above this level
  for (unsigned mm=mlevel; mm<=multistep; mm++) {
    com_mas[mm] = 0.0;
    for (unsigned k=0; k<3; k++) 
      com_lev[3*mm+k] = cov_lev[3*mm+k] = coa_lev[3*mm+k] = 0.0;
  }

  vector<thrd_pass_posn> data(nthrds);
  vector<pthread_t>      thrd(nthrds);

  if (nthrds==1) {

    data[0].id         = 0;
    data[0].c          = this;
    data[0].consp      = consp;
    data[0].tidal      = tidal;
    data[0].com_system = com_system;
    data[0].mlevel     = mlevel;

    data[0].com  = vector<double>(3*(multistep+1), 0.0);
    data[0].cov  = vector<double>(3*(multistep+1), 0.0);
    data[0].coa  = vector<double>(3*(multistep+1), 0.0);
    data[0].mtot = vector<double>(multistep+1, 0.0);

    if (consp && com_system) {
      data[0].comE  = vector<double>(3*(multistep+1), 0.0);
      data[0].covE  = vector<double>(3*(multistep+1), 0.0);
      data[0].mtotE = vector<double>(multistep+1, 0.0);
    }

    fix_positions_thread(&data[0]);

    for (unsigned mm=mlevel; mm<=multistep; mm++) {
      for (unsigned k=0; k<3; k++) {
	com_lev[3*mm + k] += data[0].com[3*mm + k];
	cov_lev[3*mm + k] += data[0].cov[3*mm + k];
	coa_lev[3*mm + k] += data[0].coa[3*mm + k];
      }
      com_mas[mm] += data[0].mtot[mm];

      if (consp && com_system) {
	for (unsigned k=0; k<3; k++) {
	  comE_lev[3*mm + k] += data[0].comE[3*mm + k];
	  covE_lev[3*mm + k] += data[0].covE[3*mm + k];
	}
	comE_mas[mm] += data[0].mtotE[mm];
      }
    }

  } else {

    int errcode;
    void *retval;
  
    for (int i=0; i<nthrds; i++) {

      data[i].id         = i;
      data[i].c          = this;
      data[i].consp      = consp;
      data[i].tidal      = tidal;
      data[i].com_system = com_system;
      data[i].mlevel     = mlevel;

      data[i].com  = vector<double>(3*(multistep+1), 0.0);
      data[i].cov  = vector<double>(3*(multistep+1), 0.0);
      data[i].coa  = vector<double>(3*(multistep+1), 0.0);
      data[i].mtot = vector<double>(multistep+1, 0.0);

      if (consp && com_system) {
	data[i].comE  = vector<double>(3*(multistep+1), 0.0);
	data[i].covE  = vector<double>(3*(multistep+1), 0.0);
	data[i].mtotE = vector<double>(multistep+1, 0.0);
      }

      errcode =  pthread_create(&thrd[i], 0, fix_positions_thread, &data[i]);

      if (errcode) {
	std::ostringstream sout;
	sout << "Process " << myid
	     << " Component::fix_positions: cannot make thread " << i
	     << ", errcode=" << errcode;
	throw GenericError(sout.str(), __FILE__, __LINE__);
      }
    }
    
    //
    // Collapse the threads
    //
    for (int i=0; i<nthrds; i++) {
      if ((errcode=pthread_join(thrd[i], &retval))) {
	std::ostringstream sout;
	sout << "Process " << myid
	     << " Component::fix_positions: thread join " << i
	     << " failed, errcode=" << errcode;
	throw GenericError(sout.str(), __FILE__, __LINE__);
      }

      for (unsigned mm=mlevel; mm<=multistep; mm++) {
	for (unsigned k=0; k<3; k++) {
	  com_lev[3*mm + k] += data[i].com[3*mm + k];
	  cov_lev[3*mm + k] += data[i].cov[3*mm + k];
	  coa_lev[3*mm + k] += data[i].coa[3*mm + k];
	}
	com_mas[mm] += data[i].mtot[mm];
      }

      if (consp && com_system) {
	for (unsigned mm=mlevel; mm<=multistep; mm++) {
	  for (unsigned k=0; k<3; k++) {
	    comE_lev[3*mm + k] += data[i].comE[3*mm + k];
	    covE_lev[3*mm + k] += data[i].covE[3*mm + k];
	  }
	  comE_mas[mm] += data[i].mtotE[mm];
	}
      }
      
    }
  }

  //
  // Sum levels
  //
  vector<double> com1(3, 0.0), cov1(3, 0.0), coa1(3, 0.0);
  double         mtot1 = 0.0;

  for (unsigned mm=0; mm<=multistep; mm++) {
    for (int k=0; k<3; k++) {
      com1[k] += com_lev[3*mm + k];
      cov1[k] += cov_lev[3*mm + k];
      coa1[k] += coa_lev[3*mm + k];
    }
    mtot1 += com_mas[mm];
  }

  MPI_Allreduce(&mtot1, &mtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&com1[0], com, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&cov1[0], cov, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&coa1[0], coa, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
  if (VERBOSE>5) {
				// Check for NaN
    bool com_nan = false, cov_nan = false, coa_nan = false;
    for (int k=0; k<3; k++)
      if (std::isnan(com[k])) com_nan = true;
    for (int k=0; k<3; k++)
      if (std::isnan(cov[k])) cov_nan = true;
    for (int k=0; k<3; k++)
      if (std::isnan(coa[k])) coa_nan = true;
    if (coa_nan && myid==0)
      cerr << "Component [" << name << "] com has a NaN" << endl;
    if (cov_nan && myid==0)
      cerr << "Component [" << name << "] cov has a NaN" << endl;
    if (coa_nan && myid==0)
      cerr << "Component [" << name << "] coa has a NaN" << endl;
  }

  if (consp && com_system) {
    
    vector<double> comE(3), covE(3);
    double         mtotE;
    
    mtot1 = 0.0;
    for (int k=0; k<3; k++) com1[k] = cov1[k] = 0.0;

    for (unsigned mm=mlevel; mm<=multistep; mm++) {
      for (int k=0; k<3; k++) {
	com1[k] += comE_lev[3*mm + k];
	cov1[k] += covE_lev[3*mm + k];
      }
      mtot1 += comE_mas[mm];
    }

    MPI_Allreduce(&mtot1,   &mtotE,   1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&com1[0], &comE[0], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&cov1[0], &covE[0], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    for (int i=0; i<3; i++) {
      com0[i] = (mtot0*com0[i] - comE[i])/(mtot0 - mtotE);
      cov0[i] = (mtot0*cov0[i] - covE[i])/(mtot0 - mtotE);
    }
    mtot0 -= mtotE;
  }
				// Compute component center of mass and
				// center of velocity, and center of accel

  if (mtot > 0.0) {
    for (int k=0; k<dim; k++) com[k]  /= mtot;
    for (int k=0; k<dim; k++) cov[k]  /= mtot;
    for (int k=0; k<dim; k++) coa[k]  /= mtot;
  }

  if (com_system and not consp) {
    for (int k=0; k<dim; k++) com0[k] = com[k];
    for (int k=0; k<dim; k++) cov0[k] = cov[k];
  }

  if (com_system) {	   // Use local center of accel for com update
    for (int k=0; k<dim; k++) acc0[k]  = coa[k];
  } else {			// No mass, no acceleration?
    for (int k=0; k<dim; k++) acc0[k]  = 0.0;
  }

  if ((EJ & Orient::CENTER) && !EJdryrun) {
    auto ctr = orient->currentCenter();
    bool ok    = true;
    for (int i=0; i<3; i++) {
      if (std::isnan(ctr[i])) ok = false;
    } 
    if (ok) {
      for (int i=0; i<3; i++) center[i] += ctr[i];
    } else if (myid==0) {
      cout << "Orient: center failure, T=" << tnow 
	   << ", adjustment skipped" << endl;
    }
  }

}


void Component::update_accel(void)
{
  if (myid==0 && com_log) {
				// Open output stream for writing
    ofstream out(comfile.c_str(), ios::out | ios::app);
    if (!out) {
      cerr << "Component: error opening <" << comfile << "> for append\n";
      return;
    }

    out << setw(15) << tnow;
    for (int k=0; k<3; k++) out << setw(15) << com0[k];
    for (int k=0; k<3; k++) out << setw(15) << cov0[k];
    for (int k=0; k<3; k++) out << setw(15) << acc0[k];
    for (int k=0; k<3; k++) out << setw(15) << center[k];
    out << endl;

  }

}


struct thrd_pass_angmom
{
  //! Thread counter id
  int id;

  //! Angular momentum for all levels for this thread
  vector<double> angm1;

  //! Current multistep level
  unsigned mlevel;

  //! Component
  Component *c;
};



void * get_angmom_thread(void *ptr)
{
  //
  // Thread ID
  //
  int id = static_cast<thrd_pass_angmom*>(ptr)->id;
  //
  // Ang mom vector
  //
  double *angm1 = &(static_cast<thrd_pass_angmom*>(ptr)->angm1[0]);
  //
  // Component
  //
  Component *c = static_cast<thrd_pass_angmom*>(ptr)->c;
  //
  // Level
  //
  unsigned mlevel = static_cast<thrd_pass_angmom*>(ptr)->mlevel;


  for (unsigned mm=mlevel; mm<=multistep; mm++) {

    unsigned ntot = c->levlist[mm].size();
    int nbeg = ntot*(id  )/nthrds;
    int nend = ntot*(id+1)/nthrds;
    double mass, *pos, *vel;
  
    //
    // Particle loop
    //
    for (int q=nbeg; q<nend; q++) {
      
      unsigned long n = c->levlist[mm][q];
      Particle     *p = c->Part(n);

      if (c->freeze(n)) continue;

      mass = p->mass;
      pos  = p->pos;
      vel  = p->vel;
    
      angm1[3*mm + 0] += mass*(pos[1]*vel[2] - pos[2]*vel[1]);

      angm1[3*mm + 1] += mass*(pos[2]*vel[0] - pos[0]*vel[2]);
      
      angm1[3*mm + 2] += mass*(pos[0]*vel[1] - pos[1]*vel[0]);
    }
  }
  
  return (NULL);
}


void Component::get_angmom(unsigned mlevel)
{
  
  //
  // Zero variables
  //
  for (unsigned mm=mlevel; mm<=multistep; mm++) {
    for (int i=0; i<3; i++) angmom_lev[3*mm+i] = 0.0;
  }

  //
  // Make the <nthrds> threads
  //
  int errcode;
  void *retval;

  vector<thrd_pass_angmom> data(nthrds);
  vector<pthread_t>        thrd(nthrds);

  if (nthrds==1) {

    data[0].id     = 0;
    data[0].c      = this;
    data[0].mlevel = mlevel;
    data[0].angm1  = vector<double>(3*(multistep+3), 0);
    
    get_angmom_thread(&data[0]);

    for (unsigned mm=mlevel; mm<=multistep; mm++) {
      for (unsigned k=0; k<3; k++) 
	angmom_lev[3*mm + k] += data[0].angm1[3*mm + k];
    }

  } else {

    for (int i=0; i<nthrds; i++) {

      data[i].id     = i;
      data[i].c      = this;
      data[i].mlevel = mlevel;
      data[i].angm1  = vector<double>(3*(multistep+3), 0);

      errcode =  pthread_create(&thrd[i], 0, get_angmom_thread, &data[i]);

      if (errcode) {
	std::ostringstream sout;
	sout << "Process " << myid
	     << " Component::get_angmom: cannot make thread " << i
	     << ", errcode=" << errcode;
	throw GenericError(sout.str(), __FILE__, __LINE__);
      }
    }
    
    //
    // Collapse the threads
    //
    for (int i=0; i<nthrds; i++) {
      if ((errcode=pthread_join(thrd[i], &retval))) {
	std::ostringstream sout;
	sout << "Process " << myid
	     << " Component::get_angmom: thread join " << i
	     << " failed, errcode=" << errcode;
	throw GenericError(sout.str(), __FILE__, __LINE__);
      }
      for (unsigned mm=mlevel; mm<=multistep; mm++) {
	for (unsigned k=0; k<3; k++) 
	  angmom_lev[3*mm + k] += data[i].angm1[3*mm + k];
      }
    }
  }


  //
  // Sum up over all levels
  //
  vector<double> angm1(3, 0);
  for (unsigned mm=0; mm<=multistep; mm++) {
    for (unsigned k=0; k<3; k++) angm1[k] += angmom_lev[3*mm + k];
  }

  MPI_Allreduce(&angm1[0], angmom, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

}



int Component::round_up(double dnumb)
{
  unsigned numb = (unsigned)(dnumb + 1.0);
  return numb;
}


void Component::setup_distribution(void)
{
				// Needed for both master and slaves
  nbodies_index = vector<unsigned int>(numprocs);
  nbodies_table = vector<unsigned int>(numprocs);

  if (myid == 0) {

    orates = vector<double>(numprocs);
    trates = vector<double>(numprocs);

    for (int n=0; n<numprocs; n++) {

      if (n == 0)
	nbodies_table[n] = nbodies_index[n] = 
	  max<int>(1, min<int>((int)(comp->rates[n] * nbodies_tot), nbodies_tot));
      else {
	if (n < numprocs-1)
	  nbodies_index[n] = (int)(comp->rates[n] * nbodies_tot) + 
	    nbodies_index[n-1];
	else
	  nbodies_index[n] = nbodies_tot;
      
	nbodies_table[n] = nbodies_index[n] - nbodies_index[n-1];
      }

    }

    std::string outrates = outdir + "current.processor.rates." + runtag;

    std::ofstream out(outrates, ios::out | ios::app);

    if (out.good()) {
      out << "# " << endl;
      out << "# Time=" << tnow << " Component=" << name << endl;
      out << "# " 
	  << setw(15) << "Norm rate"
	  << setw(15) << "Delta rate"
	  << setw(15) << "Index"
	  << setw(15) << "Current #"
	  << endl
	  << "# "
	  << setw(15) << "---------"
	  << setw(15) << "----------"
	  << setw(15) << "--------"
	  << setw(15) << "---------"
	  << endl;
      
      for (int n=0; n<numprocs; n++)
	out << "  "
	    << setw(15) << comp->rates[n]
	    << setw(15) << 1.0 - comp->rates[n]*nbodies_tot/nbodies_table[n]
	    << setw(15) << nbodies_index[n]
	    << setw(15) << nbodies_table[n]
	    << endl;

      out.close();
    }

  }


  MPI_Bcast(&nbodies_index[0], numprocs, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nbodies_table[0], numprocs, MPI_INT, 0, MPI_COMM_WORLD);

}

void Component::update_indices(void)
{
  // Make sure arrays have correct size
  //
  nbodies_index.resize(numprocs);
  nbodies_table.resize(numprocs);

  // Gather current size from all processes
  //
  nbodies = particles.size();
  MPI_Allgather(&nbodies, 1, MPI_UNSIGNED, nbodies_table.data(), 1, MPI_UNSIGNED,
		MPI_COMM_WORLD);

  // Cumulate
  //
  nbodies_index[0] = nbodies_table[0];
  for (int n=0; n<numprocs; n++)
    nbodies_index[n] = nbodies_index[n-1] + nbodies_table[n];

}

void Component::load_balance(void)
{
  MPI_Status status;
  vector<unsigned int> nbodies_index1(numprocs);
  vector<unsigned int> nbodies_table1(numprocs);
  std::ofstream out, log;

  update_indices();		// Refresh particle counts

  if (myid == 0) {

    std::vector<double> orates1(numprocs);
    std::vector<double> trates1(numprocs);

    for (int n=0; n<numprocs; n++) {

      if (n == 0)
	nbodies_table1[n] = nbodies_index1[n] = 
	  std::max<int>(1, min<int>((int)(comp->rates[n] * nbodies_tot), nbodies_tot));
      else {
	if (n < numprocs-1)
	  nbodies_index1[n] = (int)(comp->rates[n] * nbodies_tot) +  nbodies_index1[n-1];
	else
	  nbodies_index1[n] = nbodies_tot;
      
	nbodies_table1[n] = nbodies_index1[n] - nbodies_index1[n-1];
      }

    }

    std::string outrates =
      outdir + "current.processor.rates." + name + "." + runtag;

    std::string rateslog =
      outdir + "current.processor.rates.log." + name + "." + runtag;

    out.open(outrates, ios::out | ios::app);
    log.open(rateslog, ios::out | ios::app);

    if (not out.good()) {
      std::cout << "*** ERROR: Component::load_balance error opening <"
		<< outrates << ">" << std::endl;
    }

    if (not log.good()) {
      std::cout << "*** ERROR: Component::load_balance error opening <"
		<< rateslog << ">" << std::endl;
    }

    if (out) {
      out << "# " << endl;
      out << "# Time=" << tnow << " Component=" << name << endl;
      out << "# " 
	  << setw(15) << "Norm rate"
	  << setw(15) << "Delta rate"
	  << setw(15) << "Index"
	  << setw(15) << "Current #"
	  << setw(15) << "Old Index"
	  << setw(15) << "Previous #"
	  << endl
	  << "# "
	  << setw(15) << "--------"
	  << setw(15) << "----------"
	  << setw(15) << "--------"
	  << setw(15) << "---------"
	  << setw(15) << "---------"
	  << setw(15) << "---------"
	  << endl;
      
      for (int n=0; n<numprocs; n++)
	out << "  "
	    << setw(15) << comp->rates[n]
	    << setw(15) << 1.0 - comp->rates[n]*nbodies_tot/nbodies_table1[n]
	    << setw(15) << nbodies_index1[n]
	    << setw(15) << nbodies_table1[n]
	    << setw(15) << nbodies_index[n]
	    << setw(15) << nbodies_table[n]
	    << endl;
    }

  }

  MPI_Bcast(&nbodies_index1[0], numprocs, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nbodies_table1[0], numprocs, MPI_INT, 0, MPI_COMM_WORLD);

				// Compute index
  loadb.clear();
  loadb_datum datum0, datum1;
  datum0.s = 0;
  datum1.s = 1;
  for (int i=0; i<numprocs; i++) {
    datum0.top = nbodies_index[i];
    datum1.top = nbodies_index1[i];
    datum0.indx = datum1.indx = i;
    loadb.push_back(datum0);
    loadb.push_back(datum1);
  }

  sort(loadb.begin(), loadb.end(), less_loadb);

  if (myid==0 && log.good()) 
    {
      log << std::setw(72) << std::setfill('.') << ".\n" << std::setfill(' ');
      log << "Time=" << tnow << " Component=" << name << std::endl;
      log << endl;
      log << "List:\n";
      log.setf(ios::left);
      log << setw(4) << "N"
	  << setw(6) << "Index"
	  << setw(10) << "Old"
	  << setw(10) << "New"
	  << endl;

      char c = log.fill('-');
      log << setw(4) << "|"
	  << setw(6) << "|"
	  << setw(10) << "|"
	  << setw(10) << "|"
	  << endl;
      log.fill(c);
      
      for (int i=0; i<2*numprocs; i++) {
	
	log << setw(4) << i
	    << setw(6) << loadb[i].indx;
	if (loadb[i].s)
	  log << setw(10) << " " << setw(10) << loadb[i].top;
	else 
	  log << setw(10) << loadb[i].top;
	log << endl;
      }

    }


  if (myid==0 && log.good()) 
    {
      log << "\nAnalysis:\n";
      log.setf(ios::left);
      log << setw(10) << "Interval"
	  << setw(10) << "Number"
	  << setw(10) << "Old"
	  << setw(10) << "New"
	  << setw(10) << "Action"
	  << endl;
      
      char c = log.fill('-');
      log << setw(10) << "|"
	  << setw(10) << "|"
	  << setw(10) << "|"
	  << setw(10) << "|"
	  << setw(10) << "|"
	  << endl;
      log.fill(c);
    }


  int iold=0, inew=0;
  
				// Offset will keep track of position in
				// original vector
  int nump;
  std::vector<int> loc(numprocs, 0);

  std::vector<PartPtr> nlist;

  for (int i=0; i<2*numprocs-2; i++) {

				// Assign new interval
    if (loadb[i].s) inew = loadb[i].indx+1;
    else            iold = loadb[i].indx+1;
    
    if (myid==0 && log.good())
      log << setw(10) << i
	  << setw(10) << loadb[i+1].top - loadb[i].top
	  << setw(10) << iold
	  << setw(10) << inew;
    
				// Number of particles to be shifted
    nump = loadb[i+1].top - loadb[i].top;
    
    std::ostringstream msg;
    
    if (inew==iold || nump==0) 
      msg << "Do nothing";
    else if (inew>iold) {
      msg << "Add " << nump << " from #" << iold << " to #" << inew;
      
      nlist.clear();

      PartMap::iterator it = particles.begin();
      for (int n=0; n<nump; n++) {
	nlist.push_back(it->second);
	it++;
      }
      
      add_particles(iold, inew, nlist);
      
    } else if (iold>inew) {
      msg << "Add " << nump << " from #" << iold << " to #" << inew;

      nlist.clear();

      PartMapItr it = particles.begin();
      for (int n=0; n<nump; n++) {
	nlist.push_back(it->second);
	it++;
      }

      add_particles(iold, inew, nlist);

    }

    if (myid==0 && log.good()) log << setw(10) << msg.str() << endl;
  }

  
				// update indices
  nbodies = nbodies_table1[myid];
  nbodies_index = nbodies_index1;
  nbodies_table = nbodies_table1;
  
  if (myid==0) {
    try {
      out.close();
      log.close();
    }
    catch (const ofstream::failure& e) {
      std::cout << "Component: exception closing out and log files: "
		<< e.what() << std::endl;
    }
  }

}


template< typename T >
typename std::vector<T>::iterator 
insert_sorted( std::vector<T> & vec, T const& item )
{
  return vec.insert
    ( 
     std::upper_bound( vec.begin(), vec.end(), item ),
     item 
      );
}

void Component::add_particles(int from, int to, std::vector<PartPtr>& plist)
{
  unsigned number = plist.size();
  std::vector<PartPtr>::iterator it=plist.begin();

  unsigned icount, counter=0;

  pf->ShipParticles(to, from, number);

  if (myid == from) {
    
    while (counter < number) {

      icount = 0;
      while (icount < PFbufsz && counter < number) {

	pf->SendParticle(*it);

	// Remove particle from lev list
	//
	bool success = false;	// Sanity check
	for (auto & v : levlist) {
	  auto jt = std::find(v.begin(), v.end(), (*it)->indx);
	  if (jt != v.end()) {
	    v.erase(jt);
	    success = true;
	    break;
	  }
	}

	// Levlist sanity check
	//
	if (not success) {
	  std::cout << "***ERROR*** "
		    << "Component::add_particles: could not find indx="
		    << (*it)->indx << " in levlist in any of "
		    << multistep+1 << " levels" << std::endl;
	}

	particles.erase((*it)->indx);
      
	icount++;
	counter++;
	it++;
      }
      
#ifdef DEBUG
      cout << "Process " << myid 
	   << ": sent " << icount << " particles to Process " << to
	   << " for append, counter value=" << counter
	   << endl << flush;
#endif    
    }

  }

  if (myid == to) {
  
    while (counter < number) {

      while (PartPtr temp=pf->RecvParticle()) {
	particles[temp->indx] = temp;
	insert_sorted<int>(levlist[temp->level], temp->indx);
	counter++;
      }

#ifdef DEBUG
      cout << "Process " << myid 
	   << ": received " << icount << " particles from Process " << from
	   << " for append" << endl << flush;
#endif    

    }

  }

}


bool Component::freeze(unsigned indx)
{
  double r2 = 0.0;
  for (int i=0; i<3; i++) r2 += 
			    (particles[indx]->pos[i] - com0[i] - center[i])*
			    (particles[indx]->pos[i] - com0[i] - center[i]);
  if (r2 > rtrunc*rtrunc) return true;
  else return false;
}

bool Component::escape_com(const Particle& p)
{
  double r2 = 0.0;
  for (int i=0; i<3; i++) r2 += 
			    (p.pos[i] - com0[i] - center[i])*
			    (p.pos[i] - com0[i] - center[i]);
  if (r2 > rcom*rcom) return true;
  else return false;
}

double Component::Adiabatic()
{
  if (!adiabatic) return 1.0;
  return 0.25*
    ( 1.0 + erf((tnow - ton )/twid) ) *
    ( 1.0 + erf((toff - tnow)/twid) ) ;
}


void Component::redistributeByList(vector<int>& redist)
{
  // Initialize the particle ferry instance with dynamic attribute sizes
  if (not pf) pf = ParticleFerryPtr(new ParticleFerry(niattrib, ndattrib));


  vector<int>::iterator it = redist.begin();
  vector<unsigned> tlist;

  PartPtr part;
  unsigned int icount;
  int indx, curnode, tonode, lastnode, M;

  while (it != redist.end()) {
    curnode = *(it++);		// Current owner
    M       = *(it++);		// Number to transfer to another node
    if (M) {
      indx   = *(it++);		// Index
      tonode = *(it++);		// Destination
      icount = 0;		// Number transferred to destination


      // Do the first particle
      //
      tlist.clear();
      tlist.push_back(indx);
      icount++;

      lastnode = tonode;

      // Do the remaining particles
      //
      for (int m=1; m<M; m++) {
	indx   = *(it++);
	tonode = *(it++);
				// Next destination?
	if (tonode != lastnode && icount) {
	  pf->ShipParticles(tonode, curnode, icount);

	  if (myid==curnode) {
	    for (unsigned i=0; i<icount; i++) {
	      pf->SendParticle(particles[tlist[i]]);
	      particles.erase(tlist[i]);
	    }
	  }
	  if (myid==lastnode) {
	    while (part=pf->RecvParticle())
	      particles[part->indx] = part;
	  }
	  tlist.clear();
	  icount = 0;
	}

				// Add the particle
	tlist.push_back(indx);
	icount++;
      }
	
    } // End of particles on this node
    
  } // Next stanza

}


Particle* Component::GetNewPart()
{
  // Create new particle
  //
  PartPtr newp = std::make_shared<Particle>(niattrib, ndattrib);

  // Denote unsequenced particle
  //
  newp->indx = -1;
  
  // Add to new particle list
  //
  new_particles.push_back(newp);

  modified++;

  return newp.get();
}

void Component::seq_new_particles()
{
  // Update total number of bodies
  //
  nbodies = particles.size();
  MPI_Allreduce(&nbodies, &nbodies_tot, 1, MPI_UNSIGNED, MPI_SUM,
		MPI_COMM_WORLD);

  // Are there new particles to sequence?
  //
  unsigned newTot, newCur = new_particles.size();
  MPI_Allreduce(&newCur, &newTot, 1, MPI_UNSIGNED, MPI_SUM,
		MPI_COMM_WORLD);

  if (newTot==0) {
    modified = 0;
    return;
  }

  // Begin sequence numbering
  //
  for (int n=0; n<numprocs; n++) {
    // Run through the unsequenced loop
    //
    if (myid==n) {
      for (auto p : new_particles) {
	p->indx  = ++top_seq;
	p->level = multistep;
	particles[p->indx] = p;
				// Add to level list:
				// last element by construction
	levlist[p->level].push_back(p->indx); 
      }
    }

    // Share top_seq with all processes
    //
    MPI_Bcast(&top_seq, 1, MPI_UNSIGNED_LONG, n, MPI_COMM_WORLD);
  }

  // Erase the new particle list
  //
  new_particles.clear();

  // Update total number of bodies
  //
  nbodies = particles.size();
  MPI_Allreduce(&nbodies, &nbodies_tot, 1, MPI_UNSIGNED, MPI_SUM,
		MPI_COMM_WORLD);

#ifdef DEBUG
  std::cout << "seq_new_particles [Node " << myid << "] # part="
	    << particles.size() << ", nbodies=" << nbodies << std::endl;
#endif

  // Reset the change counter
  //
  modified = 0;
}


void Component::DestroyPart(PartPtr p)
{
  particles.erase(p->indx);

  // Remove from level list
  //
  bool success = false;		// For sanity check . . .
  for (auto & v : levlist) {
    auto it = std::find(v.begin(), v.end(), p->indx);
    if (it != v.end()) {
      v.erase(it);
      success = true;
      break;
    }
  }

  // Levlist sanity check
  //
  if (not success) {
    std::cout << "***ERROR*** "
	      << "Component::DestroyPart: could not find indx=" << p->indx
	      << " in levlist in any of " << multistep+1 << " levels"
	      << std::endl;
  }

  // Update counters
  //
  nbodies--;
  modified++;
}
