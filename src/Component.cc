#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
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

#include "expand.h"

// For sort algorithm below
bool less_loadb(const loadb_datum& one, const loadb_datum& two)
{
  return (one.top < two.top);
}

// Constructor
Component::Component(string NAME, string ID, string CPARAM, string PFILE, 
		     string FPARAM) : 
  name(NAME), id(ID), cparam(CPARAM), pfile(PFILE), fparam(FPARAM)
{
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
  tidal       = 0;

  com_system  = false;
  com_log     = false;

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
  comI        = 0;
  covI        = 0;

  seq_check   = false;
  indexing    = false;
  aindex      = false;
  umagic      = true;

  nlevel      = -1;

  read_bodies_and_distribute_ascii();

  mdt_ctr = vector< vector<unsigned> > (multistep+1);
  for (unsigned n=0; n<=multistep; n++) mdt_ctr[n] = vector<unsigned>(mdtDim, 0);

  angmom_lev  = vector<double>(3*(multistep+1), 0);
  com_lev     = vector<double>(3*(multistep+1), 0);
  cov_lev     = vector<double>(3*(multistep+1), 0);
  coa_lev     = vector<double>(3*(multistep+1), 0);
  com_mas     = vector<double>(multistep+1, 0);

  reset_level_lists();

  time_so_far.Microseconds();

  tree = 0;

  pbuf = new Particle [PFbufsz];

  // Already done in read_bodies_and_distribut_ascii()
  // initialize();
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

  for (int n=0; n<nbeg; n++) it++;
  for (int n=nbeg; n<nend; n++) {
    (static_cast<thrd_pass_reset*>(ptr)->newlist)[it->second.level].push_back(it->first);
    it++;
  }
  
  return (NULL);
}

void Component::reset_level_lists()
{
  if (td.size()==0) {
    td = vector<thrd_pass_reset>(nthrds);

    t  = new pthread_t [nthrds];

    if (!t) {
      cerr << "Process " << myid
	   << ": reset_level_lists: error allocating memory for thread\n";
      exit(18);
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
	cerr << "Process " << myid
	     << " reset_level_lists: cannot make thread " << i
	     << ", errcode=" << errcode << endl;
	exit(19);
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
	cerr << "Process " << myid
	     << " reset_level_lists: thread join " << i
	     << " failed, errcode=" << errcode << endl;
	exit(20);
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

    vector< vector<unsigned> > cntr(multistep+1);
    for (unsigned n=0; n<=multistep; n++) {
      cntr[n] = vector<unsigned>(mdtDim, 0);
      MPI_Reduce(&mdt_ctr[n][0], &cntr[n][0], mdtDim, MPI_UNSIGNED,
		 MPI_SUM, 0, MPI_COMM_WORLD);
      
      for (int k=0; k<mdtDim; k++) mdt_ctr[n][k] = 0;
    }
    
    if (myid==0) {

      unsigned tot=0;
      for (unsigned n=0; n<=multistep; n++) tot += cntr[n][mdtDim-1];

      if (!tot && myid==0) cout << "print_level_lists [" << name 
				<< ", T=" << tnow << "]: tot=" << tot << endl;

      if (tot) {

	ostringstream ofil;
	ofil << runtag << ".levels";
	ofstream out(ofil.str().c_str(), ios::app);

	unsigned curn, sum=0;
	out << setw(80) << setfill('-') << '-' << endl;
	ostringstream sout;
	sout << "--- Component <" << name 
	     << ", " << id  << ">, T=" << T;
	out << setw(80) << left << sout.str().c_str() << endl;
	out << setw(80) << '-' << endl << setfill(' ');
	out << setw(3)  << "L" 
	    << setw(10) << "Number" 
	    << setw(10) << "dN/dL" 
	    << setw(10) << "N(<=L)";
	if (DTold) {
	  out << setw(10) << "f(r/v)"
	      << setw(10) << "f(s/v)"
	      << setw(10) << "f(v/a)"
	      << setw(10) << "f(r/a)";
	} else {
	  out << setw(10) << "f(v/a)"
	      << setw(10) << "f(s/v)"
	      << setw(10) << "f(p/v*a)" 
	      << setw(10) << "f([r/a])";
	}
	out << setw(10) << "f(int)" << endl;
	out << setw(80) << setfill('-') << '-' << endl << setfill(' ');
	for (unsigned n=0; n<=multistep; n++) {
	  curn = cntr[n][mdtDim-1];
	  sum += curn;
	  out << setw(3)  << n 
	      << setw(10) << curn << setprecision(3) << fixed
	      << setw(10) << static_cast<double>(curn)/tot
	      << setw(10) << static_cast<double>(sum) /tot;
	  if (curn)		// If there are counts at this level:
	    out << setw(10) << static_cast<double>(cntr[n][0])/curn
		<< setw(10) << static_cast<double>(cntr[n][1])/curn
		<< setw(10) << static_cast<double>(cntr[n][2])/curn
		<< setw(10) << static_cast<double>(cntr[n][3])/curn
		<< setw(10) << static_cast<double>(cntr[n][4])/curn;
	  else			// No counts at this level:
	    out << setw(10) << "*"
		<< setw(10) << "*"
		<< setw(10) << "*"
		<< setw(10) << "*"
		<< setw(10) << "*";
	  out << endl;
	}
	out << endl << setw(3) << "T" << setw(10) << tot << endl << endl 
	    << right;
      }
    }
  }
}

Component::Component(istream *in)
{
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
  tidal       = 0;

  com_system  = false;
  com_log     = false;
  com_restart = 0;

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
  comI        = 0;
  covI        = 0;

  seq_check   = false;
  indexing    = false;
  aindex      = false;
  umagic      = true;

  nlevel      = -1;

  read_bodies_and_distribute_binary(in);

  mdt_ctr = vector< vector<unsigned> > (multistep+1);
  for (unsigned n=0; n<=multistep; n++) mdt_ctr[n] = vector<unsigned>(mdtDim, 0);

  angmom_lev  = vector<double>(3*(multistep+1), 0);
  com_lev     = vector<double>(3*(multistep+1), 0);
  cov_lev     = vector<double>(3*(multistep+1), 0);
  coa_lev     = vector<double>(3*(multistep+1), 0);
  com_mas     = vector<double>(multistep+1, 0);

  reset_level_lists();

  tree = 0;

  pbuf = new Particle [PFbufsz];
}


void Component::initialize(void)
{
				// Parse the parameters
  StringTok<string> tokens(cparam);
  pair<string, string> datum;

  string token = tokens(",");	// Comma separated tokens

  while (token.size()) {
    StringTok<string> parse(token);
    datum.first  = trimLeft(trimRight(parse("=")));
    datum.second = trimLeft(trimRight(parse("=")));

    if (!datum.first.compare("com"))      com_system = atoi(datum.second) ? true : false;

    if (!datum.first.compare("comlog"))   com_log = atoi(datum.second) ? true : false;

    if (!datum.first.compare("tidal"))    {tidal = atoi(datum.second); consp=true;}

    if (!datum.first.compare("EJ"))       EJ = atoi(datum.second.c_str());
    
    if (!datum.first.compare("eEJ0"))     {if (myid==0) cout << "Component: eEJ0 is no longer used, Ecurr is computed from the bodies using the expansion directly" << endl;}

    if (!datum.first.compare("nEJkeep"))  nEJkeep = atoi(datum.second.c_str());

    if (!datum.first.compare("nEJwant"))  nEJwant = atoi(datum.second.c_str());

    if (!datum.first.compare("EJx0"))     EJx0 = atof(datum.second.c_str());

    if (!datum.first.compare("EJy0"))     EJy0 = atof(datum.second.c_str());

    if (!datum.first.compare("EJz0"))     EJz0 = atof(datum.second.c_str());

    if (!datum.first.compare("EJu0"))     EJu0 = atof(datum.second.c_str());

    if (!datum.first.compare("EJv0"))     EJv0 = atof(datum.second.c_str());

    if (!datum.first.compare("EJw0"))     EJw0 = atof(datum.second.c_str());

    if (!datum.first.compare("EJdT"))     EJdT = atof(datum.second.c_str());

    if (!datum.first.compare("EJkinE"))   EJkinE = atoi(datum.second.c_str()) ? true : false;

    if (!datum.first.compare("EJext"))    EJext =  atoi(datum.second.c_str()) ? true : false;

    if (!datum.first.compare("EJdiag"))   EJdiag = atoi(datum.second.c_str()) ? true : false;

    if (!datum.first.compare("EJdryrun")) EJdryrun = atoi(datum.second.c_str()) ? true : false;

    if (!datum.first.compare("EJlinear")) EJlinear = atoi(datum.second.c_str()) ? true : false;

    if (!datum.first.compare("EJdamp"))   EJdamp = atof(datum.second.c_str());

    if (!datum.first.compare("rmax"))     rmax = atof(datum.second.c_str());

    if (!datum.first.compare("ton"))      {ton = atof(datum.second.c_str()); adiabatic = true;}

    if (!datum.first.compare("toff"))     {toff= atof(datum.second.c_str()); adiabatic = true;}

    if (!datum.first.compare("twid"))     {twid = atof(datum.second.c_str()); adiabatic = true;}

    if (!datum.first.compare("rtrunc"))   {rtrunc = atof(datum.second.c_str());}

    if (!datum.first.compare("rcom"))     {rcom = atof(datum.second.c_str());}
    
    if (!datum.first.compare("scheck"))   seq_check = atoi(datum.second.c_str()) ? true : false;

    if (!datum.first.compare("magic"))    umagic = atoi(datum.second.c_str()) ? true : false;

    if (!datum.first.compare("indexing")) indexing = atoi(datum.second.c_str()) ? true : false;

    if (!datum.first.compare("aindex"))   aindex = atoi(datum.second.c_str()) ? true : false;

    if (!datum.first.compare("nlevel"))   nlevel = atoi(datum.second.c_str());



				// Next parameter
    token = tokens(",");
  }


				// Instantiate the force ("reflection" by hand)

  if ( !id.compare("bessel") ) {
    force = new Bessel(fparam);
  }
  else if ( !id.compare("c_brock") ) {
    force = new CBrock(fparam);
  }
  else if ( !id.compare("c_brock_disk") ) {
    force = new CBrockDisk(fparam);
  }
  else if ( !id.compare("hernq") ) {
    force = new Hernquist(fparam);
  }
  else if ( !id.compare("sphereSL") ) {
    force = new Sphere(fparam);
  }
  else if ( !id.compare("EJcom") ) {
    force = new EJcom(fparam);
  }
  else if ( !id.compare("cube") ) {
    force = new Cube(fparam);
  }
  else if ( !id.compare("slab") ) {
    force = new Slab(fparam);
  }
  else if ( !id.compare("slabSL") ) {
    force = new SlabSL(fparam);
  }
  else if ( !id.compare("cylinder") ) {
    force = new Cylinder(fparam);
  }
  else if ( !id.compare("direct") ) {
    force = new Direct(fparam);
  }
  else if ( !id.compare("shells") ) {
    force = new Shells(fparam);
  }
  else if ( !id.compare("noforce") ) {
    force = new NoForce(fparam);
  }
  else {
    string msg("I don't know about the force: ");
    msg += id;
    bomb(msg);
  }

  dim = force->dof;

  force->RegisterComponent(this);

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
  comI   = new double[3];
  covI   = new double[3];

  for (int k=0; k<3; k++) 
    com[k] = center[k] = cov[k] = coa[k] = 
      com0[k] = cov0[k] = acc0[k] = 
      comI[k] = covI[k] = angmom[k] = 0.0;
  

  if (com_system) {

    if (consp) {
      comE_lev = vector<double>(3*(multistep+1), 0);
      covE_lev = vector<double>(3*(multistep+1), 0);
      comE_mas = vector<double>(multistep+1, 0);
    }

    initialize_com_system();

    if (myid==0) {
      cout << " Component <" <<  name 
	   << ">: center of mass system is *ON*, rtrunc=" << rtrunc;
      if (consp) cout << ", conserving com momentum [iattr #=" << tidal << "]";
      cout << ", computed COM system:";
      cout << endl << "\t\t(x, y, z)=("
	   << setw(15) << comI[0] << ", "
	   << setw(15) << comI[1] << ", "
	   << setw(15) << comI[2] << ") "
	   << endl << "\t\t"
	   << "(u, v, w)=("
	   << setw(15) << covI[0] << ", "
	   << setw(15) << covI[1] << ", "
	   << setw(15) << covI[2] << ") "
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
		  for (int k=0; k<3; k++) istr >> comI[k];
		  for (int k=0; k<3; k++) istr >> covI[k];
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

		  cout << "\t\tInitial com at T=" << ttim0 << " is:";
		  cout << endl << "\t\t(x, y, z)=("
		       << setw(15) << comI[0] << ", "
		       << setw(15) << comI[1] << ", "
		       << setw(15) << comI[2] << ") "
		       << endl << "\t\t"
		       << "(u, v, w)=("
		       << setw(15) << covI[0] << ", "
		       << setw(15) << covI[1] << ", "
		       << setw(15) << covI[2] << ") "
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
	cout << "Component <" << name << ">: EJ center finding is *ON*";
	cout << " with damping=" << EJdamp;
	if (EJkinE)   cout << ", using particle kinetic energy";
	if (EJext)    cout << ", using external potential";
	if (EJdryrun) cout << ", dryrun";
	cout << endl;
      }
      if (EJ & Orient::AXIS) {
	cout << "Component <" << name << ">: AXIS orientation is *ON*";
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
      for (int i=0; i<3; i++) center[i] = (orient->currentCenter())[i+1];
    } else {
      orient -> set_center(EJx0, EJy0, EJz0);
      orient -> set_cenvel(EJu0, EJv0, EJw0);
      if (EJlinear) orient -> set_linear();
      center[0] = EJx0;
      center[1] = EJy0;
      center[2] = EJz0;
    }

    if (EJdiag) cout << "Process " << myid << ": Orient successful\n";
  }

  if (myid == 0) {		// Center status
    cout << "Component <" << name;
    if (restart)
      cout << ">: current center on restart: x, y, z: " 
	   << center[0] << ", " 
	   << center[1] << ", " 
	   << center[2];
    else
      cout << ">: user specified initial center: x, y, z: " 
	   << EJx0 << ", " 
	   << EJy0 << ", " 
	   << EJz0;

    cout << "Component <" << name << ">: ";

    if (nlevel<0)
      cout << "no multistep level reporting";
    else
      cout << "multistep level reporting every " << nlevel << " steps";

    cout << endl << endl;
  }
  
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
  delete [] comI;
  delete [] covI;

  delete [] pbuf;

  delete tree;
}

void Component::bomb(const string& msg)
{
  cerr << "Component <" << name << ", " << id << ">: " << msg << endl;
  exit(-1);
}

void Component::read_bodies_and_distribute_ascii(void)
{
				// Open file

  ifstream *fin;
  const int nline = 2048;
  char line[nline];
  
  if (myid == 0) {
    fin = new ifstream(pfile.c_str());

    if (!*fin) {
      cerr << "Couldn't open " << pfile << " . . . quitting\n";
      MPI_Abort(MPI_COMM_WORLD, 11);
      exit(-1);
    }

    fin->getline(line, nline);
    istringstream ins(line);
    
    ins >> nbodies_tot;		
    if (!ins) {
      cerr << "Error reading nbodies_tot . . . quitting\n";
      MPI_Abort(MPI_COMM_WORLD, 11);
      exit(-1);
    }
    ins >> niattrib;
    if (!ins) {
      cerr << "Error reading integer attribute # . . . quitting\n";
      MPI_Abort(MPI_COMM_WORLD, 11);
      exit(-1);
    }
    ins >> ndattrib;
    if (!ins) {
      cerr << "Error reading double attribute # . . . quitting\n";
      MPI_Abort(MPI_COMM_WORLD, 11);
      exit(-1);
    }

				// Sanity check
    if (niattrib > nimax) {
      cerr << "Component: niattrib=" << niattrib << " >  nimax="
	   << nimax << ".  Change constant in Component.H and recompile\n";
      MPI_Abort(MPI_COMM_WORLD, 11);
      exit(-1);
    }

    if (ndattrib > ndmax) {
      cerr << "Component: ndattrib=" << ndattrib << " >  ndmax="
	   << ndmax << ".  Change constant in Component.H and recompile\n";
      MPI_Abort(MPI_COMM_WORLD, 11);
      exit(-1);
    }

  }
				// Broadcast attributes for this
				// phase-space component
  MPI_Bcast(&nbodies_tot, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&niattrib,    1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ndattrib,    1, MPI_INT, 0, MPI_COMM_WORLD);

  double rmax1=0.0, r2;

  if (nbodies_tot > nbodmax*numprocs) {
    if (myid==0) {
      cerr << "Not enough space on all processors to hold phase space\n";
      cerr << "nbodmax is currently " << nbodmax*numprocs
	   << " but should be at least "
	   << (int)( (double)nbodies_tot/numprocs + 1) << endl;
    }
    MPI_Finalize();
    exit(-1);
  }

  is_init = 1;
  setup_distribution();
  is_init = 0;
  initialize();

  Particle part(niattrib, ndattrib);

  if (myid==0) {
				// Read in Node 0's particles
    for (unsigned i=1; i<=nbodies_table[0]; i++) {

      part.readAscii(aindex, i, fin);
				// Get the radius
      double r2 = 0.0;
      for (int j=0; j<3; j++) r2 += part.pos[j]*part.pos[j];
      rmax1 = max<double>(r2, rmax1);
      
				// Load the particle
      particles[part.indx] = part;
    }

    nbodies = nbodies_table[0];

    unsigned icount, ibufcount;
    for (int n=1; n<numprocs; n++) {

      pf.ShipParticles(n, 0, nbodies_table[n]);

      icount = 0;
      ibufcount = 0;
      while (icount < nbodies_table[n]) {

	int i = nbodies_index[n-1] + 1 + icount;
	part.readAscii(aindex, i, fin);

	r2 = 0.0;
	for (int k=0; k<3; k++) r2 += part.pos[k]*part.pos[k];
	rmax1 = max<double>(r2, rmax1);

	pf.SendParticle(part);
	icount++;

      }

    }

  } else {

    pf.ShipParticles(myid, 0, nbodies);
      
#ifdef DEBUG
    int icount = 0;
#endif
    while (pf.RecvParticle(part)) {
      particles[part.indx] = part;
#ifdef DEBUG
      if (icount<5) {
	cout << "Process " << myid << ": received ";
	cout << setw(14) << part.mass;
	for (int k=0; k<3; k++) cout << setw(14) << part.pos[k];
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

				// COM HERE?
  if (myid==0) delete fin;

#ifdef DEBUG
  if (particles.size()) {
    cout << "read_bodies_and_distribute_ascii: process " << myid 
	 << " name=" << name << " bodies ["
	 << particles.begin() ->second.indx << ", "
	 << particles.rbegin()->second.indx << "], ["
	 << particles.begin() ->first << ", "
	 << particles.rbegin()->first << "]"
	 << " #=" << particles.size() << endl;
  } else {
    cout << "read_bodies_and_distribute_ascii: process " << myid 
	 << " name=" << name
	 << " #=" << particles.size() << endl;
  }
#endif
}

void Component::read_bodies_and_distribute_binary(istream *in)
{
				// Get component header
  ComponentHeader header;
				// Node local parameter buffer
  int ninfochar;
  boost::shared_array<char> info;
  
  if (myid == 0) {

    rsize = sizeof(double);

    if (umagic) {
      unsigned long cmagic;
      in->read((char*)&cmagic, sizeof(unsigned long));
      if ( (cmagic & nmask) != magic ) {
	cerr << "Error identifying new PSP.  Is this an old PSP?\n";
	MPI_Abort(MPI_COMM_WORLD, 11);
	exit(-1);
      }
      rsize = cmagic & mmask;
    }

    if(!header.read(in)) {
      cerr << "Error reading component header\n";
      MPI_Abort(MPI_COMM_WORLD, 12);
      exit(-1);
    }

    nbodies_tot = header.nbod;
    niattrib    = header.niatr;
    ndattrib    = header.ndatr;
    ninfochar   = header.ninfochar;

    info = boost::shared_array<char>(new char [ninfochar+1]);
    memcpy(info.get(), header.info.get(), ninfochar);
    info[ninfochar] = '\0';
  }

  if (umagic)
    MPI_Bcast(&rsize, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

				// Broadcast attributes for this
				// phase-space component
  MPI_Bcast(&nbodies_tot, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&niattrib,    1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ndattrib,    1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ninfochar,   1, MPI_INT, 0, MPI_COMM_WORLD);
  if (myid) info = boost::shared_array<char>(new char [ninfochar+1]);
  MPI_Bcast(info.get(), ninfochar, MPI_CHAR, 0, MPI_COMM_WORLD);


  if (myid==0) {
				// Assign particle structure buffer
				// Sanity check
    if (niattrib > nimax) {
      cerr << "Component: niattrib=" << niattrib << " >  nimax="
	   << nimax << ".  Change constant in Component.H and recompile\n";
      MPI_Abort(MPI_COMM_WORLD, 11);
      exit(-1);
    }

    if (ndattrib > ndmax) {
      cerr << "Component: ndattrib=" << ndattrib << " >  ndmax="
	   << ndmax << ".  Change constant in Component.H and recompile\n";
      MPI_Abort(MPI_COMM_WORLD, 11);
      exit(-1);
    }

  }

				// Parse info field to get 
				// id and parameter strings
  StringTok<string> tokens(info.get());
  name   = trimLeft(trimRight(tokens(":")));
  id     = trimLeft(trimRight(tokens(":")));
  cparam = trimLeft(trimRight(tokens(":")));
  fparam = trimLeft(trimRight(tokens(":")));

				// Informational output
  if (myid==0)
    cout << setw(60) << setfill('-') << "-" << endl << setfill(' ')
	 << "--- New Component" << endl
	 << setw(20) << " name   :: " << name           << endl
	 << setw(20) << " id     :: " << id             << endl
	 << setw(20) << " cparam :: " << cparam         << endl
	 << setw(20) << " fparam :: " << fparam         << endl
	 << setw(60) << setfill('-') << "-" << endl << setfill(' ');

  double rmax1=0.0, r2;

  if (nbodies_tot > nbodmax*numprocs) {
    if (myid==0) {
      cerr << "Not enough space on all processors to hold phase space\n";
      cerr << "nbodmax is currently " << nbodmax*numprocs
	   << " but should be at least "
	   << (int)( (double)nbodies_tot/numprocs + 1) << endl;
    }
    MPI_Finalize();
    exit(-1);
  }

  is_init = 1;
  setup_distribution();
  is_init = 0;
  initialize();

				// Form cumulative and differential bodies list
  
  Particle part(niattrib, ndattrib);

  unsigned int ipart=0;

  if (myid==0) {
				// Read root node particles
    seq_cur = 0;

    rmax1 = 0.0;
    for (unsigned i=1; i<=nbodies_table[0]; i++)
    {
      part.readBinary(rsize, indexing, ++seq_cur, in);

      r2 = 0.0;
      for (int j=0; j<3; j++) r2 += part.pos[j]*part.pos[j];
      rmax1 = max<double>(r2, rmax1);

				// Load the particle
      particles[part.indx] = part;
    }

    nbodies = nbodies_table[0];


				// Now load the other nodes
    unsigned icount;
    for (int n=1; n<numprocs; n++) {

      cout << "Component [" << name << "]: loading node <" << n << ">\n";

      pf.ShipParticles(n, 0, nbodies_table[n]);

      icount = 0;
      while (icount < nbodies_table[n]) {

	part.readBinary(rsize, indexing, ++seq_cur, in);

	r2 = 0.0;
	for (int k=0; k<3; k++) 
	  r2 += part.pos[k]*part.pos[k];

	rmax1 = max<double>(r2, rmax1);

	icount++;
	pf.SendParticle(part);
      }

    }

  } else {

    pf.ShipParticles(myid, 0, nbodies);
      
    int icount = 0;
    while (pf.RecvParticle(part)) {
      particles[part.indx] = part;
      icount++;
    }
  }


				// Default: set to max radius
  rmax = sqrt(fabs(rmax1));
  MPI_Bcast(&rmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#ifdef DEBUG
  cout << "read_bodies_and_distribute_binary: process " << myid 
       << " name=<" << name << "> bodies ["
       << particles.begin()->second.indx << ", "
       << particles.rbegin()->second.indx << "], ["
       << particles.begin()->first << ", "
       << particles.rbegin()->first << "]"
       << " #=" << particles.size() << endl;
#endif
}


struct Particle * Component::get_particles(int* number)
{
  static unsigned counter = 1;	// Sequence begins at 1
  static Particle part;
  static bool seq_state_ok = true;
  
  int curcount = 0;		// Counter for this bunch
  
#ifdef DEBUG
  if (*number < 0) {
    if (particles.size())
      cout << "get_particles: process " << myid 
	   << " <name=" << name << "> bodies ["
	   << particles.begin() ->second.indx << ", "
	   << particles.rbegin()->second.indx << "], ["
	   << particles.begin() ->first << ", "
	   << particles.rbegin()->first << "]" 
	   << " #=" << particles.size() << endl;
    else
      cout << "get_particles: process " << myid 
	   << " <name=" << name << "> #=" 
	   << particles.size() << endl;
  }
#endif
				// Reset
  if (*number < 0) {
    counter = 1;
    part = Particle(niattrib, ndattrib);
    seq_state_ok = true;
  }
				// Done?
  if (counter > nbodies_tot) {
    *number = 0;
    return 0;
  }

  map<unsigned int,  Particle> tlist;
  map<unsigned int,  Particle>::iterator cur;
  PartMapItr icur, ibeg, iend;

  unsigned icount;
  int beg = counter;
  int end = counter + PFbufsz - 1;

  MPI_Bcast(&beg, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&end, 1, MPI_INT, 0, MPI_COMM_WORLD);

#ifdef DEBUG
  cout << "get_particles: process " << myid 
       << " particles=<" << name << ">"
       << " number=" << particles.size()
       << " beg=" << beg
       << " end=" << end 
       << endl;
#endif

  for (int node=0; node<numprocs; node++) {
    
    if (myid==0) {
				// Do root's particle first
      if (node==0) {
      
	ibeg = particles.lower_bound(beg);
	iend = particles.lower_bound(end+1);

	icount = 0;
	for (icur=ibeg; icur!=iend; icur++)
	  pbuf[icount++] = icur->second;
#ifdef DEBUG
	cout << "get_particles: master loaded " 
	     << icount << " of its own particles" << endl << flush;
#endif    

#ifdef SANITY
	cout << "Process " << myid << ": count=" << icount
	     << ": want [" << beg << ", " << end << "]";
	if (particles.size()) {
	  cout << ": have [";
	  if (ibeg != particles.end()) cout << ibeg->first << ", ";
	  else cout << "end, ";
	  if (iend != particles.end()) cout << iend->first << "]";
	  else cout << "end)";
	  cout << ": cnts [" << particles.begin()->first << ", " << particles.rbegin()->first << "]"
	       << endl;
	} else {
	  cout << ": have none!" << endl;
	}
#endif

      } else {
	  
	unsigned number;
	pf.ShipParticles(0, node, number);

	icount = 0;
	while (pf.RecvParticle(part)) pbuf[icount++] = part;
#ifdef DEBUG
	cout << "Process " << myid 
	     << ": received " << icount << " particles from Slave " << node
	     << ", expected " << number
	     << endl << flush;
#endif    
      }

      // Load the ordered array
      for (unsigned n=0; n<icount; n++) {
	tlist[pbuf[n].indx] = pbuf[n];
	curcount++;
	counter++;
      }
      
				// Nodes send particles to master
    } else if (myid == node) {
      
	
      ibeg = particles.lower_bound(beg);
      iend = particles.lower_bound(end+1);

      icount = 0;
      for (icur=ibeg; icur!=iend; icur++) icount++;

#ifdef SANITY
				// Sanity
      cout << "Process " << myid << ": count=" << icount
	   << ": want [" << beg << ", " << end << "]";
      if (particles.size()) {
	cout << ": have [";
	if (ibeg != particles.end()) cout << ibeg->first << ", ";
	else cout << "end, ";
	if (iend != particles.end()) cout << iend->first << "]";
	else cout << "end)";
	cout << ": cnts [" << particles.begin()->first << ", " << particles.rbegin()->first << "]"
	     << endl;
      } else {
	cout << ": have none!" << endl;
      }
#endif

      pf.ShipParticles(0, myid, icount);
#ifdef DEBUG
      icount = 0;
#endif
      for (icur=ibeg; icur!=iend; icur++) {
#ifdef DEBUG
	if (icount<2) {
	  cout << "Component [" << myid << "]: sending ";
	  cout << setw(3) << icount
	       << setw(14) << icur->second.mass
#ifdef INT128
	       << setw(18) << icur->second.key.toHex()
#else
	       << setw(18) << hex << icur->second.key << dec
#endif
	    ;
	  for (int k=0; k<3; k++) cout << setw(14) << icur->second.pos[k];
	  cout << endl;
	}
	icount++;
#endif
	pf.SendParticle(icur->second);
      }

#ifdef DEBUG
      cout << "get_particles: process " << myid 
	   << ": sent " << icount << " particles to master"
	   << ", counter value=" << counter;
      if (particles.size())
	cout << ", nbodies_index=" << nbodies_index[node]
	     << ", seq_beg=" << ibeg->second.indx
	     << ", seq_end=" << iend->second.indx
	     << ", number found =" << icount
	     << ", first=" << particles.begin()->second.indx
	     << ", last=" << particles.rbegin()->second.indx;
      cout << endl << flush;
#endif    
	
    }

    MPI_Barrier(MPI_COMM_WORLD);

  }

  MPI_Bcast(&counter, 1, MPI_INT, 0, MPI_COMM_WORLD);

				// Return values
  *number = curcount;

#ifdef DEBUG
  if (myid==0) {
    cout << "get_particles: master size of tlist=" << tlist.size() 
    	 << " current count=" << curcount << endl;
  }
#endif

  int n=0;
  for (cur=tlist.begin(); cur!=tlist.end(); cur++) {
    pbuf[n++] = cur->second;
  }

#ifdef DEBUG
  cout << "Process " << myid 
       << ": received next counter=" << counter
       << " icount=" << icount;
  if (counter > nbodies_tot) cout << " [this means we are done]";
  cout << endl << flush;
#endif    

  if (myid==0 && seq_check && seq_state_ok) {
    bool seq_ok = true;
    unsigned n = beg;
    for (cur=tlist.begin(); cur!=tlist.end(); cur++) {
      if (cur->first != n++) {
	cout << "get_particles sequence error:"
	     << " expected=" << n
	     << " found=" << cur->first
	     << endl << flush;
	unsigned n = beg;
	cout << setw(80) << setfill('-') << '-' << endl << setfill(' ');
	cout << setw(10) << "Expect" << setw(10) << "Found" << endl;
	for (cur=tlist.begin(); cur!=tlist.end(); cur++)
	  cout << setw(10) << n++ << setw(10) << cur->first << endl;
	cout << setw(80) << setfill('-') << '-' << endl << setfill(' ');
	seq_ok = false;
	break;
      }
    }

    if (!seq_ok && seq_state_ok) {
      cout << "get_particles sequence failure in [" << beg
	   << ", " << end << "]" << endl;
      seq_state_ok = false;
    }

    if (counter > nbodies_tot) {
#ifdef DEBUG
      if (seq_state_ok)
	cout << "get_particles [" << name << "]: GOOD sequence!" << endl;
#endif
      if (!seq_state_ok)
	cout << "get_particles [" << name << "]: sequence ERROR!" << endl;
    }
  }

  return pbuf;
}


void Component::write_binary(ostream* out, bool real4)
{
  ComponentHeader header;

  if (myid == 0) {

    header.nbod  = nbodies_tot;
    header.niatr = niattrib;
    header.ndatr = ndattrib;
  
    ostringstream outs;
    outs << name << " : " << id << " : " << cparam << " : " << fparam;
    strncpy(header.info.get(), outs.str().c_str(), header.ninfochar);

    if (real4) rsize = sizeof(float);
    else       rsize = sizeof(double);
    unsigned long cmagic = magic + rsize;

    out->write((const char*)&cmagic, sizeof(unsigned long));

    if (!header.write(out)) {
      cerr << "Component::write_binary: Error writing particle header\n";
      MPI_Abort(MPI_COMM_WORLD, 34);
    }
  }

				// First bunch of particles
  int number = -1;
  Particle *p = get_particles(&number);

  float tf;
  double pot0, pv;
  while (p) {

    if (myid == 0) {
      for (int k=0; k<number; k++) {
	p[k].writeBinary(rsize, com0, comI, cov0, covI, indexing, out);
      }
    }
				// Next bunch of particles
    p = get_particles(&number);

  }
    
}


void Component::write_ascii(ostream* out, bool accel)
{
  int number = -1;
  Particle *p = get_particles(&number);

  while (p) {
    if (myid == 0) {
      for (int k=0; k<number; k++) {
	p[k].writeAscii(com0, comI, cov0, covI, indexing, accel, out);
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
    
    mtot1 += p->second.mass;

    for (int k=0; k<dim; k++) com1[k] += p->second.mass*p->second.pos[k];
    for (int k=0; k<dim; k++) cov1[k] += p->second.mass*p->second.vel[k];
    
  }
  
  MPI_Allreduce(&mtot1, &mtot0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(com1, com0, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(cov1, cov0, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  if (mtot0 > 0.0) {
    for (int k=0; k<dim; k++) com0[k] /= mtot0;
    for (int k=0; k<dim; k++) cov0[k] /= mtot0;
  }

  for (int k=0; k<dim; k++) {
    comI[k]   = com0[k];
    covI[k]   = cov0[k];
    center[k] = 0.0;
  }

  delete [] com1;
  delete [] cov1;
}

void Component::restart_com_system()
{
  MPI_Bcast(&com_restart, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (com_restart) {
    MPI_Bcast(&comI[0],   3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&covI[0],   3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&com0[0],   3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cov0[0],   3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&acc0[0],   3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&center[0], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

				// Particle loop
    PartMapItr p, pend = particles.end();
    for (p=particles.begin(); p != pend; p++) {

      for (int i=0; i<3; i++) {
	p->second.pos[i] -= com0[i] - comI[i];
	p->second.vel[i] -= cov0[i] - covI[i];
      }
    }

  }

}


struct thrd_pass_posn
{
  int id;
  Component *c;
  bool consp;
  bool tidal;
  bool com_system;
  unsigned mlevel;
  vector<double> com,  cov,  coa,  mtot;
  vector<double> comE, covE, mtotE;
};



void * fix_positions_thread(void *ptr)
{
  int id          =   static_cast<thrd_pass_posn*>(ptr)->id;
  Component *c    =   static_cast<thrd_pass_posn*>(ptr)->c;

  bool consp      =   static_cast<thrd_pass_posn*>(ptr)->consp;
  bool tidal      =   static_cast<thrd_pass_posn*>(ptr)->tidal;
  bool com_system =   static_cast<thrd_pass_posn*>(ptr)->com_system;

  unsigned mlevel =   static_cast<thrd_pass_posn*>(ptr)->mlevel;

  double *com     = &(static_cast<thrd_pass_posn*>(ptr)->com[0]);
  double *cov     = &(static_cast<thrd_pass_posn*>(ptr)->cov[0]);
  double *coa     = &(static_cast<thrd_pass_posn*>(ptr)->coa[0]);
  double *mtot    = &(static_cast<thrd_pass_posn*>(ptr)->mtot[0]);

  double *comE, *covE, *mtotE;

  if (consp && com_system) {
    comE          = &(static_cast<thrd_pass_posn*>(ptr)->com[0]);
    covE          = &(static_cast<thrd_pass_posn*>(ptr)->cov[0]);
    mtotE         = &(static_cast<thrd_pass_posn*>(ptr)->mtot[0]);
  }

  for (unsigned mm=mlevel; mm<=multistep; mm++) {

    int nbodies = c->levlist[mm].size();
    int nbeg    = nbodies*(id  )/nthrds;
    int nend    = nbodies*(id+1)/nthrds;

				// Particle loop
    for (int q=nbeg; q<nend; q++) {
    
      unsigned long n = c->levlist[mm][q];

      if (consp) {
	if (c->escape_com(*c->Part(n)) && c->Part(n)->iattrib[tidal]==0) {
				// Set flag indicating escaped particle
	  c->Part(n)->iattrib[tidal] = 1;

	  if (com_system) {	// Conserve momentum of center of mass
				// and compute center of acceleration
	    mtotE[mm] += c->Part(n)->mass;
	    for (unsigned k=0; k<3; k++) {
	      comE[3*mm+k] += c->Part(n)->mass*c->Part(n)->pos[k]; 
	      covE[3*mm+k] += c->Part(n)->mass*c->Part(n)->vel[k]; 
	    }
	  }
	  continue;
	}
	
	if (c->Part(n)->iattrib[tidal]==1) continue;
      }

      mtot[mm] += c->Part(n)->mass;

      // Compute new center of mass quantities
      //
      for (int k=0; k<c->dim; k++) {
	com[3*mm+k] += c->Part(n)->mass*c->Part(n)->pos[k];
	cov[3*mm+k] += c->Part(n)->mass*c->Part(n)->vel[k];
	coa[3*mm+k] += c->Part(n)->mass*c->Part(n)->acc[k];
      }
    }
  }

  return (NULL);
}
  

void Component::fix_positions(unsigned mlevel)
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
	cerr << "Process " << myid
	     << " Component::fix_positions: cannot make thread " << i
	     << ", errcode=" << errcode << endl;
	exit(19);
      }
    }
    
    //
    // Collapse the threads
    //
    for (int i=0; i<nthrds; i++) {
      if ((errcode=pthread_join(thrd[i], &retval))) {
	cerr << "Process " << myid
	     << " Component::fix_positions: thread join " << i
	     << " failed, errcode=" << errcode << endl;
	exit(20);
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
      if (isnan(com[k])) com_nan = true;
    for (int k=0; k<3; k++)
      if (isnan(cov[k])) cov_nan = true;
    for (int k=0; k<3; k++)
      if (isnan(coa[k])) coa_nan = true;
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
      comI[i] = (mtot0*comI[i] - comE[i])/(mtot0 - mtotE);
      com0[i] = (mtot0*comI[i] - comE[i])/(mtot0 - mtotE);
      covI[i] = (mtot0*covI[i] - covE[i])/(mtot0 - mtotE);
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

  if (com_system) {	   // Use local center of accel for com update
    for (int k=0; k<dim; k++) acc0[k]  = coa[k];
  } else {			// No mass, no acceleration?
    for (int k=0; k<dim; k++) acc0[k]  = 0.0;
  }

  if ((EJ & Orient::CENTER) && !EJdryrun) {
    Vector ctr = orient->currentCenter();
    bool ok    = true;
    for (int i=0; i<3; i++) {
      if (isnan(ctr[i+1])) ok = false;
    } 
    if (ok) {
      for (int i=0; i<3; i++) center[i] += ctr[i+1];
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

      if (c->freeze(n)) continue;

      mass = c->Part(n)->mass;
      pos  = c->Part(n)->pos;
      vel  = c->Part(n)->vel;
    
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
	cerr << "Process " << myid
	     << " Component::get_angmom: cannot make thread " << i
	     << ", errcode=" << errcode << endl;
	exit(19);
      }
    }
    
    //
    // Collapse the threads
    //
    for (int i=0; i<nthrds; i++) {
      if ((errcode=pthread_join(thrd[i], &retval))) {
	cerr << "Process " << myid
	     << " Component::get_angmom: thread join " << i
	     << " failed, errcode=" << errcode << endl;
	exit(20);
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
  if (numb >= nbodmax) numb = nbodmax;
  return numb;
}


void Component::setup_distribution(void)
{
  ofstream *out;
  int n;

				// Needed for both master and slaves
  nbodies_index = vector<unsigned int>(numprocs);
  nbodies_table = vector<unsigned int>(numprocs);

  if (myid == 0) {

    orates = vector<double>(numprocs);
    trates = vector<double>(numprocs);

    for (n=0; n<numprocs; n++) {

      if (n == 0)
	nbodies_table[n] = nbodies_index[n] = 
	  max<int>(1, min<int>((int)(comp.rates[n] * nbodies_tot), nbodies_tot));
      else {
	if (n < numprocs-1)
	  nbodies_index[n] = (int)(comp.rates[n] * nbodies_tot) + 
	    nbodies_index[n-1];
	else
	  nbodies_index[n] = nbodies_tot;
      
	nbodies_table[n] = nbodies_index[n] - nbodies_index[n-1];
      }

    }

    string outrates = outdir + "current.processor.rates." + runtag;

    out = new ofstream(outrates.c_str(), ios::out | ios::app);
    if (out) {
      *out << "# " << endl;
      *out << "# Time=" << tnow << " Component=" << name << endl;
      *out << "# " 
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
      
      for (n=0; n<numprocs; n++)
	*out << "  "
	    << setw(15) << comp.rates[n]
	    << setw(15) << 1.0 - comp.rates[n]*nbodies_tot/nbodies_table[n]
	    << setw(15) << nbodies_index[n]
	    << setw(15) << nbodies_table[n]
	    << endl;

      delete out;
    }

  }


  MPI_Bcast(&nbodies_index[0], numprocs, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nbodies_table[0], numprocs, MPI_INT, 0, MPI_COMM_WORLD);

}

void Component::load_balance(void)
{
  MPI_Status status;
  vector<unsigned int> nbodies_index1(numprocs);
  vector<unsigned int> nbodies_table1(numprocs);
  ofstream *out, *log;


  if (myid == 0) {

    vector<double> orates1(numprocs);
    vector<double> trates1(numprocs);

    for (int n=0; n<numprocs; n++) {

      if (n == 0)
	nbodies_table1[n] = nbodies_index1[n] = 
	  max<int>(1, min<int>((int)(comp.rates[n] * nbodies_tot), nbodies_tot));
      else {
	if (n < numprocs-1)
	  nbodies_index1[n] = (int)(comp.rates[n] * nbodies_tot) + 
	    nbodies_index1[n-1];
	else
	  nbodies_index1[n] = nbodies_tot;
      
	nbodies_table1[n] = nbodies_index1[n] - nbodies_index1[n-1];
      }

    }

    string outrates = outdir + "current.processor.rates." + runtag;
    string rateslog = outdir + "current.processor.rates.log." + runtag;

    out = new ofstream(outrates.c_str(), ios::out | ios::app);
    log = new ofstream(rateslog.c_str(), ios::out | ios::app);

    if (*out) {
      *out << "# " << endl;
      *out << "# Time=" << tnow << " Component=" << name << endl;
      *out << "# " 
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
	*out << "  "
	     << setw(15) << comp.rates[n]
	     << setw(15) << 1.0 - comp.rates[n]*nbodies_tot/nbodies_table1[n]
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
  loadb.erase(loadb.begin(), loadb.end());
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

  if (myid==0 && *log) 
    {
      *log << setw(72) << setfill('.') << ".\n" << setfill(' ');
      *log << "Time=" << tnow << " Component=" << name << endl;
      *log << endl;
      *log << "List:\n";
      log->setf(ios::left);
      *log << setw(4) << "N"
	   << setw(6) << "Index"
	   << setw(10) << "Old"
	   << setw(10) << "New"
	   << endl;

      char c = log->fill('-');
      *log << setw(4) << "|"
	   << setw(6) << "|"
	   << setw(10) << "|"
	   << setw(10) << "|"
	   << endl;
      log->fill(c);
      
      for (int i=0; i<2*numprocs; i++) {
	
	*log << setw(4) << i
	     << setw(6) << loadb[i].indx;
	if (loadb[i].s)
	  *log << setw(10) << " " << setw(10) << loadb[i].top;
	else 
	  *log << setw(10) << loadb[i].top;
	*log << endl;
      }

    }


  if (myid==0 && *log) 
    {
      *log << "\nAnalysis:\n";
      log->setf(ios::left);
      *log << setw(10) << "Interval"
	   << setw(10) << "Number"
	   << setw(10) << "Old"
	   << setw(10) << "New"
	   << setw(10) << "Action"
	   << endl;

      char c = log->fill('-');
      *log << setw(10) << "|"
	   << setw(10) << "|"
	   << setw(10) << "|"
	   << setw(10) << "|"
	   << setw(10) << "|"
	   << endl;
      log->fill(c);
    }


  int iold=0, inew=0;
  
				// Offset will keep track of position in
				// original vector
  int nump;
  vector<int> loc(numprocs, 0);

  vector<int> nlist;

  for (int i=0; i<2*numprocs-2; i++) {

				// Assign new interval
    if (loadb[i].s) inew = loadb[i].indx+1;
    else            iold = loadb[i].indx+1;
    
    if (myid==0 && *log)
      *log << setw(10) << i
	   << setw(10) << loadb[i+1].top - loadb[i].top
	   << setw(10) << iold
	   << setw(10) << inew;
    
				// Number of particles to be shifted
    nump = loadb[i+1].top - loadb[i].top;
    
    ostringstream msg;
    
    if (inew==iold || nump==0) 
      msg << "Do nothing";
    else if (inew>iold) {
      msg << "Add " << nump << " from #" << iold << " to #" << inew;
      
      nlist.erase(nlist.begin(), nlist.end());

      map<unsigned long, Particle>::reverse_iterator it = particles.rbegin();
      for (int n=0; n<nump; n++) {
	nlist.push_back(it->first);
	it++;
      }
      
      add_particles(iold, inew, nlist);
      
#ifdef DEBUG
      if (myid==iold) {
	if (particles.size())
	  cout << "Process " << myid << ": new ends :"
	       << "  beg seq=" << particles.begin() ->second.indx
	       << "  end seq=" << particles.rbegin()->second.indx
	       << endl;
	else
	  cout << "Process " << myid << ": no particles!"
	       << endl;
      }
#endif
    } else if (iold>inew) {
      msg << "Add " << nump << " from #" << iold << " to #" << inew;

      nlist.erase(nlist.begin(), nlist.end());

      PartMapItr it = particles.begin();
      for (int n=0; n<nump; n++) {
	nlist.push_back(it->first);
	it++;
      }

      add_particles(iold, inew, nlist);

#ifdef DEBUG
      if (myid==iold) {
	if (particles.size())
	  cout << "Process " << myid << ": new ends :"
	       << "  beg seq=" << particles.begin() ->second.indx
	       << "  end seq=" << particles.rbegin()->second.indx
	       << endl;
	else
	  cout << "Process " << myid << ": no particles!"
	       << endl;
      }
#endif

    }

    if (myid==0 && *log) *log << setw(10) << msg.str() << endl;
  }

  
				// update indices
  nbodies = nbodies_table1[myid];
  nbodies_index = nbodies_index1;
  nbodies_table = nbodies_table1;
  
  if (seq_check) {
    
    char msgbuf[200];		// Only need 31 characters . . .

    if (myid==0) *log << endl << "Post load-balance sequence check:" 
		      << endl << endl;

    for (int i=0; i<numprocs; i++) {
      if (myid==i) {
	ostringstream msg;
	msg << "Process " << setw(4) << myid << ":"
	    << setw(9) << particles[0].indx
	    << setw(9) << particles[nbodies-1].indx;
	strcpy(msgbuf, msg.str().c_str());
	if (myid!=0) 
	  MPI_Send(msgbuf, 200, MPI_CHAR, 0, 81, MPI_COMM_WORLD);
      }
      
      if (myid==0) {
	if (myid!=i)
	  MPI_Recv(msgbuf, 200, MPI_CHAR, i, 81, MPI_COMM_WORLD, &status);

	*log << msgbuf << endl;
      }
      

    }

    if (myid==0) seq_beg = 1;
    else         seq_beg = nbodies_index[myid-1]+1;
    
				// Explicit check
    int nbad1 = 0, nbad=0;
    for (unsigned i=0; i<nbodies; i++) {
      if (particles[i].iattrib[0] != static_cast<int>(seq_beg+i)) {
	cout << "Process " << myid << ": sequence error on load balance,"
	     << " component=" << name
	     << " i=" << i
	     << " seq=" << particles[i].iattrib[0]
	     << " expected=" << seq_beg+i
	     << endl << flush;
	nbad1++;
      }
    }

    MPI_Allreduce(&nbad1, &nbad, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (nbad) {
      if (myid==0) cout << nbad << " bad states\n";
      MPI_Finalize();
      exit(-1);
    }
    
    if (myid==0) *log << "\nSequence check ok!\n";
  }


  if (myid==0) {
    out->close();
    delete out;
    log->close();
    delete log;
  }

}


void Component::add_particles(int from, int to, vector<int>& plist)
{
  unsigned number = plist.size();
  vector<int>::iterator it=plist.begin();

  unsigned icount, counter=0;

  pf.ShipParticles(to, from, number);

  if (myid == from) {
    
    while (counter < number) {

      icount = 0;
      while (icount < PFbufsz && counter < number) {

	pf.SendParticle(particles[*it]);
	particles.erase(*it);
      
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
  
    Particle temp;

    while (counter < number) {

      while (pf.RecvParticle(temp)) {
	particles[temp.indx] = temp;
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
			    (particles[indx].pos[i] - comI[i] - center[i])*
			    (particles[indx].pos[i] - comI[i] - center[i]);
  if (r2 > rtrunc*rtrunc) return true;
  else return false;
}

bool Component::escape_com(const Particle& p)
{
  double r2 = 0.0;
  for (int i=0; i<3; i++) r2 += 
			    (p.pos[i] - comI[i] - center[i])*
			    (p.pos[i] - comI[i] - center[i]);
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
  Particle part(niattrib, ndattrib);

  vector<int>::iterator it = redist.begin();
  vector<unsigned> tlist;

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
      tlist.erase(tlist.begin(), tlist.end());
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
	  pf.ShipParticles(tonode, curnode, icount);
	  if (myid==curnode) {
	    for (unsigned i=0; i<icount; i++) {
	      pf.SendParticle(particles[tlist[i]]);
	      particles.erase(tlist[i]);
	    }
	  }
	  if (myid==lastnode) {
	    while (pf.RecvParticle(part))
	      particles[part.indx] = part;
	  }
	  tlist.erase(tlist.begin(), tlist.end());
	  icount = 0;
	}

				// Add the particle
	tlist.push_back(indx);
	icount++;
      }
	
    } // End of particles on this node
    
  } // Next stanza

}
