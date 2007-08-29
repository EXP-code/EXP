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
#include <SphereEJCOM.H>
#include <Cylinder.H>
#include <Cube.H>
#include <Slab.H>
#include <SlabSL.H>
#include <Direct.H>
#include <NoForce.H>
#include <Orient.H>

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
  EJ = 0;
  nEJkeep = 100;
  nEJwant = 500;
  eEJ0 = -0.5;
  EJkinE = true;
  EJext = false;
  EJdiag = false;
  EJdryrun = false;
  EJx0 = 0.0;
  EJy0 = 0.0;
  EJz0 = 0.0;
  EJu0 = 0.0;
  EJv0 = 0.0;
  EJw0 = 0.0;
  EJlinear = false;

  binary = false;
  npart = false;

  adiabatic = false;
  ton  = -1.0e20;
  toff =  1.0e20;
  twid = 0.1;

  rtrunc = 1.0e20;
  rcom = 1.0e20;
  consp = false;
  tidal = 0;

  com_system = false;
  com_log = false;

  force   = 0;			// Null out pointers
  orient  = 0;
  buf     = 0;

  com     = 0;
  cov     = 0;
  center  = 0;
  EJcen   = 0;
  angmom  = 0;
  ps      = 0;

  com0    = 0;
  cov0    = 0;
  acc0    = 0;
  comI    = 0;
  covI    = 0;

  ordered = true;
  seq_check = false;

  read_bodies_and_distribute_ascii();
}


Component::Component(istream *in)
{
  EJ = 0;
  nEJkeep = 100;
  nEJwant = 500;
  eEJ0 = -0.5;
  EJkinE = true;
  EJext = false;
  EJdiag = false;
  EJdryrun = false;
  EJx0 = 0.0;
  EJy0 = 0.0;
  EJz0 = 0.0;
  EJu0 = 0.0;
  EJv0 = 0.0;
  EJw0 = 0.0;
  EJlinear = false;

  binary = true;
  npart = false;

  adiabatic = false;
  ton  = -1.0e20;
  toff =  1.0e20;
  twid = 0.1;

  rtrunc = 1.0e20;
  rcom = 1.0e20;
  consp = false;
  tidal = 0;

  com_system = false;
  com_log = false;
  com_restart = 0;

  force   = 0;			// Null out pointers
  orient  = 0;
  buf     = 0;

  com     = 0;
  cov     = 0;
  center  = 0;
  EJcen   = 0;
  angmom  = 0;
  ps      = 0;

  com0    = 0;
  cov0    = 0;
  acc0    = 0;
  comI    = 0;
  covI    = 0;

  ordered = true;

  read_bodies_and_distribute_binary(in);
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
    
    if (!datum.first.compare("nEJkeep"))  nEJkeep = atoi(datum.second.c_str());

    if (!datum.first.compare("nEJwant"))  nEJwant = atoi(datum.second.c_str());

    if (!datum.first.compare("eEJ0"))     eEJ0 = atof(datum.second.c_str());

    if (!datum.first.compare("EJx0"))     EJx0 = atof(datum.second.c_str());

    if (!datum.first.compare("EJy0"))     EJy0 = atof(datum.second.c_str());

    if (!datum.first.compare("EJz0"))     EJz0 = atof(datum.second.c_str());

    if (!datum.first.compare("EJu0"))     EJu0 = atof(datum.second.c_str());

    if (!datum.first.compare("EJv0"))     EJv0 = atof(datum.second.c_str());

    if (!datum.first.compare("EJw0"))     EJw0 = atof(datum.second.c_str());

    if (!datum.first.compare("EJkinE"))   EJkinE = atoi(datum.second.c_str()) ? true : false;

    if (!datum.first.compare("EJext"))    EJext =  atoi(datum.second.c_str()) ? true : false;

    if (!datum.first.compare("EJdiag"))   EJdiag = atoi(datum.second.c_str()) ? true : false;

    if (!datum.first.compare("EJdryrun")) EJdryrun = atoi(datum.second.c_str()) ? true : false;

    if (!datum.first.compare("EJlinear")) EJlinear = atoi(datum.second.c_str()) ? true : false;

    if (!datum.first.compare("rmax"))     rmax = atof(datum.second.c_str());

    if (!datum.first.compare("ton"))      {ton = atof(datum.second.c_str()); adiabatic = true;}

    if (!datum.first.compare("toff"))     {toff= atof(datum.second.c_str()); adiabatic = true;}

    if (!datum.first.compare("twid"))     {twid = atof(datum.second.c_str()); adiabatic = true;}

    if (!datum.first.compare("rtrunc"))   {rtrunc = atof(datum.second.c_str());}

    if (!datum.first.compare("rcom"))     {rcom = atof(datum.second.c_str());}
    
    if (!datum.first.compare("scheck"))   seq_check = atoi(datum.second.c_str()) ? true : false;

				// Next parameter
    token = tokens(",");
  }


				// Instantiate the force ("reflection" by hand)

  if ( !id.compare("bessel") ) {
    force = new Bessel(fparam);
    dim = 3;
  }
  else if ( !id.compare("c_brock") ) {
    force = new CBrock(fparam);
    dim = 3;
  }
  else if ( !id.compare("c_brock_disk") ) {
    force = new CBrockDisk(fparam);
    dim = 2;
  }
  else if ( !id.compare("hernq") ) {
    force = new Hernquist(fparam);
    dim = 3;
  }
  else if ( !id.compare("sphereSL") ) {
    force = new Sphere(fparam);
    dim = 3;
  }
  else if ( !id.compare("sphereEJCOM") ) {
    force = new SphereEJCOM(fparam);
    dim = 3;
  }
  else if ( !id.compare("cube") ) {
    force = new Cube(fparam);
    dim = 3;
  }
  else if ( !id.compare("slab") ) {
    force = new Slab(fparam);
    dim = 3;
  }
  else if ( !id.compare("slabSL") ) {
    force = new SlabSL(fparam);
    dim = 3;
  }
  else if ( !id.compare("cylinder") ) {
    force = new Cylinder(fparam);
    dim = 3;
  }
  else if ( !id.compare("direct") ) {
    force = new Direct(fparam);
    dim = 3;
  }
  else if ( !id.compare("noforce") ) {
    force = new NoForce(fparam);
    dim = 3;
  }
  else {
    string msg("I don't know about the force: ");
    msg += id;
    bomb(msg);
  }

  force->RegisterComponent(this);

  com    = new double [3];
  center = new double [3];
  cov    = new double [3];
  angmom = new double [3];
  ps     = new double [6];
  EJcen  = new double [3];

				// For COM system
  com0   = new double[3];
  cov0   = new double[3];
  acc0   = new double[3];
  comI   = new double[3];
  covI   = new double[3];

  for (int k=0; k<3; k++) 
    com[k] = center[k] = cov[k] = com0[k] = cov0[k] = acc0[k] = 
      comI[k] = covI[k] = angmom[k] = EJcen[k] = 0.0;
  

  if (com_system) {

    initialize_com_system();

    if (myid==0) {
      cout << name << ": center of mass system is *ON*, rtrunc=" << rtrunc;
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

		if ( fabs(tnow - ttim) < 1.0e-8 ) {
		  istringstream istr(line);

		  istr >> ttim;

		  for (int k=0; k<3; k++) istr >> com0[k];
		  for (int k=0; k<3; k++) istr >> cov0[k];
		  for (int k=0; k<3; k++) istr >> acc0[k];
		  for (int k=0; k<3; k++) istr >> center[k];
	    
		  cout << "\t\tRead com log for " << name 
		       << " at T=" << ttim << ", using:";

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
		     << " nkeep=" << nEJkeep
		     << " nwant=" << nEJwant
		     << " eEJ=" << eEJ0
		     << " EJkinE=" << EJkinE
		     << " EJext=" << EJext;

    else if (myid==0) {
      cout << name << ": EJ centering *ON*";
      if (EJkinE) cout << ", using particle kinetic energy";
      if (EJext)  cout << ", using external potential";
      if (EJdryrun) cout << ", dryrun";
      cout << "\n";
    }
      
    
    string EJlogfile = outdir + name + ".orient." + runtag;

    unsigned EJctl = 0;
    if (EJdiag)		EJctl |= Orient::DIAG;
    if (EJkinE)		EJctl |= Orient::KE;
    if (EJext)		EJctl |= Orient::EXTERNAL;

    orient = new Orient(nEJkeep, nEJwant, eEJ0, EJ, EJctl, EJlogfile);

    if (restart && (EJ & Orient::CENTER)) {
      for (int i=0; i<3; i++) EJcen[i] = (orient->currentCenter())[i+1];
    } else {
      orient -> set_center(EJx0, EJy0, EJz0);
      orient -> set_cenvel(EJu0, EJv0, EJw0);
      if (EJlinear) orient -> set_linear();
      EJcen[0] = EJx0;
      EJcen[1] = EJy0;
      EJcen[2] = EJz0;
    }

    if (EJdiag) cout << "Process " << myid << ": Orient successful\n";
  }

  if (myid == 0) {		// Flag messages for diagnostics
    
    if (EJ & Orient::AXIS)
      cout << name << ": AXIS orientation is *ON*\n";

    if (EJ & Orient::CENTER) {
      cout << name 
	   << ": CENTER finding is *ON*";

      if (restart)
	cout << ", current center on restart: x, y, z: " 
	     << EJcen[0] << ", " 
	     << EJcen[1] << ", " 
	     << EJcen[2] << "\n";
      else
	cout << ", user specified initial center: x, y, z: " 
	     << EJx0 << ", " 
	     << EJy0 << ", " 
	     << EJz0 << "\n";
    }

  }
    
}


Component::~Component(void)
{
  delete force;

  delete orient;
  delete [] buf;

  delete [] com;
  delete [] center;
  delete [] cov;
  delete [] angmom;
  delete [] ps;
  delete [] EJcen;

  delete [] com0;
  delete [] cov0;
  delete [] acc0;
  delete [] comI;
  delete [] covI;
}

void Component::bomb(const string& msg)
{
  cerr << "Component [" << name << ", " << id << "]: " << msg << endl;
  exit(-1);
}

void Component::part_to_Particle(Partstruct& str, Particle& cls)
{
  cls.mass = str.mass;
  for (int j=0; j<3; j++) {
      cls.pos[j] = str.pos[j];
      cls.vel[j] = str.vel[j];
      cls.acc[j] = str.acc[j];
  }
  cls.pot = str.pot;
  cls.potext = str.potext;

  cls.level = str.level;
  cls.indx = str.indx;

  cls.iattrib = vector<int>(niattrib);
  for (int j=0; j<niattrib; j++) cls.iattrib[j] = str.iatr[j];

  cls.dattrib = vector<double>(ndattrib);
  for (int j=0; j<ndattrib; j++) cls.dattrib[j] = str.datr[j];

}

void Component::Particle_to_part(Partstruct& str, Particle& cls)
{
  str.mass = cls.mass;
  for (int j=0; j<3; j++) {
      str.pos[j] = cls.pos[j];
      str.vel[j] = cls.vel[j];
      str.acc[j] = cls.acc[j];
  }
  str.pot = cls.pot;
  str.potext = cls.potext;

  str.level = cls.level;
  str.indx = cls.indx;

  for (int j=0; j<niattrib; j++) str.iatr[j] = cls.iattrib[j];

  for (int j=0; j<ndattrib; j++) str.datr[j] = cls.dattrib[j];
}

void Component::read_bodies_and_distribute_ascii(void)
{
				// MPI buffers
  MPI_Status status0;
				// Open file

  ifstream *fin;
  const int nline = 256;
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
  MPI_Bcast(&niattrib, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ndattrib, 1, MPI_INT, 0, MPI_COMM_WORLD);

				// Assign particle structure buffer
  buf = new Partstruct [nbuf];

				// Make MPI datatype
  
  MPI_Datatype type[10] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
			   MPI_DOUBLE, MPI_DOUBLE, 
			   MPI_UNSIGNED, MPI_UNSIGNED_LONG,
			   MPI_INT,    MPI_DOUBLE};

				// Get displacements
  MPI_Aint disp[9];
  MPI_Get_address(&buf[0].mass,		&disp[0]);
  MPI_Get_address(&buf[0].pos,		&disp[1]);
  MPI_Get_address(&buf[0].vel,		&disp[2]);
  MPI_Get_address(&buf[0].acc,		&disp[3]);
  MPI_Get_address(&buf[0].pot,		&disp[4]);
  MPI_Get_address(&buf[0].potext,	&disp[5]);
  MPI_Get_address(&buf[0].potext,	&disp[5]);
  MPI_Get_address(&buf[0].level,	&disp[6]);
  MPI_Get_address(&buf[0].indx,		&disp[7]);
  MPI_Get_address(&buf[0].iatr,		&disp[8]);
  MPI_Get_address(&buf[0].datr,		&disp[9]);

  for (int i=9; i>=0; i--) disp[i] -= disp[0];
  
				// Block offsets
  int blocklen[10] = {1, 3, 3, 3, 1, 1, 1, 1, nimax, ndmax};
  
  MPI_Type_create_struct(10, blocklen, disp, type, &Particletype);
  MPI_Type_commit(&Particletype);


  double rmax1=0.0, r2;

  if (nbodies_tot > nbodmax*numprocs) {
    if (myid==0) {
      cerr << "Not enough space on all processors to hold phase space\n";
      cerr << "nbodmax should be at least "
	   << (int)( (double)nbodies_tot/numprocs + 1) << endl;
    }
    MPI_Finalize();
    exit(-1);
  }

  is_init = 1;
  setup_distribution();
  is_init = 0;

				// Form cumulative and differential bodies list
  
  plist = vector<int>(numprocs);
  for (int i=0; i<numprocs; i++) plist[i] = nbodies_index[i];
  ncount = plist;
  for (int i=1; i<numprocs; i++) ncount[i] -= plist[i-1];

  if (myid==0) seq_beg = 1;
  else seq_beg = plist[myid-1]+1;
  seq_end = plist[myid];

  Particle part;
  part.iattrib = vector<int>(niattrib);
  part.dattrib = vector<double>(ndattrib);

  if (myid==0) {
				// Read root node particles

				// First line
    {
      fin->getline(line, nline);
      istringstream ins(line);

      ins >> part.mass;
      for (int j=0; j<3; j++) ins >> part.pos[j];
      for (int j=0; j<3; j++) ins >> part.vel[j];
      for (int j=0; j<3; j++) part.acc[j] = 0.0;
      part.pot = part.potext = 0.0;

      part.level = 0;
      part.indx = 1;

      for (int j=0; j<niattrib; j++) {
	ins >> part.iattrib[j];
	if (!ins) part.iattrib[j] = 0;
      }

      for (int j=0; j<ndattrib; j++) {
	ins >> part.dattrib[j];
	if (!ins) part.dattrib[j] = 0;
      }

      rmax1 = 0.0;
      for (int j=0; j<3; j++) rmax1 += part.pos[j]*part.pos[j];
    }

    particles.push_back(part);

				// Remainder of Node 0's particles

    for (int i=2; i<=ncount[0]; i++) {
      
      fin->getline(line, nline);
      istringstream ins(line);

      ins >> part.mass;
      for (int j=0; j<3; j++) ins >> part.pos[j];
      for (int j=0; j<3; j++) ins >> part.vel[j];
      part.pot = part.potext = 0.0;

      part.level = 0;
      part.indx = i;

      for (int j=0; j<niattrib; j++) {
	ins >> part.iattrib[j];
	if (!ins) part.iattrib[j] = 0;
      }

      for (int j=0; j<ndattrib; j++) {
	ins >> part.dattrib[j];
	if (!ins) part.dattrib[j] = 0;
      }

      r2 = 0.0;
      for (int j=0; j<3; j++) r2 += part.pos[j]*part.pos[j];
      rmax1 = max<double>(r2, rmax1);

				// Load the particle
      particles.push_back(part);
    }

    nbodies = ncount[0];

    int icount, ibufcount;
    for (int n=1; n<numprocs; n++) {

      MPI_Send(&ncount[n], 1, MPI_INT, n, 1, MPI_COMM_WORLD);

      icount = 0;
      ibufcount = 0;
      while (icount < ncount[n]) {

	fin->getline(line, nline);
	istringstream ins(line);

	ins >> buf[ibufcount].mass;
	for (int j=0; j<3; j++) ins >> buf[ibufcount].pos[j];
	for (int j=0; j<3; j++) ins >> buf[ibufcount].vel[j];
	for (int j=0; j<3; j++) buf[ibufcount].acc[j] = 0.0;
	buf[ibufcount].pot = buf[ibufcount].potext = 0.0;

	buf[ibufcount].level = 0;
	buf[ibufcount].indx = nbodies_index[n-1] + 1 + icount;

	for (int j=0; j<niattrib; j++) {
	  ins >> buf[ibufcount].iatr[j];
	  if (!ins) buf[ibufcount].iatr[j] = 0;
	}

	for (int j=0; j<ndattrib; j++) {
	  ins >> buf[ibufcount].datr[j];
	  if (!ins) buf[ibufcount].datr[j] = 0;
	}

	r2 = 0.0;
	for (int k=0; k<3; k++) 
	  r2 += buf[ibufcount].pos[k]*buf[ibufcount].pos[k];

	rmax1 = max<double>(r2, rmax1);

	ibufcount++;
	icount++;

	if (ibufcount == nbuf || icount == ncount[n]) {
	  MPI_Send(&ibufcount, 1, MPI_INT, n, 2, MPI_COMM_WORLD);
	  MPI_Send(buf, ibufcount, Particletype, n, 3, MPI_COMM_WORLD);
	  ibufcount = 0;
	}

      }

    }

  } else {

    MPI_Recv(&nbodies, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status0);
      
    int icount = 0, isent;
    while (icount < nbodies) {
      MPI_Recv(&isent, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &status0);
      MPI_Recv(buf, isent, Particletype, 0, 3, MPI_COMM_WORLD, &status0);
      for (int i=0; i<isent; i++) {
	part_to_Particle(buf[i], part);
	particles.push_back(part);
	icount++;
      }
    }
  }
  
				// Sanity check
  if (seq_check) {
    if (particles.size()) {
      if (seq_beg != particles[0].indx || 
	  seq_end != particles[nbodies-1].indx) {
	cout << "Process " << myid << ": sequence error on init,"
	     << " seq_beg=" << seq_beg
	     << " seq_end=" << seq_end
	     << " level[1]=" << particles[0].level
	     << " level[N]=" << particles[nbodies-1].level
	     << " seq[1]=" << particles[0].level
	     << " seq[1]=" << particles[0].indx
	     << " seq[N]=" << particles[nbodies-1].indx
	     << " nbodies=" << nbodies
	     << endl << flush;
	MPI_Abort(MPI_COMM_WORLD, -1);
      }
    }
  }
				// Default: set to max radius
				// can be overriden by parameter

  rmax = sqrt(fabs(rmax1));
  MPI_Bcast(&rmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

				// COM HERE?
  if (myid==0) delete fin;

}

void Component::get_next_particle_from_file(Partstruct *onepart, istream *in)
{
  in->read((char *)&(onepart->mass), sizeof(double));
  for (int i=0; i<3; i++) in->read((char *)&(onepart->pos[i]), sizeof(double));
  for (int i=0; i<3; i++) in->read((char *)&(onepart->vel[i]), sizeof(double));
  in->read((char *)&(onepart->pot), sizeof(double));
  onepart->potext = 0.0;
  onepart->level = 0;
  onepart->indx = ++seq_cur;
  for (int i=0; i<niattrib; i++) 
    in->read((char *)&(onepart->iatr[i]), sizeof(int));
  for (int i=0; i<ndattrib; i++) 
    in->read((char *)&(onepart->datr[i]), sizeof(double));
}


int Component::get_next_particle(Partstruct *onepart)
{
				// Initialize read
  if (!npart) {
    npart_tot = -1;
    npart_p = get_particles(&npart_tot);
    if (npart_p == NULL) {
      cerr << "Component::get_next_particle: no particles found!?!\n";
      return 0;
    }
    npart = true;
    npart_cur = 0;
  }

				// Get a new bunch
  if (npart_cur == npart_tot) {
    npart_p = get_particles(&npart_tot);
    if (npart_p == NULL) {
      npart = false;
      return 0;
    }
    npart_cur = 0;
  }

  if (myid == 0) 
    onepart = &(npart_p[npart_cur++]);
  else
    onepart = NULL;


  return 1;
}


void Component::read_bodies_and_distribute_binary(istream *in)
{
				// MPI buffers
  MPI_Status status0;

				// Get component header
  ComponentHeader header;

				// Node local parameter buffer
  int ninfochar;
  char *info;
  
  if (myid == 0) {

    if(!header.read(in)) {
      cerr << "Error reading component header\n";
      MPI_Abort(MPI_COMM_WORLD, 12);
      exit(-1);
    }

    nbodies_tot = header.nbod;
    niattrib = header.niatr;
    ndattrib = header.ndatr;
    ninfochar = header.ninfochar;
    info = new char [ninfochar+1];
    memcpy(info, header.info, ninfochar);
    info[ninfochar] = '\0';
  }

				// Broadcast attributes for this
				// phase-space component
  MPI_Bcast(&nbodies_tot, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&niattrib, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ndattrib, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ninfochar, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (myid) info = new char [ninfochar+1];
  MPI_Bcast(info, ninfochar, MPI_CHAR, 0, MPI_COMM_WORLD);


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

  buf = new Partstruct [nbuf];

				// Parse info field to get 
				// id and parameter strings
  StringTok<string> tokens(info);
  name = trimLeft(trimRight(tokens(":")));
  id = trimLeft(trimRight(tokens(":")));
  cparam = trimLeft(trimRight(tokens(":")));
  fparam = trimLeft(trimRight(tokens(":")));

				// Backward compatibility
  if (fparam.size() == 0) {
    fparam = cparam;
    cparam = "";
  }

  delete [] info;

				// Informational output
  if (myid==0)
    cout << setw(60) << setfill('-') << "-" << endl << setfill(' ')
	 << "--- New Component" << endl
	 << setw(20) << " name   :: " << name           << endl
	 << setw(20) << " id     :: " << id             << endl
	 << setw(20) << " cparam :: " << cparam         << endl
	 << setw(20) << " fparam :: " << fparam         << endl
	 << setw(60) << setfill('-') << "-" << endl << setfill(' ');

				// Make MPI datatype
  
  MPI_Datatype	type[10] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
			    MPI_DOUBLE, MPI_DOUBLE, 
			    MPI_UNSIGNED, MPI_UNSIGNED_LONG,
			    MPI_INT,    MPI_DOUBLE};
  
				// Get displacements
  MPI_Aint disp[10];
  MPI_Get_address(&buf[0].mass,		&disp[0]);
  MPI_Get_address(&buf[0].pos,		&disp[1]);
  MPI_Get_address(&buf[0].vel,		&disp[2]);
  MPI_Get_address(&buf[0].acc,		&disp[3]);
  MPI_Get_address(&buf[0].pot,		&disp[4]);
  MPI_Get_address(&buf[0].potext,	&disp[5]);
  MPI_Get_address(&buf[0].level,	&disp[6]);
  MPI_Get_address(&buf[0].indx,		&disp[7]);
  MPI_Get_address(&buf[0].iatr,		&disp[8]);
  MPI_Get_address(&buf[0].datr,		&disp[9]);

  for (int i=9; i>=0; i--) disp[i] -= disp[0];
  
				// Block offsets
  int blocklen[10] = {1, 3, 3, 3, 1, 1, 1, 1, nimax, ndmax};
  
  MPI_Type_create_struct(10, blocklen, disp, type, &Particletype);
  MPI_Type_commit(&Particletype);

  double rmax1, r2;

  if (nbodies_tot > nbodmax*numprocs) {
    if (myid==0) {
      cerr << "Not enough space on all processors to hold phase space\n";
      cerr << "nbodmax should be at least "
	   << (int)( (double)nbodies_tot/numprocs + 1) << endl;
    }
    MPI_Finalize();
    exit(-1);
  }

  is_init = 1;
  setup_distribution();
  is_init = 0;

				// Form cumulative and differential bodies list
  
  plist = vector<int>(numprocs);
  for (int i=0; i<numprocs; i++) plist[i] = nbodies_index[i];
  ncount = plist;
  for (int i=1; i<numprocs; i++) ncount[i] -= plist[i-1];

  if (myid==0) seq_beg = 1;
  else seq_beg = plist[myid-1]+1;
  seq_end = plist[myid];

  Partstruct onepart;
  Particle part;
  part.iattrib = vector<int>(niattrib);
  part.dattrib = vector<double>(ndattrib);


  if (myid==0) {
				// Read root node particles

    cout << "Count=" << ncount[0] << endl;
    seq_cur = 0;

    rmax1 = 0.0;
    for (int i=1; i<=ncount[0]; i++)
    {
      get_next_particle_from_file(&onepart, in);

      part_to_Particle(onepart, part);

      r2 = 0.0;
      for (int j=0; j<3; j++) r2 += part.pos[j]*part.pos[j];
      rmax1 = max<double>(r2, rmax1);

				// Load the particle
      particles.push_back(part);
    }

    nbodies = ncount[0];


				// Now load the other nodes
    int icount, ibufcount;
    for (int n=1; n<numprocs; n++) {

      cout << "Loading node <" << n << ">\n";

      MPI_Send(&ncount[n], 1, MPI_INT, n, 1, MPI_COMM_WORLD);

      icount = 0;
      ibufcount = 0;
      while (icount < ncount[n]) {

	get_next_particle_from_file(&buf[ibufcount], in);

	r2 = 0.0;
	for (int k=0; k<3; k++) 
	  r2 += buf[ibufcount].pos[k]*buf[ibufcount].pos[k];

	rmax1 = max<double>(r2, rmax1);

	/*
	cout << "#=" << icount << " r2=" << r2 << " rmax1=" << rmax1 
	     << " mass=" << buf[ibufcount].mass << endl;
	*/

	ibufcount++;
	icount++;

	if (ibufcount == nbuf || icount == ncount[n]) {
	  MPI_Send(&ibufcount, 1, MPI_INT, n, 2, MPI_COMM_WORLD);
	  MPI_Send(buf, ibufcount, Particletype, n, 3, MPI_COMM_WORLD);
	  ibufcount = 0;
	}

      }

    }

  } else {

    MPI_Recv(&nbodies, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status0);
      
    int icount = 0, isent;
    while (icount < nbodies) {
      MPI_Recv(&isent, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &status0);
      MPI_Recv(buf, isent, Particletype, 0, 3, MPI_COMM_WORLD, &status0);
      for (int i=0; i<isent; i++) {
	part_to_Particle(buf[i], part);
	particles.push_back(part);
	icount++;
      }
    }
  }


				// Sanity check
  if (seq_check) {
    if (seq_beg != particles[0].indx ||
	seq_end != particles[nbodies-1].indx) {
      cout << "Process " << myid << ": sequence error on init,"
	   << " seq_beg=" << seq_beg
	   << " seq_end=" << seq_end
	   << " seq[1]=" << particles[0].indx
	   << " seq[N]=" << particles[nbodies-1].indx
	   << endl << flush;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

				// Default: set to max radius
  rmax = sqrt(fabs(rmax1));
  MPI_Bcast(&rmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

}


struct Component::Partstruct * Component::get_particles(int* number)
{
  if (ordered) return get_particles_ordered(number);
  else         return get_particles_unordered(number);
}

struct Component::Partstruct * Component::get_particles_ordered(int* number)
{
  MPI_Status status;

  static int counter = 0;
  static int node = 0;

  if (*number < 0) {
    counter = 0;
    node = 0;
  }

  if (node==numprocs) {
    *number = 0;
    return 0;
  }

  int icount, indx, curnode = node;
  if (node == 0) {
    
    if (myid == 0) {
      
      icount = 0;
      while (icount < nbuf && counter < nbodies_index[node]) {

	Particle_to_part(buf[icount], particles[counter]);

	icount++;
	counter++;
      }

      if (counter == nbodies_index[node]) node++;

#ifdef DEBUG
    cout << "Process " << myid 
	 << ": buffered " << icount << " particles from master"
	 << endl << flush;
#endif    

    }

  } else {

    if (myid == 0) {

      MPI_Recv(&icount, 1, MPI_INT, node, 52, MPI_COMM_WORLD, &status);
      MPI_Recv(buf, icount, Particletype, node, 53, MPI_COMM_WORLD, &status);

#ifdef DEBUG
    cout << "Process " << myid 
	 << ": received " << icount << " particles from slave"
	 << endl << flush;
#endif    

    } else if (myid == node) {
      
      icount = 0;
      while (icount < nbuf && counter < nbodies_index[node]) {

	indx = counter - nbodies_index[node-1];
	
				// Pack structure
	Particle_to_part(buf[icount], particles[indx]);

	icount++;
	counter++;
      }

      MPI_Send(&icount, 1, MPI_INT, 0, 52, MPI_COMM_WORLD);
      MPI_Send(buf, icount, Particletype, 0, 53, MPI_COMM_WORLD);
#ifdef DEBUG
    cout << "Process " << myid 
	 << ": sent " << icount << " particles to master"
	 << ", counter value=" << counter
	 << ", nbodies_index=" << nbodies_index[node]
	 << endl << flush;
#endif    

      if (counter == nbodies_index[node]) node++;
      
    }

  }
  
  MPI_Bcast(&counter, 1, MPI_INT, curnode, MPI_COMM_WORLD);
  MPI_Bcast(&node, 1, MPI_INT, curnode, MPI_COMM_WORLD);
  MPI_Bcast(&icount, 1, MPI_INT, curnode, MPI_COMM_WORLD);

				// Return values
  *number = icount;
  return buf;

#ifdef DEBUG
  cout << "Process " << myid 
       << ": received counter=" << counter << " node=" << node
       << " number=" << *number
       << " address=" << buf
       << endl << flush;
#endif    

}


struct Component::Partstruct * Component::get_particles_unordered(int* number)
{
  MPI_Status status;

  static int counter;
  static map<unsigned long, int> mapping;

  if (*number < 0) {
    counter = 1;
				// Remake the particle index map
    mapping.erase(mapping.begin(), mapping.end());
    for (int i=0; i<particles.size(); i++)  mapping[particles[i].indx] = i;
  }

				// Done?
  if (counter==nbodies_tot) {
    *number = 0;
    return 0;
  }

  if (myid == 0) {
				// Owned by the root node?
    if (mapping.find(counter) != mapping.end()) {
      // Pack structure
      Particle_to_part(buf[0], particles[mapping[counter]]);
      // No communication necessary . . .
    } else {			// Receive the particle from a node
      MPI_Recv(buf, 1, Particletype, MPI_ANY_SOURCE, 55, MPI_COMM_WORLD, &status);
    }
  } else {
				// Owned by me?
    if (mapping.find(counter) != mapping.end()) {
				// Pack structure
      Particle_to_part(buf[0], particles[mapping[counter]]);
				// Send to master
      MPI_Send(buf, 1, Particletype, 0, 55, MPI_COMM_WORLD);
    }
  }

				// Next counter . . . 
  counter++;

				// Return values
  *number = 1;
  return buf;
}


void Component::write_binary(ostream* out)
{
  ComponentHeader header;

  if (myid == 0) {

    header.nbod = nbodies_tot;
    header.niatr = niattrib;
    header.ndatr = ndattrib;
  
    ostringstream outs;
    outs << name << " : " << id << " : " << cparam << " : " << fparam;
    strncpy(header.info, outs.str().c_str(), header.ninfochar);

    if (!header.write(out)) {
      cerr << "Component::write_binary: Error writing particle header\n";
      MPI_Abort(MPI_COMM_WORLD, 34);
    }
  }

  bool first = true;
  double mass0, pot0, pv;

  int number = -1;
  Partstruct *p = get_particles(&number);

  while (p) {

    if (myid == 0) {

      for (int k=0; k<number; k++) {
	out->write((char *)&(p[k].mass), sizeof(double));
	if (first) {
	  mass0 = p[k].mass;
	  first = false;
	}

	for (int i=0; i<3; i++) {
	  pv = p[k].pos[i] + com0[i] - comI[i];
	  out->write((char *)&pv, sizeof(double));
	}
	for (int i=0; i<3; i++) {
	  pv = p[k].vel[i] + cov0[i] - covI[i];
	  out->write((char *)&pv, sizeof(double));
	}

	pot0 = p[k].pot + p[k].potext;
	out->write((char *)&pot0, sizeof(double));

	for (int i=0; i<header.niatr; i++) 
	  out->write((char *)&(p[k].iatr[i]), sizeof(int));
	for (int i=0; i<header.ndatr; i++) 
	  out->write((char *)&(p[k].datr[i]), sizeof(double));
      }
    }

    
    p = get_particles(&number);

  }
    
}


void Component::write_ascii(ostream* out, bool accel)
{
  int number = -1;
  Partstruct *p = get_particles(&number);

  while (p) {

    if (myid == 0) {

      for (int k=0; k<number; k++) {
	*out << setw(18) << p[k].mass;
	for (int i=0; i<3; i++) *out << setw(18) << p[k].pos[i]+com0[i]-comI[i];
	for (int i=0; i<3; i++) *out << setw(18) << p[k].vel[i]+cov0[i]-covI[i];
	if (accel)
	  for (int i=0; i<3; i++) *out << setw(18) << p[k].acc[i];

	*out << setw(18) << p[k].pot;
	*out << setw(18) << p[k].potext;
	  
	for (int i=0; i<niattrib; i++) 
	  *out << setw(10) << p[k].iatr[i];
	for (int i=0; i<ndattrib; i++) 
	  *out << setw(18) << p[k].datr[i];

	*out << endl;
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
  
  
  vector<Particle>::iterator p, pend;

				// Zero stuff out
  mtot0 = mtot1 = 0.0;
  for (int k=0; k<dim; k++) com1[k] = cov1[k] = 0.0;

				// Particle loop
  pend = particles.end();
  for (p=particles.begin(); p != pend; p++) {
    
    mtot1 += p->mass;

    for (int k=0; k<dim; k++) com1[k] += p->mass*p->pos[k];
    for (int k=0; k<dim; k++) cov1[k] += p->mass*p->vel[k];
    
  }
  
  MPI_Allreduce(&mtot1, &mtot0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(com1, com0, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(cov1, cov0, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  if (mtot0 > 0.0) {
    for (int k=0; k<dim; k++) com0[k] /= mtot0;
    for (int k=0; k<dim; k++) cov0[k] /= mtot0;
  }

  for (int k=0; k<dim; k++) {
    comI[k] = com0[k];
    covI[k] = cov0[k];
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
    vector<Particle>::iterator p, pend = particles.end();
    for (p=particles.begin(); p != pend; p++) {

      for (int i=0; i<3; i++) {
	p->pos[i] -= com0[i] - comI[i];
	p->vel[i] -= cov0[i] - covI[i];
      }
    }

  }

}

void Component::fix_positions(void)
{
  double mtot1;
  double *com1 = new double [3];
  double *cov1 = new double [3];

  
  vector<Particle>::iterator p, pend;

				// Zero stuff out
  mtot = mtot1 = 0.0;
  for (int k=0; k<dim; k++)
    com[k] = cov[k] = com1[k] = cov1[k] = 0.0;

  for (int i=0; i<3; i++) center[i] = 0.0;


				// Particle loop
  pend = particles.end();
  for (p=particles.begin(); p != pend; p++) {
    
    if (escape_com(*p)) {
      if (consp) {		// Set flag indicating escaped particle
	if (p->iattrib[tidal]==0) {
	  p->iattrib[tidal] = 1;
	  if (com_system) {	// Conserve momentum of center of mass
	    for (int i=0; i<3; i++) {
	      covI[i] = (mtot0*covI[i] - p->mass*p->vel[i])/(mtot0 - p->mass);
	      cov0[i] = (mtot0*cov0[i] - p->mass*p->vel[i])/(mtot0 - p->mass);
	    }
	    mtot0 -= p->mass;
	  }
	}
      }
      
    } else {
    
      mtot1 += p->mass;

      for (int k=0; k<dim; k++) com1[k] += p->mass*p->pos[k];
      for (int k=0; k<dim; k++) cov1[k] += p->mass*p->vel[k];
    }
    
  }
  
  MPI_Allreduce(&mtot1, &mtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(com1, com, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(cov1, cov, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    

				// Compute component center of mass and
				// center of velocity, and center of accel

  if (mtot > 0.0) {
    for (int k=0; k<dim; k++) com[k] /= mtot;
    for (int k=0; k<dim; k++) cov[k] /= mtot;
  }

  delete [] com1;
  delete [] cov1;

  if ((EJ & Orient::CENTER) && !EJdryrun) {
    Vector ctr = orient->currentCenter();
    for (int i=0; i<3; i++) center[i] += ctr[i+1];
  }

}


void Component::update_accel(void)
{
  double *acc1 = new double [3];

  
  vector<Particle>::iterator p;

  for (int k=0; k<dim; k++) acc0[k] = acc1[k] = 0.0;

				// Particle loop
  for (p=particles.begin(); p != particles.end(); p++) {
    
    if (escape_com(*p)) continue;
    
    for (int k=0; k<dim; k++) acc1[k] += p->mass*p->acc[k];
  }
  
  MPI_Allreduce(acc1, acc0, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    

				// Compute component center of accel
  if (mtot > 0.0) {
    for (int k=0; k<dim; k++) acc0[k] /= mtot;
  }
  

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

  delete [] acc1;

}


void Component::get_angmom(void)
{
  double *angm1 = new double [3];

  
  vector<Particle>::iterator p, pend;

				// Zero stuff out
  for (int k=0; k<3; k++)
    angmom[k] = angm1[k] = 0.0;
  
				// Particle loop
  pend = particles.end();
  for (p=particles.begin(); p != pend; p++) {
    
    if (freeze(*p)) continue;
    
    angm1[0] += p->mass*(p->pos[1]*p->vel[2]-p->pos[2]*p->vel[1]);
    angm1[1] += p->mass*(p->pos[2]*p->vel[0]-p->pos[0]*p->vel[2]);
    angm1[2] += p->mass*(p->pos[0]*p->vel[1]-p->pos[1]*p->vel[0]);
        
  }
  
  MPI_Allreduce(angm1, angmom, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  delete [] angm1;
}


int Component::round_up(double dnumb)
{
  int numb = (int)(dnumb + 1.0);
  if (numb >= nbodmax) numb = nbodmax;
  return numb;
}


void Component::setup_distribution(void)
{
  ofstream *out;
  ifstream *in;
  int n;

				/* Needed for both master and slaves */

  nbodies_index = vector<int>(numprocs);
  nbodies_table = vector<int>(numprocs);

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
  vector<int> nbodies_index1(numprocs);
  vector<int> nbodies_table1(numprocs);
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
  int bot, ptr, nump, offset=0;
  vector<int> loc(numprocs, 0);

				// Beginning particle index
  if (myid) bot = nbodies_index[myid-1];
  else      bot = 0;

  for (int i=0; i<2*numprocs-2; i++) {

				// Assign new interval
    if (loadb[i].s) inew = loadb[i].indx+1;
    else            iold = loadb[i].indx+1;
    
    if (myid==0 && *log)
      *log << setw(10) << i
	   << setw(10) << loadb[i+1].top - loadb[i].top
	   << setw(10) << iold
	   << setw(10) << inew;
    
				// Relative pointer index of particles
				// to be shifted
    if (myid==iold) {
      ptr = loadb[i].top - bot + offset;
				// Sanity check
      if (ptr<0)
	cerr << "oops in load_balance: iold=" << iold 
	     << "  beg=" << ptr << "\n";
    }

				// Number of particles to be shifted
    nump = loadb[i+1].top - loadb[i].top;

    ostringstream msg;
    
    if (inew==iold || nump==0) 
      msg << "Do nothing";
    else if (inew>iold) {
      msg << "Insert " << nump << " from #" << iold << " to #" << inew;
      insert_particles(iold, inew, ptr, nump, loc[inew]);
      if (myid==iold) {
#ifdef DEBUG
	cout << "Process " << myid << ":"
	     << "  beg seq=" << particles[ptr].iattrib[0] 
	     << "  end seq=" << particles[ptr+nump-1].iattrib[0]
	     << "  erasing=[" << ptr << ", " << ptr+nump-1 << "]" << endl;
#endif
	if (ptr+nump==nbodies)
	  particles.erase(particles.begin()+ptr, particles.end());
	else
	  particles.erase(particles.begin()+ptr, particles.begin()+ptr+nump);
	offset -= nump;
#ifdef DEBUG
	if (particles.size())
	  cout << "Process " << myid << ": new ends :"
	       << "  beg seq=" << particles[0].iattrib[0] 
	       << "  end seq=" << particles[particles.size()-1].iattrib[0] 
	       << endl;
	else
	  cout << "Process " << myid << ": no particles!"
	       << endl;
#endif
      }
				// Update counters for shifted particles
      if (myid==inew) {
	offset += nump;
	loc[inew] += nump;
      }
    }
    else if (iold>inew) {
      msg << "Append " << nump << " from #" << iold << " to #" << inew;
      append_particles(iold, inew, ptr, nump);
      if (myid==iold) {
#ifdef DEBUG
	cout << "Process " << myid << ":"
	     << "  beg seq=" << particles[ptr].iattrib[0] 
	     << "  end seq=" << particles[ptr+nump-1].iattrib[0] 
	     << "  erasing=[" << ptr << ", " << ptr+nump-1 << "]" << endl;
#endif
	if (ptr+nump==nbodies) {
	  particles.erase(particles.begin()+ptr, particles.end());
	}
	else {
	  particles.erase(particles.begin()+ptr, particles.begin()+ptr+nump);
	}
	offset -= nump;
#ifdef DEBUG
	if (particles.size())
	  cout << "Process " << myid << ": new ends :"
	       << "  beg seq=" << particles[0].iattrib[0] 
	       << "  end seq=" << particles[particles.size()-1].iattrib[0] 
	       << endl;
	else
	  cout << "Process " << myid << ": no particles!"
	       << endl;
#endif
      }
    }
    
    if (myid==0 && *log) *log << setw(10) << msg.str() << endl;

  }

  
				// Update indices
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
    for (int i=0; i<nbodies; i++) {
      if (particles[i].iattrib[0] != seq_beg+i) {
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


void Component::append_particles(int from, int to, int begin, int number)
{
  MPI_Status status;

  int icount, indx, counter=0;

  if (myid == from) {
    
    while (counter < number) {

      icount = 0;
      while (icount < nbuf && counter < number) {

	Particle_to_part(buf[icount], particles[begin+counter]);
      
	icount++;
	counter++;
      }
      
      MPI_Send(&icount, 1, MPI_INT, to, 72, MPI_COMM_WORLD);
      MPI_Send(buf, icount, Particletype, to, 73, MPI_COMM_WORLD);
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

      MPI_Recv(&icount, 1, MPI_INT, from, 72, MPI_COMM_WORLD, &status);
      MPI_Recv(buf, icount, Particletype, from, 73, MPI_COMM_WORLD, &status);

#ifdef DEBUG
      cout << "Process " << myid 
	   << ": received " << icount << " particles from Process " << from
	   << " for append" << endl << flush;
#endif    

      for (int j=0; j<icount; j++) {

	part_to_Particle(buf[j], temp);
	particles.push_back(temp);

	counter++;
      }
      
    }

  }

}



void Component::insert_particles(int from, int to, int begin, int number, 
				 int loc)
{
  MPI_Status status;

  int icount, indx, counter=0;

  if (myid == from) {
    
    while (counter < number) {

      icount = 0;
      while (icount < nbuf && counter < number) {

	Particle_to_part(buf[icount], particles[begin+counter]);
      
	icount++;
	counter++;
      }
      
      MPI_Send(&icount, 1, MPI_INT, to, 72, MPI_COMM_WORLD);
      MPI_Send(buf, icount, Particletype, to, 73, MPI_COMM_WORLD);
#ifdef DEBUG
      cout << "Process " << myid 
	   << ": sent " << icount << " particles to Process " << to
	   << " for insert, counter value=" << counter
	   << endl << flush;
#endif    
    }

  }

  if (myid == to) {
  
    Particle temp;
    vector<Particle> newp;

    while (counter < number) {

      MPI_Recv(&icount, 1, MPI_INT, from, 72, MPI_COMM_WORLD, &status);
      MPI_Recv(buf, icount, Particletype, from, 73, MPI_COMM_WORLD, &status);

#ifdef DEBUG
      cout << "Process " << myid 
	   << ": received " << icount << " particles from Process " << from
	   << " for insert" << endl << flush;
#endif    

      for (int j=0; j<icount; j++) {

	part_to_Particle(buf[j], temp);
	newp.push_back(temp);

	counter++;
      }
      
    }

    for (int i=0; i<number; i++) 
      particles.insert(particles.begin()+loc+i, newp[i]);

  }

}


bool Component::freeze(const Particle& p)
{
  double r2 = 0.0;
  for (int i=0; i<3; i++) r2 += 
			    (p.pos[i] - comI[i] - center[i])*
			    (p.pos[i] - comI[i] - center[i]);
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
  MPI_Status status;
  Particle part;
  vector<int>::iterator it = redist.begin();

				// Keep list of particles to delete
  vector<int> myDelete;		// from current node

  int indx, curnode, tonode, lastnode, M, icount;

  while (it != redist.end()) {
    curnode = *(it++);		// Current owner
    M       = *(it++);		// Number to transfer to another node
    if (M) {
      indx   = *(it++);		// Index on owner's list
      tonode = *(it++);		// Destination
      icount = 0;		// Number in transfer buffer so far

      // Do the first particle
      //
      if (myid==curnode) {
	Particle_to_part(buf[icount], particles[indx]);
	myDelete.push_back(indx);
      }
      icount++;
      lastnode = tonode;

      // Do the remaining particles
      //
      for (int m=1; m<M; m++) {
	indx   = *(it++);
	tonode = *(it++);
				// Send the particles
	if (tonode != lastnode && icount) {
	  if (myid==curnode) 
	    MPI_Send(buf, icount, Particletype, lastnode, 53, MPI_COMM_WORLD);
	  if (myid==lastnode) {
	    MPI_Recv(buf, icount, Particletype, curnode , 53, MPI_COMM_WORLD, &status);
	    for (int i=0; i<icount; i++) {
	      part_to_Particle(buf[i], part);
	      particles.push_back(part);
	    }
	  }
	}

				// Buffer this particle
	if (myid==curnode) {
	  Particle_to_part(buf[icount], particles[indx]);
	  myDelete.push_back(indx);
	}
	icount++;
	lastnode = tonode;
	
				// Buffer is full?
	if (icount == nbuf) {
	  if (myid==curnode) 
	    MPI_Send(buf, icount, Particletype, tonode,  53, MPI_COMM_WORLD);
	  if (myid==tonode) {
	    MPI_Recv(buf, icount, Particletype, curnode, 53, MPI_COMM_WORLD, &status);
	    for (int i=0; i<icount; i++) {
	      part_to_Particle(buf[i], part);
	      particles.push_back(part);
	    }
	  }
	  icount = 0;
	}
      } // End of particles on this node
    }
    
  } // Next stanza

  // Delete particles that have been transferred to a new node
  
  for (unsigned int n=myDelete.size()-1; n>=0; n--) {
    particles.erase(particles.begin()+myDelete[n]);
  }

}
