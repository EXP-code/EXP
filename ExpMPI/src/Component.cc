#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>

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
#include <Orient.H>

#include "expand.h"

bool less_loadb(const loadb_datum& one, const loadb_datum& two)
{
  return (one.top < two.top);
}

Component::Component(string NAME, string ID, string CPARAM, string PFILE, 
		     string FPARAM) : 
  name(NAME), id(ID), cparam(CPARAM), pfile(PFILE), fparam(FPARAM)
{
  use_com = false;
  use_cov = false;

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

  binary = false;
  npart = false;
  buf = new Partstruct [nbuf];

  adiabatic = false;
  ton  = -1.0e20;
  toff =  1.0e20;
  twid = 0.1;

  rtrunc = 1.0e20;

  read_bodies_and_distribute_ascii();
}


Component::Component(istream *in)
{
  use_com = false;
  use_cov = false;

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

  binary = true;
  npart = false;
  buf = new Partstruct [nbuf];


  adiabatic = false;
  ton  = -1.0e20;
  toff =  1.0e20;
  twid = 0.1;

  rtrunc = 1.0e20;

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

    if (!datum.first.compare("use_com"))  use_com = true;

    if (!datum.first.compare("use_cov"))  use_cov = true;

    if (!datum.first.compare("com_tie"))  get_com_component(datum.second);

    if (!datum.first.compare("EJ"))       EJ = atoi(datum.second.c_str());
    
    if (!datum.first.compare("nEJkeep"))  nEJkeep = atoi(datum.second.c_str());

    if (!datum.first.compare("nEJwant"))  nEJwant = atoi(datum.second.c_str());

    if (!datum.first.compare("eEJ0"))     eEJ0 = atof(datum.second.c_str());

    if (!datum.first.compare("EJx0"))     EJx0 = atof(datum.second.c_str());

    if (!datum.first.compare("EJy0"))     EJy0 = atof(datum.second.c_str());

    if (!datum.first.compare("EJz0"))     EJz0 = atof(datum.second.c_str());

    if (!datum.first.compare("EJkinE"))   EJkinE = true ? atoi(datum.second.c_str()) : false;

    if (!datum.first.compare("EJext"))    EJext = true ? atoi(datum.second.c_str()) : false;

    if (!datum.first.compare("EJdiag"))   EJdiag = true ? atoi(datum.second.c_str()) : false;

    if (!datum.first.compare("EJdryrun")) EJdryrun = true ? atoi(datum.second.c_str()) : false;


    if (!datum.first.compare("rmax"))     rmax = atof(datum.second.c_str());

    if (!datum.first.compare("ton"))      {ton = atof(datum.second.c_str()); adiabatic = true;}

    if (!datum.first.compare("toff"))     {toff= atof(datum.second.c_str()); adiabatic = true;}

    if (!datum.first.compare("twid"))     {twid = atof(datum.second.c_str()); adiabatic = true;}

    if (!datum.first.compare("rtrunc"))   {rtrunc = atof(datum.second.c_str()); adiabatic = true;}

				// Next parameter
    token = tokens(",");
  }


				// DEBUG
  if (myid==0) {
    if (use_com) cout << name << ": self com is *ON*\n";
    list<Component*>::iterator i;
    for (i=com_tie.begin(); i != com_tie.end(); i++) {
      cout << name << ": use <" << (*i)->name << "> for com\n";
    }
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
  else {
    string msg("I don't know about the force: ");
    msg += id;
    bomb(msg);
  }

  force->RegisterComponent(this);

  com = new double [3];
  center = new double [3];
  cov = new double [3];
  angmom = new double [3];
  ps = new double [6];

  for (int k=0; k<3; k++) com[k] = center[k] = cov[k] = angmom[k] = 0.0;

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
      
    
    string EJlogfile = name + ".orient." + runtag;

    unsigned EJctl = 0;
    if (EJdiag)		EJctl |= Orient::DIAG;
    if (EJkinE)		EJctl |= Orient::KE;
    if (EJext)		EJctl |= Orient::EXTERNAL;

    orient = new Orient(nEJkeep, nEJwant, eEJ0, EJ, EJctl, EJlogfile);
    orient -> set_center(EJx0, EJy0, EJz0);
    center[0] = EJx0;
    center[1] = EJy0;
    center[2] = EJz0;

    if (EJdiag) cout << "Process " << myid << ": Orient successful\n";
  }

  if (myid == 0) {		// Flag messages for diagnostics
    
    if (EJ & Orient::AXIS)
      cout << name << ": AXIS orientation is *ON*\n";

    if (EJ & Orient::CENTER) {
      if (use_com) 
	cout << name 
	     << ": CENTER finding is *ON* and will supercede COM centering,";
      else
	cout << name 
	     << ": CENTER finding is *ON*,";

      cout << ", user specified initial center: x, y, z: " 
	   << EJx0 << ", " 
	   << EJy0 << ", " 
	   << EJz0 << "\n";
    }

  }
    
}


void Component::get_com_component(string name)
{
				// Loop through components to find com target
  list<Component*>::iterator cc;
  Component *c;
  bool found = false;
  
  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    c = *cc;

				// Is this the one?
    if (c->name.compare(name) == 0) {
      com_tie.push_back(&(*c));
      found = true;
      break;
    }
  }

  if (!found) {
    string msg = "Could not find target component named <" + name +  "> for COM\n";
    bomb(msg);
  }

}


Component::~Component(void)
{
  delete force;
  delete [] buf;
  delete [] com;
  delete [] cov;
  delete [] center;
  delete [] angmom;
  delete [] ps;
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
  }
#ifdef SEQCHECK
  niattrib++;
#endif
				// Check buffer dimensions
  if (myid==0) {

    if (niattrib>nimax) {
      cerr << "Too many integer attributes: redimension nimax and recompile\n";
      MPI_Abort(MPI_COMM_WORLD, 12);
      exit(-1);
    }

    if (ndattrib>ndmax) {
      cerr << "Too many real # attributes: redimension nimax and recompile\n";
      MPI_Abort(MPI_COMM_WORLD, 13);
      exit(-1);
    }

  }
				// Broadcast attributes for this
				// phase-space component
  MPI_Bcast(&nbodies_tot, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&niattrib, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ndattrib, 1, MPI_INT, 0, MPI_COMM_WORLD);


				// Make MPI datatype
  
  MPI_Datatype type[8] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
			  MPI_DOUBLE, MPI_DOUBLE, MPI_INT,    MPI_DOUBLE};

				// Get displacements
  MPI_Aint disp[8];
  MPI_Get_address(&buf[0].mass,		&disp[0]);
  MPI_Get_address(&buf[0].pos,		&disp[1]);
  MPI_Get_address(&buf[0].vel,		&disp[2]);
  MPI_Get_address(&buf[0].acc,		&disp[3]);
  MPI_Get_address(&buf[0].pot,		&disp[4]);
  MPI_Get_address(&buf[0].potext,	&disp[5]);
  MPI_Get_address(&buf[0].iatr,		&disp[6]);
  MPI_Get_address(&buf[0].datr,		&disp[7]);

  for (int i=7; i>=0; i--) disp[i] -= disp[0];
  
				// Block offsets
  int blocklen[8] = {1, 3, 3, 3, 1, 1, nimax, ndmax};
  
  MPI_Type_create_struct(8, blocklen, disp, type, &Particletype);
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
      part.pot = part.potext = 0.0;

#ifdef SEQCHECK
      for (int j=0; j<niattrib-1; j++) ins >> part.iattrib[j];
      part.iattrib[niattrib-1] = 1;
#else
      for (int j=0; j<niattrib; j++) ins >> part.iattrib[j];
#endif
      for (int j=0; j<ndattrib; j++) ins >> part.dattrib[j];

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

#ifdef SEQCHECK
      for (int j=0; j<niattrib-1; j++) ins >> part.iattrib[j];
      part.iattrib[niattrib-1] = i;
#else
      for (int j=0; j<niattrib; j++) ins >> part.iattrib[j];
#endif
      for (int j=0; j<ndattrib; j++) ins >> part.dattrib[j];

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
	buf[ibufcount].pot = buf[ibufcount].potext = 0.0;

#ifdef SEQCHECK
	for (int j=0; j<niattrib-1; j++) ins >> buf[ibufcount].iatr[j];
	buf[ibufcount].iatr[niattrib-1] = nbodies_index[n-1] + 1 + icount;
#else
	for (int j=0; j<niattrib; j++) ins >> buf[ibufcount].iatr[j];
#endif
	for (int j=0; j<ndattrib; j++) ins >> buf[ibufcount].datr[j];

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
  
#ifdef SEQCHECK			// Sanity check
  if (particles.size()) {
    if (seq_beg != particles[0].iattrib[0] || 
	seq_end != particles[nbodies-1].iattrib[0]) {
      cout << "Process " << myid << ": sequence error on init,"
	   << " seq_beg=" << seq_beg
	   << " seq_end=" << seq_end
	   << " seq[1]=" << particles[0].iattrib[0]
	   << " seq[N]=" << particles[nbodies-1].iattrib[0]
	   << " nbodies=" << nbodies
	   << endl << flush;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }
#endif

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
				// Check buffer dimensions

    if (header.niatr>nimax) {
      cerr << "Too many integer attributes: " << header.niatr 
	   << ", redimension nimax and recompile\n";
      MPI_Abort(MPI_COMM_WORLD, 12);
      exit(-1);
    }

    if (header.ndatr>ndmax) {
      cerr << "Too many real # attributes: " << header.ndatr
	   << ", redimension ndmax and recompile\n";
      MPI_Abort(MPI_COMM_WORLD, 13);
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
  
  MPI_Datatype	type[8] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
			   MPI_DOUBLE, MPI_DOUBLE, MPI_INT,    MPI_DOUBLE};

				// Get displacements
  MPI_Aint disp[8];
  MPI_Get_address(&buf[0].mass,		&disp[0]);
  MPI_Get_address(&buf[0].pos,		&disp[1]);
  MPI_Get_address(&buf[0].vel,		&disp[2]);
  MPI_Get_address(&buf[0].acc,		&disp[3]);
  MPI_Get_address(&buf[0].pot,		&disp[4]);
  MPI_Get_address(&buf[0].potext,	&disp[5]);
  MPI_Get_address(&buf[0].iatr,		&disp[6]);
  MPI_Get_address(&buf[0].datr,		&disp[7]);

  for (int i=7; i>=0; i--) disp[i] -= disp[0];
  
				// Block offsets
  int blocklen[8] = {1, 3, 3, 3, 1, 1, nimax, ndmax};
  
  MPI_Type_create_struct(8, blocklen, disp, type, &Particletype);
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


#ifdef SEQCHECK			// Sanity check
  if (seq_beg != particles[0].iattrib[0] || 
      seq_end != particles[nbodies-1].iattrib[0]) {
    cout << "Process " << myid << ": sequence error on init,"
	 << " seq_beg=" << seq_beg
	 << " seq_end=" << seq_end
	 << " seq[1]=" << particles[0].iattrib[0]
	 << " seq[N]=" << particles[nbodies-1].iattrib[0]
	 << endl << flush;
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
#endif
				// Default: set to max radius
  rmax = sqrt(fabs(rmax1));
  MPI_Bcast(&rmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

}

struct Component::Partstruct * Component::get_particles(int* number)
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


void Component::write_binary(ostream* out)
{
  ComponentHeader header;

  if (myid == 0) {

    header.nbod = nbodies_tot;
    header.niatr = niattrib;
    header.ndatr = ndattrib;
  
    ostringstream outs;
    outs << name << " : " << id << " : " << cparam << " : " << fparam << '\0';
    strncpy(header.info, outs.str().c_str(), header.ninfochar);

    if (!header.write(out)) {
      cerr << "Component::write_binary: Error writing particle header\n";
      MPI_Abort(MPI_COMM_WORLD, 34);
    }
  }

  bool first = true;
  double mass0, pot0;

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
	for (int i=0; i<3; i++) out->write((char *)&(p[k].pos[i]), sizeof(double));
	for (int i=0; i<3; i++) out->write((char *)&(p[k].vel[i]), sizeof(double));

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
	for (int i=0; i<3; i++) *out << setw(18) << p[k].pos[i];
	for (int i=0; i<3; i++) *out << setw(18) << p[k].vel[i];
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


void Component::fix_positions(void)
{
  double mtot1;
  double *com1 = new double [3];
  double *cov1 = new double [3];

  
  vector<Particle>::iterator p, pend;

				// Zero stuff out
  mtot = mtot1 = 0.0;
  for (int k=0; k<dim; k++)
    com[k] = center[k] = cov[k] = com1[k] = cov1[k] = 0.0;

				// Particle loop
  pend = particles.end();
  for (p=particles.begin(); p != pend; p++) {
    
    if (freeze(*p)) continue;
    
    mtot1 += p->mass;

    for (int k=0; k<dim; k++) com1[k] += p->mass*p->pos[k];
    for (int k=0; k<dim; k++) cov1[k] += p->mass*p->vel[k];
    
  }
  
  MPI_Allreduce(&mtot1, &mtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(com1, com, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(cov1, cov, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
				// Add com from other specified components
				// to center

  double mtot0 = 0.0;

  if (use_com) {
    for (int k=0; k<dim; k++) center[k] = com[k];
    mtot0 += mtot;
  }
  
  list<Component*>::iterator i;
  for (i=com_tie.begin(); i != com_tie.end(); i++) {
    for (int k=0; k<dim; k++) center[k] += (*i)->com[k]*(*i)->mtot;
    mtot0 += (*i)->mtot;
  }

				// Normalize center
  if (mtot0>0.0)
    for (int k=0; k<dim; k++) center[k] /= mtot0;


				// Compute component center of mass and
				// center of velocity

  if (mtot > 0.0) {
    for (int k=0; k<dim; k++) com[k] /= mtot;
    for (int k=0; k<dim; k++) cov[k] /= mtot;
  }

  if (use_cov) {

    pend = particles.end();
    for (p=particles.begin(); p != pend; p++) {

      if (freeze(*p)) continue;
	
      for (int k=0; k<dim; k++) p->vel[k] -= cov[k];
    }
  }

  delete [] com1;
  delete [] cov1;

  Vector ctr;
  if ((EJ & Orient::CENTER) && !EJdryrun) {
    ctr = orient->currentCenter();
    for (int i=0; i<3; i++) center[i] = ctr[i+1];
  }
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

    string outrates = "current.processor.rates." + runtag;

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

    string outrates = "current.processor.rates." + runtag;
    string rateslog = "current.processor.rates.log." + runtag;

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
  

#ifdef SEQCHECK
  char msgbuf[200];		// Only need 31 characters . . .

  if (myid==0) *log << endl << "Post load-balance sequence check:" 
		    << endl << endl;

  for (int i=0; i<numprocs; i++) {
    if (myid==i) {
      ostringstream msg;
      msg << "Process " << setw(4) << myid << ":"
	  << setw(9) << particles[0].iattrib[0]
	  << setw(9) << particles[nbodies-1].iattrib[0];
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

#endif

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

  for (int i=0; i<3; i++) r2 += (p.pos[i] - center[i])*(p.pos[i] - center[i]);
  if (r2 > rtrunc*rtrunc) return true;
  else return false;
}

double Component::Adiabatic()
{
  if (!adiabatic) return 1.0;
  return 0.25*
    ( 1.0 + erf((tpos - ton )/twid) ) *
    ( 1.0 + erf((toff - tpos)/twid) ) ;
}
