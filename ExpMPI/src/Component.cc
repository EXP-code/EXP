#include <iostream>
#include <strstream>
#include <string>
#include <algorithm>

#include <Component.H>
#include <Bessel.H>
#include <CBrock.H>
#include <CBrockDisk.H>
#include <Hernquist.H>
#include <Sphere.H>
#include <Cylinder.H>
#include <Cube.H>
#include <Slab.H>
#include <SlabSL.H>
#include <Direct.H>
#include <Orient.H>

#include "expand.h"

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

  binary = false;
  npart = false;
  buf = new Partstruct [nbuf];

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

  binary = true;
  npart = false;
  buf = new Partstruct [nbuf];


  read_bodies_and_distribute_binary(in);

}


void Component::initialize(void)
{
				// Parse the parameters
  StringTok<string> tokens(cparam);
  pair<string, string> datum;

  string token = tokens("|");	// Bar separated tokens

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

    if (!datum.first.compare("rmax"))     rmax = atof(datum.second.c_str());

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
    cout << "Process " << myid << ": about to create Orient with"
	 << " nkeep=" << nEJkeep
	 << " nwant=" << nEJwant
	 << " eEJ=" << eEJ0 << endl;
    orient = new Orient(nEJkeep, nEJwant, eEJ0);
    cout << "Process " << myid << ": Orient successful\n";
  }

  if (myid == 0) {		// Flag messages for diagnostics
    
    if (EJ & Orient::AXIS)
      cout << name << ": AXIS orientation is *ON*\n";

    if (EJ & Orient::CENTER) {
      if (use_com) 
	cout << name 
	     << ": CENTER finding is *ON* and will supercede COM centering\n";
      else
	cout << name 
	     << ": CENTER finding is *ON*\n";
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
    istrstream ins(line);
    
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
  
  MPI_Datatype	type[7] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
			   MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_DOUBLE};

				// Get displacements
  MPI_Aint	disp[7];
  MPI_Get_address(&buf[0].mass,		&disp[0]);
  MPI_Get_address(&buf[0].pos,		&disp[1]);
  MPI_Get_address(&buf[0].vel,		&disp[2]);
  MPI_Get_address(&buf[0].pot,		&disp[3]);
  MPI_Get_address(&buf[0].potext,	&disp[4]);
  MPI_Get_address(&buf[0].iatr,		&disp[5]);
  MPI_Get_address(&buf[0].datr,		&disp[6]);

  for (int i=6; i>=0; i--) disp[i] -= disp[0];
  
				// Block offsets
  int		blocklen[7] = {1, 3, 3, 1, 1, nimax, ndmax};
  
  MPI_Type_create_struct(7, blocklen, disp, type, &Particletype);
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
      istrstream ins(line);

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
      istrstream ins(line);

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
	istrstream ins(line);

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
  in->read(&(onepart->mass), sizeof(double));
  for (int i=0; i<3; i++) in->read(&(onepart->pos[i]), sizeof(double));
  for (int i=0; i<3; i++) in->read(&(onepart->vel[i]), sizeof(double));
  in->read(&(onepart->pot), sizeof(double));
  onepart->potext = 0.0;
  for (int i=0; i<niattrib; i++) 
    in->read(&(onepart->iatr[i]), sizeof(int));
  for (int i=0; i<ndattrib; i++) 
    in->read(&(onepart->datr[i]), sizeof(double));
}


int Component::get_next_particle(Partstruct *onepart)
{
				// Initialize read
  if (!npart) {
    npart_tot = -1;
    npart_p = get_particles(&npart_tot);
    cerr << "Component::get_next_particle: no particles found!?!\n";
    if (!npart_tot) return 0;
    npart = true;
    npart_cur = 0;
  }

				// Get a new bunch
  if (npart_cur == npart_tot) {
    npart_p = get_particles(&npart_tot);
    if (!npart_tot) {
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
  fparam = trimLeft(trimRight(tokens(":")));

  delete [] info;

				// Informational output
  if (myid==0)
    cout << setw(60) << setfill('-') << "-" << endl << setfill(' ')
	 << "--- New Component" << endl
	 << setw(20) << " name :: "  << name          << endl
	 << setw(20) << " id :: "    << id            << endl
	 << setw(20) << " param :: " << fparam         << endl
	 << setw(60) << setfill('-') << "-" << endl << setfill(' ');

				// Make MPI datatype
  
  MPI_Datatype	type[7] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
			   MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_DOUBLE};

				// Get displacements
  MPI_Aint	disp[7];
  MPI_Get_address(&buf[0].mass,		&disp[0]);
  MPI_Get_address(&buf[0].pos,		&disp[1]);
  MPI_Get_address(&buf[0].vel,		&disp[2]);
  MPI_Get_address(&buf[0].pot,		&disp[3]);
  MPI_Get_address(&buf[0].potext,	&disp[4]);
  MPI_Get_address(&buf[0].iatr,		&disp[5]);
  MPI_Get_address(&buf[0].datr,		&disp[6]);

  for (int i=6; i>=0; i--) disp[i] -= disp[0];
  
				// Block offsets
  int		blocklen[7] = {1, 3, 3, 1, 1, nimax, ndmax};
  
  MPI_Type_create_struct(7, blocklen, disp, type, &Particletype);
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
  
    ostrstream outs;
    outs << name << " : " << id << " : " << fparam << '\0';
    strncpy(header.info, outs.str(), header.ninfochar);

    if (!header.write(out)) {
      cerr << "Component::write_binary: Error writing particle header\n";
      MPI_Abort(MPI_COMM_WORLD, 34);
    }
  }

  bool first = true;
  double mass0, pot0;

  int number = -1;
  Partstruct *p = get_particles(&number);

  while (number) {

    if (myid == 0) {

      for (int k=0; k<number; k++) {
	out->write(&(p[k].mass), sizeof(double));
	if (first) {
	  mass0 = p[k].mass;
	  first = false;
	}
	for (int i=0; i<3; i++) out->write(&(p[k].pos[i]), sizeof(double));
	for (int i=0; i<3; i++) out->write(&(p[k].vel[i]), sizeof(double));

	pot0 = p[k].pot + p[k].potext;
	out->write(&pot0, sizeof(double));

	for (int i=0; i<header.niatr; i++) 
	  out->write(&(p[k].iatr[i]), sizeof(int));
	for (int i=0; i<header.ndatr; i++) 
	  out->write(&(p[k].datr[i]), sizeof(double));
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
    
    if (p->freeze()) continue;
    
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

      if (p->freeze()) continue;
	
      for (int k=0; k<dim; k++) p->vel[k] -= cov[k];
    }
  }

  delete [] com1;
  delete [] cov1;

  Vector ctr;
  if (EJ & Orient::CENTER) {
    ctr = orient->currentCenter();
    for (int i=0; i<3; i++) center[i] = ctr[i+1];
  }
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
  double norm=0.0;

				/* Needed for both master and slaves */

  nbodies_index = vector<int>(numprocs);
  nbodies_table = vector<int>(numprocs);

  if (myid == 0) {
    in = new ifstream("processor.rates");
    if (!*in) {
      cerr << "setup: Error opening <processor.rates> . . . will assume homogeneous cluster\n";
    }

    rates =  vector<double>(numprocs);
    orates = vector<double>(numprocs);
    trates = vector<double>(numprocs);

    if (*in) {			// We are reading from a file
      for (n=0; n<numprocs; n++) {
	*in >> rates[n];
	if (!*in) {
	  cerr << "Error reading <processor.rates>\n";
	  MPI_Abort(MPI_COMM_WORLD, 33);
	  exit(0);
	}
	norm += rates[n];
      }

      delete in;
    
    } else {			// Assign equal rates to all nodes
      for (n=0; n<numprocs; n++) {
	rates[n] = 1.0;
	norm += rates[n];
      }
    }

    for (n=0; n<numprocs; n++) {

      rates[n] /= norm;
      
      if (n == 0)
	nbodies_table[n] = nbodies_index[n] = 
	  max<int>(1, min<int>((int)(rates[n] * nbodies_tot), nbodies_tot));
      else {
	if (n < numprocs-1)
	  nbodies_index[n] = (int)(rates[n] * nbodies_tot) + 
	    nbodies_index[n-1];
	else
	  nbodies_index[n] = nbodies_tot;
      
	nbodies_table[n] = nbodies_index[n] - nbodies_index[n-1];
      }

    }

    out = new ofstream("current.processor.rates", ios::out | ios::app);
    if (out) {
      *out << "# " << endl;
      *out << "# Time=" << tnow << " Component=" << name << endl;
      *out << "# " 
	  << setw(15) << "Norm rate"
	  << setw(15) << "Raw rate"
	  << setw(15) << "Delta rate"
	  << setw(15) << "Index"
	  << setw(15) << "Current #"
	  << endl
	  << "# "
	  << setw(15) << "---------"
	  << setw(15) << "--------"
	  << setw(15) << "----------"
	  << setw(15) << "--------"
	  << setw(15) << "---------"
	  << endl;
      
      for (n=0; n<numprocs; n++)
	*out << "  "
	    << setw(15) << rates[n]
	    << setw(15) << rates[n]*norm
	    << setw(15) << 1.0 - rates[n]*nbodies_tot/nbodies_table[n]
	    << setw(15) << nbodies_index[n]
	    << setw(15) << nbodies_table[n]
	    << endl;

      delete out;
    }

  }


  MPI_Bcast(&nbodies_index[0], numprocs, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nbodies_table[0], numprocs, MPI_INT, 0, MPI_COMM_WORLD);

}

