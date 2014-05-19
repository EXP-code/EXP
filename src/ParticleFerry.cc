#include <iostream>
#include <iomanip>
#include "global.H"
#include "ParticleFerry.H"
#include "pHOT.H"

//#define DEBUG 1

Partstruct::Partstruct()
{
  // Start with all fields initialized
  mass = 0.0;
  for (int i=0; i<3; i++) pos[i] = vel[i] = acc[i] = 0.0;
  pot = potext = 0.0;
  dtreq = scale = effort = 0.0;
  level = 0;
  indx  = 0;
  tree  = 0u;
  key   = 0u;
  nicnt = 0;
  ndcnt = 0;
  for (int i=0; i<nimax; i++) iatr[i] = 0;
  for (int i=0; i<ndmax; i++) datr[i] = 0.0;
}


ParticleFerry::ParticleFerry()
{
				// Assign particle structure buffer
  buf = new Partstruct [PFbufsz];
  ibufcount = 0;

  const int nf = 17;		// Number of fields

				// Make MPI datatype
#ifdef I128
  MPI_Datatype type[nf] = 
    {
      MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 3 (1)
      MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 3 (6)
      MPI_FLOAT,  MPI_FLOAT,  MPI_FLOAT,  // 3 (9)
      MPI_UNSIGNED, MPI_UNSIGNED_LONG,    // 2 (11)
      MPI_UNSIGNED, MPI_EXP_KEYTYPE,      // 2 (13)
      MPI_UNSIGNED, MPI_UNSIGNED,         // 2 (15) 
      MPI_INT, MPI_DOUBLE                 // 2 (17)
  };
#else
  MPI_Datatype type[nf] = 
    {
      MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 3 (1)
      MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, // 3 (6)
      MPI_FLOAT,  MPI_FLOAT,  MPI_FLOAT,  // 3 (9)
      MPI_UNSIGNED, MPI_UNSIGNED_LONG,    // 2 (11)
      MPI_UNSIGNED, MPI_UNSIGNED_LONG,    // 2 (13)
      MPI_UNSIGNED, MPI_UNSIGNED,         // 2 (15) 
      MPI_INT, MPI_DOUBLE                 // 2 (17)
    };
#endif

				// Get displacements
  MPI_Aint disp[nf];
  MPI_Get_address(&buf[0].mass,		&disp[0]);
  MPI_Get_address(&buf[0].pos,		&disp[1]);
  MPI_Get_address(&buf[0].vel,		&disp[2]);
  MPI_Get_address(&buf[0].acc,		&disp[3]);
  MPI_Get_address(&buf[0].pot,		&disp[4]);
  MPI_Get_address(&buf[0].potext,	&disp[5]);
  MPI_Get_address(&buf[0].dtreq,	&disp[6]);
  MPI_Get_address(&buf[0].scale,	&disp[7]);
  MPI_Get_address(&buf[0].effort,	&disp[8]);  
  MPI_Get_address(&buf[0].level,	&disp[9]);
  MPI_Get_address(&buf[0].indx,		&disp[10]);
  MPI_Get_address(&buf[0].tree,		&disp[11]);
  MPI_Get_address(&buf[0].key,		&disp[12]);
  MPI_Get_address(&buf[0].nicnt,	&disp[13]);
  MPI_Get_address(&buf[0].ndcnt,	&disp[14]);
  MPI_Get_address(&buf[0].iatr,		&disp[15]);
  MPI_Get_address(&buf[0].datr,		&disp[16]);
  
  for (int i=nf-1; i>=0; i--) disp[i] -= disp[0];
  
				// Block offsets
				// 
  int blocklen[nf] = {1, 3, 3, 3, 1, 1, // Doubles
		      1, 1, 1,	        // Floats
		      1, 1,	        // Uint, Ulong
		      1, 1,	        // Uint, U2long
		      1, 1,	        // Uint, Uint
		      nimax, ndmax};	// Uint, Double
  
  MPI_Type_create_struct(nf, blocklen, disp, type, &Particletype);
  MPI_Type_commit(&Particletype);


  pk_lo = 1u << (3*pkbits);
  pk_hi = 1u << (3*pkbits+1);

  key_lo = 1u;
  key_lo <<= (3*nbits);
  key_hi = 1u;
  key_hi <<= (3*nbits+1);

}

ParticleFerry::~ParticleFerry()
{
  delete [] buf;
}

void ParticleFerry::part_to_Particle(Partstruct& str, Particle& cls)
{
  cls.mass = str.mass;
  if(str.mass == 0 or cls.mass == 0) {
	cout << "Error in ptP" << endl;
  }
  for (int j=0; j<3; j++) {
      cls.pos[j] = str.pos[j];
      cls.vel[j] = str.vel[j];
      cls.acc[j] = str.acc[j];
  }
  cls.pot    = str.pot;
  cls.potext = str.potext;

  cls.dtreq  = str.dtreq;
  cls.scale  = str.scale;
  cls.effort = str.effort;
  cls.level  = str.level;
  cls.indx   = str.indx;
  cls.tree   = str.tree;
  cls.key    = str.key;

  cls.iattrib = vector<int>(str.nicnt);
  for (unsigned j=0; j<str.nicnt; j++) cls.iattrib[j] = str.iatr[j];
  
  cls.dattrib = vector<double>(str.ndcnt);
  for (unsigned j=0; j<str.ndcnt; j++) cls.dattrib[j] = str.datr[j];

  if (cls.tree > 0) {
    if ( (cls.tree < pk_lo) || (cls.tree >= pk_hi) ) {
      cout << "Error!! [4], id=" << myid 
	   << ": tree_0=" << str.tree << " tree_1=" << cls.tree
	   << " seq=" << cls.indx
	   << " (x, y, z)={" << cls.pos[0] << ", " << cls.pos[1]
	   << ", " << cls.pos[2]
	   << endl;
    }
  }
}

void ParticleFerry::Particle_to_part(Partstruct& str, Particle& cls)
{

  str.mass = cls.mass;
  if(cls.mass == 0 or cls.indx == 0 or isnan(cls.mass)) {
	cout << "Error in Ptp" << endl;
  }
  for (int j=0; j<3; j++) {
      str.pos[j] = cls.pos[j];
      str.vel[j] = cls.vel[j];
      str.acc[j] = cls.acc[j];
  }
  str.pot    = cls.pot;
  str.potext = cls.potext;

  str.dtreq  = cls.dtreq;
  str.scale  = cls.scale;
  str.effort = cls.effort;
  str.level  = cls.level;
  str.indx   = cls.indx;
  str.tree   = cls.tree;
  str.key    = cls.key;

  str.nicnt  = min<int>(nimax, cls.iattrib.size());
  str.ndcnt  = min<int>(ndmax, cls.dattrib.size());

  for (unsigned j=0; j<str.nicnt; j++) str.iatr[j] = cls.iattrib[j];

  for (unsigned j=0; j<str.ndcnt; j++) str.datr[j] = cls.dattrib[j];

  if (str.tree > 0) {
    if ( (str.tree < pk_lo) || (str.tree >= pk_hi) ) {
      cout << "Error!! [5], id=" << myid 
	   << ": tree_0=" << cls.tree << " tree_1=" << str.tree
	   << " seq=" << str.indx
	   << " (x, y, z)={" << cls.pos[0] << ", " << cls.pos[1]
	   << ", " << cls.pos[2]
	   << endl;
    }
  }

}


void ParticleFerry::ShipParticles(unsigned to, unsigned from, unsigned& total)
{
  MPI_Status status;

  _to    = to;
  _from  = from;
  _total = total;

  if (_from == myid) {
    MPI_Send(&_total, 1, MPI_UNSIGNED, _to, 29, MPI_COMM_WORLD);
    ibufcount = 0;
    itotcount = 0;
  }
  
  if (_to == myid) {
    MPI_Recv(&_total, 1, MPI_UNSIGNED, _from, 29, MPI_COMM_WORLD, &status);
    ibufcount = 0;
    itotcount = 0;
    total = _total;
  }
}

#ifdef I128
void ParticleFerry::SendParticle(Particle& ptc, unsigned seq, uint128 key)
#else
void ParticleFerry::SendParticle(Particle& ptc, unsigned seq, unsigned long key)
#endif
{
  // Add particle to buffer
  //
  Particle_to_part(buf[ibufcount], ptc);
  buf[ibufcount].indx = seq;
  buf[ibufcount].key  = key;

  // If buffer is full, send the buffer and reset
  //
  ibufcount++;
  itotcount++;
  if (ibufcount == PFbufsz || itotcount == _total) BufferSend();
}

void ParticleFerry::SendParticle(Particle& part)
{
  // Add particle to buffer
  //
  Particle_to_part(buf[ibufcount], part);
  ibufcount++;
  itotcount++;

  // If buffer is full, send the buffer and reset
  //
  if (ibufcount == PFbufsz || itotcount == _total) BufferSend();
}

#ifdef I128
bool ParticleFerry::RecvParticle(Particle& ptc, unsigned& seq, uint128& key)
#else
bool ParticleFerry::RecvParticle(Particle& ptc, unsigned& seq, unsigned long& key)
#endif
{
  if (itotcount++ == _total) return false;

  if (ibufcount==0) BufferRecv();

  part_to_Particle(buf[ibufcount], ptc);
  seq = buf[ibufcount].indx;
  key = buf[ibufcount].key;

  ibufcount--;

  return true;
}

bool ParticleFerry::RecvParticle(Particle& part)
{
  if (itotcount++ == _total) return false;
  if (ibufcount==0) BufferRecv();
  part_to_Particle(buf[--ibufcount], part);
  if (part.indx==0 || part.mass<=0.0 || isnan(part.mass)) {
	cout << "BAD MASS!" << endl;
  }
#ifdef DEBUG
  if (part.indx==0) {
    cout << "ParticleFerry: process " << myid << " error in sequence" << endl;
  }
#endif
  return true;
}

void ParticleFerry::BufferSend()
{
  MPI_Send(&ibufcount, 1,         MPI_INT,      _to, 2, MPI_COMM_WORLD);
  MPI_Send(buf,        ibufcount, Particletype, _to, 3, MPI_COMM_WORLD);
#ifdef DEBUG
  cout << "ParticleFerry: process " << myid  << " send, tot=" << itotcount << endl;
  bufferKeyCheck();
#endif
  ibufcount = 0;
}

void ParticleFerry::BufferRecv()
{
  MPI_Status s;

  MPI_Recv(&ibufcount, 1,         MPI_INT,      _from, 2, MPI_COMM_WORLD, &s);
  MPI_Recv(buf,        ibufcount, Particletype, _from, 3, MPI_COMM_WORLD, &s);
#ifdef DEBUG
  cout << "ParticleFerry: process " << myid  << " recv, tot=" << itotcount-1+ibufcount << endl;
  bufferKeyCheck();
#endif
}


void ParticleFerry::bufferKeyCheck()
{
				// Sanity check for pHOT keys
#ifdef I128
  uint128 minkey = 1u;
  minkey <<= 128 - 1;
  uint128 maxkey = 0u;
  unsigned err0 = 0;
  for (unsigned n=0; n<ibufcount; n++) {
    minkey = min<uint128>(minkey, buf[n].key);
    maxkey = max<uint128>(maxkey, buf[n].key);
    if ((buf[n].key < key_lo || buf[n].key >= key_hi) && buf[n].key > 0u) err0++;
  }
#else
  unsigned long minkey = 1u;
  minkey <<= sizeof(unsigned long)*8 - 1;
  unsigned long maxkey = 0u;
  unsigned err0 = 0;
  for (unsigned n=0; n<ibufcount; n++) {
    minkey = min<unsigned long>(minkey, buf[n].key);
    maxkey = max<unsigned long>(maxkey, buf[n].key);
    if ((buf[n].key < key_lo || buf[n].key >= key_hi) && buf[n].key > 0u) err0++;
  }
#endif

  unsigned maxpexp = 1u;
  maxpexp <<= (3*pkbits);
  unsigned minpkey = 1u;
  minpkey <<= (32 - 1);
  unsigned maxpkey = 0u;
  unsigned err1 = 0;
  for (unsigned n=0; n<ibufcount; n++) {
    minpkey = min<unsigned>(minpkey, buf[n].tree);
    maxpkey = max<unsigned>(maxpkey, buf[n].tree);
    if ((buf[n].tree < pk_lo || buf[n].tree >= pk_hi) && buf[n].tree > 0u) err1++;
  }

  unsigned wid = 3*nbits/4 + 3;

  if (err0)
    cerr << "ParticleFerry: Key err=" << err0 << endl << hex
	 << "ParticleFerry: min key=" << right << setw(wid) << minkey << endl 
	 << "ParticleFerry: max key=" << right << setw(wid) << maxkey << endl;
  if (err1)
    cout << "ParticleFerry: Cel err=" << dec << err1 << endl
	 << "ParticleFerry: min cel=" << right << setw(12) << minpkey << endl 
	 << "ParticleFerry: max cel=" << right << setw(12) << maxpkey << endl;

  return;
}
