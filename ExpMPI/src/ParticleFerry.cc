#include <iostream>
#include <iomanip>
#include "ParticleFerry.H"

unsigned ParticleFerry::nbuf = 2000;
MPI_Datatype ParticleFerry::Particletype;

ParticleFerry::ParticleFerry()
{
				// Assign particle structure buffer
  buf = new Partstruct [nbuf];
  ibufcount = 0;
				// Make MPI datatype
  
  MPI_Datatype type[14] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
			   MPI_DOUBLE, MPI_DOUBLE, MPI_FLOAT,
			   MPI_UNSIGNED, MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG_LONG, 
			   MPI_UNSIGNED, MPI_UNSIGNED,
			   MPI_INT, MPI_DOUBLE};

				// Get displacements
  MPI_Aint disp[14];
  MPI_Get_address(&buf[0].mass,		&disp[0]);
  MPI_Get_address(&buf[0].pos,		&disp[1]);
  MPI_Get_address(&buf[0].vel,		&disp[2]);
  MPI_Get_address(&buf[0].acc,		&disp[3]);
  MPI_Get_address(&buf[0].pot,		&disp[4]);
  MPI_Get_address(&buf[0].potext,	&disp[5]);
  MPI_Get_address(&buf[0].dtreq,	&disp[6]);
  MPI_Get_address(&buf[0].level,	&disp[7]);
  MPI_Get_address(&buf[0].indx,		&disp[8]);
  MPI_Get_address(&buf[0].key,		&disp[9]);
  MPI_Get_address(&buf[0].nicnt,	&disp[10]);
  MPI_Get_address(&buf[0].ndcnt,	&disp[11]);
  MPI_Get_address(&buf[0].iatr,		&disp[12]);
  MPI_Get_address(&buf[0].datr,		&disp[13]);

  for (int i=13; i>=0; i--) disp[i] -= disp[0];
  
				// Block offsets
  int blocklen[14] = {1, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, nimax, ndmax};
  
  MPI_Type_create_struct(14, blocklen, disp, type, &Particletype);
  MPI_Type_commit(&Particletype);
}


void ParticleFerry::part_to_Particle(Partstruct& str, Particle& cls)
{
  cls.mass = str.mass;
  for (int j=0; j<3; j++) {
      cls.pos[j] = str.pos[j];
      cls.vel[j] = str.vel[j];
      cls.acc[j] = str.acc[j];
  }
  cls.pot    = str.pot;
  cls.potext = str.potext;

  cls.dtreq  = str.dtreq;
  cls.level  = str.level;
  cls.indx   = str.indx;
  cls.key    = str.key;

  cls.iattrib = vector<int>(str.nicnt);
  for (int j=0; j<str.nicnt; j++) cls.iattrib[j] = str.iatr[j];

  cls.dattrib = vector<double>(str.ndcnt);
  for (int j=0; j<str.ndcnt; j++) cls.dattrib[j] = str.datr[j];

}

void ParticleFerry::Particle_to_part(Partstruct& str, Particle& cls)
{
  str.mass = cls.mass;
  for (int j=0; j<3; j++) {
      str.pos[j] = cls.pos[j];
      str.vel[j] = cls.vel[j];
      str.acc[j] = cls.acc[j];
  }
  str.pot = cls.pot;
  str.potext = cls.potext;

  str.dtreq = cls.dtreq;
  str.level = cls.level;
  str.indx  = cls.indx;
  str.key   = cls.key;

  str.nicnt = min<int>(nimax, cls.iattrib.size());
  str.ndcnt = min<int>(ndmax, cls.dattrib.size());

  for (int j=0; j<str.nicnt; j++) str.iatr[j] = cls.iattrib[j];

  for (int j=0; j<str.ndcnt; j++) str.datr[j] = cls.dattrib[j];
}


void ParticleFerry::ShipParticles(int to, int from, unsigned& total)
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

void ParticleFerry::SendParticle(Particle& ptc, 
				 unsigned seq, unsigned long long key)
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
  if (ibufcount == nbuf || itotcount == _total) BufferSend();
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
  if (ibufcount == nbuf || itotcount == _total) BufferSend();
}

bool ParticleFerry::RecvParticle(Particle& ptc,
				  unsigned& seq, unsigned long long& key)
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
  // bufferKeyCheck();
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
  // bufferKeyCheck();
#endif
}

void ParticleFerry::bufferKeyCheck()
{
				// Sanity check
  unsigned long long maxexp = 1;
  maxexp = maxexp << 48;
  unsigned long long minkey = 1;
  minkey = minkey << 63;
  unsigned long long maxkey = 0;
  for (unsigned n=0; n<ibufcount; n++) {
    minkey = min<unsigned long long>(minkey, buf[n].key);
    maxkey = max<unsigned long long>(maxkey, buf[n].key);
    if (buf[n].key < maxexp) {
      cerr << "oops" << endl;
    }
  }
  cout << "min key=" << hex << right << setw(12) << minkey << endl 
       << "max key=" << hex << right << setw(12) << maxkey << endl
       << "max exp=" << hex << right << setw(12) << maxexp << endl
       << dec;

  return;
}
