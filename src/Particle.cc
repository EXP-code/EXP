#include "expand.h"

#include <Particle.H>

float Particle::effort_default = 1.0e-12;
const speciesKey Particle::defaultKey {-1, -1};


Particle::Particle()
{
  //
  // Initialize basic dynamical fields
  //
  mass = pot = potext = 0.0;
  for (int k=0; k<3; k++) pos[k] = vel[k] = acc[k] = 0.0;

  //
  // Time stepping and partitioning info
  //
  level   = 0;
  dtreq   = -1;
  scale   = -1;
  effort  = effort_default;
  indx    = 0;
  tree    = 0u;
  key     = 0u;
  skey    = defaultKey;
}

Particle::Particle(unsigned niatr, unsigned ndatr)
{
  //
  // Initialize basic fields
  //
  mass = pot = potext = 0.0;
  for (int k=0; k<3; k++) pos[k] = vel[k] = acc[k] = 0.0;
  level   = 0;
  dtreq   = -1;
  scale   = -1;
  effort  = effort_default;
  indx    = 0;
  tree    = 0u;
  key     = 0u;
  iattrib = vector<int   >(niatr, 0);
  dattrib = vector<double>(ndatr, 0);
  skey    = defaultKey;
}

Particle::Particle(const Particle &p)
{
  mass = p.mass;
  for (int k=0; k<3; k++) {
    pos[k] = p.pos[k];
    vel[k] = p.vel[k];
    acc[k] = p.acc[k];
  }
  pot     = p.pot;
  potext  = p.potext;
  iattrib = p.iattrib;
  dattrib = p.dattrib;
  level   = p.level;
  dtreq   = p.dtreq;
  scale   = p.scale;
  effort  = p.effort;
  indx    = p.indx;
  tree    = p.tree;
  key     = p.key;
  skey    = p.skey;
}


void Particle::readBinary(unsigned rsize, bool indexing, int seq, 
			  std::istream *in)
{
  //
  // Read index value if this field is recorded
  //
  if (indexing) 
    in->read((char *)&(indx), sizeof(unsigned long));
  else
    indx = seq;
  
  //
  // Floating (4-byte version)
  //
  if (rsize == sizeof(float)) {

    float tf;			// Temporary float value

    in->read((char *)&tf, sizeof(float));
    mass = tf;

    for (int i=0; i<3; i++) {
      in->read((char *)&tf, sizeof(float));
      pos[i] = tf;
    }

    for (int i=0; i<3; i++) {
      in->read((char *)&tf, sizeof(float));
      vel[i] = tf;
    }

    in->read((char *)&tf, sizeof(float));
    pot    = tf;
    potext = 0.0;

    level = multistep;

    for (auto &it : iattrib)
      in->read((char *)&it, sizeof(int));

    for (auto &jt : dattrib) {
      in->read((char *)&tf, sizeof(float));
      jt = tf;
    }
    
  } else {
    //
    // Floating (8-byte version)
    //
    in->read((char *)&mass, sizeof(double));

    for (int i=0; i<3; i++) in->read((char *)&(pos[i]), sizeof(double));

    for (int i=0; i<3; i++) in->read((char *)&(vel[i]), sizeof(double));

    in->read((char *)&pot, sizeof(double));
    potext = 0.0;

    level = multistep;

    for (auto& it : iattrib)
      in->read((char *)&it, sizeof(int));

    for (auto& jt : dattrib)
      in->read((char *)&jt, sizeof(double));

  }
}


void Particle::writeBinary(unsigned rsize, 
			   double* com0, double* comI,
			   double* cov0, double* covI,
			   bool indexing, std::ostream *out)
{
  // Working variable
  float tf;

  if (indexing) 		// Cache index if desired
    out->write((const char *)&(indx), sizeof(unsigned long));

  if (rsize == sizeof(float)) {
    tf = static_cast<float>(mass);
    out->write((const char *)&tf, sizeof(float));
  }
  else
    out->write((const char *)&(mass), sizeof(double));
  
  for (int i=0; i<3; i++) {
    double pv = pos[i] + com0[i] - comI[i];
    if (rsize == sizeof(float)) {
      tf = static_cast<float>(pv);
      out->write((const char *)&tf, sizeof(float));
    }
    else
      out->write((const char *)&pv, sizeof(double));
  }
  
  for (int i=0; i<3; i++) {
    double pv = vel[i] + cov0[i] - covI[i];
    if (rsize == sizeof(float)) {
      tf = static_cast<float>(pv);
      out->write((const char *)&tf, sizeof(float));
    }
    else
      out->write((const char *)&pv, sizeof(double));
  }

  double pot0 = pot + potext;
  if (rsize == sizeof(float)) {
    tf = static_cast<float>(pot0);
    out->write((const char *)&tf, sizeof(float));
  }
  else
    out->write((const char *)&pot0, sizeof(double));

  for (auto it : iattrib)
    out->write((const char *)&it, sizeof(int));
  
  for (auto jt: dattrib) {
    if (rsize == sizeof(float)) {
      tf = static_cast<float>(jt);
      out->write((const char *)&tf, sizeof(float));
    }
    else
      out->write((const char *)&jt, sizeof(double));
  }
}


void Particle::readAscii(bool indexing, int seq, std::istream* fin)
{
  //
  // Character array for file reading
  //
  const int nline = 2048;
  char line[nline];

  //
  // Read the line
  //
  fin->getline(line, nline);
  istringstream ins(line);

  if (indexing)
    ins >> indx;
  else
    indx = seq;

  ins >> mass;
  for (int j=0; j<3; j++) ins >> pos[j];
  for (int j=0; j<3; j++) ins >> vel[j];
  for (int j=0; j<3; j++) acc[j] = 0.0;
  pot = potext = 0.0;
  
  level = multistep;

  for (auto &it : iattrib) {
    ins >> it;
    if (!ins) it = 0;
  }

  for (auto &jt : dattrib) {
    ins >> jt;
    if (!ins) jt = 0;
  }
}

void Particle::writeAscii(double* com0, double* comI, 
			  double* cov0, double* covI, 
			  bool indexing, bool accel, std::ostream* out)
{
  if (indexing) *out << std::setw(12) << indx;
  *out << std::setw(18) << mass;
  for (int i=0; i<3; i++) *out << std::setw(18) << pos[i]+com0[i]-comI[i];
  for (int i=0; i<3; i++) *out << std::setw(18) << vel[i]+cov0[i]-covI[i];
  if (accel)
    for (int i=0; i<3; i++) *out << std::setw(18) << acc[i];
    
  *out << std::setw(18) << pot;
  *out << std::setw(18) << potext;

  for (auto it : iattrib)
    *out << std::setw(10) << it;

  for (auto jt : dattrib)
    *out << std::setw(18) << jt;
  
  *out << std::endl;
}

// For debugging . . . 
std::ostream& operator<< (std::ostream& os, const PMapType& p)
{
  std::streamsize sp = os.precision();
  os.precision(6);
  os << std::setw(10) << p.second->indx
     << std::setw( 4) << p.second->level
     << std::setw(16) << p.second->mass;
  for (int k=0; k<3; k++) os << std::setw(16) << p.second->pos[k];
  // for (int k=0; k<3; k++) os << std::setw(16) << p->second.vel[k];
  for (int k=0; k<3; k++) os << std::setw(16) << p.second->acc[k];
  // os << std::setw(16) << p.second.pot << std::setw(16) << p.second->potext;
  os.precision(sp);
  return os;
}
