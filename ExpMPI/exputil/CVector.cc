
#include <cmath>
#include <cstdlib>
#include <string>
#include <sstream>

#include <kevin_complex.h>
#include <CVector.h>
#include <Vector.h>

#include <pthread.h>
#include <mpi.h>

using namespace std;

extern char threading_on;
extern pthread_mutex_t mem_lock;

/*
  Default constructor; make a null vector.
*/

CVector::CVector()
{
  low  = 0;
  high = 0;
  size = 0;
}

CVector::CVector(int l, int h)
{
  low  = 0;
  high = 0;
  size = 0;
  setsize(l, h);
}

CVector::CVector(int l, int h, double *v)
{
  low  = 0;
  high = 0;
  size = 0;
  setsize(l, h);
  
  for (int i=0; i<size; i++) pelement[i] = v[i+low];
}

CVector::CVector(int l, int h, KComplex *v)
{
  low  = 0;
  high = 0;
  size = 0;
  setsize(l, h);
  
  for (int i=0; i<size; i++) pelement[i] = v[i+low];
}

/*
  Copy constructor; create a new vector which is a copy of another.
*/

CVector::CVector(const CVector &v)
{
  low  = 0;
  high = 0;
  size = 0;

  if (v.elements.size()==0) return;

  setsize(v.low, v.high);
  
  elements = v.elements;
  pelement = &elements[0];
}

CVector::CVector(const Vector &v)
{
  low  = 0;
  high = 0;
  size = 0;

  if (v.elements.size()==0) return;

  setsize(v.low, v.high);
  
  for (int i=0; i<size; i++) pelement[i] = v.pelement[i];
}


/*
  Assignment operator for CVector; must be defined as a reference
  so that it can be used on lhs. Compatibility checking is performed;
  the destination vector is allocated if its elements are undefined.
*/

CVector &CVector::operator=(const CVector &v)
{
  int i;
  
  // Allow assignment of null vector
  if (v.elements.size() == 0) {
    low = high = 0;
    size = 0;
    return *this;
  }
  
  low      = v.low;
  high     = v.high;
  size     = high - low + 1;
  elements = v.elements;
  pelement = &elements[0];
  
  return *this;
}


KComplex &CVector::operator[](int i) const
{
  if (i<low || i>high) {
    bomb_CVector("subscript out of range");
  }
  
  return pelement[i-low];
}

void CVector::setsize(int l, int h)
{
  // do we need to resize at all?
  
  if (l==low && h==high) return;
  
  // is the requested size positive?
  
  if (h<l) {
    ostringstream mesg;
    mesg << "invalid size <0, l=" << l << " h=" << h;
    bomb_CVector(mesg.str());
  }
  
  // set the new size
  
  low  = l;
  high = h;
  size = high - low + 1;
  
  // allocate the new elements, and make a pointer
  
  elements = vector<KComplex>(size, 0.0);
  pelement = &elements[0];
}


void CVector::zero(void)
{
  for (int i=0; i<size; i++) pelement[i] = 0.0;
}

KComplex *CVector::array(int l, int h)
{
  return pelement - l;
}



CVector CVector::operator-(void)
{
  CVector tmp(low, high);
  
  for (int i=0; i<size; i++) tmp.pelement[i] = -pelement[i];
  
  return tmp;
}

CVector operator+(const CVector &v1, const CVector &v2)
{
  if (v1.low != v2.low || v1.high != v2.high) bomb_CVector_operation("+");
  
  CVector tmp(v1.low, v1.high);
  
  for (int i=0; i<v1.size; i++)
    tmp.pelement[i] = v1.pelement[i] + v2.pelement[i];
  
  return tmp;
}

CVector operator-(const CVector &v1, const CVector &v2)
{
  if (v1.low != v2.low || v1.high != v2.high) bomb_CVector_operation("-");
  
  CVector tmp(v1.low, v1.high);
  
  for (int i=0; i<v1.size; i++)
    tmp.pelement[i] = v1.pelement[i] - v2.pelement[i];
  
  return tmp;
}

CVector operator*(const KComplex &a, const CVector &v)
{
  CVector tmp(v.low, v.high);
  for (int i=0; i<v.size; i++) tmp.pelement[i] = a * v.pelement[i];
  
  return tmp;
}

CVector operator*(const CVector &v, const KComplex &a)
{
  CVector tmp(v.low, v.high);
  for (int i=0; i<v.size; i++) tmp.pelement[i] = a * v.pelement[i];
  
  return tmp;
}

CVector operator/(const CVector &v, const KComplex &a)
{
  CVector tmp(v.low, v.high);
  for (int i=0; i<v.size; i++) tmp.pelement[i] = v.pelement[i]/a;
  
  return tmp;
}


CVector operator&(const CVector &v1, const Vector &v2)
{
  if (v1.high != v2.gethigh() || v1.low != v2.getlow()) {
    bomb_CVector_operation("&");
  }
  
  CVector tmp(v1.low, v1.high);
  for (int i=v1.low; i<=v1.high; i++)
    tmp[i] = v1[i]*v2[i];
  
  return tmp;
}


CVector operator&(const Vector &v1, const CVector &v2)
{
  if (v1.gethigh() != v2.high || v1.getlow() != v2.low) {
    bomb_CVector_operation("&");
  }
  
  CVector tmp(v2.low, v2.high);
  for (int i=v2.low; i<=v2.high; i++)
    tmp[i] = v1[i] * v2[i];
  
  return tmp;
}


CVector operator*(double a, const CVector &v)
{
  CVector tmp(v.low, v.high);
  for (int i=0; i<v.size; i++) tmp.pelement[i] = a * v.pelement[i];
  
  return tmp;
}

CVector operator*(const CVector &v, double a)
{
  CVector tmp(v.low, v.high);
  for (int i=0; i<v.size; i++) tmp.pelement[i] = a * v.pelement[i];
  
  return tmp;
}

CVector operator/(const CVector &v, double a)
{
  CVector tmp(v.low, v.high);
  for (int i=0; i<v.size; i++) tmp.pelement[i] = v.pelement[i]/a;
  
  return tmp;
}


CVector &CVector::operator+=(const CVector &v)
{
  if (low != v.low || high != v.high) bomb_CVector_operation("+=");
  
  for (int i=0; i<size; i++) pelement[i] += v.pelement[i];
  
  return *this;
}	

CVector &CVector::operator-=(const CVector &v)
{
  if (low != v.low || high != v.high) bomb_CVector_operation("-=");
  
  for (int i=0; i<size; i++) pelement[i] -= v.pelement[i];
  
  return *this;
}	

CVector &CVector::operator*=(double a)
{
  for (int i=0; i<size; i++) pelement[i] *= a;
  
  return *this;
}

CVector &CVector::operator/=(double a)
{
  for (int i=0; i<size; i++) pelement[i] /= a;
  
  return *this;
}

CVector &CVector::operator*=(const KComplex &a)
{
  for (int i=0; i<size; i++) pelement[i] *= a;
  
  return *this;
}

CVector &CVector::operator/=(const KComplex &a)
{
  for (int i=0; i<size; i++) pelement[i] /= a;
  
  return *this;
}

KComplex operator*(const CVector &v1, const CVector &v2)
{
  if (v1.low != v2.low || v1.high != v2.high) bomb_CVector_operation("*");
  
  KComplex tmp = 0.0;
  for (int i=0; i<v1.size; i++) {
    tmp += v1.pelement[i] * v2.pelement[i];
  }
  return tmp;
}

void CVector::print(ostream& out)
{
  for (int i=0; i<size; i++)
    out << "(" 
	<< pelement[i].real() << " + " 
	<< pelement[i].imag() << " I) ";
  out << endl;
}

void CVector::binwrite(ostream& out)
{
  out.write((char *)&low,  sizeof(int));
  out.write((char *)&high, sizeof(int));
  out.write((char *)pelement, (high-low+1)*sizeof(KComplex));
}

CVector CVector_binread(istream& in)
{
  int low, high;
  
  in.read((char *)&low,  sizeof(int));
  in.read((char *)&high, sizeof(int));
  
  CVector tmp(low, high);
  
  in.read((char *)tmp.pelement, (high-low+1)*sizeof(KComplex));
  
  return tmp;
}


void bomb_CVector(const string& msg)
{
  cerr << "CVECTOR ERROR: " << msg << endl;
  exit(0);
}

void bomb_CVector_operation(const string& op)
{
  string msg("incompatible lengths in CVector operation ");
  msg += op;
  bomb_CVector(msg);
}

Vector CVector::Re(void)
{
  Vector tmp(low, high);
  
  for (int i=0; i<size; i++) tmp[i+low] = pelement[i].r;
  
  return tmp;
}	

Vector CVector::Im(void)
{
  Vector tmp(low, high);
  
  for (int i=0; i<size; i++) tmp[i+low] = pelement[i].imag();
  
  return tmp;
}	


CVector CVector::Conjg(void)
{
  static KComplex I(0.0, 1.0);
  CVector tmp(low, high);
  
  for (int i=0; i<size; i++) {
    tmp.pelement[i] = pelement[i].real() - I * pelement[i].imag();
  }
  
  return tmp;
}	


CVector Conjg(const CVector &c)
{
  CVector t(c);
  return t.Conjg();
}

Vector Re(const CVector &c)
{
  CVector t(c);
  return t.Re();
}

Vector Im(const CVector &c)
{
  CVector t(c);
  return t.Im();
}


CMatrix::CMatrix()
{
  rlow  = 0;
  rhigh = 0;
  clow  = 0;
  chigh = 0;
}

CMatrix::CMatrix(int rl, int rh, int cl, int ch)
{
  rlow  = 0;
  rhigh = 0;
  clow  = 0;
  chigh = 0;
  setsize(rl, rh, cl, ch);
}


CMatrix::CMatrix(int rl, int rh, int cl, int ch, KComplex **array)
{
  rlow  = 0;
  rhigh = 0;
  clow  = 0;
  chigh = 0;
  setsize(rl, rh, cl, ch);
  
  for (int i=0; i<rsize; i++) {
    for (int j=0; j<csize; j++) {
      rows[i].pelement[j] = array[i+rlow][j+clow];
    }
  }
}

CMatrix::CMatrix(int rl, int rh, int cl, int ch, double **array)
{
  rlow  = 0;
  rhigh = 0;
  clow  = 0;
  chigh = 0;
  setsize(rl, rh, cl, ch);
  
  for (int i=0; i<rsize; i++) {
    for (int j=0; j<csize; j++) {
      rows[i].pelement[j] = array[i+rlow][j+clow];
    }
  }
}

CMatrix::CMatrix(const CMatrix &m)
{
  rlow  = 0;
  rhigh = 0;
  clow  = 0;
  chigh = 0;

  if (m.rows.size() == 0) return;

  setsize(m.rlow, m.rhigh, m.clow, m.chigh);
  
  for (int i=0; i<rsize; i++) {
    for (int j=0; j<csize; j++) {
      rows[i].pelement[j] = m.rows[i].pelement[j];
    }
  }
}	


CMatrix::CMatrix(const Matrix &m)
{
  rlow  = 0;
  rhigh = 0;
  clow  = 0;
  chigh = 0;
  setsize(m.rlow, m.rhigh, m.clow, m.chigh);
  
  for (int i=0; i<rsize; i++) {
    for (int j=0; j<csize; j++) {
      rows[i].pelement[j] = m.rows[i].pelement[j];
    }
  }
}	

void CMatrix::setsize(int rl, int rh, int cl, int ch)
{
  // do we need to resize at all?
  
  if (rl==rlow && cl==clow && rh==rhigh && ch==chigh) return;
  
  // is the new size positive?
  
  if (rh<rl || ch<cl) {
    ostringstream mesg;
    mesg << "invalid size <0, rl=" << rl << " rh=" << rh
	 << " cl=" << cl << " ch" << ch;
    bomb_CMatrix(mesg.str());
  }
  
  
  // set the new size
  
  rlow  = rl;
  rhigh = rh;
  clow  = cl;
  chigh = ch;
  
  rsize = rhigh - rlow + 1;
  csize = chigh - clow + 1;
  
  // allocate the array of rows
  
  rows = vector<CVector>(rsize);
  prow = &rows[0];

  // create the individual rows
  
  for (int i=0; i<rsize; i++) rows[i].setsize(clow, chigh);
}


void bomb_CMatrix(const string& msg)
{
  cerr << "CMATRIX ERROR: " << msg << endl;
  exit(0);
}

void bomb_CMatrix_operation(const string& op)
{
  string str("incompatible sizes in CMatrix operator ");
  str += op;
  bomb_Matrix(str);
}


CVector CMatrix::fastcol(int j) const
{
  CVector tmp(rlow, rhigh);
  for (int i=rlow; i<=rhigh; i++) tmp[i] = rows[i-rlow].pelement[j-clow];
  
  return tmp;
}

CVector CMatrix::col(int j)
{
  if (j<clow || j>chigh) bomb_CMatrix("column subscript out of range");
  
  CVector tmp = fastcol(j);
  
  return tmp;
}



CVector &CMatrix::row(int i)
{
  if (i<rlow || i>rhigh) bomb_CMatrix("row subscript out of range");
  
  return rows[i-rlow];
}

CVector &CMatrix::fastrow(int i) const
{
  return prow[i-rlow];
}

void CMatrix::fastsetrow(int i, const CVector &v)
{
  rows[i-rlow] = v;
}

void CMatrix::setrow(int i, const CVector &v)
{
  if (v.low != clow || v.high != chigh) 
    bomb_CMatrix("row-vector size mismatch");

  if (i<rlow || i>rhigh) bomb_CMatrix("row subscript out of range");

  rows[i-rlow] = v;
}

void CMatrix::fastsetcol(int j, const CVector &v)
{
  for (int i=0; i<rsize; i++)
    rows[i].pelement[j-clow] = v.pelement[i];
}

void CMatrix::setcol(int j, const CVector &v)
{
  if (v.low != rlow || v.high != rhigh) 
    bomb_Matrix("column-vector size mismatch");
  if (j<clow || j>chigh) bomb_CMatrix("column subscript out of range");
  
  fastsetcol(j, v);
}

void CMatrix::zero(void)
{
  for (int i=0; i<rsize; i++) rows[i].zero();
}


CVector &CMatrix::operator[](int i) const
{
  if (i<rlow || i>rhigh) bomb_CMatrix("row subscript out of range");
  return prow[i-rlow];
}

CMatrix &CMatrix::operator=(const CMatrix &m)
{
  // Allow assignment of null matrix
  if (m.rows.size() == 0) {
    rlow  = rhigh = clow = chigh = 0;
    rsize = csize = 0;
    return *this;
  }
  
  setsize(m.rlow, m.rhigh, m.clow, m.chigh);
  
  if (m.rlow!=rlow || m.rhigh!=rhigh || m.clow!=clow || m.chigh!=chigh) 
    bomb_CMatrix_operation("=");
  
  for (int i=rlow; i<=rhigh; i++) fastsetrow(i, m.rows[i-rlow]);

  return *this;
}


CMatrix CMatrix::operator-(void)
{
  CMatrix tmp(rlow, rhigh, clow, chigh);
  
  for (int i=0; i<rsize; i++) {
    for (int j=0; j<csize; j++) {
      tmp.rows[i].pelement[j] = -rows[i].pelement[j];
    }
  }
  return tmp;
}


CMatrix operator+(const CMatrix &m1, const CMatrix &m2)
{
  if (m1.rlow!=m2.rlow || m1.rhigh!=m2.rhigh || 
      m1.clow!=m2.clow || m1.chigh!=m2.chigh) 
    bomb_CMatrix_operation("+");
  
  CMatrix tmp(m1.rlow, m1.rhigh, m1.clow, m1.chigh);
  
  for (int i=0; i<m1.rsize; i++) {
    for (int j=0; j<m1.csize; j++) {
      tmp.rows[i].pelement[j] = 
	m1.rows[i].pelement[j] + m2.rows[i].pelement[j];
    }
  }
  
  return tmp;
}

CMatrix operator-(const CMatrix &m1, const CMatrix &m2)
{
  if (m1.rlow!=m2.rlow || m1.rhigh!=m2.rhigh || 
      m1.clow!=m2.clow || m1.chigh!=m2.chigh) 
    bomb_CMatrix_operation("-");
  
  CMatrix tmp(m1.rlow, m1.rhigh, m1.clow, m1.chigh);
  
  for (int i=0; i<m1.rsize; i++) {
    for (int j=0; j<m1.csize; j++) {
      tmp.rows[i].pelement[j] = 
	m1.rows[i].pelement[j] - m2.rows[i].pelement[j];
    }
  }
  
  return tmp;
}



CMatrix operator*(const CMatrix &m, double a)
{
  CMatrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);
  
  for (int i=0; i<m.rsize; i++) {
    for (int j=0; j<m.csize; j++) {
      tmp.rows[i].pelement[j] = m.rows[i].pelement[j] * a;
    }
  }
  return tmp;
}


CMatrix operator*(double a, const CMatrix &m)
{
  CMatrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);
  
  for (int i=0; i<m.rsize; i++) {
    for (int j=0; j<m.csize; j++) {
      tmp.rows[i].pelement[j] = m.rows[i].pelement[j] * a;
    }
  }
  return tmp;
}


CMatrix operator+(const CMatrix &m, double a)
{
  CMatrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);
  
  for (int i=0; i<m.rsize; i++) {
    for (int j=0; j<m.csize; j++) {
      tmp.rows[i][j] += a;
    }
  }
  return tmp;
}


CMatrix operator+(double a, const CMatrix &m)
{
  CMatrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);
  
  for (int i=0; i<m.rsize; i++) {
    for (int j=0; j<m.csize; j++) {
      tmp.rows[i][j] += a;
    }
  }
  return tmp;
}


CMatrix operator-(const CMatrix &m, double a)
{
  CMatrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);
  
  for (int i=0; i<m.rsize; i++) {
    for (int j=0; j<m.csize; j++) {
      tmp.rows[i][j] -= a;
    }
  }
  return tmp;
}


CMatrix operator-(double a, const CMatrix &m)
{
  CMatrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);
  
  for (int i=0; i<m.rsize; i++) {
    for (int j=0; j<m.csize; j++) {
      tmp[i][j] = a - m[i][j];
    }
  }
  return tmp;
}



CMatrix operator/(const CMatrix &m, double a)
{
  CMatrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);
  
  for (int i=0; i<m.rsize; i++) {
    for (int j=0; j<m.csize; j++) {
      tmp.rows[i].pelement[j] = m.rows[i].pelement[j] / a;
    }
  }
  return tmp;
}


CMatrix &CMatrix::operator*=(double a)
{
  for (int i=0; i<rsize; i++) {
    for (int j=0; j<csize; j++) {
      rows[i].pelement[j] *= a;
    }
  }
  return *this;
}

CMatrix &CMatrix::operator/=(double a)
{
  for (int i=0; i<rsize; i++) {
    for (int j=0; j<csize; j++){
      rows[i].pelement[j] /= a;
    }
  }
  return *this;
}

CMatrix operator*(const CMatrix &m, const KComplex &a)
{
  CMatrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);
  
  for (int i=0; i<m.rsize; i++) {
    for (int j=0; j<m.csize; j++) {
      tmp.rows[i].pelement[j] = m.rows[i].pelement[j] * a;
    }
  }
  return tmp;
}




CMatrix operator*(const KComplex &a, const CMatrix &m)
{
  CMatrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);
  
  for (int i=0; i<m.rsize; i++) {
    for (int j=0; j<m.csize; j++) {
      tmp.rows[i].pelement[j] = m.rows[i].pelement[j] * a;
    }
  }
  return tmp;
}




CMatrix operator/(const CMatrix &m, const KComplex &a)
{
  CMatrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);
  
  for (int i=0; i<m.rsize; i++) {
    for (int j=0; j<m.csize; j++) {
      tmp.rows[i].pelement[j] = m.rows[i].pelement[j] / a;
    }
  }
  return tmp;
}


CMatrix &CMatrix::operator*=(const KComplex &a)
{
  for (int i=0; i<rsize; i++) {
    for (int j=0; j<csize; j++) {
      rows[i].pelement[j] *= a;
    }
  }
  return *this;
}

CMatrix &CMatrix::operator/=(const KComplex &a)
{
  for (int i=0; i<rsize; i++) {
    for (int j=0; j<csize; j++) {
      rows[i].pelement[j] /= a;
    }
  }
  return *this;
}


CMatrix &CMatrix::operator+=(const CMatrix &m)
{
  if (m.rlow!=rlow || m.rhigh!=rhigh || m.clow!=clow || m.chigh!=chigh) 
    bomb_CMatrix_operation("+=");
  
  for (int i=0; i<rsize; i++) {
    for (int j=0; j<csize; j++) {
      rows[i].pelement[j] += m.rows[i].pelement[j];
    }
  }
  return *this;
}


CMatrix &CMatrix::operator-=(const CMatrix &m)
{
  if (m.rlow!=rlow || m.rhigh!=rhigh || m.clow!=clow || m.chigh!=chigh) 
    bomb_CMatrix_operation("+=");
  
  for (int i=0; i<rsize; i++) {
    for (int j=0; j<csize; j++) {
      rows[i].pelement[j] -= m.rows[i].pelement[j];
    }
  }
  return *this;
}


CVector operator*(const CMatrix &m, const CVector &v)
{
  if (m.clow != v.low || m.chigh != v.high) bomb_CMatrix_operation("M*v");
  
  CVector tmp(m.rlow, m.rhigh);
  
  for (int i=0; i<m.rsize; i++) {
    tmp.pelement[i] = 0.0;
    for (int j=0; j<m.csize; j++) {
      tmp.pelement[i] += 
	m.rows[i].pelement[j]*v.pelement[j];
    }
  }
  return tmp;
}

CVector operator*(const CMatrix &m, const Vector &v)
{
  if (m.clow != v.low || m.chigh != v.high) bomb_CMatrix_operation("M*v");
  
  CVector tmp(m.rlow, m.rhigh);
  
  for (int i=0; i<m.rsize; i++) {
    tmp.pelement[i] = 0.0;
    for (int j=0; j<m.csize; j++) {
      tmp.pelement[i] += 
	m.rows[i].pelement[j]*v.pelement[j];
    }
  }
  return tmp;
}

CVector operator*(const CVector &v, const CMatrix &m)
{
  if (m.rlow != v.low || m.rhigh != v.high) bomb_CMatrix_operation("v*M");
  
  CVector tmp(m.clow, m.chigh);
  
  for (int j=0; j<m.csize; j++) {
    tmp.pelement[j] = 0.0;
    for (int i=0; i<m.rsize; i++) {
      tmp.pelement[j] +=
	m.rows[i].pelement[j]*v.pelement[i];
    }
  }
  return tmp;
}


CMatrix operator*(const CMatrix &m1, const CMatrix &m2)
{
  if (m1.clow != m2.rlow || m1.chigh != m2.rhigh)
    bomb_CMatrix_operation("M*M");
  
  CMatrix tmp(m1.rlow, m1.rhigh, m2.clow, m2.chigh);
  
  for (int i=0; i<m1.rsize; i++) {
    for (int j=0; j<m2.csize; j++) {
      tmp.rows[i].pelement[j] = 0.0;
      for (int k=0; k<m1.csize; k++){
	tmp.rows[i].pelement[j] += 
	  m1.rows[i].pelement[k] * m2.rows[k].pelement[j];
      }
    }
  }
  return tmp;
  
}


CMatrix CMatrix::Conjg(void) const
{
  CMatrix tmp(rlow, rhigh, clow, chigh);
  
  for (int i=0; i<rsize; i++) {
    for (int j=clow; j<=chigh; j++) {
      tmp.rows[i][j].r =  rows[i][j].r;
      tmp.rows[i][j].r = -rows[i][j].i;
    }
  }      
  
  return tmp;
}





CMatrix CMatrix::Transpose(void) const
{
  CMatrix t(clow, chigh, rlow, rhigh);
  
  t.rlow  = clow;
  t.rhigh = chigh;
  t.chigh = rhigh;
  t.clow  = rlow;
  
  for (int i=t.rlow; i<=t.rhigh; i++) t.fastsetrow(i, fastcol(i));
  
  return t;
}

KComplex CMatrix::Trace(void) const
{
  if (rlow != clow || rhigh != chigh) 
    bomb_CVector("cannot take trace of non-square matrix");
  
  
  KComplex t = 0.0;
  for (int i=0; i<rsize; i++) t += rows[i][i];
  
  return t;
}



Matrix CMatrix::Re(void) const
{
  Matrix tmp(rlow, rhigh, clow, chigh);
  
  for (int i=0; i<rsize; i++) {
    for (int j=0; j<csize; j++) 
      tmp.rows[i].pelement[j] = rows[i].pelement[j].r;
  }
  return tmp;
}	


Matrix CMatrix::Im(void) const
{
  Matrix tmp(rlow, rhigh, clow, chigh);
  
  for (int i=0; i<rsize; i++) {
    for (int j=0; j<csize; j++) 
      tmp.rows[i].pelement[j] = rows[i].pelement[j].i;
  }
  
  return tmp;
}	


Matrix Re(const CMatrix &c)
{
  return c.Re();
}

Matrix Im(const CMatrix &c)
{
  return c.Im();
}
CMatrix Conjg(const CMatrix &c)
{
  return c.Conjg();
}

CMatrix Transpose(const CMatrix &c)
{
  return c.Transpose();
}

CMatrix Adjoint(const CMatrix &c)
{
  return c.Conjg().Transpose();
}

KComplex Trace(const CMatrix &c)
{
  return c.Trace();
}

void CMatrix::print(ostream& out)
{
  for (int i=0; i<rsize; i++) rows[i].print(out);
}



void CMatrix::binwrite(ostream& out)
{
  out.write((char *)&rlow,  sizeof(int));
  out.write((char *)&rhigh, sizeof(int));
  out.write((char *)&clow,  sizeof(int));
  out.write((char *)&chigh, sizeof(int));
  
  for (int i=0; i<rsize; i++) rows[i].binwrite(out);
}


CMatrix CMatrix_binread(istream& in)
{
  
  int rlow, rhigh, clow, chigh;
  
  in.read((char *)&rlow,  sizeof(int));
  in.read((char *)&rhigh, sizeof(int));
  in.read((char *)&clow,  sizeof(int));
  in.read((char *)&chigh, sizeof(int));
  
  CMatrix tmp(rlow, rhigh, clow, chigh);
  
  for (int i=0; i<tmp.rsize; i++)
    tmp.rows[i] = CVector_binread(in);
  
  return tmp;
}

void ComplexSynchronize(KComplex& c, int id)
{
  MPI_Bcast(&c.real(), 1, MPI_DOUBLE, id, MPI_COMM_WORLD);
  MPI_Bcast(&c.imag(), 1, MPI_DOUBLE, id, MPI_COMM_WORLD);
}


void CVectorSynchronize(CVector& cv, int id)
{
  int bb[2], myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  if (myid == id) {
    bb[0] = cv.getlow();
    bb[1] = cv.gethigh();
  }

  MPI_Bcast(bb, 2, MPI_INT, id, MPI_COMM_WORLD);
  
  if (bb[0]==0 && bb[1]==0) return;
  
  if (myid!=id) cv.setsize(bb[0], bb[1]);
  
  int sz = bb[1] - bb[0] + 1;
  vector<double> tmp(sz);

  // Real part
  
  if (myid==id) {
    for (int i=bb[0]; i<=bb[1]; i++) tmp[i-bb[0]] = cv[i].real();
  }

  MPI_Bcast(&tmp[0], sz, MPI_DOUBLE, id, MPI_COMM_WORLD);

  if (myid!=id) {
    for (int i=bb[0]; i<=bb[1]; i++) cv[i].real() = tmp[i-bb[0]];
  }

  // Imag part
    
  if (myid==id) {
    for (int i=bb[0]; i<=bb[1]; i++) tmp[i-bb[0]] = cv[i].imag();
  }

  MPI_Bcast(&tmp[0], sz, MPI_DOUBLE, id, MPI_COMM_WORLD);

  if (myid!=id) {
    for (int i=bb[0]; i<=bb[1]; i++) cv[i].imag() = tmp[i-bb[0]];
  }
}

void CMatrixSynchronize(CMatrix& mat, int id)
{
  int bb[4], myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  if (myid == id) {
    bb[0] = mat.getrlow();
    bb[1] = mat.getrhigh();
    bb[2] = mat.getclow();
    bb[3] = mat.getchigh();
  }

  MPI_Bcast(bb, 4, MPI_INT, id, MPI_COMM_WORLD);

  if (bb[0]==0 && bb[1]==0 && bb[2]==0 && bb[3]==0) return;

  if (myid!=id) mat.setsize(bb[0], bb[1], bb[2], bb[3]);
  
  int sz = bb[3] - bb[2] + 1;
  vector<double> tmp(sz);
  for (int j=bb[0]; j<=bb[1]; j++) {

    // Real part of column

    if (myid==id) {
      for (int i=bb[2]; i<=bb[3]; i++) 
	tmp[i-bb[2]] = mat[j][i].real();
    }

    MPI_Bcast(&tmp[0], sz, MPI_DOUBLE, id, MPI_COMM_WORLD);

    if (myid!=id) {
      for (int i=bb[2]; i<=bb[3]; i++) 
	mat[j][i].real() = tmp[i-bb[2]];
    }

    // Imag part of column
    
    if (myid==id) {
      for (int i=bb[2]; i<=bb[3]; i++) 
	tmp[i-bb[2]] = mat[j][i].imag();
    }

    MPI_Bcast(&tmp[0], sz, MPI_DOUBLE, id, MPI_COMM_WORLD);

    if (myid!=id) {
      for (int i=bb[2]; i<=bb[3]; i++) 
	mat[j][i].imag() = tmp[i-bb[2]];
    }
  }
}

