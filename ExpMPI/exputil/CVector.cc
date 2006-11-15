
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <string>
#include <sstream>

#include <kevin_complex.h>
#include <CVector.h>
#include <Vector.h>

#include <pthread.h>

using namespace std;

extern char threading_on;
extern pthread_mutex_t mem_lock;

/*
  Default constructor; make a null vector.
*/

CVector::CVector()
{
  low=0;
  high=0;
  elements = NULL;
}


CVector::CVector(int l, int h)
{
  low=0;
  high=0;
  elements = NULL;
  setsize(l, h);
}

CVector::CVector(int l, int h, double *v)
{
  int i;
  
  low=0;
  high=0;
  elements = NULL;
  setsize(l, h);
  
  for (i=low; i<=high; i++) elements[i] = v[i];
}

CVector::CVector(int l, int h, Complex *v)
{
  int i;
  
  low=0;
  high=0;
  elements = NULL;
  setsize(l, h);
  
  for (i=low; i<=high; i++) elements[i] = v[i];
}


/*
  Copy constructor; create a new vector which is a copy of another.
*/

CVector::CVector(const CVector &v)
{
  int i;
  
  low=0;
  high=0;
  elements = NULL;
  setsize(v.low, v.high);
  
  for (i=low; i<=high; i++) elements[i] = v.elements[i];
}

CVector::CVector(const Vector &v)
{
  int i;
  
  low=0;
  high=0;
  elements = NULL;
  setsize(v.low, v.high);
  
  for (i=low; i<=high; i++) elements[i] = v.elements[i];
}




/*
  Destructor. Free elements if it exists.
*/

CVector::~CVector()
{
  if (threading_on) pthread_mutex_lock(&mem_lock);
  delete [] (elements+low);
  if (threading_on) pthread_mutex_unlock(&mem_lock);
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
  if (v.elements==NULL) {
    elements = NULL;
    low = high = 0;
    return *this;
  }
  
  if (low!=v.low && high!=v.high 
      && (high != 0 || (elements+low) != NULL))
    {
      bomb_CVector_operation("=");
    }
  
  if (high == 0 && low==0 && elements == NULL)
    {
      setsize(v.low, v.high);
    }
  
  for (i=low; i<=high; i++) elements[i] = v.elements[i];
  
  return *this;
}



Complex &CVector::operator[](int i) const
{
  if (i<low || i>high)
    {
      bomb_CVector("subscript out of range");
    }
  
  return elements[i];
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
  
  
  // delete the old elements if they already exist
  
  if (threading_on) pthread_mutex_lock(&mem_lock);
  delete [] (elements + low);
  if (threading_on) pthread_mutex_unlock(&mem_lock);
  
  
  
  // set the new size
  
  low = l;
  high = h;
  
  
  
  // allocate the new elements, and offset the pointer
  
  if (threading_on) pthread_mutex_lock(&mem_lock);
  elements = new Complex[high-low+1];
  if (threading_on) pthread_mutex_unlock(&mem_lock);
  if (elements == NULL) bomb_CVector("could not allocate");
  elements -= low;
}


void CVector::zero(void)
{
  int i;
  
  for (i=low; i<=high; i++) elements[i] = 0.0;
}

Complex *CVector::array(int l, int h)
{
  return elements+low-l;
}



CVector CVector::operator-(void)
{
  int i;
  
  CVector tmp(low, high);
  
  for (i=low; i<=high; i++) tmp.elements[i] = -elements[i];
  
  return tmp;
}




CVector operator+(const CVector &v1, const CVector &v2)
{
  int i;
  
  if (v1.low != v2.low || v1.high != v2.high) bomb_CVector_operation("+");
  
  CVector tmp(v1.low, v1.high);
  
  for (i=v1.low; i<=v1.high; i++)
    {
      tmp.elements[i] = v1.elements[i] + v2.elements[i];
    }
  
  return tmp;
}

CVector operator-(const CVector &v1, const CVector &v2)
{
  int i;
  
  if (v1.low != v2.low || v1.high != v2.high) bomb_CVector_operation("-");
  
  CVector tmp(v1.low, v1.high);
  
  for (i=v1.low; i<=v1.high; i++)
    {
      tmp.elements[i] = v1.elements[i] - v2.elements[i];
    }
  
  return tmp;
}

CVector operator*(const Complex &a, const CVector &v)
{
  int i;
  
  CVector tmp(v.low, v.high);
  for (i=v.low; i<=v.high; i++) tmp.elements[i] = a * v.elements[i];
  
  return tmp;
}

CVector operator*(const CVector &v, const Complex &a)
{
  int i;
  
  CVector tmp(v.low, v.high);
  for (i=v.low; i<=v.high; i++) tmp.elements[i] = a * v.elements[i];
  
  return tmp;
}

CVector operator/(const CVector &v, const Complex &a)
{
  int i;
  
  CVector tmp(v.low, v.high);
  for (i=v.low; i<=v.high; i++) tmp.elements[i] = v.elements[i]/a;
  
  return tmp;
}


CVector operator&(const CVector &v1, const Vector &v2)
{
  int i;
  
  if (v1.high != v2.gethigh() || v1.low != v2.getlow())
    {
      bomb_CVector_operation("&");
    }
  
  CVector tmp(v1.low, v1.high);
  for (i=v1.low; i<=v1.high; i++)
    tmp[i] = v1[i]*v2[i];
  
  return tmp;
}


CVector operator&(const Vector &v1, const CVector &v2)
{
  int i;
  
  if (v1.gethigh() != v2.high || v1.getlow() != v2.low)
    {
      bomb_CVector_operation("&");
    }
  
  CVector tmp(v2.low, v2.high);
  for (i=v2.low; i<=v2.high; i++)
    tmp[i] = v1[i]*v2[i];
  
  return tmp;
}


CVector operator*(double a, const CVector &v)
{
  int i;
  
  CVector tmp(v.low, v.high);
  for (i=v.low; i<=v.high; i++) tmp.elements[i] = a * v.elements[i];
  
  return tmp;
}

CVector operator*(const CVector &v, double a)
{
  int i;
  
  CVector tmp(v.low, v.high);
  for (i=v.low; i<=v.high; i++) tmp.elements[i] = a * v.elements[i];
  
  return tmp;
}

CVector operator/(const CVector &v, double a)
{
  int i;
  
  CVector tmp(v.low, v.high);
  for (i=v.low; i<=v.high; i++) tmp.elements[i] = v.elements[i]/a;
  
  return tmp;
}


CVector &CVector::operator+=(const CVector &v)
{
  int i;
  
  if (low != v.low || high != v.high) bomb_CVector_operation("+=");
  
  for (i=low; i<=high; i++) elements[i] += v[i];
  
  
  return *this;
}	

CVector &CVector::operator-=(const CVector &v)
{
  int i;
  
  if (low != v.low || high != v.high) bomb_CVector_operation("-=");
  
  for (i=low; i<=high; i++) elements[i] -= v[i];
  
  
  return *this;
}	

CVector &CVector::operator*=(double a)
{
  int i;
  
  for (i=low; i<=high; i++) elements[i] *= a;
  
  
  return *this;
}

CVector &CVector::operator/=(double a)
{
  int i;
  
  for (i=low; i<=high; i++) elements[i] /= a;
  
  
  return *this;
}

CVector &CVector::operator*=(const Complex &a)
{
  int i;
  
  for (i=low; i<=high; i++) elements[i] *= a;
  
  
  return *this;
}

CVector &CVector::operator/=(const Complex &a)
{
  int i;
  
  for (i=low; i<=high; i++) elements[i] /= a;
  
  return *this;
}



Complex operator*(const CVector &v1, const CVector &v2)
{
  int i;
  Complex tmp;
  
  if (v1.low != v2.low || v1.high != v2.high) bomb_CVector_operation("*");
  
  tmp = 0.0;
  for (i=v1.low; i<=v1.high; i++)
    {
      tmp += v1.elements[i] * v2.elements[i];
    }
  return tmp;
}



void CVector::print(ostream& out)
{
  int i;
  
  for (i=low; i<=high; i++)
    out << "(" 
	<< elements[i].real() << " + " 
	<< elements[i].imag() << " I) ";
  out << endl;
}

void CVector::binwrite(ostream& out)
{
  out.write((char *)&low, sizeof(int));
  out.write((char *)&high, sizeof(int));
  out.write((char *)(elements+low), (high-low+1)*sizeof(Complex));
}

CVector CVector_binread(istream& in)
{
  int low, high;
  
  in.read((char *)&low, sizeof(int));
  in.read((char *)&high, sizeof(int));
  
  CVector tmp(low, high);
  
  in.read((char *)(tmp.elements+low), (high-low+1)*sizeof(Complex));
  
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
  int i;
  Vector tmp(low, high);
  
  for (i=low; i<=high; i++) tmp[i] = elements[i].real();
  
  return tmp;
}	

Vector CVector::Im(void)
{
  int i;
  Vector tmp(low, high);
  
  for (i=low; i<=high; i++) tmp[i] = elements[i].imag();
  
  return tmp;
}	


CVector CVector::Conjg(void)
{
  int i;
  static Complex I(0.0, 1.0);
  CVector tmp(low, high);
  
  for (i=low; i<=high; i++) 
    {
      tmp.elements[i] = elements[i].real() - I * elements[i].imag();
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
  rlow=0;
  rhigh=0;
  clow=0;
  chigh=0;
  rows = NULL;
}


CMatrix::CMatrix(int rl, int rh, int cl, int ch)
{
  rlow=0;
  rhigh=0;
  clow=0;
  chigh=0;
  rows = NULL;
  setsize(rl, rh, cl, ch);
}


CMatrix::CMatrix(int rl, int rh, int cl, int ch, Complex **array)
{
  int i, j;
  
  rlow=0;
  rhigh=0;
  clow=0;
  chigh=0;
  rows = NULL;
  setsize(rl, rh, cl, ch);
  
  for (i=rlow; i<=rhigh; i++)
    {
      for (j=clow; j<=chigh; j++)
	{
	  rows[i].elements[j] = array[i][j];
	}
    }
}

CMatrix::CMatrix(int rl, int rh, int cl, int ch, double **array)
{
  int i, j;
  
  rlow=0;
  rhigh=0;
  clow=0;
  chigh=0;
  rows = NULL;
  setsize(rl, rh, cl, ch);
  
  for (i=rlow; i<=rhigh; i++)
    {
      for (j=clow; j<=chigh; j++)
	{
	  rows[i].elements[j] = array[i][j];
	}
    }
}

CMatrix::CMatrix(const CMatrix &m)
{
  int i, j;
  
  rlow=0;
  rhigh=0;
  clow=0;
  chigh=0;
  rows = NULL;
  setsize(m.rlow, m.rhigh, m.clow, m.chigh);
  
  for (i=rlow; i<=rhigh; i++)
    {
      for (j=clow; j<=chigh; j++)
	{
	  rows[i].elements[j] = m.rows[i].elements[j];
	}
    }
}	


CMatrix::CMatrix(const Matrix &m)
{
  int i, j;
  
  rlow=0;
  rhigh=0;
  clow=0;
  chigh=0;
  rows = NULL;
  setsize(m.rlow, m.rhigh, m.clow, m.chigh);
  
  for (i=rlow; i<=rhigh; i++)
    {
      for (j=clow; j<=chigh; j++)
	{
	  rows[i].elements[j] = m.rows[i].elements[j];
	}
    }
}	


CMatrix::~CMatrix(void)
{
  rows += rlow;
  if (threading_on) pthread_mutex_lock(&mem_lock);
  delete [] rows;
  if (threading_on) pthread_mutex_unlock(&mem_lock);
}


void CMatrix::setsize(int rl, int rh, int cl, int ch)
{
  int i;
  
  // do we need to resize at all?
  
  if (rl==rlow && cl==clow && rh==rhigh && ch==chigh) return;
  
  
  
  // is the new size positive?
  
  if (rh<rl || ch<cl) {
    ostringstream mesg;
    mesg << "invalid size <0, rl=" << rl << " rh=" << rh
	 << " cl=" << cl << " ch" << ch;
    bomb_CMatrix(mesg.str());
  }
  
  
  // delete old storage if it exists
  
  if (threading_on) pthread_mutex_lock(&mem_lock);
  delete [] (rows+rlow);
  if (threading_on) pthread_mutex_unlock(&mem_lock);
  
  
  
  // set the new size
  
  rlow = rl;
  rhigh = rh;
  clow = cl;
  chigh = ch;
  
  
  
  // allocate the array of rows
  
  if (threading_on) pthread_mutex_lock(&mem_lock);
  rows = new CVector[rhigh+1-rlow];
  if (threading_on) pthread_mutex_unlock(&mem_lock);
  if (rows==NULL) bomb_CMatrix("could not allocate rows");
  rows -= rlow;
  
  
  
  // create the individual rows
  
  for (i=rlow; i<=rhigh; i++) rows[i].setsize(clow, chigh);
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
  /*
    Complex *array;
    int i;
    
    if (threading_on) pthread_mutex_lock(&mem_lock);
    array = new Complex[rhigh-rlow+1];
    if (threading_on) pthread_mutex_unlock(&mem_lock);
    if (array == NULL) bomb_CMatrix("could not allocate in fastcol");
    array -= rlow;
    
    for (i=rlow; i<=rhigh; i++) array[i] = rows[i].elements[j];
  */
  
  CVector tmp(rlow, rhigh);
  for (int i=rlow; i<=rhigh; i++) tmp[i] = rows[i].elements[j];
  
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
  
  return rows[i];
}

CVector &CMatrix::fastrow(int i) const
{
  return rows[i];
}



void CMatrix::fastsetrow(int i, const CVector &v)
{
  rows[i] = v;
}

void CMatrix::setrow(int i, const CVector &v)
{
  if (v.low != clow || v.high != chigh) 
    bomb_CMatrix("row-vector size mismatch");
  if (i<rlow || i>rhigh) bomb_CMatrix("row subscript out of range");
  rows[i] = v;
}

void CMatrix::fastsetcol(int j, const CVector &v)
{
  int i;
  
  for (i=rlow; i<=rhigh; i++)
    {
      rows[i].elements[j] = v.elements[i];
    }
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
  int i;
  
  for (i=rlow; i<=rhigh; i++) rows[i].zero();
}


CVector &CMatrix::operator[](int i) const
{
  if (i<rlow || i>rhigh) bomb_CMatrix("row subscript out of range");
  return rows[i];
}

CMatrix &CMatrix::operator=(const CMatrix &m)
{
  int i;
  
  // Allow assignment of null matrix
  if (m.rows==NULL) {
    rows = NULL;
    rlow = rhigh = clow = chigh = 0;
    return *this;
  }
  
  if ((rows+rlow)==NULL) 
    {
      setsize(m.rlow, m.rhigh, m.clow, m.chigh);
    }
  
  if (m.rlow!=rlow || m.rhigh!=rhigh || m.clow!=clow || m.chigh!=chigh) 
    bomb_CMatrix_operation("=");
  
  for (i=rlow; i<=rhigh; i++) 
    {
      fastsetrow(i, m.rows[i]);
    }
  return *this;
}


CMatrix CMatrix::operator-(void)
{
  int i, j;
  CMatrix tmp(rlow, rhigh, clow, chigh);
  
  for (i=rlow; i<=rhigh; i++)
    {
      for (j=clow; j<=chigh; j++) 
	tmp.rows[i].elements[j] = -rows[i].elements[j];
    }
  return tmp;
}


CMatrix operator+(const CMatrix &m1, const CMatrix &m2)
{
  int i, j;
  
  if (m1.rlow!=m2.rlow || m1.rhigh!=m2.rhigh || 
      m1.clow!=m2.clow || m1.chigh!=m2.chigh) 
    bomb_CMatrix_operation("+");
  
  CMatrix tmp(m1.rlow, m1.rhigh, m1.clow, m1.chigh);
  
  for (i=m1.rlow; i<=m1.rhigh; i++) 
    {
      for (j=m1.clow; j<=m1.chigh; j++)
	{
	  tmp.rows[i].elements[j] = m1.rows[i].elements[j]
	    + m2.rows[i].elements[j];
	}
    }
  
  return tmp;
}

CMatrix operator-(const CMatrix &m1, const CMatrix &m2)
{
  int i, j;
  
  if (m1.rlow!=m2.rlow || m1.rhigh!=m2.rhigh || 
      m1.clow!=m2.clow || m1.chigh!=m2.chigh) 
    bomb_CMatrix_operation("-");
  
  CMatrix tmp(m1.rlow, m1.rhigh, m1.clow, m1.chigh);
  
  for (i=m1.rlow; i<=m1.rhigh; i++) 
    {
      for (j=m1.clow; j<=m1.chigh; j++)
	{
	  tmp.rows[i].elements[j] = m1.rows[i].elements[j]
	    - m2.rows[i].elements[j];
	}
    }
  
  return tmp;
}



CMatrix operator*(const CMatrix &m, double a)
{
  int i, j;
  
  CMatrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);
  
  for (i=m.rlow; i<=m.rhigh; i++) 
    {
      for (j=m.clow; j<=m.chigh; j++)
	{
	  tmp.rows[i].elements[j] = m.rows[i].elements[j] * a;
	}
    }
  return tmp;
}


CMatrix operator*(double a, const CMatrix &m)
{
  int i, j;
  
  CMatrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);
  
  for (i=m.rlow; i<=m.rhigh; i++) 
    {
      for (j=m.clow; j<=m.chigh; j++)
	{
	  tmp.rows[i].elements[j] = m.rows[i].elements[j] * a;
	}
    }
  return tmp;
}


CMatrix operator+(const CMatrix &m, double a)
{
  int i, j;
  
  CMatrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);
  
  for (i=m.rlow; i<=m.rhigh; i++) 
    {
      for (j=m.clow; j<=m.chigh; j++)
	{
	  tmp.rows[i][j] += a;
	}
    }
  return tmp;
}


CMatrix operator+(double a, const CMatrix &m)
{
  int i, j;
  
  CMatrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);
  
  for (i=m.rlow; i<=m.rhigh; i++) 
    {
      for (j=m.clow; j<=m.chigh; j++)
	{
	  tmp.rows[i][j] += a;
	}
    }
  return tmp;
}


CMatrix operator-(const CMatrix &m, double a)
{
  int i, j;
  
  CMatrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);
  
  for (i=m.rlow; i<=m.rhigh; i++) 
    {
      for (j=m.clow; j<=m.chigh; j++)
	{
	  tmp.rows[i][j] -= a;
	}
    }
  return tmp;
}


CMatrix operator-(double a, const CMatrix &m)
{
  int i, j;
  
  CMatrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);
  
  for (i=m.rlow; i<=m.rhigh; i++) 
    {
      for (j=m.clow; j<=m.chigh; j++)
	{
	  tmp[i][j] = a - m[i][j];
	}
    }
  return tmp;
}



CMatrix operator/(const CMatrix &m, double a)
{
  int i, j;
  
  CMatrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);
  
  for (i=m.rlow; i<=m.rhigh; i++) 
    {
      for (j=m.clow; j<=m.chigh; j++)
	{
	  tmp.rows[i].elements[j] = m.rows[i].elements[j] / a;
	}
    }
  return tmp;
}


CMatrix &CMatrix::operator*=(double a)
{
  int i, j;
  
  for (i=rlow; i<=rhigh; i++) 
    {
      for (j=clow; j<=chigh; j++)
	{
	  rows[i].elements[j] *= a;
	}
    }
  return *this;
}

CMatrix &CMatrix::operator/=(double a)
{
  int i, j;
  
  for (i=rlow; i<=rhigh; i++) 
    {
      for (j=clow; j<=chigh; j++)
	{
	  rows[i].elements[j] /= a;
	}
    }
  return *this;
}



CMatrix operator*(const CMatrix &m, const Complex &a)
{
  int i, j;
  
  CMatrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);
  
  for (i=m.rlow; i<=m.rhigh; i++) 
    {
      for (j=m.clow; j<=m.chigh; j++)
	{
	  tmp.rows[i].elements[j] = m.rows[i].elements[j] * a;
	}
    }
  return tmp;
}




CMatrix operator*(const Complex &a, const CMatrix &m)
{
  int i, j;
  
  CMatrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);
  
  for (i=m.rlow; i<=m.rhigh; i++) 
    {
      for (j=m.clow; j<=m.chigh; j++)
	{
	  tmp.rows[i].elements[j] = m.rows[i].elements[j] * a;
	}
    }
  return tmp;
}




CMatrix operator/(const CMatrix &m, const Complex &a)
{
  int i, j;
  
  CMatrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);
  
  for (i=m.rlow; i<=m.rhigh; i++) 
    {
      for (j=m.clow; j<=m.chigh; j++)
	{
	  tmp.rows[i].elements[j] = m.rows[i].elements[j] / a;
	}
    }
  return tmp;
}


CMatrix &CMatrix::operator*=(const Complex &a)
{
  int i, j;
  
  for (i=rlow; i<=rhigh; i++) 
    {
      for (j=clow; j<=chigh; j++)
	{
	  rows[i].elements[j] *= a;
	}
    }
  return *this;
}

CMatrix &CMatrix::operator/=(const Complex &a)
{
  int i, j;
  
  for (i=rlow; i<=rhigh; i++) 
    {
      for (j=clow; j<=chigh; j++)
	{
	  rows[i].elements[j] /= a;
	}
    }
  return *this;
}


CMatrix &CMatrix::operator+=(const CMatrix &m)
{
  int i, j;
  
  if (m.rlow!=rlow || m.rhigh!=rhigh || m.clow!=clow || m.chigh!=chigh) 
    bomb_CMatrix_operation("+=");
  
  
  for (i=rlow; i<=rhigh; i++) 
    {
      for (j=clow; j<=chigh; j++)
	{
	  rows[i].elements[j] += m.rows[i].elements[j];
	}
    }
  return *this;
}


CMatrix &CMatrix::operator-=(const CMatrix &m)
{
  int i, j;
  
  if (m.rlow!=rlow || m.rhigh!=rhigh || m.clow!=clow || m.chigh!=chigh) 
    bomb_CMatrix_operation("+=");
  
  
  for (i=rlow; i<=rhigh; i++) 
    {
      for (j=clow; j<=chigh; j++)
	{
	  rows[i].elements[j] -= m.rows[i].elements[j];
	}
    }
  return *this;
}


CVector operator*(const CMatrix &m, const CVector &v)
{
  int i, j;
  
  if (m.clow != v.low || m.chigh != v.high) bomb_CMatrix_operation("M*v");
  
  CVector tmp(m.rlow, m.rhigh);
  
  for (i=m.rlow; i<=m.rhigh; i++)
    {
      tmp.elements[i] = 0.0;
      for (j=m.clow; j<=m.chigh; j++)
	{
	  tmp.elements[i] += 
	    m.rows[i].elements[j]*v.elements[j];
	}
    }
  return tmp;
}

CVector operator*(const CMatrix &m, const Vector &v)
{
  int i, j;
  
  if (m.clow != v.low || m.chigh != v.high) bomb_CMatrix_operation("M*v");
  
  CVector tmp(m.rlow, m.rhigh);
  
  for (i=m.rlow; i<=m.rhigh; i++)
    {
      tmp.elements[i] = 0.0;
      for (j=m.clow; j<=m.chigh; j++)
	{
	  tmp.elements[i] += 
	    m.rows[i].elements[j]*v.elements[j];
	}
    }
  return tmp;
}

CVector operator*(const CVector &v, const CMatrix &m)
{
  int i, j;
  
  if (m.rlow != v.low || m.rhigh != v.high) bomb_CMatrix_operation("v*M");
  
  CVector tmp(m.clow, m.chigh);
  
  for (j=m.clow; j<=m.chigh; j++)
    {
      tmp.elements[j] = 0.0;
      for (i=m.rlow; i<=m.rhigh; i++) 
	{
	  tmp.elements[j] +=
	    m.rows[i].elements[j]*v.elements[i];
	}
    }
  return tmp;
}


CMatrix operator*(const CMatrix &m1, const CMatrix &m2)
{
  int i, j, k;
  
  if (m1.clow != m2.rlow || m1.chigh != m2.rhigh)
    bomb_CMatrix_operation("M*M");
  
  CMatrix tmp(m1.rlow, m1.rhigh, m2.clow, m2.chigh);
  
  for (i=m1.rlow; i<=m1.rhigh; i++)
    {
      for (j=m2.clow; j<=m2.chigh; j++)
	{
	  tmp.rows[i].elements[j] = 0.0;
	  for (k=m1.clow; k<=m1.chigh; k++) 
	    {
	      tmp.rows[i].elements[j] += m1.rows[i].elements[k]
		*m2.rows[k].elements[j];
	    }
	}
    }
  return tmp;
  
}


CMatrix CMatrix::Conjg(void) const
{
  int i;
  CMatrix tmp(rlow, rhigh, clow, chigh);
  
  for (i=rlow; i<=rhigh; i++) 
    {
      tmp.rows[i] = rows[i].Conjg();
    }
  
  return tmp;
}





CMatrix CMatrix::Transpose(void) const
{
  int i;
  CMatrix t(clow, chigh, rlow, rhigh);
  
  t.rlow = clow;
  t.rhigh = chigh;
  t.chigh = rhigh;
  t.clow = rlow;
  
  for (i=t.rlow; i<=t.rhigh; i++)
    {
      t.fastsetrow(i, fastcol(i));
    }
  
  return t;
}

Complex CMatrix::Trace(void) const
{
  Complex t;
  int i;
  
  if (rlow != clow || rhigh != chigh) 
    bomb_CVector("cannot take trace of non-square matrix");
  
  
  t = 0.0;
  for (i=rlow; i<=rhigh; i++) t += rows[i][i];
  
  return t;
}



Matrix CMatrix::Re(void) const
{
  int i, j;
  Matrix tmp(rlow, rhigh, clow, chigh);
  
  for (i=rlow; i<=rhigh; i++) 
    {
      for (j=clow; j<=chigh; j++) tmp.rows[i].elements[j] = rows[i].elements[j].real();
    }
  return tmp;
}	


Matrix CMatrix::Im(void) const
{
  int i, j;
  Matrix tmp(rlow, rhigh, clow, chigh);
  
  for (i=rlow; i<=rhigh; i++) 
    {
      for (j=clow; j<=chigh; j++) tmp.rows[i].elements[j] = rows[i].elements[j].imag();
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

Complex Trace(const CMatrix &c)
{
  return c.Trace();
}







void CMatrix::print(ostream& out)
{
  int i;
  
  for (i=rlow; i<=rhigh; i++) rows[i].print(out);
}



void CMatrix::binwrite(ostream& out)
{
  int i;
  
  out.write((char *)&rlow, sizeof(int));
  out.write((char *)&rhigh, sizeof(int));
  out.write((char *)&clow, sizeof(int));
  out.write((char *)&chigh, sizeof(int));
  
  for (i=rlow; i<=rhigh; i++)
    {
      rows[i].binwrite(out);
    }
}


CMatrix CMatrix_binread(istream& in)
{
  
  int rlow, rhigh, clow, chigh;
  int i;
  
  in.read((char *)&rlow, sizeof(int));
  in.read((char *)&rhigh, sizeof(int));
  in.read((char *)&clow, sizeof(int));
  in.read((char *)&chigh, sizeof(int));
  
  CMatrix tmp(rlow, rhigh, clow, chigh);
  
  for (i=rlow; i<=rhigh; i++)
    {
      tmp.rows[i] = CVector_binread(in);
    }
  
  return tmp;
}

