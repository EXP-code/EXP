#include <assert.h>
#include <stdlib.h>
#include <unistd.h>

#include <mpi.h>

#include <Vector.h>
#include <kevin_complex.h>

#include <sstream>
#include <string>

#include <pthread.h>

using namespace std;

extern char threading_on;
extern pthread_mutex_t mem_lock;

extern int myid;

/*
  Default constructor; make a null vector.
*/

Vector::Vector()
{
  low  = 0;
  high = 0;
  size = 0;
}


Vector::Vector(int l, int h)
{
  low  = 0;
  high = 0;
  setsize(l, h);
}

Vector::Vector(int l, int h, double *v)
{
  low  = 0;
  high = 0;
  setsize(l, h);
  
  for (int i=0; i<size; i++) pelement[i] = v[i+low];
}

/*
  Copy constructor; create a new vector which is a copy of another.
*/

Vector::Vector(const Vector &v)
{
  // Allow consruction of null vector
  if (v.elements.size() == 0) {
    low = high = size = 0;
    pelement = 0;			// Null pointer
  } else {
    low      = v.low;
    high     = v.high;
    size     = high - low + 1;
    elements = v.elements;
    pelement = &elements[0];
  }
}

/* 
   conversion operator (Three_Vector)
*/

Vector::operator Three_Vector(void)
{
  if (low!=1 || high != 3) {
    bomb_Vector("Vector->3Vector conversion");
  }
  
  for (int i=0; i<3; i++) (*this)[i+1] = pelement[i];
  
  return *this;
}


/*
  Assignment operator for Vector; must be defined as a reference
  so that it can be used on lhs. Compatibility checking is performed;
  the destination vector is allocated if its elements are undefined.
*/


Vector &Vector::operator=(const Vector &v)
{
  // Allow assignment of null vector
  if (v.elements.size() == 0) {
    low = high = size = 0;
    return *this;
  }
  
  if (low!=v.low && high!=v.high 
      && (high != 0 || elements.size() != 0)) {
    bomb_Vector_operation("=");
  }
  
  low      = v.low;
  high     = v.high;
  size     = high - low + 1;
  elements = v.elements;
  pelement = &elements[0];
  
  return *this;
}


double &Vector::operator[](int i) const
{
  if (i<low || i>high) {
    cerr << "VECTOR ERROR: " 
	 << low << " < " << i << " < " << high << "\n";
    bomb_Vector("subscript out of range");
  }

  return pelement[i-low];
}


void Vector::setsize(int l, int h)
{
  // do we need to resize at all?
  
  if (l==low && h==high) return;
  
  
  // is the requested size positive?
  
  if (h<l) {
    ostringstream mesg;
    mesg << "invalid size <0, l=" << l << " h=" << h;
    bomb_Vector(mesg.str());
  }
  
  // set the new size
  
  low  = l;
  high = h;
  size = high - low + 1;
  
  // allocate the new elements, and make a pointer
  
  elements = vector<double>(size, 0.0);
  pelement = &elements[0];
}

void Vector::load(double *array)
{
  for (int i=0; i<size; i++) pelement[i] = array[i+low];
}


void Vector::zero(void)
{
  for (int i=0; i<size; i++) pelement[i] = 0.0;
}

double *Vector::array(int l, int h)
{
  return pelement - l;
}

double *Vector::array_copy(int l, int h)
{
  if (threading_on) pthread_mutex_lock(&mem_lock);
  double *ptr = new double [h-l+1];
  if (threading_on) pthread_mutex_unlock(&mem_lock);
  assert(ptr);
  ptr -= l;
  
  for (int i=l; i<=h; i++) ptr[i] = pelement[i-l];
  return ptr;
}


void destroy_Vector_array(double *p, int l, int h)
{
  if (threading_on) pthread_mutex_lock(&mem_lock);
  delete [] &p[l];
  if (threading_on) pthread_mutex_unlock(&mem_lock);
}


Vector Vector::operator-(void)
{
  Vector tmp(low, high);
  
  for (int i=0; i<size; i++) tmp.pelement[i] = -pelement[i];
  
  return tmp;
}




Vector operator+(const Vector &v1, const Vector &v2)
{
  if (v1.low != v2.low || v1.high != v2.high) bomb_Vector_operation("+");
  
  Vector tmp(v1.low, v1.high);
  
  for (int i=0; i<v1.size; i++)
    tmp.pelement[i] = v1.pelement[i] + v2.pelement[i];
  
  return tmp;
}

Vector operator-(const Vector &v1, const Vector &v2)
{
  if (v1.low != v2.low || v1.high != v2.high) bomb_Vector_operation("-");
  
  Vector tmp(v1.low, v1.high);
  
  for (int i=0; i<v1.size; i++)
    tmp.pelement[i] = v1.pelement[i] - v2.pelement[i];
  
  return tmp;
}

Vector operator^(const Vector &v1, const Vector &v2)
{
  if (v1.high != 3 || v2.high != 3) {
    bomb_Vector_operation("^");
  }

  Vector tmp(1, 3);
  
  tmp[1] = v1[2]*v2[3] - v1[3]*v2[2];
  tmp[2] = v1[3]*v2[1] - v1[1]*v2[3];
  tmp[3] = v1[1]*v2[2] - v1[2]*v2[1];
  
  return tmp;
}

Vector operator&(const Vector &v1, const Vector &v2)
{
  if (v1.high != v2.high || v1.low != v2.low) {
    bomb_Vector_operation("&");
  }
  
  Vector tmp(v1.low, v1.high);
  for (int i=v1.low; i<=v1.high; i++) tmp[i] = v1[i]*v2[i];
  
  return tmp;
}

Vector operator*(double a, const Vector &v)
{
  Vector tmp(v.low, v.high);
  for (int i=0; i<=v.high-v.low; i++) tmp.pelement[i] = a * v.pelement[i];
  
  return tmp;
}

Vector operator*(const Vector &v, double a)
{
  Vector tmp(v.low, v.high);
  for (int i=0; i<=v.high-v.low; i++) tmp.pelement[i] = a * v.pelement[i];
  
  return tmp;
}

Vector operator/(const Vector &v, double a)
{
  Vector tmp(v.low, v.high);
  for (int i=0; i<=v.high-v.low; i++) tmp.pelement[i] = v.pelement[i]/a;
  
  return tmp;
}

Vector &Vector::operator+=(double a)
{
  for (int i=0; i<size; i++) pelement[i] += a;
  
  return *this;
}

Vector &Vector::operator-=(double a)
{
  for (int i=0; i<size; i++) pelement[i] -= a;
  
  return *this;
}


Vector &Vector::operator+=(const Vector &v)
{
  if (low != v.low || high != v.high) bomb_Vector_operation("+=");
  
  for (int i=0; i<size; i++) pelement[i] += v.pelement[i];
  
  return *this;
}	

Vector &Vector::operator-=(const Vector &v)
{
  if (low != v.low || high != v.high) bomb_Vector_operation("-=");
  
  for (int i=0; i<size; i++) pelement[i] -= v.pelement[i];
  
  return *this;
}	

Vector &Vector::operator*=(double a)
{
  for (int i=0; i<size; i++) pelement[i] *= a;
  
  return *this;
}

Vector &Vector::operator/=(double a)
{
  for (int i=0; i<size; i++) pelement[i] /= a;
  
  return *this;
}


double operator*(const Vector &v1, const Vector &v2)
{
  if (v1.low != v2.low || v1.high != v2.high) bomb_Vector_operation("*");
  
  double tmp = 0.0;
  for (int i=0; i<v1.size; i++) tmp += v1.pelement[i] * v2.pelement[i];

  return tmp;
}



void Vector::print(ostream& out)
{
  for (int i=0; i<size; i++) out << pelement[i] << " ";
  out << endl;
}

void Vector::binwrite(ostream& out)
{
  out.write((char *)&low,  sizeof(int));
  out.write((char *)&high, sizeof(int));
  out.write((char *)pelement, size*sizeof(double));
}

Vector Vector_binread(istream& in)
{
  int low, high;
  
  in.read((char *)&low,  sizeof(int));
  in.read((char *)&high, sizeof(int));
  
  Vector tmp(low, high);
  
  in.read((char *)tmp.pelement, tmp.size*sizeof(double));
  
  return tmp;
}


void bomb_Vector(const string& msg)
{
  if (myid>=0)
    cerr << "VECTOR ERROR [mpi_id=" << myid << "]: " << msg << '\n';
  else
    cerr << "VECTOR ERROR: " << msg << '\n';
#ifdef DEBUG
  chdir("/tmp");
  abort();
#endif
  exit(0);
}


void bomb_Vector_operation(const string& op)
{
  string msg("incompatible lengths in operation ");
  msg += op;
  bomb_Vector(msg);
}



Matrix::Matrix()
{
  rlow  = 0;
  rhigh = 0;
  clow  = 0;
  chigh = 0;
  rsize = 0;
  csize = 0;
}


Matrix::Matrix(int rl, int rh, int cl, int ch)
{
  rlow  = 0;
  rhigh = 0;
  clow  = 0;
  chigh = 0;
  setsize(rl, rh, cl, ch);
}


Matrix::Matrix(int rl, int rh, int cl, int ch, double **array)
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


Matrix::Matrix(const Matrix &m)
{
  // Allow construction of null matrix
  if (m.rows.size() == 0) {
    rlow = rhigh = clow = chigh = 0;
  }
  else {
    
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
}	


void Matrix::setsize(int rl, int rh, int cl, int ch)
{
  // do we need to resize at all?
  
  if (rl==rlow && cl==clow && rh==rhigh && ch==chigh) return;
  
  // is the new size positive?
  
  if (rh<rl || ch<cl) {
    ostringstream mesg;
    mesg << "invalid size <0, rl=" << rl << " rh=" << rh
	 << " cl=" << cl << " ch" << ch;
    bomb_Matrix(mesg.str());
  }
  
  
  // set the new size
  
  rlow  = rl;
  rhigh = rh;
  clow  = cl;
  chigh = ch;

  rsize = rhigh - rlow + 1;
  csize = chigh - clow + 1;
  
  // allocate the array of rows
  
  rows = vector<Vector>(rsize);
  prow = &rows[0];
  
  // create the individual rows
  
  for (int i=0; i<rsize; i++) rows[i].setsize(clow, chigh);
}


void bomb_Matrix(const string& msg)
{
  if (myid>=0)
    cerr << "MATRIX ERROR [mpi_id=" << myid << "]: " << msg << '\n';
  else
    cerr << "MATRIX ERROR: " << msg << '\n';
#ifdef DEBUG
  abort();
#endif
  exit(0);
}

void bomb_Matrix_operation(const string& op)
{
  string str("incompatible sizes in operator ");
  str += op;
  bomb_Matrix(str);
}

Vector Matrix::fastcol(int j)
{
  Vector tmp(rlow, rhigh);
  for (int i=rlow; i<rhigh; i++) tmp[i] = rows[i-rlow].pelement[j-clow];
  
  return tmp;
}

Vector Matrix::col(int j)
{
  if (j<clow || j>chigh) bomb_Matrix("column subscript out of range");
  
  return fastcol(j);
}



Vector &Matrix::row(int i)
{
  if (i<rlow || i>rhigh) {
    ostringstream out;
    out << "row subscript out of range"
	<< ": i=" << i << " rlow=" << rlow << " rhigh=" << rhigh;
    bomb_Matrix(out.str());
  }
  
  return rows[i-rlow];
}

Vector &Matrix::fastrow(int i)
{
  return rows[i-rlow];
}



void Matrix::fastsetrow(int i, const Vector &v)
{
  rows[i-rlow] = v;
}

void Matrix::setrow(int i, Vector &v)
{
  if (v.low != clow || v.high != chigh) 
    bomb_Matrix("row-vector size mismatch");
  
  if (i<rlow || i>rhigh) {
    ostringstream out;
    out << "row subscript out of range"
	<< ": i=" << i << " rlow=" << rlow << " rhigh=" << rhigh;
    bomb_Matrix(out.str());
  }

  for (int j=0; j<csize; j++)
    rows[i-rlow].pelement[j] = v.pelement[j];
}

void Matrix::fastsetcol(int j, Vector &v)
{
  for (int i=0; i<rsize; i++) {
    rows[i].pelement[j-clow] = v.pelement[i];
  }
}

void Matrix::setcol(int j, Vector &v)
{
  if (v.low != rlow || v.high != rhigh) 
    bomb_Matrix("column-vector size mismatch");
  if (j<clow || j>chigh) bomb_Matrix("column subscript out of range");
  
  fastsetcol(j, v);
}

void Matrix::zero(void)
{
  for (int i=0; i<rsize; i++) rows[i].zero();
}


double **Matrix::array(int rl, int rh, int cl, int ch)
{
  if (threading_on) pthread_mutex_lock(&mem_lock);
  double **ptr = new double* [rh-rl+1];
  if (threading_on) pthread_mutex_unlock(&mem_lock);
  if (ptr == NULL) bomb_Matrix("could not allocate array");
  ptr -= rl;
  
  for (int i=rl; i<=rh; i++) ptr[i] = rows[i-rl].array(cl, ch);
  
  return ptr;
}

double **Matrix::array_copy(int rl, int rh, int cl, int ch)
{
  if (threading_on) pthread_mutex_lock(&mem_lock);
  double **ptr = new double* [rh-rl+1];
  if (threading_on) pthread_mutex_unlock(&mem_lock);
  if (ptr == NULL) bomb_Matrix("could not allocate array");
  ptr -= rl;
  
  for (int i=rl; i<=rh; i++) ptr[i] = rows[i-rl].array_copy(cl, ch);
  
  return ptr;
}

void destroy_Matrix_array(double **ptr, int rl, int rh, int cl, int ch)
{
  for (int i=rl; i<=rh; i++) destroy_Vector_array(ptr[i], cl, ch);
  
  if (threading_on) pthread_mutex_lock(&mem_lock);
  delete [] &ptr[rl];
  if (threading_on) pthread_mutex_unlock(&mem_lock);
}



Vector &Matrix::operator[](int i) const
{
  if (i<rlow || i>rhigh) {
    ostringstream out;
    out << "row subscript out of range"
	<< ": i=" << i << " rlow=" << rlow << " rhigh=" << rhigh;
    bomb_Matrix(out.str());
  }
  return prow[i-rlow];
}

Matrix &Matrix::operator=(const Matrix &m)
{
  // Allow assignment of null matrix
  if (m.rows.size()==0) {
    rlow = rhigh = clow = chigh = 0;
    rsize = csize = 0;
    return *this;
  }
  
  setsize(m.rlow, m.rhigh, m.clow, m.chigh);
  
  if (m.rlow!=rlow || m.rhigh!=rhigh || m.clow!=clow || m.chigh!=chigh)  {
    ostringstream sout;
    sout << "=" 
	 << " rlow[" << m.rlow << "|" << rlow << "]"
	 << " rhigh[" << m.rhigh << "|" << rhigh << "]"
	 << " clow[" << m.clow << "|" << clow << "]"
	 << " chigh[" << m.chigh << "|" << chigh << "]";
    bomb_Matrix_operation(sout.str());
  }
  
  for (int i=rlow; i<=rhigh; i++)  fastsetrow(i, m.rows[i-rlow]);

  return *this;
}

Matrix Matrix::operator-(void)
{
  for (int i=0; i<rsize; i++) {
    for (int j=0; j<csize; j++) {
      rows[i].pelement[j] = -rows[i].pelement[j];
    }
  }
  
  return *this;
}

Matrix operator+(const Matrix &m1, const Matrix &m2)
{
  if (m1.rlow!=m2.rlow || m1.rhigh!=m2.rhigh || 
      m1.clow!=m2.clow || m1.chigh!=m2.chigh) 
    bomb_Matrix_operation("+");
  
  Matrix tmp(m1.rlow, m1.rhigh, m1.clow, m1.chigh);
  
  for (int i=0; i<m1.rsize; i++) {
    for (int j=0; j<m1.csize; j++) {
      tmp.rows[i].pelement[j] = 
	m1.rows[i].pelement[j] + m2.rows[i].pelement[j];
    }
  }
  
  return tmp;
}

Matrix operator-(const Matrix &m1, const Matrix &m2)
{
  if (m1.rlow!=m2.rlow || m1.rhigh!=m2.rhigh || 
      m1.clow!=m2.clow || m1.chigh!=m2.chigh) 
    bomb_Matrix_operation("-");
  
  Matrix tmp(m1.rlow, m1.rhigh, m1.clow, m1.chigh);
  
  for (int i=0;i<m1.rsize; i++) {
    for (int j=0; j<m1.csize; j++) {
      tmp.rows[i].pelement[j] = 
	m1.rows[i].pelement[j] - m2.rows[i].pelement[j];
    }
  }
  
  return tmp;
}


Matrix operator*(const Matrix &m, double a)
{
  Matrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);
  
  for (int i=0; i<m.rsize; i++) {
    for (int j=0; j<m.csize; j++) {
      tmp.rows[i].pelement[j] = m.rows[i].pelement[j] * a;
    }
  }
  return tmp;
}


Matrix operator*(double a, const Matrix &m)
{
  Matrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);
  
  for (int i=0; i<m.rsize; i++) {
    for (int j=0; j<m.csize; j++) {
      tmp.rows[i].pelement[j] = m.rows[i].pelement[j] * a;
    }
  }
  return tmp;
}


Matrix operator/(const Matrix &m, double a)
{
  Matrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);
  
  for (int i=0; i<m.rsize; i++) {
    for (int j=0; j<m.csize; j++) {
      tmp.rows[i].pelement[j] = m.rows[i].pelement[j] / a;
    }
  }
  return tmp;
}


Matrix operator+(const Matrix &m, double a)
{
  Matrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);
  
  for (int i=0; i<m.rsize; i++) {
    for (int j=0; j<m.csize; j++) {
      tmp.rows[i][j] += a;
    }
  }
  return tmp;
}

Matrix operator-(const Matrix &m, double a)
{
  Matrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);
  
  for (int i=0; i<m.rsize; i++) {
    for (int j=0; j<m.csize; j++) {
      tmp.rows[i][j] -= a;
    }
  }
  return tmp;
}


Matrix operator+(double a, const Matrix &m)
{
  Matrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);
  
  for (int i=0; i<m.rsize; i++) {
    for (int j=0; j<m.csize; j++) {
      tmp.rows[i][j] += a;
    }
  }
  return tmp;
}

Matrix operator-(double a, const Matrix &m)
{
  Matrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);
  
  for (int i=0; i<m.rsize; i++) {
    for (int j=0; j<m.csize; j++) {
      tmp.rows[i][j] -= a;
    }
  }
  return tmp;
}


Matrix &Matrix::operator*=(double a)
{
  for (int i=0; i<rsize; i++) {
    for (int j=0; j<csize; j++) {
      rows[i].pelement[j] *= a;
    }
  }
  return *this;
}

Matrix &Matrix::operator/=(double a)
{
  for (int i=0; i<rsize; i++) {
    for (int j=0; j<csize; j++) {
      rows[i].pelement[j] /= a;
    }
  }
  return *this;
}

Matrix &Matrix::operator+=(double a)
{
  for (int i=0; i<rsize; i++) {
    for (int j=0; j<csize; j++) {
      rows[i].pelement[j] += a;
    }
  }
  return *this;
}

Matrix &Matrix::operator-=(double a)
{
  for (int i=0; i<rsize; i++) {
    for (int j=0; j<csize; j++) {
      rows[i].pelement[j] -= a;
    }
  }
  return *this;
}


Matrix &Matrix::operator+=(const Matrix &m)
{
  if (m.rlow!=rlow || m.rhigh!=rhigh || m.clow!=clow || m.chigh!=chigh) 
    bomb_Matrix_operation("+=");
  
  
  for (int i=0; i<rsize; i++) {
    for (int j=0; j<csize; j++) {
      rows[i].pelement[j] += m.rows[i].pelement[j];
    }
  }
  return *this;
}


Matrix &Matrix::operator-=(const Matrix &m)
{
  if (m.rlow!=rlow || m.rhigh!=rhigh || m.clow!=clow || m.chigh!=chigh) 
    bomb_Matrix_operation("+=");
  
  
  for (int i=0; i<rsize; i++) {
    for (int j=0; j<csize; j++) {
      rows[i].pelement[j] -= m.rows[i].pelement[j];
    }
  }
  return *this;
}


double operator^(const Matrix &m1, const Matrix &m2)
{
  if (m1.getrlow()!=m2.getrlow() || m1.getrhigh()!=m2.getrhigh() || 
      m1.getclow()!=m2.getclow() || m1.getchigh()!=m2.getchigh()) 
    bomb_Matrix_operation("^");
  
  double tmp = 0.0;
  
  for (int i=0; i<m1.rsize; i++) {
    for (int j=0; j<m1.csize; j++) {
      tmp +=  m1[i][j]*m2[i][j];
    }
  }
  
  return tmp;
}

Three_Vector operator*(const Three_Vector &v, const Matrix &m)
{
  if (m.clow != 1 || m.chigh != 3 || m.rlow != 1 || m.rhigh != 3) 
    bomb_Matrix_operation("3v*M");
  
  Three_Vector tmp;
  
  for (int j=0; j<3; j++) {
    tmp.x[j] = 0.0;
    for (int i=0; i<3; i++) {
      tmp.x[j] += m.rows[i].pelement[j] * v.x[j];
    }
  }
  return tmp;
}


Three_Vector operator*(const Matrix &m, const Three_Vector &v)
{
  if (m.clow != 1 || m.chigh != 3 || m.rlow != 1 || m.rhigh != 3) 
    bomb_Matrix_operation("3v*M");
  
  Three_Vector tmp;
  
  for (int i=0; i<3; i++) {
    tmp.x[i] = 0.0;
    for (int j=0; j<3; j++) {
      tmp.x[i] += m.rows[i].pelement[j] * v.x[j];
    } 
  }
  return tmp;
}

Vector operator*(const Matrix &m, const Vector &v)
{
  if (m.clow != v.low || m.chigh != v.high) bomb_Matrix_operation("M*v");
  
  Vector tmp(m.rlow, m.rhigh);
  
  for (int i=0; i<m.rsize; i++) {
    tmp.pelement[i] = 0.0;
    for (int j=0; j<m.csize; j++) {
      tmp.pelement[i] += m.rows[i].pelement[j] * v.pelement[j];
    }
  }
  return tmp;
}

Vector operator*(const Vector &v, const Matrix &m)
{
  if (m.rlow != v.low || m.rhigh != v.high) bomb_Matrix_operation("v*M");
  
  Vector tmp(m.clow, m.chigh);
  
  for (int j=0; j<m.csize; j++) {
    tmp.pelement[j] = 0.0;
    for (int i=0; i<m.rsize; i++)  {
      tmp.pelement[j] += m.rows[i].pelement[j] * v.pelement[i];
    }
  }
  return tmp;
}


Matrix operator*(const Matrix &m1, const Matrix &m2)
{
  if (m1.clow != m2.rlow || m1.chigh != m2.rhigh)
    bomb_Matrix_operation("M*M");
  
  Matrix tmp(m1.rlow, m1.rhigh, m2.clow, m2.chigh);
  
  for (int i=0; i<m1.rsize; i++) {
    for (int j=0; j<m2.csize; j++) {
      tmp.rows[i].pelement[j] = 0.0;
      for (int k=0; k<m1.csize; k++) {
	tmp.rows[i].pelement[j] +=
	  m1.rows[i].pelement[k] * m2.rows[k].pelement[j];
      }
    }
  }
  return tmp;
  
}

Matrix Matrix::Transpose(void)
{
  Matrix t(clow, chigh, rlow, rhigh);
  
  
  for (int i=t.rlow; i<=t.rhigh; i++) {
    t.fastsetrow(i, fastcol(i));
  }
  
  return t;
}

double Matrix::Trace(void)
{
  double t = 0.0;
  for (int i=0; i<rsize; i++) t += rows[i][i];
  return t;
}

double Trace(Matrix &m)
{
  return m.Trace();
}

Matrix Transpose(Matrix &m)
{
  return m.Transpose();
}


void Matrix::print(ostream& out)
{
  for (int i=0; i<rsize; i++) rows[i].print(out);
}

void Matrix::binwrite(ostream& out)
{
  out.write((char *)&rlow,  sizeof(int));
  out.write((char *)&rhigh, sizeof(int));
  out.write((char *)&clow,  sizeof(int));
  out.write((char *)&chigh, sizeof(int));
  
  for (int i=0; i<rsize; i++) rows[i].binwrite(out);
}


Matrix Matrix_binread(istream& in)
{
  
  int rlow, rhigh, clow, chigh;
  
  in.read((char *)&rlow,  sizeof(int));
  in.read((char *)&rhigh, sizeof(int));
  in.read((char *)&clow,  sizeof(int));
  in.read((char *)&chigh, sizeof(int));
  
  Matrix tmp(rlow, rhigh, clow, chigh);
  
  for (int i=0; i<tmp.rsize; i++)
    tmp.rows[i] = Vector_binread(in);
  
  return tmp;
}


inline double fabs(Vector& v) {return sqrt(v*v);}

