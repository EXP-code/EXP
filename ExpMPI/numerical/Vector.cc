#pragma implementation

#include <assert.h>
#include <stdlib.h>
#include <unistd.h>
#include <Vector.h>
#include <kevin_complex.h>
#include <string>

#include <pthread.h>

extern char threading_on;
extern pthread_mutex_t mem_lock;

#include "localmpi.h"

/*
	Default constructor; make a null vector.
*/

Vector::Vector()
{
	low=0;
	high=0;
	elements = NULL;

}


Vector::Vector(int l, int h)
{
	low=0;
	high=0;
	elements = NULL;
	setsize(l, h);
}

Vector::Vector(int l, int h, double *v)
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

Vector::Vector(const Vector &v)
{
	int i;

				// Allow consruction of null vector
	if (v.elements==NULL) {
	  elements = NULL;
	  low = high = 0;
	}
	else {

	  low=0;
	  high=0;
	  elements = NULL;
	  setsize(v.low, v.high);
	
	  for (i=low; i<=high; i++) elements[i] = v.elements[i];
	}
}




/*
	Destructor. Free elements if it exists.
*/

Vector::~Vector()
{
	elements += low;
	//	if (elements != NULL) delete [] elements;
	if (threading_on) pthread_mutex_lock(&mem_lock);
	delete [] elements;
	if (threading_on) pthread_mutex_unlock(&mem_lock);
	//	else 
	//	{
	//		puts("WARNING: destructor called for NULL Vector elements");
	//	}
}


/* 
	conversion operator (Three_Vector)
*/


Vector::operator Three_Vector(void)
{
	if (low!=1 || high != 3)
	{
		bomb_Vector("Vector->3Vector conversion");
	}

//	static Three_Vector tmp;
	int i;

//	for (i=1; i<=3; i++) tmp[i] = elements[i];
	for (i=1; i<=3; i++) (*this)[i] = elements[i];

	return *this;
}


/*
	conversion operator (CVector)
*/

/*
Vector::operator CVector(void)
{
	int i;
	CVector tmp;
	tmp.setsize(low, high);

	for (i=low; i<=high; i++) tmp[i] = elements[i];

	return tmp;
}
*/


/*
	Assignment operator for Vector; must be defined as a reference
so that it can be used on lhs. Compatibility checking is performed;
the destination vector is allocated if its elements are undefined.
*/


Vector &Vector::operator=(const Vector &v)
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
		bomb_Vector_operation("=");
	}

	if (high == 0 && low==0 && elements == NULL)
	{
		setsize(v.low, v.high);
	}

	for (i=low; i<=high; i++) elements[i] = v.elements[i];

	return *this;
}

	
double &Vector::operator[](int i) const
{
	if (i<low || i>high)
	{
	        cerr << "VECTOR ERROR: " 
		<< low << " < " << i << " < " << high << "\n";
		bomb_Vector("subscript out of range");
	}

	return elements[i];
}


void Vector::setsize(int l, int h)
{
	/* do we need to resize at all? */

	if (l==low && h==high) return;



	
	/* is the requested size positive? */

	if (h<l) bomb_Vector("invalid size <0");



	/* delete the old elements if they already exist */

	elements += low;
	// if (elements != NULL) delete [] elements;
	if (threading_on) pthread_mutex_lock(&mem_lock);
	delete [] elements;
	if (threading_on) pthread_mutex_unlock(&mem_lock);



	/* set the new size */

	low = l;
	high = h;



	/* allocate the new elements, and offset the pointer */
	
	if (threading_on) pthread_mutex_lock(&mem_lock);
	elements = new double[high-low+1];
	if (threading_on) pthread_mutex_unlock(&mem_lock);
	if (elements == NULL) bomb_Vector("could not allocate");
	elements -= low;
}

void Vector::load(double *array)
{
	int i;

	for (i=low; i<=high; i++) elements[i] = array[i];
}


void Vector::zero(void)
{
	int i;

	for (i=low; i<=high; i++) elements[i] = 0.0;
}

double *Vector::array(int l, int h)
{
	return elements+low-l;
}
	
double *Vector::array_copy(int l, int h)
{
	int i;
	double *ptr;

	if (threading_on) pthread_mutex_lock(&mem_lock);
	ptr = new double[h-l+1];
	if (threading_on) pthread_mutex_unlock(&mem_lock);
	assert(ptr);
	ptr-=l;

	for (i=low; i<=high; i++) ptr[i] = elements[i];
	return ptr;
}


void destroy_Vector_array(double *p, int l, int h)
{
	p += l;
	if (threading_on) pthread_mutex_lock(&mem_lock);
	delete [] p;
	if (threading_on) pthread_mutex_unlock(&mem_lock);
}


Vector Vector::operator-(void)
{
	int i;
	
	Vector tmp(low, high);
	
	for (i=low; i<=high; i++) tmp.elements[i] = -elements[i];

	return tmp;
}
	
	


Vector operator+(const Vector &v1, const Vector &v2)
{
	int i;
	
	if (v1.low != v2.low || v1.high != v2.high) bomb_Vector_operation("+");

	Vector tmp(v1.low, v1.high);

	for (i=v1.low; i<=v1.high; i++)
	{
		tmp.elements[i] = v1.elements[i] + v2.elements[i];
	}

	return tmp;
}

Vector operator-(const Vector &v1, const Vector &v2)
{
	int i;
	
	if (v1.low != v2.low || v1.high != v2.high) bomb_Vector_operation("-");

	Vector tmp(v1.low, v1.high);

	for (i=v1.low; i<=v1.high; i++)
	{
		tmp.elements[i] = v1.elements[i] - v2.elements[i];
	}

	return tmp;
}

Vector operator^(const Vector &v1, const Vector &v2)
{
	if (v1.high != 3 || v2.high != 3)
	{
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
        int i;

	if (v1.high != v2.high || v1.low != v2.low)
	{
		bomb_Vector_operation("&");
	}

	Vector tmp(v1.low, v1.high);
	for (i=v1.low; i<=v1.high; i++)
	  tmp[i] = v1[i]*v2[i];

	return tmp;
}

Vector operator*(double a, const Vector &v)
{
	int i;

	Vector tmp(v.low, v.high);
	for (i=v.low; i<=v.high; i++) tmp.elements[i] = a * v.elements[i];

	return tmp;
}

Vector operator*(const Vector &v, double a)
{
	int i;

	Vector tmp(v.low, v.high);
	for (i=v.low; i<=v.high; i++) tmp.elements[i] = a * v.elements[i];

	return tmp;
}

Vector operator/(const Vector &v, double a)
{
	int i;

	Vector tmp(v.low, v.high);
	for (i=v.low; i<=v.high; i++) tmp.elements[i] = v.elements[i]/a;

	return tmp;
}

Vector &Vector::operator+=(double a)
{
	for (int i=low; i<=high; i++) elements[i] += a;

	return *this;
}

Vector &Vector::operator-=(double a)
{
	for (int i=low; i<=high; i++) elements[i] -= a;

	return *this;
}


Vector &Vector::operator+=(const Vector &v)
{
	int i;
	
	if (low != v.low || high != v.high) bomb_Vector_operation("+=");

	for (i=low; i<=high; i++) elements[i] += v[i];

	return *this;
}	

Vector &Vector::operator-=(const Vector &v)
{
	int i;
	
	if (low != v.low || high != v.high) bomb_Vector_operation("-=");

	for (i=low; i<=high; i++) elements[i] -= v[i];

	return *this;
}	

Vector &Vector::operator*=(double a)
{
	int i;

	for (i=low; i<=high; i++) elements[i] *= a;

	return *this;
}

Vector &Vector::operator/=(double a)
{
	int i;

	for (i=low; i<=high; i++) elements[i] /= a;

	return *this;
}



double operator*(const Vector &v1, const Vector &v2)
{
	int i;
	double tmp;
	
	if (v1.low != v2.low || v1.high != v2.high) bomb_Vector_operation("*");

	tmp = 0.0;
	for (i=v1.low; i<=v1.high; i++)
	{
		tmp += v1.elements[i] * v2.elements[i];
	}
	return tmp;
}



void Vector::print(ostream& out)
{
	int i;
	
	for (i=low; i<=high; i++) out << elements[i] << " ";
	out << endl;
}

void Vector::binwrite(ostream& out)
{
	out.write((char *)&low, sizeof(int));
	out.write((char *)&high, sizeof(int));
	out.write((char *)(elements+low), (high-low+1)*sizeof(double));
}

Vector Vector_binread(istream& in)
{
	int low, high;

	in.read((char *)&low, sizeof(int));
	in.read((char *)&high, sizeof(int));

	Vector tmp(low, high);

	in.read((char *)(tmp.elements+low), (high-low+1)*sizeof(double));

	return tmp;
}
	



void bomb_Vector(const char *msg)
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

void bomb_Vector_operation(const char *op)
{
        string msg("incompatible lengths in operation ");
	msg += op;
	bomb_Vector(msg.c_str());
}
			
	
		

	
	




Matrix::Matrix()
{
	rlow=0;
	rhigh=0;
	clow=0;
	chigh=0;
	rows = NULL;
}


Matrix::Matrix(int rl, int rh, int cl, int ch)
{
	rlow=0;
	rhigh=0;
	clow=0;
	chigh=0;
	rows = NULL;
	setsize(rl, rh, cl, ch);
}


Matrix::Matrix(int rl, int rh, int cl, int ch, double **array)
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


Matrix::Matrix(const Matrix &m)
{
	int i, j;

				// Allow construction of null matrix
	if (m.rows==NULL) {
	  rows = NULL;
	  rlow = rhigh = clow = chigh = 0;
	}
	else {

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
}	


Matrix::~Matrix(void)
{
	rows += rlow;
	// if (rows != NULL) delete [] rows;
	if (threading_on) pthread_mutex_lock(&mem_lock);
	delete [] rows;
	if (threading_on) pthread_mutex_unlock(&mem_lock);
}


/*
CMatrix Matrix::operator CMatrix(void)
{
	int i, j;
	CMatrix tmp(rlow, rhigh, clow, chigh);

	for (i=rlow; i<=rhigh; i++) 
	{	
		for (j=clow; j<=chigh; j++) 
			tmp.rows[i].elements[j] = rows[i].elements[j];
	}
}
*/

void Matrix::setsize(int rl, int rh, int cl, int ch)
{
	int i;

	/* do we need to resize at all? */

	if (rl==rlow && cl==clow && rh==rhigh && ch==chigh) return;



	/* is the new size positive? */

	if (rh<rl || ch<cl) bomb_Matrix("invalid size<=0");



	/* delete old storage if it exists */

	rows += rlow;
	// if (rows != NULL) delete [] rows;
	if (threading_on) pthread_mutex_lock(&mem_lock);
	delete [] rows;
	if (threading_on) pthread_mutex_unlock(&mem_lock);



	/* set the new size */

	rlow = rl;
	rhigh = rh;
	clow = cl;
	chigh = ch;



	/* allocate the array of rows */

	if (threading_on) pthread_mutex_lock(&mem_lock);
	rows = new Vector[rhigh+1-rlow];
	if (threading_on) pthread_mutex_unlock(&mem_lock);
	if (rows==NULL) bomb_Matrix("could not allocate rows");
	rows -= rlow;



	/* create the individual rows */

	for (i=rlow; i<=rhigh; i++) rows[i].setsize(clow, chigh);
}


void bomb_Matrix(const char *msg)
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

void bomb_Matrix_operation(const char *op)
{
	string str("incompatible sizes in operator ");
	str += op;
	bomb_Matrix(str.c_str());
}






Vector Matrix::fastcol(int j)
{
/*
	double *array;
	int i;

	if (threading_on) pthread_mutex_lock(&mem_lock);
	array = new double[rhigh-rlow+1];
	if (threading_on) pthread_mutex_unlock(&mem_lock);
	if (array == NULL) bomb_Matrix("could not allocate in fastcol");
	array -= rlow;

	for (i=rlow; i<=rhigh; i++) array[i] = rows[i].elements[j];
*/

	Vector tmp(rlow, rhigh);
	for (int i=rlow; i<=rhigh; i++) tmp[i] = rows[i].elements[j];

	return tmp;
}

Vector Matrix::col(int j)
{
	if (j<clow || j>chigh) bomb_Matrix("column subscript out of range");

	Vector tmp = fastcol(j);

	return tmp;
}



Vector &Matrix::row(int i)
{
	if (i<rlow || i>rhigh) bomb_Matrix("row subscript out of range");

	return rows[i];
}

Vector &Matrix::fastrow(int i)
{
	return rows[i];
}



void Matrix::fastsetrow(int i, const Vector &v)
{
	rows[i] = v;
}

void Matrix::setrow(int i, Vector &v)
{
	if (v.low != clow || v.high != chigh) 
		bomb_Matrix("row-vector size mismatch");
	if (i<rlow || i>rhigh) bomb_Matrix("row subscript out of range");
	rows[i] = v;
}

void Matrix::fastsetcol(int j, Vector &v)
{
	int i;

	for (i=rlow; i<=rhigh; i++)
	{
		rows[i].elements[j] = v.elements[i];
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
	int i;

	for (i=rlow; i<=rhigh; i++) rows[i].zero();
}


double **Matrix::array(int rl, int rh, int cl, int ch)
{
	int i;
	double **ptr;

	if (threading_on) pthread_mutex_lock(&mem_lock);
	ptr = new double*[rh-rl+1];
	if (threading_on) pthread_mutex_unlock(&mem_lock);
	if (ptr == NULL) bomb_Matrix("could not allocate array");
	ptr -= rl;

	for (i=rl; i<=rh; i++) ptr[i] = rows[i].array(cl, ch);

	return ptr;
}

double **Matrix::array_copy(int rl, int rh, int cl, int ch)
{
	int i;
	double **ptr;

	if (threading_on) pthread_mutex_lock(&mem_lock);
	ptr = new double*[rh-rl+1];
	if (threading_on) pthread_mutex_unlock(&mem_lock);
	if (ptr == NULL) bomb_Matrix("could not allocate array");
	ptr -= rl;

	for (i=rl; i<=rh; i++) ptr[i] = rows[i].array_copy(cl, ch);

	return ptr;
}

void destroy_Matrix_array(double **ptr, int rl, int rh, int cl, int ch)
{
	int i;

	for (i=rl; i<=rh; i++) destroy_Vector_array(ptr[i], cl, ch);

	ptr+=rl;
	if (threading_on) pthread_mutex_lock(&mem_lock);
	delete [] ptr;
	if (threading_on) pthread_mutex_unlock(&mem_lock);
}
	


Vector &Matrix::operator[](int i) const
{
	if (i<rlow || i>rhigh) bomb_Matrix("row subscript out of range");
	return rows[i];
}

Matrix &Matrix::operator=(const Matrix &m)
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
		bomb_Matrix_operation("=");

	for (i=rlow; i<=rhigh; i++) 
	{
		fastsetrow(i, m.rows[i]);
	}
	return *this;
}

Matrix Matrix::operator-(void)
{
	int i, j;

	for (i=rlow; i<=rhigh; i++)
	{
		for (j=clow; j<=chigh; j++) 
			rows[i].elements[j] = -rows[i].elements[j];
	}

	return *this;
}

Matrix operator+(const Matrix &m1, const Matrix &m2)
{
	int i, j;

	if (m1.rlow!=m2.rlow || m1.rhigh!=m2.rhigh || 
		m1.clow!=m2.clow || m1.chigh!=m2.chigh) 
		bomb_Matrix_operation("+");

	Matrix tmp(m1.rlow, m1.rhigh, m1.clow, m1.chigh);

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

Matrix operator-(const Matrix &m1, const Matrix &m2)
{
	int i, j;

	if (m1.rlow!=m2.rlow || m1.rhigh!=m2.rhigh || 
		m1.clow!=m2.clow || m1.chigh!=m2.chigh) 
		bomb_Matrix_operation("-");

	Matrix tmp(m1.rlow, m1.rhigh, m1.clow, m1.chigh);

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



Matrix operator*(const Matrix &m, double a)
{
	int i, j;

	Matrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);

	for (i=m.rlow; i<=m.rhigh; i++) 
	{
		for (j=m.clow; j<=m.chigh; j++)
		{
			tmp.rows[i].elements[j] = m.rows[i].elements[j] * a;
		}
	}
	return tmp;
}


Matrix operator*(double a, const Matrix &m)
{
	int i, j;

	Matrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);

	for (i=m.rlow; i<=m.rhigh; i++) 
	{
		for (j=m.clow; j<=m.chigh; j++)
		{
			tmp.rows[i].elements[j] = m.rows[i].elements[j] * a;
		}
	}
	return tmp;
}

	


Matrix operator/(const Matrix &m, double a)
{
	int i, j;

	Matrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);

	for (i=m.rlow; i<=m.rhigh; i++) 
	{
		for (j=m.clow; j<=m.chigh; j++)
		{
			tmp.rows[i].elements[j] = m.rows[i].elements[j] / a;
		}
	}
	return tmp;
}


Matrix operator+(const Matrix &m, double a)
{
	int i, j;

	Matrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);

	for (i=m.rlow; i<=m.rhigh; i++) 
	{
		for (j=m.clow; j<=m.chigh; j++)
		{
			tmp.rows[i][j] += a;
		}
	}
	return tmp;
}

Matrix operator-(const Matrix &m, double a)
{
	int i, j;

	Matrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);

	for (i=m.rlow; i<=m.rhigh; i++) 
	{
		for (j=m.clow; j<=m.chigh; j++)
		{
			tmp.rows[i][j] -= a;
		}
	}
	return tmp;
}


Matrix operator+(double a, const Matrix &m)
{
	int i, j;

	Matrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);

	for (i=m.rlow; i<=m.rhigh; i++) 
	{
		for (j=m.clow; j<=m.chigh; j++)
		{
			tmp.rows[i][j] += a;
		}
	}
	return tmp;
}

Matrix operator-(double a, const Matrix &m)
{
	int i, j;

	Matrix tmp(m.rlow, m.rhigh, m.clow, m.chigh);

	for (i=m.rlow; i<=m.rhigh; i++) 
	{
		for (j=m.clow; j<=m.chigh; j++)
		{
			tmp.rows[i][j] -= a;
		}
	}
	return tmp;
}


Matrix &Matrix::operator*=(double a)
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

Matrix &Matrix::operator/=(double a)
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

Matrix &Matrix::operator+=(double a)
{
	int i, j;

	for (i=rlow; i<=rhigh; i++) 
	{
		for (j=clow; j<=chigh; j++)
		{
			rows[i].elements[j] += a;
		}
	}
	return *this;
}

Matrix &Matrix::operator-=(double a)
{
	int i, j;

	for (i=rlow; i<=rhigh; i++) 
	{
		for (j=clow; j<=chigh; j++)
		{
			rows[i].elements[j] -= a;
		}
	}
	return *this;
}


Matrix &Matrix::operator+=(const Matrix &m)
{
	int i, j;

	if (m.rlow!=rlow || m.rhigh!=rhigh || m.clow!=clow || m.chigh!=chigh) 
		bomb_Matrix_operation("+=");


	for (i=rlow; i<=rhigh; i++) 
	{
		for (j=clow; j<=chigh; j++)
		{
			rows[i].elements[j] += m.rows[i].elements[j];
		}
	}
	return *this;
}

	
Matrix &Matrix::operator-=(const Matrix &m)
{
	int i, j;

	if (m.rlow!=rlow || m.rhigh!=rhigh || m.clow!=clow || m.chigh!=chigh) 
		bomb_Matrix_operation("+=");


	for (i=rlow; i<=rhigh; i++) 
	{
		for (j=clow; j<=chigh; j++)
		{
			rows[i].elements[j] -= m.rows[i].elements[j];
		}
	}
	return *this;
}


double operator^(const Matrix &m1, const Matrix &m2)
{
	int i, j;

	if (m1.getrlow()!=m2.getrlow() || m1.getrhigh()!=m2.getrhigh() || 
	    m1.getclow()!=m2.getclow() || m1.getchigh()!=m2.getchigh()) 
	  bomb_Matrix_operation("^");

	double tmp = 0.0;

	for (i=m1.rlow; i<=m1.rhigh; i++) 
	{
		for (j=m1.clow; j<=m1.chigh; j++)
		{
		  tmp +=  m1[i][j]*m2[i][j];
		}
	}

	return tmp;
}

Three_Vector operator*(const Three_Vector &v, const Matrix &m)
{
	int i, j;

	if (m.clow != 1 || m.chigh != 3 || m.rlow != 1 || m.rhigh != 3) 
		bomb_Matrix_operation("3v*M");

	Three_Vector tmp;

	for (j=1; j<=3; j++)
	{
		tmp.x[j-1] = 0.0;
		for (i=1; i<=3; i++)
		{
			tmp.x[j-1] += m.rows[i].elements[j]*v.x[j-1];
		}
	}
	return tmp;
}


Three_Vector operator*(const Matrix &m, const Three_Vector &v)
{
	int i, j;

	if (m.clow != 1 || m.chigh != 3 || m.rlow != 1 || m.rhigh != 3) 
		bomb_Matrix_operation("3v*M");

	Three_Vector tmp;

	for (i=1; i<=3; i++)
	{
		tmp.x[i-1] = 0.0;
		for (j=1; j<=3; j++)
		{
			tmp.x[i-1] += m.rows[i].elements[j]*v.x[j-1];
		}
	}
	return tmp;
}

Vector operator*(const Matrix &m, const Vector &v)
{
	int i, j;

	if (m.clow != v.low || m.chigh != v.high) bomb_Matrix_operation("M*v");

	Vector tmp(m.rlow, m.rhigh);

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

Vector operator*(const Vector &v, const Matrix &m)
{
	int i, j;

	if (m.rlow != v.low || m.rhigh != v.high) bomb_Matrix_operation("v*M");

	Vector tmp(m.clow, m.chigh);

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


Matrix operator*(const Matrix &m1, const Matrix &m2)
{
        int i, j, k;

	if (m1.clow != m2.rlow || m1.chigh != m2.rhigh)
		bomb_Matrix_operation("M*M");

	Matrix tmp(m1.rlow, m1.rhigh, m2.clow, m2.chigh);

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

Matrix Matrix::Transpose(void)
{
	int i;
	Matrix t(clow, chigh, rlow, rhigh);


	for (i=t.rlow; i<=t.rhigh; i++)
	{
		t.fastsetrow(i, fastcol(i));
	}

	return t;
}

double Matrix::Trace(void)
{
	double t;
	int i;

	for (i=rlow; i<=rhigh; i++) t += rows[i][i];
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
	int i;

	for (i=rlow; i<=rhigh; i++) rows[i].print(out);
}



void Matrix::binwrite(ostream& out)
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


Matrix Matrix_binread(istream& in)
{

	int rlow, rhigh, clow, chigh;
	int i;

	in.read((char *)&rlow, sizeof(int));
	in.read((char *)&rhigh, sizeof(int));
	in.read((char *)&clow, sizeof(int));
	in.read((char *)&chigh, sizeof(int));

	Matrix tmp(rlow, rhigh, clow, chigh);

	for (i=rlow; i<=rhigh; i++)
	{
		tmp.rows[i] = Vector_binread(in);
	}

	return tmp;
}


inline double fabs(Vector& v) {return sqrt(v*v);}

	
		










	




	




	

	




	




