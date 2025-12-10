
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <string>

#include "cpoly.H"

/*
	Default constructor; make a null vector.
*/

CPoly::CPoly(void) : Eigen::VectorXcd()
{
  order = 0;
}


CPoly::CPoly(int n): Eigen::VectorXcd(n+1)
{
  order = n;
  setZero();
}


CPoly::CPoly(int n, double * vec) : Eigen::VectorXcd(n+1)
{
  for (int i=0; i<=n; i++) (*this)[i] = vec[i];
  order = n;
}



CPoly::CPoly(const Eigen::VectorXcd& vec) : Eigen::VectorXcd(vec)
{
  order = vec.size()+1;
}



/*
	Conversion constructor
*/

CPoly::CPoly(const Poly &p) : Eigen::VectorXcd(static_cast<Eigen::VectorXcd>(p))
{
  order = p.getorder();
}

/*
	Copy constructor; create a new CPoly which is a copy of another.
*/

CPoly::CPoly(const CPoly &p) : Eigen::VectorXcd(static_cast<Eigen::VectorXcd>(p))
{
  order = p.order;
}




/*
	Destructor. Free elements if it exists.
	[Does nothing, yet]
*/

// CPoly::~CPoly() {}


/*
	Function to reduce order if there are leading zero coefficients
*/


void CPoly::reduce_order(void)
{
  while ((*this)[order].real() == 0.0 && (*this)[order].imag() == 0.0 && 
	 order>0) order--;
}


/*
	Assignment operator for CPoly; must be defined as a reference
so that it can be used on lhs. Compatibility checking is performed;
the destination vector is allocated if its elements are undefined.
*/


CPoly &CPoly::operator=(const CPoly &v)
{
  if (v.size()<1) {
    bomb_CPoly_operation("=");
  }

  resize(v.order+1);
  for (int i=0; i<=v.order; i++) (*this)[i] = v[i];
  order = v.order;

  return *this;
}


CPoly CPoly::operator-(void)
{
  int i;
  for (i=0; i<=order; i++) (*this)[i] = -(*this)[i];
  return *this;
}


CPoly &CPoly::operator+=(const CPoly &p2)
{
  int n2 = p2.order;
  Eigen::VectorXcd tmp;

  if (order <= n2) {
    tmp = static_cast<Eigen::VectorXcd>(*this);
    resize(n2+1);
    setZero();
    for (int i=0; i<=order; i++) (*this)[i] = tmp[i];
    order = n2;
  }
  for (int i=0; i<=n2; i++) (*this)[i] += p2[i];

  reduce_order();
  return *this;
}
	
CPoly &CPoly::operator-=(const CPoly &p2)
{
  int i;
  int n2 = p2.order;
  Eigen::VectorXcd tmp;

  if (order <= n2) {
    tmp = static_cast<Eigen::VectorXcd>(*this);
    resize(n2+1);
    setZero();
    for (int i=0; i<=order; i++) (*this)[i] = tmp[i];
    order = n2;
  }
  for (int i=0; i<=n2; i++) (*this)[i] -= p2[i];

  reduce_order();
  return *this;
}


CPoly operator+(const CPoly &p1, const CPoly &p2)
{
  int i;
  CPoly tmp;
  int n1 = p1.order;
  int n2 = p2.order;

  if (n1 <= n2) {
    tmp = p2;
    for (i=0; i<=n1; i++) tmp[i] += p1[i];
  }
  else {
    tmp = p1;
    for (i=0; i<=n2; i++) tmp[i] += p2[i];
  }

  tmp.reduce_order();
  return tmp;
}
	

CPoly operator-(const CPoly &p1, const CPoly &p2)
{
  int i;
  CPoly tmp;
  int n1 = p1.order;
  int n2 = p2.order;

  if (n1 <= n2) {
    tmp = CPoly(n2);
    tmp.setZero();
    for (i=0; i<=n1; i++) tmp[i] = p1[i] - p2[i];
    for (i=n1+1; i<=n2; i++) tmp[i] = - p2[i];
  }
  else {
    tmp = p1;
    for (i=0; i<=n2; i++) tmp[i] -= p2[i];
  }

  tmp.reduce_order();
  return tmp;
}
	

				// Cauchy product
CPoly operator&(const CPoly &p1, const CPoly &p2)
{
  int i, j;
  int n1 = p1.order;
  int n2 = p2.order;
  int neworder = n1+n2;
  CPoly tmp(neworder);

  for (i=0; i<=n1; i++) {
    for (j=0; j<=n2; j++) tmp[i+j] += p1[i]*p2[j];
  }

  tmp.reduce_order();
  return tmp;
}
	

CPoly &CPoly::operator&=(const CPoly &p2)
{
  int n2 = p2.order;
  int neworder = order + n2;
  CPoly tmp = *this;
  resize(neworder+1);
  setZero();

  for (int i=0; i<=order; i++) {
    for (int j=0; j<=n2; j++) (*this)[i+j] += tmp[i]*p2[j];
  }

  order = neworder;
  reduce_order();
  
  return *this;
}
	

				// Euclidian division for power series
CPoly operator%(const CPoly &p1, const CPoly &p2)
{

  int k,j;
  int n1 = p1.order;
  int n2 = p2.order;

  CPoly quotient(n1);
  CPoly remainder = p1;

  for (k=0; k<=n1; k++) {
    quotient[k] = remainder[k]/p2[0];
    for (j=k+1; j<=n1 && j-k<=n2; j++)
      remainder[j] -= quotient[k]*p2[j-k];
  }
  
  quotient.reduce_order();
  return quotient;
}
	
CPoly &CPoly::operator%=(const CPoly &p2)
{

  int k,j;
  int n1 = order;
  int n2 = p2.order;

  for (k=0; k<=n1; k++) {
    (*this)[k] /= p2[0];
    for (j=k+1; j<=n1 && j-k<=n2; j++)
      (*this)[j] -= (*this)[k]*p2[j-k];
  }
  
  reduce_order();

  return *this;
}

/*
CPoly &CPoly::operator%=(const CPoly &p2)
{

  int k,j;
  int n1 = order;
  int n2 = p2.getorder();

  CPoly quotient = CPoly(n1);
  CPoly remainder = (*this);

  for (k=0; k<=n1; k++) {
    quotient[k] = remainder[k]/p2[0];
    for (j=k+1; j<=n1 && j-k<=n2; j++)
      remainder[j] -= quotient[k]*p2[j-k];
  }
  
  quotient.reduce_order();
  *this = quotient;

  return *this;
}
*/

std::complex<double> CPoly::eval(std::complex<double> z)
{
  int j;
  std::complex<double> p = (*this)[j=order];
  while (j>0) p = p*z + (*this)[--j];

  return p;
}
	
std::complex<double> CPoly::deriv(std::complex<double> z)
{
  int j;
  std::complex<double> p = (*this)[j=order];
  std::complex<double> dp = 0.0;

  while (j>0) {
    dp = dp*z + p;
    p = p*z + (*this)[--j];
  }
  return dp;
}
	
void CPoly::print(ostream& out)
{
  int i;
	
  out << "[" << order << "]: ";
  for (i=0; i<=order; i++)
    out << "(" << (*this)[i].real() << " + " << (*this)[i].imag() << ") ";
  out << endl;

}

void bomb_CPoly(const char *msg)
{
  cerr << "CPoly ERROR: " << msg << '\n';
  exit(0);
}

void bomb_CPoly_operation(const char *op)
{
  string msg("incompatible arguments in operation ");
  msg += op;
  bomb_CPoly(msg.c_str());
}
