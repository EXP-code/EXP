#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <string>

#include <Vector.h>

#include <poly.h>

using namespace std;

/*
	Default constructor; make a null vector.
*/

Poly::Poly(void) : Vector()
{
  order = 0;
}


Poly::Poly(int n): Vector(0, n)
{
  order = n;
  (*this).zero();
}


Poly::Poly(int n, double * vec) : Vector(0, n, vec)
{
  order = n;
  reduce_order();
}



Poly::Poly(const Vector& vec) : Vector((Vector&)vec)
{
  if (vec.getlow() != 0) bomb_Poly("Error constructing Poly with Vector");
  order = vec.gethigh();
  reduce_order();
}



/*
	Copy constructor; create a new Poly which is a copy of another.
*/

Poly::Poly(const Poly &p) : Vector((Vector &)p)
{
  order = p.order;
}




/*
	Destructor. Free elements if it exists.
	[Does nothing, yet]
*/

Poly::~Poly() {}


/*
	Function to reduce order if there are leading zero coefficients
*/


void Poly::reduce_order(void)
{
  while ((*this)[order] == 0.0 && order>0) order--;
}


/*
	Assignment operator for Poly; must be defined as a reference
so that it can be used on lhs. Compatibility checking is performed;
the destination vector is allocated if its elements are undefined.
*/


Poly &Poly::operator=(Poly &v)
{
  int i;
	
  if (v.getlow()!=0 || v.gethigh()<1) {
    bomb_Poly_operation("=");
  }

  setsize(0,v.order);
  for (i=0; i<=v.order; i++) (*this)[i] = v[i];
  order = v.order;

  (*this).reduce_order();

  return *this;
}


Poly Poly::operator-(void)
{
  (Vector)*this *= -1.0;
  return *this;
}


Poly &Poly::operator+=(Poly &p2)
{
  int i;
  int n2 = p2.getorder();

  if (order <= n2) {
    Vector p1 = *this;
    (*this).setsize(0,n2);
    (*this).zero();
    for (i=0; i<=order; i++) (*this)[i] = p1[i];
  }
  for (i=0; i<=n2; i++) (*this)[i] += p2[i];

  (*this).reduce_order();
  return *this;
}
	
Poly &Poly::operator-=(Poly &p2)
{
  int i;
  int n2 = p2.getorder();

  if (order <= n2) {
    Vector p1 = *this;
    (*this).setsize(0,n2);
    (*this).zero();
    for (i=0; i<=order; i++) (*this)[i] = p1[i];
  }
  for (i=0; i<=n2; i++) (*this)[i] -= p2[i];

  (*this).reduce_order();
  return *this;
}


Poly operator+(Poly &p1, Poly &p2)
{
  int i;
  Vector tmp;
  int n1 = p1.getorder();
  int n2 = p2.getorder();

  if (n1 <= n2) {
    tmp = (Vector)p2;
    for (i=0; i<=n1; i++) tmp[i] += p1[i];
  }
  else {
    tmp = (Vector)p1;
    for (i=0; i<=n2; i++) tmp[i] += p2[i];
  }

  Poly tmp2 = Poly(tmp);
  tmp2.reduce_order();
  return tmp2;
}
	

Poly operator-(Poly &p1, Poly &p2)
{
  int i;
  Vector tmp;
  int n1 = p1.getorder();
  int n2 = p2.getorder();

  if (n1 <= n2) {
    tmp = Vector(0, n2);
    tmp.zero();
    for (i=0; i<=n1; i++) tmp[i] = p1[i] - p2[i];
    for (i=n1+1; i<=n2; i++) tmp[i] = - p2[i];
  }
  else {
    tmp = (Vector)p1;
    for (i=0; i<=n2; i++) tmp[i] -= p2[i];
  }

  Poly tmp2 = Poly(tmp);
  tmp2.reduce_order();
  return tmp2;
}
	

				// Cauchy product
Poly operator&(Poly &p1, Poly &p2)
{
  int i, j;
  int n1 = p1.getorder();
  int n2 = p2.getorder();
  int neworder = n1+n2;
  Poly tmp = Poly(neworder);

  for (i=0; i<=n1; i++) {
    for (j=0; j<=n2; j++) tmp[i+j] += p1[i]*p2[j];
  }

  tmp.reduce_order();
  return tmp;
}
	

Poly &Poly::operator&=(Poly &p2)
{
  int i, j;
  int n2 = p2.getorder();
  int neworder = order + n2;
  Poly tmp = Poly(neworder);

  for (i=0; i<=order; i++) {
    for (j=0; j<=n2; j++) tmp[i+j] += (*this)[i]*p2[j];
  }

  tmp.reduce_order();
  *this = tmp;
  
  return *this;
}
	

				// Euclidian division for power series
Poly operator%(Poly &p1, Poly &p2)
{

  int k,j;
  int n1 = p1.getorder();
  int n2 = p2.getorder();

  Poly quotient = Poly(n1);
  Poly remainder = p1;

  for (k=0; k<=n1; k++) {
    quotient[k] = remainder[k]/p2[0];
    for (j=k+1; j<=n1 && j-k<=n2; j++)
      remainder[j] -= quotient[k]*p2[j-k];
  }
  
  quotient.reduce_order();
  return quotient;
}
	
Poly &Poly::operator%=(Poly &p2)
{

  int k,j;
  int n1 = order;
  int n2 = p2.getorder();

  Poly quotient = Poly(n1);
  Poly remainder = (*this);

  for (k=0; k<=n1; k++) {
    quotient[k] = remainder[k]/p2[0];
    for (j=k+1; j<=n1 && j-k<=n2; j++)
      remainder[j] -= quotient[k]*p2[j-k];
  }
  
  quotient.reduce_order();
  *this = quotient;

  return *this;
}
	

double Poly::eval(double z)
{
  int j;
  double p = (*this)[j=order];
  while (j>0) p = p*z + (*this)[--j];

  return p;
}

 
double Poly::deriv(double z)
{
  int j;
  double p = (*this)[j=order];
  double dp = 0.0;

  while (j>0) {
    dp = dp*z + p;
    p = p*z + (*this)[--j];
  }
  return dp;
}
	

void Poly::print(ostream& out)
{
  int i;
	
  out << "[" << order << "]: ";
  for (i=0; i<=order; i++) out << (*this)[i] << " ";
  out << '\n';

}

void bomb_Poly(const char *msg)
{
  cerr << "POLY ERROR: " << msg << '\n';
  exit(0);
}

void bomb_Poly_operation(const char *op)
{
  string msg("incompatible arguments in operation ");
  msg += op;
  bomb_Poly(msg.c_str());
}
