#ifndef _poly_h

#define _poly_h 1

#ifdef __GNUG__
#pragma interface
#endif

#include <iostream.h>
#include <fstream.h>

class Poly : public Vector
{

private:

  int order;
  void reduce_order(void);

public:

  // constructors
		
  Poly(void);
  Poly(int);
  Poly(int, double *);
  Poly(const Vector &);
  Poly(const Poly &);
		

  // the destructor
		
  ~Poly(void);

  // access to privates

  Poly &operator=(Poly &);
  int getorder(void) const {return order;}

  // unary plus and minus

  Poly operator+() {return *this;}
  Poly operator-();

  // Vector addition and subtraction
		
  Poly &operator+=(Poly &);
  friend Poly operator+(Poly &, Poly &);

  Poly &operator-=(Poly &);
  friend Poly operator-(Poly &, Poly &);

  Poly &operator&=(Poly &);	/* Cauchy product */
  friend Poly operator&(Poly &, Poly &);

  Poly &operator%=(Poly &);	/* Synthetic division for power series */
  friend Poly operator%(Poly &, Poly &);

  // Evaluate polynomial
    
  double eval(double z);
  double deriv(double z);

  /* IO */

  void print(ostream &);
};

void bomb_Poly_operation(const char *);
void bomb_Poly(const char *);


#endif
