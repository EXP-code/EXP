#ifndef _cpoly_H
#define _cpoly_H

#include "poly.H"

class CPoly : public Eigen::VectorXcd
{

private:

  int order;
  void reduce_order(void);
  
public:
  
  // constructors
		
  CPoly(void);
  CPoly(const int);
  CPoly(int, double *);
  CPoly(const Eigen::VectorXcd &);
  CPoly(const CPoly &);
  CPoly(const Poly &);
  
  // the destructor
		
//  ~CPoly(void);

  // reduce order by n (BE CAREFUL with this one)
    
  void Pop(int i) {order -= i;}

  // access to privates

  CPoly &operator=(const CPoly &);
  const int getorder(void) {return order;}
  
  // unary plus and minus

  CPoly operator+() {return *this;}
  CPoly operator-();

  // Vector addition and subtraction
		
  CPoly &operator+=(const CPoly &);
  friend CPoly operator+(const CPoly &, const CPoly &);

  CPoly &operator-=(const CPoly &);
  friend CPoly operator-(const CPoly &, const CPoly &);

  CPoly &operator&=(const CPoly &);	/* Cauchy product */
  friend CPoly operator&(const CPoly &, const CPoly &);

  CPoly &operator%=(const CPoly &);	/* Synthetic division for power series */
  friend CPoly operator%(const CPoly &, const CPoly &);

  // Evaluate polynomial
    
  std::complex<double> eval(std::complex<double> z);
  std::complex<double> deriv(std::complex<double> z);

  /* IO */

  void print(ostream &);
};

void bomb_CPoly_operation(const char *);
void bomb_CPoly(const char *);

#endif
