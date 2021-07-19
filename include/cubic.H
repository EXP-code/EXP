#ifndef cubic_H
#define cubic_H

#include <Eigen/Eigen>

class Cubic_Table
{
protected:
  Eigen::VectorXd xgrid, dx;
  Eigen::MatrixXd a;

public:
  Cubic_Table(void);
  ~Cubic_Table(void);
  
  Cubic_Table(Cubic_Table &c)
  {
    dx = c.dx;
    xgrid = c.xgrid;
    a = c.a;
  }
  
  Cubic_Table &operator=(Cubic_Table &c)
  {
    dx=c.dx; 
    xgrid=c.xgrid; 
    a=c.a; 
    
    return *this;
  }
  
  Cubic_Table &operator+=(Cubic_Table &c)
  {
    a += c.a;
    return *this;
  }
  
  Cubic_Table &operator-=(Cubic_Table &c)
  {
    a -= c.a;
    return *this;
  }
  
  Cubic_Table &operator*=(double x)
  {
    a *= x;
    return *this;
  }
  
  Cubic_Table &operator/=(double x)
  {
    a /= x;
    return *this;
  }
  
  
  
  friend Cubic_Table operator+(Cubic_Table &c1, Cubic_Table &c2); 
  friend Cubic_Table operator-(Cubic_Table &c1, Cubic_Table &c2);
  friend Cubic_Table operator*(Cubic_Table &c, double x);
  friend Cubic_Table operator/(Cubic_Table &c, double x);
  friend Cubic_Table operator*(double x, Cubic_Table &c);
  
  
  void build(Eigen::VectorXd &, Eigen::VectorXd &, Eigen::VectorXd &);
  double value(double, double *);
  
};

#endif
