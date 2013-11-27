#ifndef _Func1d_H
#define _Func1d_H

class Func1d
{
 public:
  
  virtual ~Func1d() {};

  virtual double CostFunction(double *) = 0;

};

#endif
