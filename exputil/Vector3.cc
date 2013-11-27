using namespace std;

#include <cstdlib>
#include <iostream>

#include <Vector.h>

int Three_Vector::nlive=0;


Three_Vector::operator Vector(void)
{
  Vector tmp(1, 3);
  for (int i=1; i<=3; i++) tmp[i] = x[i-1];
  
  return tmp;
}

Three_Vector operator+(const Three_Vector &v1, const Three_Vector &v2)
{
  Three_Vector tmp;

  tmp.x[0] = v1.x[0] + v2.x[0];
  tmp.x[1] = v1.x[1] + v2.x[1];
  tmp.x[2] = v1.x[2] + v2.x[2];
  
  return tmp;
}


Three_Vector operator-(const Three_Vector &v1, const Three_Vector &v2)
{
  Three_Vector tmp;

  tmp.x[0] = v1.x[0] - v2.x[0];
  tmp.x[1] = v1.x[1] - v2.x[1];
  tmp.x[2] = v1.x[2] - v2.x[2];
  
  return tmp;
}


Three_Vector operator^(const Three_Vector &v1, const Three_Vector &v2)
{
  Three_Vector tmp;

  tmp.x[0] = v1.x[1]*v2.x[2] - v1.x[2]*v2.x[1];
  tmp.x[1] = v1.x[2]*v2.x[0] - v1.x[0]*v2.x[2];
  tmp.x[2] = v1.x[0]*v2.x[1] - v1.x[1]*v2.x[0];
  
  return tmp;
}

Three_Vector Cross(const Three_Vector &v1, const Three_Vector &v2)
{
  Three_Vector tmp;
  
  tmp.x[0] = v1.x[1]*v2.x[2] - v1.x[2]*v2.x[1];
  tmp.x[1] = v1.x[2]*v2.x[0] - v1.x[0]*v2.x[2];
  tmp.x[2] = v1.x[0]*v2.x[1] - v1.x[1]*v2.x[0];
  
  return tmp;
}


Three_Vector operator*(const Three_Vector &v, double a)
{
  Three_Vector tmp;
  
  tmp.x[0] = v.x[0]*a;
  tmp.x[1] = v.x[1]*a;
  tmp.x[2] = v.x[2]*a;

  return tmp;
}

Three_Vector operator*(double a, const Three_Vector &v)
{
  Three_Vector tmp;

  tmp.x[0] = v.x[0]*a;
  tmp.x[1] = v.x[1]*a;
  tmp.x[2] = v.x[2]*a;
  
  return tmp;
}


Three_Vector operator/(const Three_Vector &v, double a)
{
  Three_Vector tmp;

  tmp.x[0] = v.x[0]/a;
  tmp.x[1] = v.x[1]/a;
  tmp.x[2] = v.x[2]/a;
  
  return tmp;
}




void Three_Vector::print(ostream& out)
{
  out << "(" << x[0] << ", " << x[1] << ", " << x[2] << "\n";
}

void Three_Vector::binwrite(ostream& out)
{
  out.write((char *)x, 3*sizeof(double));
}

void Three_Vector::binread(istream& in)
{
        in.read((char *)x, 3*sizeof(double));
}

double operator*(const Three_Vector &v1, const Three_Vector &v2)
		{return v1.x[0]*v2.x[0] + v1.x[1]*v2.x[1] + v1.x[2]*v2.x[2];}


