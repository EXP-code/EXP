
/* Interpolation routines operating on vectors */

#ifndef _interp_h
#define _interp_h 1
#include <Vector.h>

void Spline(const Vector &x, const Vector &y, double yp1, double ypn, Vector &y2);

void Splint1(const Vector &xa, const Vector &ya, const Vector &y2a, double x, 
	     double &y, int even=0);

void Splint2(const Vector &xa, const Vector &ya, const Vector &y2a, double x, 
	     double &y, double &dy, int even=0);

void Splint3(const Vector &xa, const Vector &ya, const Vector &y2a, double x, 
	     double &y, double &dy, double &dyy, int even=0);

double Splsum(const Vector& x, const Vector& y);
void Splsum(const Vector& x, const Vector& y, Vector& z);

double odd2(double x, const Vector &xtab, const Vector &ftab, int even=0);

double drv2(double x, const Vector &xtab, const Vector &ftab, int even=0);

int Vlocate(double x, const Vector xtab);

#endif
