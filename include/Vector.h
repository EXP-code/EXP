// This is really -*- C++ -*-

#ifndef _Vector_h

#define _Vector_h 1

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

class CVector;
class Three_Vector;
class Matrix;
class CMatrix;

void bomb_Vector_operation(const string&);
void bomb_Vector(const string&);
void bomb_Matrix_operation(const string&);
void bomb_Matrix(const string&);

class Vector
{
  friend class CMatrix;
  friend class CVector;
  friend class Matrix;

 protected:
  int low, high, size;

  // memory management uses STL vector
  double *pelement;
  vector<double> elements;

 public:
  static int mpi_id;

  // constructors

  Vector(void);
  Vector(int, int);
  Vector(int, int, double *);
  Vector(const Vector &);


  // conversion

  operator Three_Vector(void);


  // assignment, access, resizing

  void setsize(int, int);

  double &operator[](int) const;	// safe access

  Vector &operator=(const Vector &);
  void load(double *);
  void zero(void);
  int getlength(void) const {return high-low+1;}
  int getlow(void) const {return low;}
  int gethigh(void) const {return high;}


  // make an array

  double *array(int l, int h);
  double *array_copy(int l, int h);
  friend void destroy_Vector_array(double *, int, int);




  // unary plus and minus

  Vector operator+() {return *this;}
  Vector operator-();



  // Vector addition and subtraction

  Vector &operator+=(const Vector &);
  Vector &operator-=(const Vector &);
  friend Vector operator+(const Vector &, const Vector &);
  friend Vector operator-(const Vector &, const Vector &);



  // dot product

  friend double operator*(const Vector &, const Vector &);
  friend Vector operator^(const Vector &, const Vector &);

  // combine product
  friend Vector operator&(const Vector&, const Vector &);

  // operations with scalars

  Vector &operator*=(double);
  Vector &operator/=(double);
  Vector &operator+=(double);
  Vector &operator-=(double);
  friend Vector operator*(const Vector &, double);
  friend Vector operator*(double, const Vector &);
  friend Vector operator/(const Vector &, double);


  // magnitude of vector

  friend double fabs(Vector& v);

  // operations with matrices

  friend Vector operator*(const Matrix &, const Vector &);
  friend Vector operator*(const Vector &, const Matrix &);
  friend Matrix operator+(const Matrix &, const Matrix &);
  friend Matrix operator-(const Matrix &, const Matrix &);
  friend Matrix operator*(const Matrix &, const Matrix &);
  friend Matrix operator*(const Matrix &, double);
  friend Matrix operator/(const Matrix &, double);
  friend Matrix operator*(double, const Matrix &);
  friend CVector operator*(const CMatrix &a, const Vector &v);
  friend Three_Vector operator*(const Matrix &, const Three_Vector &);
  friend Three_Vector operator*(const Three_Vector &, const Matrix &);


  // IO

  void print(ostream&);
  void binwrite(ostream&);
  friend Vector Vector_binread(istream&);



  void Sort(void);


};



class Matrix
{
  friend class Vector;
  friend class Three_Vector;
  friend class CMatrix;
  friend class CVector;

 protected:
  int rlow,  rhigh;	// low and high row indices
  int clow,  chigh;	// low and high column indices
  int rsize, csize;	// dimension sizes
    
  // memory management uses STL vector
  Vector *prow;
  vector<Vector> rows;

  // fast (but unsafe) access functions

  Vector fastcol(int);
  Vector &fastrow(int);
  void fastsetrow(int, const Vector &);
  void fastsetcol(int, Vector &);

 public:

  static int mpi_id;

  // constructors and the destructor

  Matrix(void);
  Matrix(int, int, int, int);
  Matrix(int, int, int, int, double **);
  Matrix(const Matrix &);


  // assignment and access

  Matrix &operator=(const Matrix &);
  Vector &operator[](int) const;

  int getnrows(void) const {return rhigh-rlow+1;}
  int getncols(void) const {return chigh-clow+1;}
  int getrlow(void) const {return rlow;}
  int getclow(void) const {return clow;}
  int getrhigh(void) const {return rhigh;}
  int getchigh(void) const {return chigh;}



  // sizing

  void setsize(int, int, int, int);
  void zero(void);
  void load(double **);



  // row and column access

  Vector col(int);
  Vector &row(int);
  void setrow(int, Vector &);
  void setcol(int, Vector &);



  // make a double pointer

  double **array(int, int, int, int);
  double **array_copy(int, int, int, int);
  friend void destroy_Matrix_array(double **, int, int, int, int);



  // operations with matrices

  friend Matrix operator+(const Matrix &, const Matrix &);
  friend Matrix operator-(const Matrix &, const Matrix &);
  friend Matrix operator*(const Matrix &, const Matrix &);
  friend double operator^(const Matrix &, const Matrix &);
  Matrix &operator+=(const Matrix &);
  Matrix &operator-=(const Matrix &);

  // unary plus and minus

  Matrix operator+(void) {return *this;}
  Matrix operator-(void);



  // operations with vectors

  friend Vector operator*(const Matrix &, const Vector &);
  friend Vector operator*(const Vector &, const Matrix &);
  friend Three_Vector operator*(const Matrix &, const Three_Vector &);
  friend Three_Vector operator*(const Three_Vector &, const Matrix &);



  // operations with scalars


  friend Matrix operator*(const Matrix &, double);
  friend Matrix operator*(double, const Matrix &);
  friend Matrix operator/(const Matrix &, double);
  friend Matrix operator+(const Matrix &, double);
  friend Matrix operator+(double, const Matrix &);
  friend Matrix operator-(const Matrix &, double);
  friend Matrix operator-(double, const Matrix &);
  Matrix &operator*=(double);
  Matrix &operator/=(double);
  Matrix &operator+=(double);
  Matrix &operator-=(double);



  // IO

  void print(ostream&);
  void binwrite(ostream&);
  friend Matrix Matrix_binread(istream&);

  Matrix Transpose(void);	
  double Trace(void);
  Vector Symmetric_Eigenvalues(Matrix &);
  Vector Symmetric_Eigenvalues_GHQL(Matrix &);
  Matrix Inverse(void);

};





class Three_Vector;
class Vector;

class Three_Vector
{
  friend class Vector;
  friend class Matrix;
 protected:
  double x[3];
  static int nlive;
 public:
  Three_Vector(void)
    {nlive++;}
  Three_Vector(double xx, double yy, double zz)
    {x[0] = xx; x[1] = yy; x[2] = zz;}
  Three_Vector(const Three_Vector &v)
    {nlive++; x[0]=v.x[0]; x[1]=v.x[1]; x[2]=v.x[2];}
  ~Three_Vector(void)
    {nlive--;}


  // conversion

  operator Vector(void);




  // access

  double &operator[](int i) {return x[i-1];}



  // assignment operators

  Three_Vector &operator=(const Three_Vector &v)	
    {x[0]=v.x[0]; x[1]=v.x[1]; x[2]=v.x[2]; return *this;}



  // reflexive arithmetic operators

  Three_Vector &operator+=(const Three_Vector &v)
    {x[0]+=v.x[0]; x[1]+=v.x[1]; x[2]+=v.x[2]; return *this;}
  Three_Vector &operator-=(const Three_Vector &v)
    {x[0]-=v.x[0]; x[1]-=v.x[1]; x[2]-=v.x[2]; return *this;}
  Three_Vector &operator*=(double a)
    {x[0]*=a; x[1]*=a; x[2]*=a; return *this;}
  Three_Vector &operator/=(double a)
    {x[0]/=a; x[1]/=a; x[2]/=a; return *this;}

  // unary plus and minus

  Three_Vector operator+(void) {return *this;}
  Three_Vector operator-(void)
    {
      x[0] = -x[0];
      x[1] = -x[1];
      x[2] = -x[2];
      return *this;
    }




  // binary operations with Three_Vectors

  friend Three_Vector operator+(const Three_Vector &, const Three_Vector &);
  friend Three_Vector operator-(const Three_Vector &, const Three_Vector &);
  friend Three_Vector operator^(const Three_Vector &, const Three_Vector &);
  friend Three_Vector Cross(const Three_Vector &, const Three_Vector &);


  // binary operations with doubles

  friend Three_Vector operator*(const Three_Vector &, double);
  friend Three_Vector operator*(double, const Three_Vector &);
  friend Three_Vector operator/(const Three_Vector &, double);



  // binary operations with Vectors

  friend Three_Vector operator+(const Vector &, const Three_Vector &);
  friend Three_Vector operator-(const Vector &, const Three_Vector &);
  friend Three_Vector operator^(const Vector &, const Three_Vector &);
  friend Three_Vector operator+(const Three_Vector &, const Vector &);
  friend Three_Vector operator-(const Three_Vector &, const Vector &);
  friend Three_Vector operator^(const Three_Vector &, const Vector &);

  friend Three_Vector operator*(const Matrix &, const Three_Vector &);
  friend Three_Vector operator*(const Three_Vector &, const Matrix &);


  // scalar product

  friend double operator*(const Three_Vector &v1, const Three_Vector &v2);

  void zero(void)		
    {x[0]=0.0; x[1]=0.0; x[2]=0.0;}


  void binwrite(ostream&);
  void binread(istream&);
  void print(ostream&);
  static void count(void)
    {cerr << "number of live three_vectors " << nlive << endl;}

};


Matrix Transpose(Matrix &m);
double Trace(Matrix &m);
Vector Symmetric_Eigenvalues(Matrix &m, Matrix &ev);
Vector Symmetric_Eigenvalues_GHQL(Matrix &m, Matrix &ev);
void jacobi(double **, int, double *, double **, int *);
void eigsrt(double *, double **, int);
int SVD(Matrix &A, Matrix &U, Matrix &V, Vector &Z);

void VectorSynchronize(Vector &c,  int id);
void MatrixSynchronize(Matrix &c,  int id);

#endif
