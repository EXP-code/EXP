#ifndef _CVector_h
#ifdef __GNUG__
#pragma interface
#endif
#define _CVector_h 1
#include <iostream.h>
#include <fstream.h>

#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

#include <kevin_complex.h>

class Vector;
class Matrix;
class CVector;
class CMatrix;
class Complex;

class CVector
{
	friend Matrix;
	friend Vector;
	friend CMatrix;
//	private:
        protected:
		int low, high;
		Complex *elements;
	public:
		/* constructors */
		
		CVector(void);
		CVector(int, int);
		CVector(int, int, double *);
		CVector(int, int, Complex *);
		CVector(const CVector &);
		CVector(const Vector &);
		

		/* the destructor */
		
		~CVector(void);




		/* assignment, access, resizing */
		
	        void setsize(int, int);
		
		Complex &operator[](int) const;	/* safe access */
		CVector &operator=(const CVector &);

		void zero(void);
		int getlength(void) const {return high-low+1;}
		int getlow(void) const {return low;}
		int gethigh(void) const {return high;}


		/* make an array */

		Complex *array(int l, int h);
		
		


		/* unary plus and minus */

		CVector operator+(void) {return *this;}
		CVector operator-(void);

		

		/* Vector addition and subtraction */
		
		CVector &operator+=(const CVector &);
		CVector &operator-=(const CVector &);
		CVector &operator+=(const Vector &v)
		        {*this += v; return *this;}
//			{CVector c(v); *this += c; return *this;}
		CVector &operator-=(const Vector &v)
			{*this -= v; return *this;}
//			{CVector c(v); *this -= c; return *this;}
		friend CVector operator+(const CVector &, const CVector &);
		friend CVector operator-(const CVector &, const CVector &);


		

		/* dot product */
		
		friend Complex operator*(const CVector &, const CVector &);

	        /* combine product */

                friend CVector operator&(const CVector&, const Vector &);
                friend CVector operator&(const Vector&, const CVector &);


		/* operations with scalars */
		
		CVector &operator*=(const Complex &);
		CVector &operator/=(const Complex &);
		friend CVector operator*(const CVector &, const Complex &);
		friend CVector operator*(const Complex &, const CVector &);
		friend CVector operator/(const CVector &, const Complex &);

		CVector &operator*=(double);
		CVector &operator/=(double);
		friend CVector operator*(const CVector &, double);
		friend CVector operator*(double, const CVector &);
		friend CVector operator/(const CVector &, double);


		/* operations with matrices */
		
		friend CVector operator*(const CMatrix &, const CVector &);
		friend CVector operator*(const CVector &, const CMatrix &);
	        friend CVector operator*(const CMatrix &, const Vector &);

		/* matrix-matrix operations */

		friend CMatrix operator+(const CMatrix &, const CMatrix &);
		friend CMatrix operator-(const CMatrix &, const CMatrix &);
		friend CMatrix operator*(const CMatrix &, const CMatrix &);

		/* scalar-matrix operations */

		friend CMatrix operator*(const CMatrix &, const Complex &);
		friend CMatrix operator/(const CMatrix &, const Complex &);
		friend CMatrix operator*(const Complex &, const CMatrix &);
		friend CMatrix operator*(const CMatrix &, double);
		friend CMatrix operator/(const CMatrix &, double);
		friend CMatrix operator*(double, const CMatrix &);

	        /* Mixed */

//	        friend CVector operator+(const Vector &v, const CVector &c);
	        friend CVector operator+(const CVector &c, const Vector &v);
	        friend CVector operator-(const Vector &v, const CVector &c);
	        friend CVector operator-(const CVector &c, const Vector &v);

//	        friend CVector operator*(const CMatrix &c, const Vector &v);
	        friend CVector operator*(const Matrix &m, const CVector &c);
	        friend CVector operator*(const CVector &c, const Matrix &m);
	        friend CVector operator*(const Vector &v, const CMatrix &c);

                friend CVector operator/(const Vector &v, const Complex &c);
                friend CVector operator*(const Vector &v, const Complex &c);
                friend CVector operator*(const Complex &c, const Vector &v);


		/* miscellaneous */

		Vector Re(void);
		Vector Im(void);
		CVector Conjg(void);



		/* IO */

		void print(ostream&);
		void binwrite(ostream&);
		friend CVector CVector_binread(istream&);



};

void bomb_CVector_operation(const char *);
void bomb_CVector(const char *);
void bomb_CMatrix_operation(const char *);
void bomb_CMatrix(const char *);




class CMatrix
{
        friend CVector;
        friend Matrix;
	friend Vector;
//	private:
        protected:
		int rlow, rhigh;	/* low and high row indices */
		int clow, chigh;	/* low and high column indices */
		CVector *rows;

		/* fast (but unsafe) access functions */

		CVector fastcol(int) const;
		CVector &fastrow(int) const;
		void fastsetrow(int, const CVector &);
		void fastsetcol(int, const CVector &);

	public:

		/* constructors and the destructor */

		CMatrix(void);
		CMatrix(int, int, int, int);
		CMatrix(int, int, int, int, double **);
		CMatrix(int, int, int, int, Complex **);
		CMatrix(const CMatrix &);
		CMatrix(const Matrix &);
		~CMatrix(void);


		
		/* assignment and access */

		CMatrix &operator=(const CMatrix &);
		CVector &operator[](int) const;

		int getnrows(void) {return rhigh-rlow+1;}
		int getncols(void) {return chigh-clow+1;}
		int getrlow(void) {return rlow;}
		int getclow(void) {return clow;}
		int getrhigh(void) {return rhigh;}
		int getchigh(void) {return chigh;}



		/* sizing */

		void setsize(int, int, int, int);
		void zero(void);



		/* row and column access */

		CVector col(int);
		CVector &row(int);
		void setrow(int, const CVector &);
		void setcol(int, const CVector &);


		/* unary plus and minus */

		CMatrix operator-(void);
		CMatrix operator+(void) {return *this;}



		
		
		/* operations with matrices */
		
		friend CMatrix operator+(const CMatrix &, const CMatrix &);
		friend CMatrix operator-(const CMatrix &, const CMatrix &);
		friend CMatrix operator*(const CMatrix &, const CMatrix &);
		CMatrix &operator+=(const CMatrix &);
		CMatrix &operator-=(const CMatrix &);
		CMatrix &operator+=(Matrix &m) 
			{CMatrix c(m); *this += c; return *this;}
		CMatrix &operator-=(Matrix &m) 
			{CMatrix c(m); *this -= c; return *this;}
			


		
		/* operations with vectors */
		
		friend CVector operator*(const CMatrix &, const CVector &);
	        friend CVector operator*(const CMatrix &, const Vector &);
		friend CVector operator*(const CVector &, const CMatrix &);


	
		/* operations with scalars */

		
		friend CMatrix operator*(const CMatrix &, double);
		friend CMatrix operator*(double, const CMatrix &);
		friend CMatrix operator+(const CMatrix &, double);
		friend CMatrix operator+(double, const CMatrix &);
		friend CMatrix operator-(const CMatrix &, double);
		friend CMatrix operator-(double, const CMatrix &);
		friend CMatrix operator/(const CMatrix &, double);
		CMatrix &operator*=(double);
		CMatrix &operator/=(double);

		friend CMatrix operator*(const CMatrix &, const Complex &);
		friend CMatrix operator*(const Complex &, const CMatrix &);
		friend CMatrix operator/(const CMatrix &, const Complex &);
		CMatrix &operator*=(const Complex &);
		CMatrix &operator/=(const Complex &);


	        /* Mixed */

                friend CMatrix operator*(const CMatrix &c, const Matrix &m);
                friend CMatrix operator*(const Matrix &m, const CMatrix &c);

                friend CMatrix operator+(const Matrix &m, const CMatrix &c);
                friend CMatrix operator+(const CMatrix &c, const Matrix &m);
                friend CMatrix operator-(const Matrix &m, const CMatrix &c);
                friend CMatrix operator-(const CMatrix &c, const Matrix &m);

                friend CMatrix operator*(const Matrix &m, const Complex &c);
                friend CMatrix operator*(const Complex &c, const Matrix &m);
                friend CMatrix operator/(const Matrix &m, const Complex &c);



		/* IO */

		void print(ostream&);
		void binwrite(ostream&);
		friend CMatrix CMatrix_binread(istream&);



		/* miscellaneous */

		Matrix Re(void) const;
		Matrix Im(void) const;
		CMatrix Conjg(void) const;	
		CMatrix Transpose(void) const;	
		CMatrix Adjoint(void) const {return Transpose().Conjg();}
		Complex Trace(void) const;


};





// declarations for mixed-type operators


Complex operator*(const Vector &v, const CVector &c);
Complex operator*(const CVector &c, const Vector &v);


	/* miscellaneous functions */

Matrix Re(const CMatrix &c);
Matrix Im(const CMatrix &c);
CMatrix Conjg(const CMatrix &c);
CMatrix Adjoint(const CMatrix &c);
CMatrix Transpose(const CMatrix &c);
Complex Trace(const CMatrix &c);

Vector Re(const CVector &c);
Vector Im(const CVector &c);
CVector Conjg(const CVector &c);


#endif
