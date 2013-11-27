
#include <Vector.h>
#include <CVector.h>


// mixed-type operations. I do these by casting to the more
// general type and calling that class's operator.


// dot product

KComplex operator*(const Vector &v, const CVector &c)
	{CVector tmp(v); return tmp*c;}
KComplex operator*(const CVector &c, const Vector &v)
	{CVector tmp(v); return tmp*c;}


// vector addition


CVector operator+(const Vector &v, const CVector &c)
	{CVector tmp(v); return tmp+c;}

CVector operator+(const CVector &c, const Vector &v)
	{CVector tmp(v); return tmp+c;}

CVector operator-(const Vector &v, const CVector &c)
	{CVector tmp(v); return tmp - c;}

CVector operator-(const CVector &c, const Vector &v)
	{CVector tmp(v); return c - tmp;}


// matrix-vector multiplication

/*
CVector operator*(const CMatrix &c, const Vector &v)
	{CVector tmp(v); return c*tmp;}
	*/

CVector operator*(const Matrix &m, const CVector &c)
	{CMatrix tmp(m); return tmp*c;}

CVector operator*(const CVector &c, const Matrix &m)
	{CMatrix tmp(m); return c*tmp;}

CVector operator*(const Vector &v, const CMatrix &c)
	{CVector tmp(v); return tmp*c;}


// matrix-matrix multiplication

CMatrix operator*(const CMatrix &c, const Matrix &m)
	{CMatrix tmp(m); return c*tmp;}

CMatrix operator*(const Matrix &m, const CMatrix &c)
	{CMatrix tmp(m); return tmp*c;}


// matrix addition

CMatrix operator+(const Matrix &m, const CMatrix &c)
	{CMatrix tmp(m); return tmp+c;}

CMatrix operator+(const CMatrix &c, const Matrix &m)
	{CMatrix tmp(m); return tmp+c;}

CMatrix operator-(const Matrix &m, const CMatrix &c)
	{CMatrix tmp(m); return tmp-c;}

CMatrix operator-(const CMatrix &c, const Matrix &m)
	{CMatrix tmp(m); return c-tmp;}


// multiplication by scalars


CMatrix operator*(const Matrix &m, const KComplex &c)
	{CMatrix tmp(m); return tmp*c;}

CMatrix operator*(const KComplex &c, const Matrix &m)
	{CMatrix tmp(m); return tmp*c;}

CVector operator*(const Vector &v, const KComplex &c)
	{CVector tmp(v); return tmp*c;}

CVector operator*(const KComplex &c, const Vector &v)
	{CVector tmp(v); return tmp*c;}

CMatrix operator/(const Matrix &m, const KComplex &c)
	{CMatrix tmp(m); return tmp/c;}

CVector operator/(const Vector &v, const KComplex &c)
	{CVector tmp(v); return tmp/c;}


	
