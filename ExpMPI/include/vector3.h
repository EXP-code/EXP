
#ifndef _vector_3_h

#define _vector_3_h

class Three_Vector;
class Vector;

class Three_Vector
{
protected:
	double x[3];
public:
	Three_Vector(void)
		{;}
	Three_Vector(Three_Vector &v)
		{x[0]=v.x[0]; x[1]=v.x[1]; x[2]=v.x[2];}
	Three_Vector(Vector &);

	void zero(void)		
		{x[0]=0.0; x[1]=0.0; x[2]=0.0;}
		
	double &operator[](int i) {return x[i-1];}

	Three_Vector &operator=(Vector &);

	Three_Vector &operator=(Three_Vector &v)	
		{x[0]=v.x[0]; x[1]=v.x[1]; x[2]=v.x[2]; return *this;}
	
	Three_Vector &operator+=(Three_Vector &v)
		{x[0]+=v.x[0]; x[1]+=v.x[1]; x[2]+=v.x[2]; return *this;}
	Three_Vector &operator-=(Three_Vector &v)
		{x[0]-=v.x[0]; x[1]-=v.x[1]; x[2]-=v.x[2]; return *this;}
	Three_Vector &operator*=(double a)
		{x[0]*=a; x[1]*=a; x[2]*=a; return *this;}
	Three_Vector &operator/=(double a)
		{x[0]/=a; x[1]/=a; x[2]/=a; return *this;}

	/* binary operations with Three_Vectors */

	friend Three_Vector operator+(Three_Vector &, Three_Vector &);
	friend Three_Vector operator-(Three_Vector &, Three_Vector &);
	friend Three_Vector operator^(Three_Vector &, Three_Vector &);



	/* binary operations with doubles */

	friend Three_Vector operator*(Three_Vector &, double);
	friend Three_Vector operator*(double, Three_Vector &);
	friend Three_Vector operator/(Three_Vector &, double);



	/* binary operations with Vectors */

	friend Three_Vector operator+(Vector &, Three_Vector &);
	friend Three_Vector operator-(Vector &, Three_Vector &);
	friend Three_Vector operator^(Vector &, Three_Vector &);
	friend Three_Vector operator+(Three_Vector &, Vector &);
	friend Three_Vector operator-(Three_Vector &, Vector &);
	friend Three_Vector operator^(Three_Vector &, Vector &);


	friend double operator*(Three_Vector &v1, Three_Vector &v2)
		{return v1.x[0]*v2.x[0] + v1.x[1]*v2.x[1] + v1.x[2]*v2.x[2];}

	/* magnitude of vector */

	friend double fabs(Three_Vector& v) {return sqrt(v*v);}

	void print(ostream&);
	void grab(void);

};

#endif
