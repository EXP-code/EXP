
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <malloc.h>
#include <numerical.h>
#include <Vector.h>

#include "phase.h"




double Ensemble::total_Energy(void)
{
	return total_Potential() + total_Kinetic();
}





double Ensemble::total_Kinetic(void)
{
	int i;
	double ktot=0.0;

	for (i=0; i<Nstars; i++) 
		ktot += 0.5*stars[i].Mass()*SQR(stars[i].Velocity());

	return ktot;
}





double Ensemble::total_Potential(void)
{
	int i;
	double ptot=0.0;

	/* factor 1/2 accounts for self-energy */

	for (i=0; i<Nstars; i++) 
		ptot += 0.5*stars[i].Mass()*
			(*Phase::Potential)(stars[i].Time(), 
			stars[i].Position());

	return ptot;
}





double Ensemble::Virial(void)
{
	return 2.0*total_Kinetic() + total_Potential();
}




Three_Vector Ensemble::CM_Position(void)
{
	static Three_Vector x_cm;
	double mtot;
	int i;

	x_cm.zero();
	for (i=0; i<Nstars; i++) 
	{
		x_cm += stars[i].Mass()*stars[i].Position();
		mtot += stars[i].Mass();
	}
	x_cm /= mtot;

	return x_cm;
}





Three_Vector Ensemble::CM_Velocity(void)
{
	static Three_Vector v_cm;
	double mtot;
	int i;

	v_cm.zero();
	for (i=0; i<Nstars; i++) 
	{
		v_cm += stars[i].Mass()*stars[i].Velocity();
		mtot += stars[i].Mass();
	}
	v_cm /= mtot;

	return v_cm;
}


Three_Vector Ensemble::total_Angular_Momentum(void)
{
	int i;
	static Three_Vector Jtot;

	Jtot.zero();

	for (i=0; i<Nstars; i++) Jtot += stars[i].Angular_Momentum();

	return Jtot;
}






/* 
	get second moments of all stars within a limiting 
radius rmax.

	M_{i,j} = \sum_{r<rmax} m_{i,j} x_i x_j


*/

Matrix Ensemble::Moment_Tensor(double rmax)
{
	int n, i, j;
	Matrix M(1, 3, 1, 3);
	double r, mtot;
	Three_Vector x, x0;

	M.zero();
	x0.zero();
	mtot = 0.0;



	// find CM position of stars within rmax of origin

	for (n=0; n<Nstars; n++)
	{
		x = stars[n].x;
		if (sqrt(x*x)>rmax) continue;
		x0 += stars[n].m*x;
		mtot += stars[n].m;
	}

	if (mtot==0.0) return M;
	x0 /= mtot;







	// find moment tensor of stars within rmax of CM position


	mtot = 0.0;

	for (n=0; n<Nstars; n++)
	{
		x = stars[n].x-x0;
		r = sqrt(x*x);
		if (r > rmax) continue;
		mtot += stars[n].m;
		for (i=1; i<=3; i++)
		{
			for (j=1; j<=3; j++)
			{
				M[i][j] += stars[n].m*x[i]*x[j];
			}
		}
	}


	if (mtot>0.0) M /= mtot;

	return M;
}


/*
	get the principal axes of the moment tensor
*/


Three_Vector Ensemble::Principal_Axes(double rmax, Matrix &directions)
{
	Matrix M(1, 3, 1, 3);
	static Three_Vector PA;
	Vector V(1, 3);
	V.zero();

	M = Moment_Tensor(rmax);

	if (M.Trace() == 0.0) 
	{
		PA.zero();
		return PA;
	}

	V = M.Symmetric_Eigenvalues(directions);
	PA[1] = V[1];
	PA[2] = V[2];
	PA[3] = V[3];

	return PA;


}


Matrix Ensemble::Inertia_Tensor(double rmax)
{
	double r2;
	Matrix M(1, 3, 1, 3), I(1,3,1,3);

	M = Moment_Tensor(rmax);

	r2 = M.Trace();

	I = -M;
	I[1][1] += r2;
	I[2][2] += r2;
	I[3][3] += r2;

	return I;
}



Three_Vector Ensemble::Solid_Body_Frequency(double rmax)
{
	Matrix I(1,3,1,3), Iinv(1,3,1,3);
	Three_Vector J;
	double mtot;
	int i;


	I = Inertia_Tensor(rmax);
	J.zero();
	mtot = 0.0;

	for (i=0; i<Nstars; i++)
	{
		if (SQR(stars[i].x) < rmax*rmax)
		{
			J += stars[i].m * Cross(stars[i].x, stars[i].v);
			mtot += stars[i].m;
		}
	}
	J /= mtot;

	Iinv = I.Inverse();

	return Iinv*J;
}


