#include <math.h>
#include <unistd.h>
#include <numerical.h>
#include <Vector.h>



Matrix build_cubic_table(Vector &xt, Vector &ft, Vector &dft)
{
	int i;
	int n;
	Matrix tmp;
	double dx;



	
	/* find size of input vectors */

	n= xt.gethigh();
	if (ft.gethigh() != n|| dft.gethigh() != n)
	{
		nrerror("input mismatch in Cubic_Table::build()");
	}




	/* assign sizes of everything */

	tmp.setsize(1, n-1, 0, 3);

	



	/* work out the coefficients of the cubic at each grid point */
	
	for (i=1; i<n; i++)
	{
		dx = xt[i+1] - xt[i];
		tmp[i][0] = ft[i];
		tmp[i][1] = dft[i]*dx;
		tmp[i][2] = -(2.0*dft[i] + dft[i+1])*dx 
			- 3.0*(ft[i] - ft[i+1]);
		tmp[i][3] = (dft[i] + dft[i+1])*dx + 2.0*(ft[i]-ft[i+1]);
	}

	return tmp;
}





/*
	
	Look up a function value from the cubic table. Also get its 
derivative.

*/



double cubic_value(Matrix &a, Vector &xgrid, double x, double *df)
{
	double u, f;
	int loc, n;
	double dx;


	n = xgrid.gethigh();

	locate(xgrid.array(1, n), n, x, &loc);

	if (loc==0 || loc==n)
	{
		cerr.form("point %12.5lg off interpolation grid\n", x);
		exit(0);
	}

	dx = xgrid[loc+1] - xgrid[loc];
	u = (x - xgrid[loc])/dx;

	f = a[loc][0] + u*(a[loc][1] + 
		u*(a[loc][2] + u*a[loc][3]));

	*df = a[loc][1] + u*(2.0*a[loc][2] + u*3.0*a[loc][3]);
	*df /= dx;


	return f;
}


	
	
	

	
	
	
