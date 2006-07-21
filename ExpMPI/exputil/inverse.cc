
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <numerical.h>
#include <Vector.h>

void ludcmp(double **, int, int *, double *);
void lubksb(double **, int, int *, double *);


void destroy_Matrix_array(double **ptr, int rl, int rh, int cl, int ch);

Matrix Matrix::Inverse(void)
{
	double *col, d, **a;
	int i, j, n, *indx;
	
	if (rlow != 1 || clow != 1)
	{
		puts("Inverse expects a unit-offset matrix");
		exit(0);
	}

	if (rhigh != chigh)
	{
		puts("I can't take the inverse of a non-square matrix");
		exit(0);
	}

	n = rhigh;
	Matrix ans(1, n, 1, n);

	a = array_copy(1, n, 1, n);
	col = nr_vector(1, n);
	indx = nr_ivector(1, n);

	ludcmp(a, n, indx, &d);

	for (j=1; j<=n; j++)
	{
		for (i=1; i<=n; i++) col[i]=0.0;
		col[j] = 1.0;
		lubksb(a, n, indx, col);
		for (i=1; i<=n; i++) ans[i][j] = col[i];
	}


	destroy_Matrix_array(a, 1, n, 1, n);
	free_nr_vector(col, 1, n);
	free_nr_ivector(indx, 1, n);

	return ans;
}




void ludcmp(double **a, int n,int *indx,double *d)
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;
//	double *vector();
//	void nrerror(),free_vector();

	vv=nr_vector(1,n);
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine LUDCMP");
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	free_nr_vector(vv,1,n);
}


	

void lubksb(double **a, int n, int *indx, double *b)
{
	int i,ii=0,ip,j;
	double sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}
	

	


	

	
	













