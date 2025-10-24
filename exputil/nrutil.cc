#include <cstdlib>
#include <iostream>
#include "numerical.H"

using namespace std;

void nrerror(const char *error_text)
{
  cerr << "Numerical Recipes run-time error...\n";
  cerr << error_text << endl;
  cerr << "...now exiting to system...\n";
  exit(1);
}



double *nr_vector(int nl, int nh)
{
	double *v;

	v=(double *)malloc((unsigned) ((nh-nl+1)*sizeof(double)));
	if (!v) nrerror("allocation failure in nr_vector()");
	return v-nl;
}

int *nr_ivector(int nl, int nh)
{
	int *v;

	v=(int *)malloc((unsigned) ((nh-nl+1)*sizeof(int)));
	if (!v) nrerror("allocation failure in nr_ivector()");
	return v-nl;
}

double *dvector(int nl, int nh)
{
	double *v;

	v=(double *)malloc((unsigned) ((nh-nl+1)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl;
}



double **nr_matrix(int nrl, int nrh, int ncl, int nch)
{
	int i;
	double **m;

	m=(double **) malloc((unsigned) ((nrh-nrl+1)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in nr_matrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) malloc((unsigned) ((nch-ncl+1)*sizeof(double)));
		if (!m[i]) nrerror("allocation failure 2 in nr_matrix()");
		m[i] -= ncl;
	}
	return m;
}

double **dmatrix(int nrl, int nrh, int ncl, int nch)
{
	int i;
	double **m;

	m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
	if (!m) nrerror("allocation failure 1 in dmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
		if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
}

int **imatrix( int nrl, int nrh, int ncl, int nch)
{
	int i,**m;

	m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
	if (!m) nrerror("allocation failure 1 in imatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
		if (!m[i]) nrerror("allocation failure 2 in imatrix()");
		m[i] -= ncl;
	}
	return m;
}



double **nr_submatrix(double **a, int oldrl, int oldrh, 
	int oldcl, int oldch, int newrl, int newcl)
{
	int i,j;
	double **m;

	m=(double **) malloc((unsigned) (oldrh-oldrl+1)*sizeof(double*));
	if (!m) nrerror("allocation failure in submatrix()");
	m -= newrl;

	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+oldcl-newcl;

	return m;
}



void free_nr_vector(double *v, int nl, int nh)
{
	free(v+nl);
}

void free_nr_ivector(int *v, int nl, int nh)
{
	free(v+nl);
}

void free_dvector(double *v, int nl, int nh)
{
	free(v+nl);
}



void free_nr_matrix(double **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free(m[i]+ncl);
	free(m+nrl);
}

void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free( (m[i]+ncl));
	free((m+nrl));
}

void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free( (m[i]+ncl));
	free( (m+nrl));
}



void free_nr_submatrix(double **b, int nrl, int nrh, int ncl, int nch)
{
	free( (b+nrl));
}



double **convert_nr_matrix(double *a, int nrl, int nrh, int ncl, int nch)
{
	int i,j,nrow,ncol;
	double **m;

	nrow=nrh-nrl+1;
	ncol=nch-ncl+1;
	m = (double **) malloc((unsigned) ((nrow)*sizeof(double*)));
	if (!m) nrerror("allocation failure in convert_matrix()");
	m -= nrl;
	for(i=0,j=nrl;i<=nrow-1;i++,j++) m[j]=a+ncol*i-ncl;
	return m;
}



void free_convert_nr_matrix(double **b, int nrl, int nrh, int ncl, int nch)
{
	free( (b+nrl));
}



