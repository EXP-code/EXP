
extern "C" {
#include <malloc.h>
#include <stdio.h>

void nrerror(char *error_text);
void myerror(char *error_text);
char *malloc(unsigned);
int free(char *ptr);

/* My standard usage template */

void std_usage(char *prog, char *usage_head,char usage_data[][80])
{
  char templte[10] = "%-25s%s\n";
  int i,j;

  fprintf(stderr,"Usage: %s %s\n\n",prog,usage_head);

  j = 0;
  while (*(usage_data[j]) != '\0') {
    fprintf(stderr,templte,usage_data[j],usage_data[j+1]);
    j+=2;
  }
  fprintf(stderr,"\n");
}

/* Check for num options in parsed commmand line */

void chk_arg(char **argv, char *prog, int num)
{
  int i;

  for (i=1; i<=num; i++)
    if (*(argv+i)[0] == '-') usage(prog);
}


/* Most of the alloc error checks were not working,   DFC 6/7/89 */

void myerror(char *error_text)
/* my standard error handler */
{
	char *p=NULL;

	fprintf(stderr,"Last system error was...\n");
	perror(p);
	fprintf(stderr,"\nRun-time error in user program...\n");
	fprintf(stderr,"%s\n",error_text);

	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

void nrerror(char *error_text)
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

float *vector(int nl, int nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	if(NULL == (v=(float *)malloc((unsigned) (nh-nl+1)*sizeof(float))))
	   myerror("allocation failure in vector()");
	v-=nl;
	return v;
}

int *ivector(int nl, int nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	if(NULL == (v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int))))
	   myerror("allocation failure in ivector()");
	v-=nl;
	return v;
}

double *dvector(int nl, int nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	if(NULL == (v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double))))
	   myerror("allocation failure in dvector()");
	v-=nl;
	return v;
}

float **matrix(int nrl, int nrh, int ncl, int nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	int i;
	float **m;

	/* allocate pointers to rows */
	if(NULL == (m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*))))
	   myerror("allocation failure 1 in matrix()");
	m-=nrl;

	/* allocate rows and set pointers to them */
	for(i=nrl;i<=nrh;i++) {
	  if(NULL == (m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float))))
	     myerror("allocation failure 2 in matrix()");
	  m[i]-=ncl;
	}
	/* return pointer to array of pointers to rows */
	return m;
}

double **dmatrix(int nrl, int nrh, int ncl, int nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	int i;
	double **m;

	/* allocate pointers to rows */
	if(NULL == (m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*))))
	   myerror("allocation failure 1 in dmatrix()");
	m-=nrl;

	/* allocate rows and set pointers to them */
	for(i=nrl;i<=nrh;i++) {
		if(NULL == (m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double))))
		  myerror("allocation failure 2 in dmatrix()");
		m[i]-=ncl;
	}
	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(int nrl, int nrh, int ncl, int nch)
                    
/* allocate an int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	int i,**m;

	/* allocate pointers to rows */
	if(NULL == (m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*))))
	  myerror("allocation failure 1 in imatrix()");
	m-=nrl;

	/* allocate rows and set pointers to them */
	for(i=nrl;i<=nrh;i++) {
		if(NULL == (m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int))))
		  myerror("allocation failure 2 in imatrix()");
		m[i]-=ncl;
	}
	/* return pointer to array of pointers to rows */
	return m;
}

float **submatrix(float **a, int oldrl, int oldrh, int oldcl, int oldch, int newrl, int newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
          
                                        
{
	int i,j;
	float **m;

	/* allocate array of pointers to rows */
	if(NULL == (m=(float **) malloc((unsigned) (oldrh-oldrl+1)*sizeof(float*))))
	   myerror("allocation failure in submatrix()");
	   m-=newrl;

	/* set pointers to rows */
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+oldcl-newcl;

	/* return pointer to array of pointers to rows */
	return m;
}

float **convert_matrix(float *a, int nrl, int nrh, int ncl, int nch)
         
                    
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
	int i,j,nrow,ncol;
	float **m;

	nrow=nrh-nrl+1;
	ncol=nch-ncl+1;
	/* allocate pointers to rows */
	if (NULL == (m=(float **) malloc((unsigned) (nrow)*sizeof(float*))))
		myerror("allocation failure in convert_matrix()");
	m -= nrl;

	/* set pointers to rows */
	for(i=0,j=nrl;i<=nrow-1;i++,j++) m[j]=a+ncol*i-ncl;
	/* return pointer to array of pointers to rows */
	return m;
}

void free_vector(float *v, int nl, int nh)
/* free a float vector allocated with vector() */
{
	free((char*) (v+nl));
}

void free_ivector(int *v, int nl, int nh)
             
/* free an int vector allocated with ivector() */
{
	free((char*) (v+nl));
}

void free_dvector(double *v, int nl, int nh)
          
          
/* free a double vector allocated with dvector() */
{
	free((char*) (v+nl));
}

void free_matrix(float **m, int nrl, int nrh, int ncl, int nch)
/* free a float matrix allocated by matrix() */
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
           
                    
/* free a double matrix allocated by dmatrix() */
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch)
        
                    
/* free an int matrix allocated by imatrix() */
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_submatrix(float **b, int nrl, int nrh, int ncl, int nch)
/* free a submatrix allocated by submatrix() */
          
                    
{
	free((char*) (b+nrl));
}

void free_convert_matrix(float **b, int nrl, int nrh, int ncl, int nch)
/* free a matrix allocated by convert_matrix() */
          
                    
{
	free((char*) (b+nrl));
}

}
