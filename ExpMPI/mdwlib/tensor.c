/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  These routines allocate and maintain arbitary order n-tensors of type
 *  double.  The tensor is allocated and freed by "ntensor" and "free_ntensor"
 *  respectively.  "write_ntensor" and "read_ntensor" will write and read
 *  an n-tensor from a binary file.
 *
 *
 *  Call sequence:
 *  -------------
 *  p = (double *..*)ntensor(n,low,high);
 *  void free_ntensor(p,n,low,high);
 *  int n,*low,*high;            The desired range of indicies for each 
 *                               dimension 1 --> n are passed in the integer
 *                               arrays low[1]-->low[n] and high[1]-->high[n].
 *
 *  double *..*p;
 *
 *  void write_ntensor(filename,p,n,low,high);
 *  void read_ntensor(filename,p,&n,low,high);        (p must be fully
 *                                                     allocated before use)
 *
 *
 *  double *p;
 *  int nrl,nrh;                 writes/read p[nrl] to p[nrh], nrl <= nrh.
 *                               the origin of p is assumed to be p[0]
 *
 *  void write_vector(filename,p, nrl, nrh);
 *  void read_vector (filename,p, nrl, nrh);          (p must be fully
 *                                                     allocated prior to use)
 *
 *
 *  Parameters:
 *  ----------
 *
 *  Debug routines are enabled by compiling with the flag -DDEBUG.  There
 *  are 4 other flages -D{SINGLE,DOUBLE,TRIPLE,QUADRUPLE} (only one should
 *  be used at a time) which checks allocation of that type.
 *
 *  Returns:
 *  -------
 *
 *  Value
 *
 *  Notes:
 *  -----
 *
 *  The arrays low[] and high[] must be allocated separately and must have
 *  large enough size to buffer any rank read in using read_ntensor();
 *  Compile sequence:
 *       cc -c tensor.c
 *
 *       cc -c -DDEBUG -D{SINGLE,DOUBLE,TRIPLE,QUADRUPLE}
 *
 *  By:
 *  --
 *
 *  MDW 11/13/88
 *  DFC  5/17/89
 *  DFC  6/4/89     check return codes on malloc's and free's
 *
 ***************************************************************************/

/*#include <malloc.h>*/
#include <stdio.h>
#include "cutil.h"

char *malloc(unsigned);
int free(char *);

char *ntensor(int n, int *low, int *high);
void free_ntensor(char **p, int n, int *low, int *high),write_ntensor(char *filename, char **v, int n, int *low, int *high),read_ntensor(char *filename, char **v, int *n, int *low, int *high);

/* 
  For tensor of rank n, call:

      ntensor(n,low,high);

*/

char *ntensor_r(int n, int m, int *low, int *high);
void free_ntensor_r(char **p, int n, int m, int *low, int *high),write_ntensor_r(FILE *fout, char **p, int n, int m, int *low, int *high),read_ntensor_r(FILE *fout, char **p, int n, int m, int *low, int *high);


/* top level user callable */

char *ntensor(int n, int *low, int *high)
{
  return ntensor_r(n,0,low+1,high+1);
}


/* recursing version */

char *ntensor_r(int n, int m, int *low, int *high)
{
  int i;
  char **p;
  int malloc_verify();

  if (m+1 == n) {
    if(NULL == (p = (char **)
		( (double *) malloc((unsigned) (high[m]-low[m]+1)*sizeof(double))))){
      myerror("ntensor_r: alloc failure \n");
    }
/*     malloc_verify();                       */
/*     printf("alloc double stor %d\n",p);    */
    p=(char **)( (double *)p - low[m]);
  }
  else {
    if(NULL == (p=(char **)
		malloc((unsigned) (high[m]-low[m]+1)*sizeof(char*)))){
      myerror("ntensor_r: alloc failure \n");
    }
/*    malloc_verify();                        */
/*     printf("alloc pointer stor %d\n",p);   */
    p -= low[m];
    for(i=low[m]; i<=high[m]; i++) {
      p[i] = ntensor_r(n,m+1,low,high);
    }

  }
  
  return (char *)p;

}


/* top level user callable */

void free_ntensor(char **p, int n, int *low, int *high)
{
  free_ntensor_r(p,n,0,low+1,high+1);
}

/* recursing version */

void free_ntensor_r(char **p, int n, int m, int *low, int *high)
{
  int i;

  if (m+1 == n){
    free((char*) ((double *)p+low[m]));
/*      printf("freed storage %d ",p);  */
  }
  else {
    for (i=low[m]; i<=high[m]; i++)
      free_ntensor_r((char **)p[i],n,m+1,low,high);
    free((char*) (p+low[m]));
/*      printf("freed pointer %d ",p);  */
  }
  
}



void write_ntensor(char *filename, char **v, int n, int *low, int *high)
{
  FILE *fout;

  if ( (fout=fopen(filename,"w")) == NULL ) {
    fprintf(stderr,"Couldn't open %s . . . returning\n",filename);
    return;
  }

  fwrite(&n,sizeof(int),1,fout);
  fwrite(low+1,sizeof(int),n,fout);
  fwrite(high+1,sizeof(int),n,fout);
  write_ntensor_r(fout,v,n,0,low+1,high+1);
  fclose(fout);
}

void read_ntensor(char *filename, char **v, int *n, int *low, int *high)
{
  FILE *fin;


  if ( (fin=fopen(filename,"r")) == NULL ) {
    fprintf(stderr,"Couldn't open %s . . . returning\n",filename);
    return;
  }

  fread(n,sizeof(int),1,fin);
  fread(low+1,sizeof(int),*n,fin);
  fread(high+1,sizeof(int),*n,fin);
  read_ntensor_r(fin,v,*n,0,low+1,high+1);
  fclose(fin);
}


void write_ntensor_r(FILE *fout, char **p, int n, int m, int *low, int *high)
{
  int i;

  if (m+1 == n) {
    for (i=low[m]; i<=high[m]; i++)
      fwrite(((double *)p+i),sizeof(double),1,fout);
  }
  else {
    for(i=low[m]; i<=high[m]; i++)
      write_ntensor_r(fout,(char **)p[i],n,m+1,low,high);
  }

}

void read_ntensor_r(FILE *fout, char **p, int n, int m, int *low, int *high)
{
  int i;

  if (m+1 == n) {
    for (i=low[m]; i<=high[m]; i++)
      fread(((double *)p+i),sizeof(double),1,fout);
  }
  else {
    for(i=low[m]; i<=high[m]; i++)
      read_ntensor_r(fout,(char **)p[i],n,m+1,low,high);
  }

}

void write_vector(char *filename, double *v, int *nrl, int *nrh)
{
  FILE *fout;

  if ( (fout=fopen(filename,"w")) == NULL ) {
    fprintf(stderr,"Couldn't open %s . . . returning\n",filename);
    return;
  }

  fwrite(nrl,sizeof(int),1,fout);
  fwrite(nrh,sizeof(int),1,fout);
  write_ntensor_r(fout,(char **)v,1,0,nrl,nrh);

  fclose(fout);
}

void read_vector(char *filename, double *v, int *nrl, int *nrh)
{
  FILE *fin;

  if ( (fin=fopen(filename,"r")) == NULL ) {
    fprintf(stderr,"Couldn't open %s . . . returning\n",filename);
    return;
  }

  fread(nrl,sizeof(int),1,fin);
  fread(nrh,sizeof(int),1,fin);
  read_ntensor_r(fin,(char **)v,1,0,nrl,nrh);
  fclose(fin);
}



/*              commented it all out to get rid of main() interactions 
#ifdef DEBUG

main()
{
#ifdef SINGLE
  double *v;
#endif
#ifdef DOUBLE
  double **v;
#endif
#ifdef TRIPLE
  double ***v,***w;
#endif
#ifdef QUADRUPLE
  double ****v;
#endif
  int *low,*high,i,j,k,l,n;
  
#ifdef SINGLE
  low = ivector(1,1);
  high = ivector(1,1);

  low[1] = 1;
  high[1] = 5;

  v = (double *)ntensor(1,low,high);

  for (i=low[1]; i<=high[1]; i++)
    v[i] = (double)i;

  for (i=low[1]; i<=high[1]; i++) {
      printf(" %le\n",v[i]);
  }

  free_ntensor(v,1,low,high);
#endif

#ifdef DOUBLE
  low = ivector(1,2);
  high = ivector(1,2);

  low[1] = 1;
  high[1] = 5;
  low[2] = 1;
  high[2] = 5;

  v = (double **)ntensor(2,low,high);

  for (i=low[1]; i<=high[1]; i++) {
    for (j=low[2]; j<=high[2]; j++) {
      v[i][j] = (double)j + 100.0*i;
    }
  }

  for (i=low[1]; i<=high[1]; i++) {
    for (j=low[2]; j<=high[2]; j++)
      printf(" %le",v[i][j]);
    printf("\n");
  }

  free_ntensor(v,2,low,high);

#endif

#ifdef TRIPLE
  low = ivector(1,3);
  high = ivector(1,3);

  low[1] = 1;
  high[1] = 10;
  low[2] = 1;
  high[2] = 5;
  low[3] = 1;
  high[3] = 5;

  v = (double ***)ntensor(3,low,high);
  w = (double ***)ntensor(3,low,high);

  for (i=low[1]; i<=high[1]; i++) {
    for (j=low[2]; j<=high[2]; j++) {
      for (k=low[3]; k<=high[3]; k++) {
	v[i][j][k] = (double)k + 10.0*j + 100.0*i;
      }
    }
  }

  n = 3;
  write_ntensor("/tmp/temp.tmp",v,3,low,high);
  read_ntensor("/tmp/temp.tmp",w,&n,low,high);

  for (i=low[1]; i<=high[1]; i++) {
    printf("\n\n==> %d:\n",i);
    for (j=low[2]; j<=high[2]; j++) {
      for (k=low[3]; k<=high[3]; k++)
	printf(" %le",w[i][j][k]);
      printf("\n");
    }
  }

  free_ntensor(v,3,low,high);
  free_ntensor(w,3,low,high);
#endif

#ifdef QUADRUPLE
  low = ivector(1,4);
  high = ivector(1,4);

  low[1] = 1;
  high[1] = 3;
  low[2] = 1;
  high[2] = 5;
  low[3] = 1;
  high[3] = 5;
  low[4] = 1;
  high[4] = 5;

  v = (double ****)ntensor(4,low,high);

  for (i=low[1]; i<=high[1]; i++) {
    for (j=low[2]; j<=high[2]; j++) {
      for (k=low[3]; k<=high[3]; k++) {
	for (l=low[4]; l<=high[4]; l++) {
	  v[i][j][k][l] = (double)l + 10.0*k + 100.0*j + 1000.0*i;
	}
      }
    }
  }

  printf("Will print out highest first index subspace . . . \n");

  i = high[1];
  for (j=low[2]; j<=high[2]; j++) {
    printf("\n\n==> %d:\n",j);
    for (k=low[3]; k<=high[3]; k++) {
      for (l=low[4]; l<=high[4]; l++)
	printf(" %le",v[i][j][k][l]);
      printf("\n");
    }
  }

  free_ntensor(v,4,low,high);

#endif

}

#endif DEBUG

                   end to comment out cause of main()    */
