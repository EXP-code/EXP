/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine reads in two body files and computes the total
 *  velocity dispersion
 *
 *
 *  Call sequence:
 *  -------------
 *
 *  Parameters:
 *  ----------
 *
 *
 *  Returns:
 *  -------
 *
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  MDW 04/29/92
 *
 ***************************************************************************/


#include <stdio.h>
#include <string.h>
#include "cutil.h"

int parse_args(int, char **);

struct BODIES {
  int n;
  double *mass;
  double *x,*y,*z;
  double *u,*v,*w;
  double time;
};

				/* Fct declaration */

void compute_vel_and_disp
  (struct BODIES *bod1, struct BODIES *bod2,
   double *vel1, double *vel2, double *disp, double *rate);

void read_bodies(char *infile, struct BODIES *bod);

				/* Global variables */

char parmfile[128] = "parm.file";
char INFIL1[128] = "infil1";
char INFIL2[128] = "infil2";

main(argc,argv)
char **argv;
int argc;
{
  struct BODIES bod1,bod2;
  double vel1,vel2,disp,rate;

  /*====================*/
  /* Parse command line */
  /*====================*/

  if (!parse_args(argc,argv)) return -1;

  read_bodies(INFIL1,&bod1);
  read_bodies(INFIL2,&bod2);

  compute_vel_and_disp(&bod1,&bod2,&vel1,&vel2,&disp,&rate);

  printf("%e %e %e %e %e %e %e\n",bod2.time-bod1.time,bod1.time,bod2.time,
	 vel1,vel2,disp,rate);

}


inline double
bod_vel2(struct BODIES *b, int i)
{
  return b->u[i]*b->u[i] + b->v[i]*b->v[i] + b->w[i]*b->w[i];
}


inline double
bod_dvel2(struct BODIES *b1, struct BODIES *b2, int i)
{
  double du,dv,dw;

  du = b2->u[i] - b1->u[i];
  dv = b2->v[i] - b1->v[i];
  dw = b2->w[i] - b1->w[i];

  return du*du + dv*dv + dw*dw;
}

void compute_vel_and_disp
  (struct BODIES *bod1, struct BODIES *bod2,
   double *vel1, double *vel2, double *disp, double *rate)
{
  int i;
  double v1,v2,dv;

  *vel1 = 0.0;
  *vel2 = 0.0;
  *disp = 0.0;
  *rate = 0.0;
  if (bod1->n != bod2->n) return;

  for (i=1; i<=bod1->n; i++) {
    v1 = bod_vel2(bod1,i);
    v2 = bod_vel2(bod2,i);
    dv = bod_dvel2(bod1,bod2,i);
    *vel1 += bod1->mass[i] * v1;
    *vel2 += bod2->mass[i] * v2;
    *disp += bod1->mass[i] * dv;
    if (v1+v2>0.0)
      *rate += bod1->mass[i] * 0.5*dv/(v1+v2);
  }
}

#define BUF 256

void read_bodies(char *infile, struct BODIES *bod)
{
  int i,iret;
  FILE *fin;
  char buf[BUF];

  /* Read in phase space */

  if ( (fin=fopen(infile,"r")) == NULL) {
    fprintf(stderr,"Couldn't open %s . . . quitting\n",infile);
    exit(-1);
  }

  fscanf(fin,"%d %lf\n",&bod->n,&bod->time);

  /* Initialize vectors */
  
  bod->mass = dvector(1,bod->n);
  bod->x = dvector(1,bod->n);
  bod->y = dvector(1,bod->n);
  bod->z = dvector(1,bod->n);
  bod->u = dvector(1,bod->n);
  bod->v = dvector(1,bod->n);
  bod->w = dvector(1,bod->n);
  
  for (i=1; i<=bod->n; i++) {
    fgets(buf,BUF,fin);
    iret = sscanf(buf,"%lf %lf %lf %lf %lf %lf %lf\n",
		  &bod->mass[i],
		  &bod->x[i],
		  &bod->y[i],
		  &bod->z[i],
		  &bod->u[i],
		  &bod->v[i],
		  &bod->w[i] );

    if (iret != 7) {
      fprintf(stderr,"Trouble reading %s . . . quitting\n",infile);
      exit(-1);
    }
  }
}

void usage(char *);
void set_parm(char *, char *);
void write_parm(void);
void print_parm(FILE *, char *);
void print_default(void);

int get_key_value(int, char **, char ***, char ***);
int get_key_value_from_file(char *, char ***, char ***);

int parse_args(argc,argv)
char **argv;
int argc;
{
  extern char *optarg;
  extern int optind;
  char c,*prog=argv[0];
  int iparmf=0,iret,i;
  char file[128],**word,**valu;

  while ((c = getopt(argc, argv, "f:dh")) != -1)
      switch (c) {
      case 'f': 
	iparmf=1;
	strcpy(file,optarg);
	break;
      case 'd':
	print_default();
	break;
      case '?':
      case 'h':
	usage(prog);
	break;
      }

  argc -= optind;
  if (iparmf)
    iret = get_key_value_from_file(file,&word,&valu);
  else
    {
      iret = get_key_value(argc,&argv[optind],&word,&valu);
      argc -= iret;
    }
  
  if (argc != 0)
    usage(prog);

				/* Set parameters */

  for (i=0; i<iret; i++)
    set_parm(word[i],valu[i]);
}

void set_parm(word,valu)
char *word,*valu;
{
  if (!strcmp("INFIL1",word))
    strcpy(INFIL1,valu);

  if (!strcmp("INFIL2",word))
    strcpy(INFIL2,valu);
}

void write_parm()
{
  FILE *fout;

  if ( (fout=fopen(parmfile,"w")) == NULL) {
    fprintf(stderr,"Couldn't open parameter file: %s . . . quitting\n",parmfile);
    exit(-1);
  }

  print_parm(fout,"\0");
  fclose(fout);
}

void print_default()
{
  fprintf(stderr,"\nDefaults:\n");
  fprintf(stderr,"----------------------------\n");
  print_parm(stderr,"\0");
  exit(0);
}


void print_parm(FILE *stream, char *comment)
{
  fprintf(stream,"%sINFIL1 =       %s\n",comment,INFIL1);
  fprintf(stream,"%sINFIL2 =       %s\n",comment,INFIL2);
}


void usage(prog)
char *prog;
{
  char usage_head[] = 
    "[-f file -d] [keyword=value [keyword=value] .. ]";

  char usage_data[][80] = {
    "     -f file",      "keyword/value parameter file",
    "     -d",           "print default parameters",
    "\nKeywords:",       " INFIL1,INFIL2",
    "\0"                 };



  (void)std_usage(prog,usage_head,usage_data);
  exit(-1);
}
