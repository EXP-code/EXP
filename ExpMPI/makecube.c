/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine make a unit-cube phase space with given gaussian dispersion
 *  and total mass=1
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
 *  MDW 04/28/92
 *
 ***************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "cutil.h"
				/* Fct. defs */
int parse_args(int, char **);
double rnd_01d(void);
double gasdev(void);
void print_parm(FILE *, char *comment);
int rnd_init(unsigned int);

				/* Global parameters */

int NUMBER = 2000;
double SIGMA = 1.0;
unsigned int SEED=3;

char parmfile[128] = "parm.file";

main(argc,argv)
char **argv;
int argc;
{
  int i;
  double mass,x,y,z,u,v,w;

  /*====================*/
  /* Parse command line */
  /*====================*/

  if (!parse_args(argc,argv)) return -1;

  rnd_init(SEED);

  mass = 1.0/(double)NUMBER;
  printf("%d %e\n",NUMBER,0.0);
  for (i=1; i<=NUMBER; i++) {
    x = rnd_01d();
    y = rnd_01d();
    z = rnd_01d();
    u = SIGMA*gasdev();
    v = SIGMA*gasdev();
    w = SIGMA*gasdev();
    printf("%e %e %e %e %e %e %e\n",mass,x,y,z,u,v,w);
  }

  print_parm(stdout,"! ");

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

void usage(char *);
void set_parm(char *, char *);
void write_parm(void);
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
  if (!strcmp("NUMBER",word))
    NUMBER = atoi(valu);
  if (!strcmp("SEED",word))
    SEED = atoi(valu);
  if (!strcmp("SIGMA",word))
    SIGMA = atof(valu);
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
  fprintf(stream,"%sNUMBER =       %d\n",comment,NUMBER);
  fprintf(stream,"%sSEED =         %d\n",comment,SEED);
  fprintf(stream,"%sSIGMA =        %e\n",comment,SIGMA);
}


void usage(prog)
char *prog;
{
  char usage_head[] = 
    "[-f file -d] [keyword=value [keyword=value] .. ]";

  char usage_data[][80] = {
    "     -f file",      "keyword/value parameter file",
    "     -d",           "print default parameters",
    "\nKeywords:",       " SIGMA,NUMBER,SEED",
    "\0"                 };



  (void)std_usage(prog,usage_head,usage_data);
  exit(-1);
}
