/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Analyze output of expand simulation
 *
 *  Call sequence:
 *  -------------
 *  analcd [-f parameter file -t -d -n] [ [keyword=value] [keyword=value] ... ]
 *
 *
 *  Parameters:
 *  ----------
 *  The parameter file may be used instead of the keyword/value pairs.
 *  The available keywords are: 
 *
 *  KEYWORD        EXPLANATION                            DEFAULT
 *  ----------     ----------------------------           --------
 *  infile =       input p-s file                         IN.FILE
 *  outname =      output file prefix                     OUT
 *  parmfile =     parameter file name                    PARMS.FILE
 *
 *
 *  Compile options:
 *  ---------------
 *
 *  FPETRAP        turn on floating point trapping
 *                 (for debugger only)
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  MDW 11/13/91
 *
 ***************************************************************************/
#define MAIN 1
#include "expand.h"
#define BUF 512

int nbins=20;
int nbeg=0;
int nend=100;
int tally=0;
int nprint=1;

static char rcsid[] = "$Id$";

int parse_args(int, char **);

main(int argc, char **argv)
{
  int i,indx,iret,istr;
  double **sortrad,**sortmas,**sortmet,**sortsiz;
  double tmp,dn,dm,imas,rad,mas,siz,met,mstars;
  char buf[BUF],filename[80];
  FILE *fin,*fout1,*fout2,*fout3;

  /*====================*/
  /* Parse command line */
  /*====================*/

  if (!parse_args(argc,argv)) return -1;
	
  
  /*===================-=*/
  /* Read in phase space */
  /*===================-=*/

  
  if ( (fin=fopen(infile,"r")) == NULL) {
    fprintf(stderr,"Couldn't open %s . . . quitting\n",infile);
    exit(-1);
  }

  fscanf(fin,"%d %lf\n",&nbodies,&tnow);


  /* Initialize vectors */

  mass = dvector(1,nbodies);
  x = dvector(1,nbodies);
  y = dvector(1,nbodies);
  z = dvector(1,nbodies);
  rr = dvector(1,nbodies);
  vx = dvector(1,nbodies);
  vy = dvector(1,nbodies);
  vz = dvector(1,nbodies);
  size = dvector(1,nbodies);
  metal = dvector(1,nbodies);
  

  /* Read in points */

  fgets(buf,BUF,fin);
  iret = sscanf(buf,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
      &mass[1],&x[1],&y[1],&z[1],&vx[1],&vy[1],&vz[1],&tmp,&size[1],&metal[1]);
  if (iret==9) cloud = 1;
  else if (iret==10) cloud = metallicity = 1;
  else {
    fprintf(stderr,"Trouble reading %s . . . quitting\n",infile);
    exit(-1);
  }

  for (i=2; i<=nbodies; i++) {
    if (cloud && !metallicity)
      fscanf(fin,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
	     &mass[i],&x[i],&y[i],&z[i],&vx[i],&vy[i],&vz[i],&tmp,&size[i]);
    if (cloud && metallicity)
      fscanf(fin,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	     &mass[i],&x[i],&y[i],&z[i],&vx[i],&vy[i],&vz[i],&tmp,&size[i],&metal[i]);
    else
      fscanf(fin,"%lf %lf %lf %lf %lf %lf %lf %lf",
	     &mass[i],&x[i],&y[i],&z[i],&vx[i],&vy[i],&vz[i],&tmp);
    rr[i] = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
  }



  /*============================*/
  /* Compute desired quantities */
  /*============================*/

				/* Form cumulative distributions in
				   mass, size, metallicity, radius */

  sortrad =   dmatrix(1,2,1,nbodies);
  sortsiz =   dmatrix(1,2,1,nbodies);
  sortmet =   dmatrix(1,2,1,nbodies);
  sortmas =   dmatrix(1,2,1,nbodies);
  for (i=1; i<=nbodies; i++) {
    sortrad[1][i] = sortsiz[1][i] = sortmet[1][i] = sortmas[1][i] = mass[i];
    sortrad[2][i] = rr[i];
    sortsiz[2][i] = size[i];
    sortmet[2][i] = metal[i];
    sortmas[2][i] = mass[i];
  }

  sort2(nbodies,sortrad[2],sortrad[1]);
  sort2(nbodies,sortsiz[2],sortsiz[1]);
  sort2(nbodies,sortmet[2],sortmet[1]);
  sort2(nbodies,sortmas[2],sortmas[1]);

  for (i=2; i<=nbodies; i++) {
    sortrad[1][i] += sortrad[1][i-1];
    sortsiz[1][i] += sortsiz[1][i-1];
    sortmet[1][i] += sortmet[1][i-1];
    sortmas[1][i] += sortmas[1][i-1];
  }

  if (nprint) {
    sprintf(filename,"%s.number",outname);
    if ( (fout1=fopen(filename,"w")) == NULL) {
      fprintf(stderr,"Couldn't open %s . . . quitting\n",filename);
      exit(-1);
    }

    sprintf(filename,"%s.mass",outname);
    if ( (fout2=fopen(filename,"w")) == NULL) {
      fprintf(stderr,"Couldn't open %s . . . quitting\n",filename);
      exit(-1);
    }
  }

  if (tally) {

    for (istr=1; istr<=nbodies; istr++) if (sortsiz[2][istr]>0.0) break;
    if (istr==1) 
      mstars = 0.0;
    else
      mstars = sortsiz[1][istr-1];

    sprintf(filename,"%s.tally",outname);
    if ( (fout3=fopen(filename,"a+")) == NULL) {
      fprintf(stderr,"Couldn't open %s . . . quitting\n",filename);
      exit(-1);
    }

    fprintf(fout3," %e %e %e %e %e %e %e\n",
	    tnow,mstars,sortsiz[2][istr],sortsiz[2][nbodies],
	    sortmas[2][1],sortmas[2][nbodies],sortmet[2][nbodies]);

    fclose(fout3);

    fprintf(stderr,"total mass =%e\n",sortmas[1][nbodies]);
    fprintf(stderr,"star mass  =%e\n",mstars);
    fprintf(stderr,"max mass   =%e\n",sortmas[2][nbodies]);
    fprintf(stderr,"max size   =%e\n",sortsiz[2][nbodies]);
    fprintf(stderr,"max radius =%e\n",sortrad[2][nbodies]);
    fprintf(stderr,"max metal  =%e\n",sortmet[2][nbodies]);
    if (nprint==0) return 0;

    fprintf(fout2,"! total mass =%e\n",sortmas[1][nbodies]);
    fprintf(fout2,"! star mass  =%e\n",mstars);
    fprintf(fout2,"! max mass   =%e\n",sortmas[2][nbodies]);
    fprintf(fout2,"! max size   =%e\n",sortsiz[2][nbodies]);
    fprintf(fout2,"! max radius =%e\n",sortrad[2][nbodies]);
    fprintf(fout2,"! max metal  =%e\n",sortmet[2][nbodies]);
  }

  fprintf(fout1,"! num    index     radius     mass     size   metallicity\n");
  fprintf(fout2,"! num    tot.mass  radius     mass     size   metallicity\n");

  nbeg = MAX(0,nbeg);
  nend = MIN(nend,100);
  dn = (nend-nbeg)/100.0 * (double)nbodies/nbins;
  dm = (nend-nbeg)/100.0 * sortmas[1][nbodies]/nbins;
  for (i=1; i<=nbins; i++) {
    indx = nbeg*nbodies/100 + (int)(dn*((double)i-0.5));
    fprintf(fout1," %d %d %e %e %e %e\n",
	    i,indx,sortrad[2][indx],sortmas[2][indx],sortsiz[2][indx],
	    sortmet[2][indx]);

    imas = nbeg/100.0 * sortmas[1][nbodies] + dm*((double)i-0.5);

    locate_with_guard(sortrad[1],nbodies,imas,&indx);
    rad = sortrad[2][indx];

    locate_with_guard(sortmas[1],nbodies,imas,&indx);
    mas = sortmas[2][indx];

    locate_with_guard(sortsiz[1],nbodies,imas,&indx);
    siz = sortsiz[2][indx];

    locate_with_guard(sortmet[1],nbodies,imas,&indx);
    met = sortmet[2][indx];

    fprintf(fout2," %d %e %e %e %e %e\n",i,imas,rad,mas,siz,met);

  }


}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

#include <string.h>

void usage(char *);
void set_parm(char *, char *);
void write_parm(void);
void print_parm(FILE *);
void print_default(void);

int get_key_value(int, char **, char ***, char ***);
int get_key_value_from_file(char *, char ***, char ***);

int parse_args(int argc, char **argv)
{
  char *c,*prog=argv[0];
  int iparmf=0,iret,i;
  char file[128],**word,**valu;

  while (--argc > 0 && (*++argv)[0] == '-')
    for (c = argv[0]+1; *c != '\0'; c++)
      switch (*c) {
      case 'f': 
	argc--;
	iparmf=1;
	strcpy(file,*++argv);
	break;
      case 'n':
	nprint=0;
	break;
      case 't':
	tally=1;
	break;
      case 'd':
	print_default();
	break;
      default:
	fprintf(stderr,"%s: illegal option %c\n",prog,*c);
      case 'h':
	usage(prog);
	break;
      }

  if (iparmf)
    iret = get_key_value_from_file(file,&word,&valu);
  else
    {
      iret = get_key_value(argc,argv,&word,&valu);
      argc -= iret;
    }
  
  if (argc != 0)
    usage(prog);

				/* Set parameters */

  for (i=0; i<iret; i++)
    set_parm(word[i],valu[i]);
}

void set_parm(char *word, char *valu)
{
  
  if (!strcmp("infile",word))
    strcpy(infile,valu);

  if (!strcmp("outname",word))
    strcpy(outname,valu);
  
  if (!strcmp("nbins",word))
    nbins = atoi(valu);
  
  if (!strcmp("nbeg",word))
    nbeg = atoi(valu);
  
  if (!strcmp("nend",word))
    nend = atoi(valu);
  
  if (!strcmp("tally",word))
    tally = atoi(valu);

  if (!strcmp("nprint",word))
    nprint = atoi(valu);
  
}

void write_parm(void)
{
  FILE *fout;

  if ( (fout=fopen(parmfile,"w")) == NULL) {
    fprintf(stderr,"Couldn't open parameter file: %s . . . quitting\n",parmfile);
    exit(-1);
  }

  print_parm(fout);
  fclose(fout);
}

void print_default(void)
{
  fprintf(stderr,"\nDefaults:\n");
  fprintf(stderr,"----------------------------\n");
  print_parm(stderr);
  exit(0);
}


void print_parm(FILE *stream)
{
  fprintf(stream,"infile =       %s\n",infile);
  fprintf(stream,"outname =      %s\n",outname);
  fprintf(stream,"nbins =        %d\n",nbins);
  fprintf(stream,"nbeg =         %d\n",nbeg);
  fprintf(stream,"nend =         %d\n",nend);
  fprintf(stream,"tally =        %d\n",tally);
  fprintf(stream,"nprint =       %d\n",nprint);
}


void usage(char *prog)
{
  char usage_head[] = 
    "[-f file -d -t] [keyword=value [keyword=value] .. ]";

  char usage_data[][80] = {
    "     -f file",      "keyword/value parameter file",
    "     -d",           "print default parameters",
    "     -t",           "print out maximum values",
    "\nKeywords:",       " infile,outname,nbins,nbeg,nend,tally,nprint",
    "\0"                 };



  (void)std_usage(prog,usage_head,usage_data);
  exit(-1);
}
