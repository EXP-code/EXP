/*
  Parse command line
*/

#include "expand.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif
#include <unistd.h>
#include <getopt.h>
#include <string.h>

void usage(char *);
void set_parm(char *, char *);
void write_parm(void);
void print_parm(FILE *, char*);
void print_default(void);

int get_key_value(int, char **, char ***, char ***);
int get_key_value_from_file(char *, char ***, char ***);

#define WBUFSIZE 80
#define VBUFSIZE 80

void MPL_parse_args(int argc, char** argv)
{
  extern char *optarg;
  extern int optind;
  char *prog=argv[0];
  int iparmf=0, iret, i, myid;
  char file[128],**word,**valu;
  int c, numb=argc;

  char wbuf[WBUFSIZE];
  char vbuf[WBUFSIZE];

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  if (myid==0) {

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

    if (iparmf)
      iret = get_key_value_from_file(file,&word,&valu);

    if (optind<numb)
      iret = get_key_value(numb-optind,&argv[optind],&word,&valu);

    /*
    fprintf(stderr, "Process 0: numb=%d argc=%d optind=%d iret=%d\n",
	    numb, argc, optind, iret);
   
    for (i=0; i<iret; i++)
      fprintf(stderr, "%3d:  %30s  %30s\n", i, word[i], valu[i]);
    */
      
  }

  MPI_Bcast(&iret, 1, MPI_INT, 0, MPI_COMM_WORLD);

				/* Set parameters */
  

  for (i=0; i<iret; i++) {
    if (myid==0) {
      strcpy(wbuf, word[i]);
      strcpy(vbuf, valu[i]);
      free(word[i]);
      free(valu[i]);
    }
    
				/* Send values to all processes */
    MPI_Bcast(wbuf, WBUFSIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(vbuf, VBUFSIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
    
    set_parm(wbuf, vbuf);

  }

  if (myid==0) {
    free(word);
    free(valu);
  }

				/* Synchronize needed? */
    MPI_Barrier(MPI_COMM_WORLD);
}


void set_parm(char *word, char *valu)
{
  if (!strcmp("seed",word))
    seed = atoi(valu);

  else if (!strcmp("nbodmax",word))
    nbodmax = atoi(valu);

  else if (!strcmp("lmax",word))
    lmax = atoi(valu);

  else if (!strcmp("nmax",word))
    nmax = atoi(valu);

  else if (!strcmp("mmax2",word))
    mmax2 = atoi(valu);

  else if (!strcmp("nmax2",word))
    nmax2 = atoi(valu);

  else if (!strcmp("nmaxx",word))
    nmaxx = atoi(valu);

  else if (!strcmp("nmaxy",word))
    nmaxy = atoi(valu);

  else if (!strcmp("nmaxz",word))
    nmaxz = atoi(valu);

  else if (!strcmp("nminx",word))
    nminx = atoi(valu);

  else if (!strcmp("nminy",word))
    nminy = atoi(valu);

  else if (!strcmp("nminz",word))
    nminz = atoi(valu);

  else if (!strcmp("nfft",word))
    nfft = atoi(valu);

  else if (!strcmp("ncylnzof",word))
    ncylnzof = atoi(valu);

  else if (!strcmp("ncylordz",word))
    ncylordz = atoi(valu);

  else if (!strcmp("ncylrecomp",word))
    ncylrecomp = atoi(valu);

  else if (!strcmp("ncylkeep",word))
    ncylkeep = atoi(valu);

  else if (!strcmp("ncylnx",word))
    ncylnx = atoi(valu);

  else if (!strcmp("ncylny",word))
    ncylny = atoi(valu);

  else if (!strcmp("nsteps",word))
    nsteps = atoi(valu);

  else if (!strcmp("nscale",word))
    nscale = atoi(valu);

  else if (!strcmp("nthrds",word))
    nthrds = atoi(valu);

  else if (!strcmp("dtime",word))
    dtime = atof(valu);

  else if (!strcmp("nlog",word))
    nlog = atoi(valu);

  else if (!strcmp("ncoef",word))
    ncoef = atoi(valu);

  else if (!strcmp("ncoefcyl",word))
    ncoefcyl = atoi(valu);

  else if (!strcmp("nrelx",word))
    nrelx = atoi(valu);

  else if (!strcmp("nlist",word))
    nlist = atoi(valu);

  else if (!strcmp("nscat",word))
    nscat = atoi(valu);

  else if (!strcmp("ntipsy",word))
    ntipsy = atoi(valu);

  else if (!strcmp("ndiag",word))
    ndiag = atoi(valu);

  else if (!strcmp("nchkpt",word))
    nchkpt = atoi(valu);

  else if (!strcmp("ncylamom",word))
    ncylamom = atoi(valu);

  else if (!strcmp("rmax",word))
    rmax = atof(valu);

  else if (!strcmp("zmax",word))
    zmax = atof(valu);

  else if (!strcmp("ecyl0",word))
    ecyl0 = atof(valu);

  else if (!strcmp("EJcyl",word))
    EJcyl = atoi(valu);

  else if (!strcmp("rdiag",word))
    rdiag = atof(valu);

  else if (!strcmp("scale",word))
    scale = atof(valu);

  else if (!strcmp("rcylEM",word))
    rcylEM = atof(valu);

  else if (!strcmp("rcylSL",word))
    rcylSL = atof(valu);

  else if (!strcmp("rsphSL",word))
    rsphSL = atof(valu);

  else if (!strcmp("acyl",word))
    acyl = atof(valu);

  else if (!strcmp("hcyl",word))
    hcyl = atof(valu);

  else if (!strcmp("hslab",word))
    hslab = atof(valu);

  else if (!strcmp("rmax_tidal",word))
    rmax_tidal = atof(valu);

  else if (!strcmp("hills_omega",word))
    hills_omega = atof(valu);

  else if (!strcmp("hills_p",word))
    hills_p = atof(valu);

  else if (!strcmp("bessel_sph",word))
    bessel_sph = atoi(valu);

  else if (!strcmp("c_brock",word))
    c_brock = atoi(valu);

  else if (!strcmp("c_brock_disk",word))
    c_brock_disk = atoi(valu);

  else if (!strcmp("hernq",word))
    hernq = atoi(valu);

  else if (!strcmp("sphereSL",word))
    sphereSL = atoi(valu);

  else if (!strcmp("cube",word))
    cube = atoi(valu);

  else if (!strcmp("slab",word))
    slab = atoi(valu);

  else if (!strcmp("slabSL",word))
    slabSL = atoi(valu);

  else if (!strcmp("cylinder",word))
    cylinder = atoi(valu);

  else if (!strcmp("nulltest",word))
    nulltest = atoi(valu);

  else if (!strcmp("fixpos",word))
    fixpos = atoi(valu);

  else if (!strcmp("fixvel",word))
    fixvel = atoi(valu);

  else if (!strcmp("fixacc",word))
    fixacc = atoi(valu);

  else if (!strcmp("selfgrav",word))
    selfgrav = atoi(valu);

  else if (!strcmp("NO_L1",word))
    NO_L1 = atoi(valu);

  else if (!strcmp("adjust",word))
    adjust = atoi(valu);

  else if (!strcmp("outbods",word))
    outbods = atoi(valu);

  else if (!strcmp("outengy",word))
    outengy = atoi(valu);

  else if (!strcmp("tipsy",word))
    tipsy = atoi(valu);

  else if (!strcmp("olist",word))
    olist = atoi(valu);

  else if (!strcmp("diag",word))
    diag = atoi(valu);

  else if (!strcmp("coef",word))
    coef = atoi(valu);

  else if (!strcmp("coefcyl",word))
    coefcyl = atoi(valu);

  else if (!strcmp("relx",word))
    relx = atoi(valu);

  else if (!strcmp("scatter",word))
    scatter = atoi(valu);

  else if (!strcmp("finish",word))
    finish = atoi(valu);

  else if (!strcmp("tides",word))
    tides = atoi(valu);

  else if (!strcmp("shock",word))
    shock = atoi(valu);

  else if (!strcmp("halo",word))
    halo = atoi(valu);

  else if (!strcmp("disk_on_halo",word))
    disk_on_halo = atoi(valu);

  else if (!strcmp("selector",word))
    selector = atoi(valu);

  else if (!strcmp("tk_type",word))
    tk_type = atoi(valu);

  else if (!strcmp("eigen_tk",word))
    eigen_tk = atoi(valu);

  else if (!strcmp("npca",word))
    npca = atoi(valu);

  else if (!strcmp("pcaout",word))
    pcaout = atoi(valu);

  else if (!strcmp("npcaout",word))
    npcaout = atoi(valu);

  else if (!strcmp("balance",word))
    balance = atoi(valu);

  else if (!strcmp("nbalance",word))
    nbalance = atoi(valu);

  else if (!strcmp("tbalance",word))
    tbalance = atof(valu);

  else if (!strcmp("NICE",word))
    NICE = atoi(valu);

  else if (!strcmp("zerocom",word))
    zerocom = atoi(valu);

  else if (!strcmp("zerovel",word))
    zerovel = atoi(valu);

  else if (!strcmp("L_pca",word))
    L_pca = atoi(valu);

  else if (!strcmp("M_pca",word))
    M_pca = atoi(valu);

  else if (!strcmp("tksmooth",word))
    tksmooth = atof(valu);

  else if (!strcmp("tkcum",word))
    tkcum = atof(valu);

  else if (!strcmp("tauscat",word))
    tauscat = atof(valu);

  else if (!strcmp("self_consistent",word))
    self_consistent = atoi(valu);

  else if (!strcmp("inertia",word))
    inertia = atoi(valu);

  else if (!strcmp("user",word))
    user = atoi(valu);

  else if (!strcmp("pmsoft",word))
    pmsoft = atof(valu);

  else if (!strcmp("U1",word))
    U1 = atof(valu);

  else if (!strcmp("U2",word))
    U2 = atof(valu);

  else if (!strcmp("U3",word))
    U3 = atof(valu);

  else if (!strcmp("U4",word))
    U4 = atof(valu);

  else if (!strcmp("U5",word))
    U5 = atof(valu);

  else if (!strcmp("U6",word))
    U6 = atof(valu);

  else if (!strcmp("homedir",word))
    strcpy(homedir, valu);

  else if (!strcmp("logfile",word))
    strcpy(logfile, valu);

  else if (!strcmp("infile",word))
    strcpy(infile, valu);

  else if (!strcmp("outname",word))
    strcpy(outname, valu);

  else if (!strcmp("parmfile",word))
    strcpy(parmfile, valu);

  else if (!strcmp("coeffile",word))
    strcpy(coeffile, valu);

  else if (!strcmp("coeffilecyl",word))
    strcpy(coeffilecyl, valu);

  else {
    fprintf(stderr, "No such parameter <%s> . . . check header file\n",
	    word);
    exit(-1);
  }

}

void write_parm(void)
{
  FILE *fout;

  if ( (fout=fopen(parmfile,"w")) == NULL) {
    fprintf(stderr,"Couldn't open parameter file: %s . . . quitting\n",parmfile);
    exit(-1);
  }

  print_parm(fout,"\0");
  fclose(fout);
}

void print_default(void)
{
  fprintf(stderr,"\nDefaults:\n");
  fprintf(stderr,"----------------------------\n");
  print_parm(stderr,"\0");
  exit(0);
}


void print_parm(FILE *stream, char *comment)
{
  fprintf(stream,"%s%-15s = %d\n",comment,"seed",seed);
  fprintf(stream,"%s%-15s = %d\n",comment,"nbodmax",nbodmax);
  fprintf(stream,"%s%-15s = %d\n",comment,"lmax",lmax);
  fprintf(stream,"%s%-15s = %d\n",comment,"nmax",nmax);
  fprintf(stream,"%s%-15s = %d\n",comment,"mmax2",mmax2);
  fprintf(stream,"%s%-15s = %d\n",comment,"nmax2",nmax2);
  fprintf(stream,"%s%-15s = %d\n",comment,"nmaxx",nmaxx);
  fprintf(stream,"%s%-15s = %d\n",comment,"nmaxy",nmaxy);
  fprintf(stream,"%s%-15s = %d\n",comment,"nmaxz",nmaxz);
  fprintf(stream,"%s%-15s = %d\n",comment,"nminx",nminx);
  fprintf(stream,"%s%-15s = %d\n",comment,"nminy",nminy);
  fprintf(stream,"%s%-15s = %d\n",comment,"nminz",nminz);
  fprintf(stream,"%s%-15s = %d\n",comment,"nfft",nfft);
  fprintf(stream,"%s%-15s = %d\n",comment,"ncylnzof",ncylnzof);
  fprintf(stream,"%s%-15s = %d\n",comment,"ncylordz",ncylordz);
  fprintf(stream,"%s%-15s = %d\n",comment,"ncylrecomp",ncylrecomp);
  fprintf(stream,"%s%-15s = %d\n",comment,"ncylkeep",ncylkeep);
  fprintf(stream,"%s%-15s = %d\n",comment,"ncylnx",ncylnx);
  fprintf(stream,"%s%-15s = %d\n",comment,"ncylny",ncylny);
  fprintf(stream,"%s%-15s = %d\n",comment,"nsteps",nsteps);
  fprintf(stream,"%s%-15s = %d\n",comment,"nscale",nscale);
  fprintf(stream,"%s%-15s = %d\n",comment,"nthrds",nthrds);
  fprintf(stream,"%s%-15s = %e\n",comment,"dtime",dtime);
  fprintf(stream,"%s%-15s = %d\n",comment,"nlog",nlog);
  fprintf(stream,"%s%-15s = %d\n",comment,"ncoef",ncoef);
  fprintf(stream,"%s%-15s = %d\n",comment,"ncoefcyl",ncoefcyl);
  fprintf(stream,"%s%-15s = %d\n",comment,"nlist",nlist);
  fprintf(stream,"%s%-15s = %d\n",comment,"nscat",nscat);
  fprintf(stream,"%s%-15s = %d\n",comment,"ntipsy",ntipsy);
  fprintf(stream,"%s%-15s = %d\n",comment,"ndiag",ndiag);
  fprintf(stream,"%s%-15s = %d\n",comment,"nchkpt",nchkpt);
  fprintf(stream,"%s%-15s = %d\n",comment,"ncylamom",ncylamom);
  fprintf(stream,"%s%-15s = %d\n",comment,"nrelx",nrelx);
  fprintf(stream,"%s%-15s = %e\n",comment,"rmax",rmax);
  fprintf(stream,"%s%-15s = %e\n",comment,"zmax",zmax);
  fprintf(stream,"%s%-15s = %e\n",comment,"ecyl0",ecyl0);
  fprintf(stream,"%s%-15s = %d\n",comment,"EJcyl",EJcyl);
  fprintf(stream,"%s%-15s = %e\n",comment,"rdiag",rdiag);
  fprintf(stream,"%s%-15s = %e\n",comment,"scale",scale);
  fprintf(stream,"%s%-15s = %e\n",comment,"acyl",acyl);
  fprintf(stream,"%s%-15s = %e\n",comment,"hcyl",hcyl);
  fprintf(stream,"%s%-15s = %e\n",comment,"hslab",hslab);
  fprintf(stream,"%s%-15s = %e\n",comment,"rmax_tidal",rmax_tidal);
  fprintf(stream,"%s%-15s = %e\n",comment,"hills_omega",hills_omega);
  fprintf(stream,"%s%-15s = %e\n",comment,"hills_p",hills_p);
  fprintf(stream,"%s%-15s = %d\n",comment,"bessel_sph",bessel_sph);
  fprintf(stream,"%s%-15s = %d\n",comment,"c_brock",c_brock);
  fprintf(stream,"%s%-15s = %d\n",comment,"c_brock_disk",c_brock_disk);
  fprintf(stream,"%s%-15s = %d\n",comment,"hernq",hernq);
  fprintf(stream,"%s%-15s = %d\n",comment,"sphereSL",sphereSL);
  fprintf(stream,"%s%-15s = %d\n",comment,"cube",cube);
  fprintf(stream,"%s%-15s = %d\n",comment,"slab",slab);
  fprintf(stream,"%s%-15s = %d\n",comment,"slabSL",slabSL);
  fprintf(stream,"%s%-15s = %d\n",comment,"cylinder",cylinder);
  fprintf(stream,"%s%-15s = %d\n",comment,"nulltest",nulltest);
  fprintf(stream,"%s%-15s = %d\n",comment,"fixpos",fixpos);
  fprintf(stream,"%s%-15s = %d\n",comment,"fixvel",fixvel);
  fprintf(stream,"%s%-15s = %d\n",comment,"fixacc",fixacc);
  fprintf(stream,"%s%-15s = %d\n",comment,"selfgrav",selfgrav);
  fprintf(stream,"%s%-15s = %d\n",comment,"NO_L1",NO_L1);
  fprintf(stream,"%s%-15s = %d\n",comment,"adjust",adjust);
  fprintf(stream,"%s%-15s = %d\n",comment,"outbods",outbods);
  fprintf(stream,"%s%-15s = %d\n",comment,"outengy",outengy);
  fprintf(stream,"%s%-15s = %d\n",comment,"tipsy",tipsy);
  fprintf(stream,"%s%-15s = %d\n",comment,"olist",olist);
  fprintf(stream,"%s%-15s = %d\n",comment,"diag",diag);
  fprintf(stream,"%s%-15s = %d\n",comment,"coef",coef);
  fprintf(stream,"%s%-15s = %d\n",comment,"coefcyl",coefcyl);
  fprintf(stream,"%s%-15s = %d\n",comment,"relx",relx);
  fprintf(stream,"%s%-15s = %d\n",comment,"scatter",scatter);
  fprintf(stream,"%s%-15s = %d\n",comment,"finish",finish);
  fprintf(stream,"%s%-15s = %d\n",comment,"tides",tides);
  fprintf(stream,"%s%-15s = %d\n",comment,"shock",shock);
  fprintf(stream,"%s%-15s = %d\n",comment,"halo",halo);
  fprintf(stream,"%s%-15s = %d\n",comment,"disk_on_halo",disk_on_halo);
  fprintf(stream,"%s%-15s = %d\n",comment,"selector",selector);
  fprintf(stream,"%s%-15s = %d\n",comment,"tk_type",tk_type);
  fprintf(stream,"%s%-15s = %d\n",comment,"eigen_tk",eigen_tk);
  fprintf(stream,"%s%-15s = %d\n",comment,"npca",npca);
  fprintf(stream,"%s%-15s = %d\n",comment,"pcaout",pcaout);
  fprintf(stream,"%s%-15s = %d\n",comment,"npcaout",npcaout);
  fprintf(stream,"%s%-15s = %d\n",comment,"balance",balance);
  fprintf(stream,"%s%-15s = %d\n",comment,"nbalance",nbalance);
  fprintf(stream,"%s%-15s = %f\n",comment,"tbalance",tbalance);
  fprintf(stream,"%s%-15s = %d\n",comment,"NICE",NICE);
  fprintf(stream,"%s%-15s = %d\n",comment,"zerocom",zerocom);
  fprintf(stream,"%s%-15s = %d\n",comment,"zerovel",zerovel);
  fprintf(stream,"%s%-15s = %d\n",comment,"L_pca",L_pca);
  fprintf(stream,"%s%-15s = %d\n",comment,"M_pca",M_pca);
  fprintf(stream,"%s%-15s = %f\n",comment,"tksmooth",tksmooth);
  fprintf(stream,"%s%-15s = %f\n",comment,"tkcum",tkcum);
  fprintf(stream,"%s%-15s = %f\n",comment,"tauscat",tauscat);
  fprintf(stream,"%s%-15s = %d\n",comment,"self_consistent",self_consistent);
  fprintf(stream,"%s%-15s = %d\n",comment,"user",user);
  fprintf(stream,"%s%-15s = %f\n",comment,"pmsoft",pmsoft);
  fprintf(stream,"%s%-15s = %f\n",comment,"U1",U1);
  fprintf(stream,"%s%-15s = %f\n",comment,"U2",U2);
  fprintf(stream,"%s%-15s = %f\n",comment,"U3",U3);
  fprintf(stream,"%s%-15s = %f\n",comment,"U4",U4);
  fprintf(stream,"%s%-15s = %f\n",comment,"U5",U5);
  fprintf(stream,"%s%-15s = %f\n",comment,"U6",U6);
  fprintf(stream,"%s%-15s = %d\n",comment,"inertia",inertia);
  fprintf(stream,"%s%-15s = %s\n",comment,"homedir",homedir);
  fprintf(stream,"%s%-15s = %s\n",comment,"logfile",logfile);
  fprintf(stream,"%s%-15s = %s\n",comment,"infile",infile);
  fprintf(stream,"%s%-15s = %s\n",comment,"outname",outname);
  fprintf(stream,"%s%-15s = %s\n",comment,"parmfile",parmfile);
  fprintf(stream,"%s%-15s = %s\n",comment,"coeffile",coeffile);
  fprintf(stream,"%s%-15s = %s\n",comment,"coeffilecyl",coeffilecyl);
}


void usage(char *prog)
{
  void std_usage(char *prog, char *usage_head,char usage_data[][80]);

  char usage_head[] = 
    "[-f file -d] [keyword=value [keyword=value] .. ]";

  char usage_data[][80] = {
    "     -f file",      "keyword/value parameter file",
    "     -d",           "print default parameters",
    "\nKeywords:",       "seed nbodmax lmax nmax nmaxx nmaxy nmaxz nminx,",
    " ",                 "nminy nminz nsteps nscale nthrds dtime nlog ncoef nlist,",
    " ",                 "ndiag nchkpt rmax rdiag scale nfft acyl hcyl cylinder", 
    " ",                 "hslab rmax_tidal hills_omega hills_p mmax2 nmax2,",
    " ",                 "bessel_sph c_brock hernq sphereSLcube slab slabSL cloud,", 
    " ",                 "fixpos fixvel fixacc selfgrav adjust olist outbods,",
    " ",                 "outengy diag coeffinish tides shock halo disk_on_halo,",
    " ",                 "nulltest self_consistent inertia logfile infile outname,",
    " ",                 "parmfile coeffile nsquare soft selector tk_type,",
    " ",                 "tksmooth tkcum eigen_tk npca pcaout npcaout L_pca,",
    " ",                 "M_pca relx nrelx user balance nbalance tbalance,",
    " ",                 "coefcyl ncoefcyl coeffilecyl ncylordz ncylnzof,",
    " ",                 "ncylrecomp ncylkeep ncylnx ncylny ncylmom ecyl0,",
    " ",                 "rcylEM,rcylSL,rsphSL,EJcyl tipsy ntipsy zerocom zerovel,",
    " ",                 "pmsoft scatter tauscat nscat NICE U1 U2 U3 U4 U5 NO_L1",
    " ",                 " ",
    " ",                 "NB: Command line parameters override command file",
    " ",                 "    parameters!",
    "\0"                 };



  (void)std_usage(prog,usage_head,usage_data);
  exit(-1);
}
