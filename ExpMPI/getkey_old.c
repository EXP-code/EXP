/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine gets key words on command line of the form <x=y> and
 *  returns two char vectors of keyword and value
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
 *  Number of key word--value combinations
 *
 *
 *  Notes:
 *  -----
 *
 *  There is no error checking here.  If you don't have keyword pairs 
 *  exactly of the form:
 *                          key=word
 *
 *  the results are unpredictable.  Similarly, if using the file form,
 *  there should be no white space save \n's.  (Yeah, this should be
 *  fixed).
 *
 *  By:
 *  --
 *
 *  MDW 11/13/91
 *
 ***************************************************************************/

#include <stdio.h>
#include <string.h>

#ifdef DEBUG

usage(prog)
char *prog;
{
  fprintf(stderr,"Usage: %s -f file [keyword pairs]\n",prog);
  exit(-1);
}

main(argc,argv)
char **argv;
int argc;
{
  char *c,*prog=argv[0];
  char **word,**valu,file[256];
  int iret,i,iflgf=0;
	
  while (--argc > 0 && (*++argv)[0] == '-')
    for (c = argv[0]+1; *c != '\0'; c++)
      switch (*c) {
      case 'f': 
	iflgf=1;
	argc--;
	strcpy(file,*++argv);
	break;
      default:
	fprintf(stderr,"%s: illegal option %c\n",prog,*c);
      case 'h':
	usage(prog);
	break;
      }

  if (iflgf)
    iret = get_key_value_from_file(file,&word,&valu);
  else {
    iret = get_key_value(argc,argv,&word,&valu);
		/* decrement argument counter by # keys */
    argc -= iret;
  }

  printf("%d pairs found, arguments left: %d\n",iret,argc);
  for (i=0; i<iret; i++) printf("%30s       %30s\n",word[i],valu[i]);

}
#endif DEBUG

int get_key_value(int argc, char **argv, char ***word, char ***valu)
{
  static char **keyword;
  static char **keyvalu;
  char *i1;
  int i,icnt,i2,i3;

  keyword = (char **)malloc((unsigned) argc*sizeof(char *));
  keyvalu = (char **)malloc((unsigned) argc*sizeof(char *));

  icnt = 0;
  for (i=0; i<argc; i++) {
    i1 = (char *)index(*argv,'=');
    if (i1 != NULL) {
      i2 = (int)( i1-*argv+1 );
      i3 = (int)( (*argv+strlen(*argv))-i1+1 );
      keyword[i] = (char *)malloc((unsigned) i2*sizeof(char));
      keyvalu[i] = (char *)malloc((unsigned) i3*sizeof(char));
      strncpy(keyword[i],*argv,i1-*argv);
      strcpy(keyvalu[i],++i1);
      icnt++;
      argv++;
    }
  }

  *word = keyword;
  *valu = keyvalu;
}

#define MAXREC 1024
#define MAXTOK 100

int get_key_value_from_file(char *file, char ***word, char ***valu)
{
  FILE *fin;
  char buf[MAXREC],*p,*argv[MAXTOK];
  int icnt,iret,argc=0;

  if ( (fin=fopen(file,"r")) == NULL) {
    fprintf(stderr,"Couldn't open parameter file: %s!\n",file);
    exit(-1);
  }

  p = buf;
  while ((icnt=(int)(p-buf)) < MAXREC) {
    if (fgets(p,MAXREC-icnt,fin) == NULL) break;
    argv[argc] = p;
    p += strlen(p);
    if (p[strlen(p)-1] == '\n') p[strlen(p)-1] = '\0';
    argc++;
  }

  return get_key_value(argc,argv,word,valu);
}




