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
 *  exactly of the form of the following regexp:
 *                          key[ ]*=[ ]*word[ ]*[#].*
 *
 *  the results are unpredictable.  E.g. you can have keywords, equals,
 *  value with spaces or TABS in between.  Any line without an '=' is 
 *  ignored.  Anything after a comment character '#' on a line with an '='
 *  is ignored.  Lines beginning with '#' are also ignored (and therefore
 *  thie '#' can be followed by any character).
 *  
 *
 *  By:
 *  --
 *
 *  MDW 11/13/91
 *      06/05/96 added comment field
 *
 ***************************************************************************/

#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>

#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

void strip_blanks(char *);

#ifdef TEST

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
    iret = get_key_value_from_file(file, &word, &valu);
  else {
    iret = get_key_value(argc, argv, &word, &valu);
		/* decrement argument counter by # keys */
    argc -= iret;
  }

  printf("%d pairs found, arguments left: %d\n",iret,argc);
  for (i=0; i<iret; i++) printf("%30s       %30s\n",word[i],valu[i]);

  return 0;
}
#endif /* TEST */

int get_key_value(int argc, char **argv, char ***word, char ***valu)
{
  static char **keyword;
  static char **keyvalu;
  char *i1, *iC;
  int i,icnt,i2;

  keyword = (char **)malloc((unsigned) argc*sizeof(char *));
  if (!keyword) {
    fprintf(stderr, "get_key_value: problem allocating <keyword>\n");
    exit(-1);
  }
  keyvalu = (char **)malloc((unsigned) argc*sizeof(char *));
  if (!keyvalu) {
    fprintf(stderr, "get_key_value: problem allocating <keyvalu>\n");
    exit(-1);
  }

  icnt = 0;
  for (i=0; i<argc; i++) {
				/* Eliminate any comments */
    iC = (char *)index(*argv,'#');
    if (iC != NULL) *iC = '\0';
				/* Look for an assigned value */
    i1 = (char *)index(*argv,'=');

				/* Skip anything without '=' */
    if (i1 != NULL) {
      i2 = (int)( i1-*argv+1 );
      keyword[icnt] = (char *)malloc((unsigned) i2*sizeof(char));
      if (!keyword[icnt]) {
	fprintf(stderr, "get_key_value: problem allocating <keyword[%d]>\n", icnt);
	exit(-1);
      }
      keyword[icnt][i2-1] = '\0';
      strncpy(keyword[icnt],*argv,i1-*argv);

      i2 = (int)( (*argv+strlen(*argv))-i1 );
      if (i2>1) {
	keyvalu[icnt] = (char *)malloc((unsigned) i2*sizeof(char));
	if (!keyvalu[icnt]) {
	  fprintf(stderr, "get_key_value: problem allocating <keyvalu[%d]>\n", icnt);
	  exit(-1);
	}
	strcpy(keyvalu[icnt],++i1);
      }
      else {
				/* No value is assigned zero */
	keyvalu[i] = (char *)malloc((unsigned) 2*sizeof(char));
	if (!keyvalu[i]) {
	  fprintf(stderr, "get_key_value: problem allocating <keyvalu[%d]>\n", i);
	  exit(-1);
	}

	*(keyvalu[i]) = '0';
	*(keyvalu[i]+1) = '\0';
      }
      icnt++;
    }
    argv++;
  }

  *word = keyword;
  *valu = keyvalu;

  return icnt;
}

#define MAXREC 8192
#define MAXTOK 200

int get_key_value_from_file(const char *file, char ***word, char ***valu)
{
  FILE *fin;
  char buf[MAXREC],*p,*argv[MAXTOK];
  int icnt, argc=0;

  if ( (fin=fopen(file,"r")) == NULL) {
    char dirbuf[64];

    fprintf(stderr,"get_key_value: Couldn't open parameter file: %s!\n",file);
    fprintf(stderr,"get_key_value: Working directory: %s\n", getcwd(dirbuf, 64));
    exit(-1);
  }

  p = buf;
  while ((icnt=(int)(p-buf)) < MAXREC) {
    if (fgets(p,MAXREC-icnt,fin) == NULL) break;
    argv[argc] = p;
    strip_blanks(p);
    p += strlen(p);
    if (p[strlen(p)-1] == '\n') p[strlen(p)-1] = '\0';
    argc++;
  }

  return get_key_value(argc, argv, word, valu);
}

/* Strip blanks from a null-terminated string */

void strip_blanks(char *p)
{
  char *t;

  while (*p != '\0') {
    if (*p == ' ' || *p == '\t') {
      t = p;
      do {*t = *(t+1);} while (*++t != '\0');
    }
    else
      ++p;
  }
}

