#define DEBUG 1

/*
  spit out phase-space, one particle per line
*/

#include <stdlib.h>

#include "expand.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif


void create_tipsy(void);

void out_tipsy(int n)
{
  char cmd[255];

  create_tipsy();

  if (myid==0) {

    sprintf(cmd, "cat %s %s %s >> %s.bin\0", 
	    "particle.header", "dark.particles", "star.particles",
	    outname);

    system(cmd);

  }

}


void out_chkpt(void)
{
  char cmd[255];

  create_tipsy();

  if (myid==0) {

    sprintf(cmd, "cat %s %s %s > %s.chkpt\0", 
	    "particle.header", "dark.particles", "star.particles",
	    outname);

    system(cmd);

  }
  
}
