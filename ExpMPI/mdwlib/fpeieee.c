/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine sets the fpmode to enable floating point exception trapping
 *  and provides an error handler with messages for the usual exceptions.
 *  This was designed to be used with dbx (put breaks in the fpetrap routine).
 *
 *
 *  Call sequence:
 *  -------------
 *  void fpeinit();       * call once at beginning to 
 *                                initialize error handler *
 *
 *  void fpereset();      * use to reset exception handler to default *
 *
 *  Parameters:
 *  ----------
 *
 *  none
 *
 *  Returns:
 *  -------
 *
 *  none
 *
 *  Notes:
 *  -----
 *
 *  See <signal.h> and SUN floating point ref. manual.
 *
 *
 *  By:
 *  --
 *
 *  MDW 12/01/88
 *
 ***************************************************************************/

#include <stdio.h>
void fpeinit(int iflg) {}

/*
#include <sys/signal.h>
#include <floatingpoint.h>

void fpetrap(int sig, int code, struct sigcontext *scp, char *addr);

void fpeinit(int iflg)
{
  sigfpe_handler_type hld;
  
  hld = (sigfpe_handler_type)fpetrap;

  if (iflg) {
    if (ieee_handler("set","all",hld) != 0)
      printf("ieee_handler can't set all!\n");
  }
  else {
    if (ieee_handler("set","common",hld) != 0)
      printf("ieee_handler can't set common!\n");
  }

}


void fpetrap(int sig, int code, struct sigcontext *scp, char *addr)
{
  switch (code) {
  case FPE_FLTINEX_TRAP:
    printf("===> Inexact result!\n");
    break;
  case FPE_FLTOVF_TRAP:
    printf("===> Overflow!\n");
    break;
  case FPE_FLTOPERR_TRAP:
    printf("===> Operand error!\n");
    break;
  case FPE_INTDIV_TRAP:
    printf("===> Integer divide by zero!\n");
    break;
  case FPE_FLTDIV_TRAP:
    printf("===> Floating ivide by zero!\n");
    break;
  case FPE_FLTUND_TRAP:
    printf("===> Floating underflow!\n");
    break;
  case FPE_INTOVF_TRAP:
    printf("===> Integer underflow!\n");
    break;
  case FPE_STARTSIG_TRAP:
    printf("===> Process using fp!\n");
    break;
  default:
    printf("===> Exception [code=%x]!\n",code);
    break;
  }
#ifdef ON_ABORT
  abort();
#endif
}

*/
