/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine computes the double precision log gamma function
 *
 *
 *  Call sequence:
 *  -------------
 *  x = dgammln(x);
 *
 *  double x;
 *
 *  Parameters:
 *  ----------
 *
 *  x        as above
 *
 *  Returns:
 *  -------
 *
 *  Value
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  MDW 11/13/88
 *
 ***************************************************************************/

#include <math.h>

double dgammln(double xx)
{
  double x,tmp,ser;
  static double cof[6]={76.18009173,-86.50532033,24.01409822,
			-1.231739516,0.120858003e-2,-0.536382e-5};
  int j;

  if (xx <= 0.0) return 0.0;
  x=xx-1.0;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.0;
  for (j=0;j<=5;j++) {
    x += 1.0;
    ser += cof[j]/x;
  }
  return -tmp+log(2.50662827465*ser);
}
