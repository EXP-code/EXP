/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine finds the index in provided array j such that
 *  xx[j] < xx < xx[j+1].
 *
 *
 *  Call sequence:
 *  -------------
 *  void locate(xx,n,x,j);
 *  void locate_with_guard(xx,n,x,j);
 *
 *  double xx[],x;
 *  int n,*j;
 *
 *  Parameters:
 *  ----------
 *
 *  xx       array
 *  x        value to be located
 *  n        range of array: xx[1] . . . xx[n]
 *  j        index returned
 *
 *  Returns:
 *  -------
 *
 *  sets j = 0 if x<xx[1] and j=n+1 if x > xx[n];
 *
 *  Notes:
 *  -----
 *  Simple binary search.  Locate_with_guard() always returns a grid
 *  point within the specfied range
 *
 *  By:
 *  --
 *
 *  MDW 11/13/88
 *  DFC  5/18/89
 *  MDW 01/24/90
 *
 ***************************************************************************/

#include "cutil.h"

void locate(double *xx, int n, double x, int *j)
{
  int ascnd,ju,jm,jl;
  
  jl=0;
  ju=n+1;
  ascnd=xx[n] > xx[1];
  while (ju-jl > 1) {
    jm=(ju+jl) >> 1;
    if ((x > xx[jm]) == ascnd)
      jl=jm;
    else
      ju=jm;
  }
  *j=jl;
}


void locate_with_guard(double *vector, int num, double value, int *i)
{
  int which;

  which = (vector[1] < vector[num]);

  if( ( (vector[1] < value)   == which ) &&
      ( (value < vector[num]) == which ) ) {
    locate(vector, num, value, i);
  }
  else{
    if( (value <= vector[1]  ) == which ){
      *i = 1;
    }
    else if( (value >= vector[num]) == which ){
      *i = num;
    }
    else{
      myerror("locate_with_guard: misorder data");
    }
  }
  
  return;
}


