/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine linearly interpolates on the grid defined by {ftab}[{xtab}]
 *  where {ftab} and {xtab} are arrays of dimension lda.
 *
 *
 *  Call sequence:
 *  -------------
 *  y = odd2(x,xtab,ftab,lda);
 *
 *  double x,xtab[],ftab[];
 *  int lda;
 *
 *  Parameters:
 *  ----------
 *
 *  x        value for desired evaluation
 *  xtab     array containing abcissa
 *  ftab     array containing ordinate
 *  lda      range of array: xtab[1] . . . xtab[lda]
 *  j        index returned
 *
 *  Returns:
 *  -------
 *
 *  interpolated value
 *
 *  Notes:
 *  -----
 *  Uses routine "locate" to do aimple binary search
 *
 *  By:
 *  --
 *
 *  MDW 4/10/89
 *
 ***************************************************************************/

#include <cstring>
#include <cstdlib>

#include <Vector.h>
#include <interp.h>


double odd2(double x, const Vector &xtab, const Vector &ftab, int even)
{
  int index,min,max;
  double ans;


  /*  find position in grid  */

  min = xtab.getlow();
  max = xtab.gethigh();

  if (even)
    index = (int)((x-xtab[min])/(xtab[max]-xtab[min])*(double)(max-min)) + min;
  else
    index = Vlocate(x, xtab);

  if (index < min) index=min;
  if (index >= max) index=max-1;

  ans = ( ftab[index+1]*(x-xtab[index]  ) -
	  ftab[index  ]*(x-xtab[index+1]) )
    /( xtab[index+1]-xtab[index] ) ;

  return ans;

}

double drv2(double x, const Vector &xtab, const Vector &ftab, int even)
{
  int index,min,max;
  double ans;


  /*  find position in grid  */

  min = xtab.getlow();
  max = xtab.gethigh();

  if (even)
    index = (int)((x-xtab[min])/(xtab[max]-xtab[min])*(double)(max-min)) + min;
  else
    index = Vlocate(x, xtab);

  if (index < min) index=min;
  if (index >= max) index=max-1;

  ans = ( ftab[index+1] -ftab[index] )/( xtab[index+1]-xtab[index] ) ;

  return ans;

}

double odd2(double x, const vector<double> &xtab, const vector<double> &ftab, int even)
{
  int index,min,max;
  double ans;


  /*  find position in grid  */

  min = 0;
  max = xtab.size() - 1;

  if (even)
    index = (int)((x-xtab[min])/(xtab[max]-xtab[min])*(double)(max-min)) + min;
  else
    index = Vlocate(x, xtab);

  if (index < min) index=min;
  if (index >= max) index=max-1;

  ans = ( ftab[index+1]*(x-xtab[index]  ) -
	  ftab[index  ]*(x-xtab[index+1]) )
    /( xtab[index+1]-xtab[index] ) ;

  return ans;

}

double drv2(double x, const vector<double> &xtab, const vector<double> &ftab, int even)
{
  int index,min,max;
  double ans;


  /*  find position in grid  */

  min = 0;
  max = xtab.size() - 1;

  if (even)
    index = (int)((x-xtab[min])/(xtab[max]-xtab[min])*(double)(max-min)) + min;
  else
    index = Vlocate(x, xtab);

  if (index < min) index=min;
  if (index >= max) index=max-1;

  ans = ( ftab[index+1] -ftab[index] )/( xtab[index+1]-xtab[index] ) ;

  return ans;

}

