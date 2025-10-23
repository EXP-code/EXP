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
#include <vector>
#include <deque>

#include "interp.H"

template <class V>
double odd2(double x, const V &xtab, const V &ftab, int even)
{
  // find position in grid
  //
  int min = 0;
  int max = xtab.size() - 1;
  int index;

  if (even)
    index = (int)((x-xtab[min])/(xtab[max]-xtab[min])*(double)(max-min)) + min;
  else
    index = Vlocate(x, xtab);

  if (index < min) index=min;
  if (index >= max) index=max-1;

  double ans = ( ftab[index+1]*(x-xtab[index]  ) -
		 ftab[index  ]*(x-xtab[index+1]) )
    /( xtab[index+1]-xtab[index] ) ;

  return ans;
}

template <class V>
double drv2(double x, const V &xtab, const V &ftab, int even)
{
  // find position in grid
  //
  int min = 0;
  int max = xtab.size() - 1;
  int index;

  if (even)
    index = (int)((x-xtab[min])/(xtab[max]-xtab[min])*(double)(max-min)) + min;
  else
    index = Vlocate(x, xtab);

  if (index <  min) index=min;
  if (index >= max) index=max-1;

  double ans = ( ftab[index+1] -ftab[index] )/( xtab[index+1]-xtab[index] ) ;

  return ans;
}

template double odd2(double x, const std::vector<double>& xtab, 
		     const std::vector<double>& ftab, int even=0);

template double odd2(double x, const std::deque<double>& xtab, 
		     const std::deque<double>& ftab, int even=0);

template double odd2(double x, const Eigen::VectorXd& xtab, 
		     const Eigen::VectorXd& ftab, int even=0);

template double drv2(double x, const std::vector<double>& xtab, 
		     const std::vector<double>& ftab, int even=0);

template double drv2(double x, const std::deque<double>& xtab, 
		     const std::deque<double>& ftab, int even=0);

template double drv2(double x, const Eigen::VectorXd& xtab, 
		     const Eigen::VectorXd& ftab, int even=0);

