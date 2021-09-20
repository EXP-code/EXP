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

#include <iostream>
#include <vector>
#include <deque>
#include <Eigen/Eigen>

template <class V>
int Vlocate(double x, const V& xx)
{
  int min = 0;
  int max = xx.size() - 1;
  int jl = min-1;
  int ju = max+1;
  int ascnd = xx[max] > xx[min];

  while (ju-jl > 1) {
    int jm = (ju+jl) >> 1;
    if ((x > xx[jm]) == ascnd)
      jl = jm;
    else
      ju = jm;
  }
  return jl;
}


template <class V>
int Vlocate_with_guard(double value, const V& vec)
{
  int min = 0;
  int max = vec.size() - 1;
  int which = (vec[min] < vec[max]);

  if( ( vec[min] < value == which ) &&
      ( value < vec[max] == which ) ) {
    return Vlocate(value, vec);
  }
  else {
    if ( (value <= vec[min]  ) == which ){
      return min;
    }
    else if( (value >= vec[max]) == which ){
      return max;
    }
    else{
      std::cerr << "WARNING: misordered data in locate_with_guard" << std::endl;
      return -1;
    }
  }
}


template int Vlocate(double x, const std::vector<double>& xtab);

template int Vlocate(double x, const std::deque<double>& xtab);

template int Vlocate(double x, const Eigen::VectorXd& xtab);

template int Vlocate_with_guard(double x, const std::vector<double>& xtab);

template int Vlocate_with_guard(double x, const std::deque<double>& xtab);

template int Vlocate_with_guard(double x, const Eigen::VectorXd& xtab);

