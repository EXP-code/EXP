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

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <deque>

#include <Vector.h>

int Vlocate(double x, const Vector& xx)
{
  int ascnd, ju, jm, jl, min, max;
  
  min = xx.getlow();
  max = xx.gethigh();
  jl = min-1;
  ju = max+1;
  ascnd = xx[max] > xx[min];
  while (ju-jl > 1) {
    jm = (ju+jl) >> 1;
    if ((x > xx[jm]) == ascnd)
      jl = jm;
    else
      ju = jm;
  }
  return jl;
}


int Vlocate_with_guard(double value, const Vector& vector)
{
  int which, min, max;

  min = vector.getlow();
  max = vector.gethigh();
  which = (vector[min] < vector[max]);

  if( ( vector[min] < value == which ) &&
      ( value < vector[max] == which ) ) {
    return Vlocate(value, vector);
  }
  else{
    if( (value <= vector[min]  ) == which ){
      return min;
    }
    else if( (value >= vector[max]) == which ){
      return max;
    }
    else{
      puts("WARNING: misordered data in locate_with_guard");
      return -1;
    }
  }
}

template <class V>
int Vlocate(double x, const V& xx)
{
  int ascnd, ju, jm, jl, min, max;
  
  min = 0;
  max = xx.size() - 1;
  jl = min-1;
  ju = max+1;
  ascnd = xx[max] > xx[min];
  while (ju-jl > 1) {
    jm = (ju+jl) >> 1;
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
  int which, min, max;

  min = 0;
  max = vec.size() - 1;
  which = (vec[min] < vec[max]);

  if( ( vec[min] < value == which ) &&
      ( value < vec[max] == which ) ) {
    return Vlocate(value, vec);
  }
  else{
    if( (value <= vec[min]  ) == which ){
      return min;
    }
    else if( (value >= vec[max]) == which ){
      return max;
    }
    else{
      puts("WARNING: misordered data in locate_with_guard");
      return -1;
    }
  }
}


template int Vlocate(double x, const vector<double>& xtab);

template int Vlocate(double x, const deque<double>& xtab);

template int Vlocate_with_guard(double x, const vector<double>& xtab);

template int Vlocate_with_guard(double x, const deque<double>& xtab);


