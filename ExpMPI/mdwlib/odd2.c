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

void locate(double *xx, int n, double x, int *j);

double odd2(double x, double *xtab, double *ftab, int lda)
{
  int index;
  double ans;


  /*  find position in grid  */

  locate(xtab,lda,x,&index);

  if (index == 0) index=1;
  if (index == lda) index=lda-1;

  ans = ( ftab[index+1]*(x-xtab[index]  ) -
	  ftab[index  ]*(x-xtab[index+1]) )
    /( xtab[index+1]-xtab[index] ) ;

  return ans;

}
