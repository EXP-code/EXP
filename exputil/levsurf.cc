/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine finds the level set of a function in a supplied matrix using
 *  a bilinear interpolation algorithm.  You need one call for each contour
 *  level.  A set of contiguous points representing the contour is returned.
 *
 *
 *  Call sequence:
 *  -------------
 *  nseg = level_surface(Matrix& z, double val, Vector& xret, Vector& yret,
 *                       int &num);
 *
 *  Parameters:
 *  ----------
 *
 *  z        matrix (z[m][n]) containing the function over the domain
 *  val      value of desired level set
 *  xret     x values of returned level set scaled to unit square
 *  yret     y values of returned level set scaled to unit square
 *  num      number of values in level set
 *
 *  Returns:
 *  -------
 *
 *  nseg = "number of distinct segments. Level set is assigned to passed
 *  pointers.  IMPORTANT!!!: note that values are scaled to the unit
 *  square.
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  MDW 		4/12/89 
 *  bug fix 		4/16/89
 *  converted to C++	10/25/98
 *
 ***************************************************************************/


#include <cstdlib>
#include <cmath>

#include <iostream>

#include <Eigen/Eigen>

/******* Function definitions *******/

int level_surface(Eigen::MatrixXd& z, double val, Eigen::VectorXd& xr, Eigen::VectorXd& yr, int& num);
void level_surface_sort_swap(Eigen::MatrixXd& surf,int i1, int i2);
void level_surface_exch(Eigen::MatrixXd& surf, int i);
void level_surface_add_point(Eigen::MatrixXd& surf,
			     double x0, double y0, double x1, double y1,
			     int& k);
void level_surface_dump_raw(Eigen::MatrixXd& surf, int k);
void level_surface_dump_cooked(Eigen::MatrixXd& surf, int k);
int level_surface_sort_find_min(Eigen::MatrixXd& surf, int begin, int end);
int level_surface_sort_loc(Eigen::MatrixXd& surf, int begin, int end);

#define TOL 1.0e-8


#ifdef DEBUG_MAIN
#define ALT 2

int
main(int argc, char** argv)
{
  char          file[255];
  double        x,y,dx,dy,level,exact;
  int           m,n,i,j,num,nseg;
  Matrix        z;
  Eigen::VectorXd        xret, yret;


  cout << "Output filename: ";
  cin >> file;

  ofstream out(file);

  if (!out) {
    cerr.form("Couldn't open %s . . . quitting\n",file);
    exit(-1);
  }

  cout << "num in grid, level: ";
  cin >> num;
  cin >> level;

  dx = dy = 1.0/(double)num;

  n = m = 2*num+1;
  z.setsize(1, m, 1, n);

  x = -1.0;
  for (i=1; i<=m; i++) {
    y = -1.0;
    for (j=1; j<=n; j++) {
      switch (ALT) {
      case 1:
	z[i][j] = x*y + y*y;
	break;
      case 2:
	z[i][j] = x + 2.0*y*y;
	break;
      default:
	z[i][j] = 2.0*x*x*x*y;
      }
      y += dy;
    }
    x += dx;
  }

  nseg = level_surface(z, level, xret, yret, num);
  cerr.form("Number of contiguous segments: %d\n",nseg);

  for (i=1; i<=num; i++) {
    x = -1.0 + 2.0*xret[i];
    y = -1.0 + 2.0*yret[i];

    switch (ALT) {
    case 2:
      {
	double tmp1 = -sqrt(0.5*(level - x));
	double tmp2 =  sqrt(0.5*(level - x));
	if (fabs(tmp1 - y) < fabs(tmp2 - y))
	  exact = tmp1;
	else
	  exact = tmp2;
      }
      break;
    case 1:
      {
	double tmp1 = -0.5*x + sqrt(0.25*x*x + level);
	double tmp2 = -0.5*x - sqrt(0.25*x*x + level);
	if (fabs(tmp1 - y) < fabs(tmp2 - y))
	  exact = tmp1;
	else
	  exact = tmp2;
      }
      break;
    default:
      exact = 0.5*level/(x*x*x);
    }
    out.form("  %le  %le  %le  %le\n",x,y,y-exact,(y-exact)/exact);
  }

}

#endif // DEBUG_MAIN

static int DIM=0;
static Eigen::MatrixXd surf;

int level_surface(Eigen::MatrixXd& z, double val,
		  Eigen::VectorXd& xret, Eigen::VectorXd& yret, int& num)
{
  double deltx,delty,x,y,v1,v2,v3,v4,x0,x1,y0,y1;
  int i,j,k,icase,next,nseg;

  int m = z.rows();
  int n = z.cols();
				/* Allocate matrix to contain segment list */
				/* Can't be more than 2(n-1)(m-1) long */
  DIM = 2*(n-1)*(m-1);
  surf.resize(DIM, 4);

  deltx = 1.0/(double)(m-1);
  delty = 1.0/(double)(n-1);
  y = 0.0;
  k = 0;
  for (int j=0; j<n; j++) {
    x = 0.0;
    for (int i=0; i<m; i++) {
      v1 = z(i,   j);
      v2 = z(i+1, j);
      v3 = z(i,   j+1);
      v4 = z(i+1, j+1);

      icase = 1;
      if (val>v1) icase = icase + 1;
      if (val>v2) icase = icase + 2;
      if (val>v3) icase = icase + 4;
      if (val>v4) icase = 9 - icase;

      switch (icase) {
      case 1:
	break;
      case 2:
	x0 = x + deltx*(val-v1)/(v2-v1);
	y0 = y;
	x1 = x;
	y1 = y + delty*(val-v1)/(v3-v1);
	level_surface_add_point(surf, x0, y0, x1, y1, k);
	break;
      case 3:
	x0 = x + deltx*(val-v1)/(v2-v1);
	y0 = y;
	x1 = x + deltx;
	y1 = y + delty*(val-v2)/(v4-v2);
	level_surface_add_point(surf, x0, y0, x1, y1, k);
	break;
      case 4:
	x0 = x;
	y0 = y + delty*(val-v1)/(v3-v1);
	x1 = x + deltx;
	y1 = y + delty*(val-v2)/(v4-v2);
	level_surface_add_point(surf, x0, y0, x1, y1, k);
	break;
      case 5:
	x0 = x;
	y0 = y + delty*(val-v1)/(v3-v1);
	x1 = x + deltx*(val-v3)/(v4-v3);
	y1 = y + delty;
	level_surface_add_point(surf, x0, y0, x1, y1, k);
	break;
      case 6:
	x0 = x + deltx*(val-v1)/(v2-v1);
	y0 = y;
	x1 = x + deltx*(val-v3)/(v4-v3);
	y1 = y + delty;
	level_surface_add_point(surf, x0, y0, x1, y1, k);
	break;
      case 7:
	x0 = x + deltx*(val-v1)/(v2-v1);
	y0 = y;
	x1 = x;
	y1 = y + delty*(val-v1)/(v3-v1);
	level_surface_add_point(surf, x0, y0, x1, y1, k);
      case 8:
	x0 = x + deltx*(val-v3)/(v4-v3);
	y0 = y + delty;
	x1 = x + deltx;
	y1 = y + delty*(val-v2)/(v4-v2);
	level_surface_add_point(surf, x0, y0, x1, y1, k);
      }
      x = x + deltx;
    }
    y = y + delty;
  }

#ifdef DEBUG
  level_surface_dump_raw(surf, k);
#endif

/* Sort vectors to produce contour */

  nseg = 0;
  for (i=1, next=0; i<k; i++) {
    if (next == 0) {
      nseg++;
      next = level_surface_sort_find_min(surf,i,k);
    }
    else
      level_surface_exch(surf,next);
    level_surface_sort_swap(surf,i,next);
    next = level_surface_sort_loc(surf,i,k);
  }
  if (next != 0)
    level_surface_exch(surf,next);

  if (k>0) {			// Did we find any points?
    xret.resize(k);
    yret.resize(k);
    num = k+1;
    for (int i=0; i<k-1; i++) {
      xret[i] = surf(i, 1);
      yret[i] = surf(i, 2);
    }
    xret[k-1] = surf(k-1, 3);
    yret[k-1] = surf(k-1, 4);
    
#ifdef DEBUG
  level_surface_dump_cooked(surf, k);
#endif
  }
  else
    num = 0;			// None found

  return nseg;

}

/*
   Add segment to segment list
*/

void level_surface_add_point(Eigen::MatrixXd& surf,
			     double x0, double y0, double x1, double y1,
			     int& k)
{
  
  if (fabs(x0-x1) < 0.9*TOL*fabs(x0) &&	     /* Don't add if zero length  */
      fabs(y0-y1) < 0.9*TOL*fabs(y0) ) {
    return;
  }

  k++;				             /* Add */
  surf(k, 1) = x0;	surf(k, 2) = y0;
  surf(k, 3) = x1;	surf(k, 4) = y1;

  }


/*
   Find segment with minimum x value; returns index of minimum segment;
   positive values means tail has minimum value, negative values means head
   has minimum value
*/

int level_surface_sort_find_min(Eigen::MatrixXd& surf, int begin, int end)
{
  double ymin;
  int indx;
				/* Find minimum "head-tail".  Begin with   */
				/* minimum y and all else equal, minimum x */
  ymin = 1.0e30;		
  indx = 0;

  for (int i=begin; i<=end; i++) {
    if (surf(i, 2) <= ymin) {
      if ( surf(i, 2) == ymin ) {
	if ( surf(i, 1) < surf(abs(indx), 1) )
	  indx = i;
      }
      else {
	indx = i;
	ymin = surf(i, 2);
      }
    }
    else if (surf(i, 4) <= ymin) {
      if ( surf(i, 4) == ymin ) {
	if ( surf(i, 3) < surf(abs(indx), 3) )
	  indx = -i;
      }
      else {
	indx = -i;
	ymin = surf(i, 4);
      }
    }
  }

				// Is it on the x edge?  Exchange . . .
  if ( fabs(surf(abs(indx), 3) - 1.0) < 1.0e-8 ) {
    double xtmp,ytmp;

    xtmp = surf(abs(indx), 1);
    ytmp = surf(abs(indx), 2);
    surf(abs(indx), 1) = surf(abs(indx), 3);
    surf(abs(indx), 2) = surf(abs(indx), 4);
    surf(abs(indx), 3) = xtmp;
    surf(abs(indx), 4) = ytmp;
  }

  return indx;
}


/* 
   Search segments from index=begin through index=end to find adjacent
   grid block containing level surface.  

   Returns segment index if found; positive values means tail matches head,
   negative values mean tail matches tail, 0 if no segment found.
*/

int level_surface_sort_loc(Eigen::MatrixXd& surf, int begin, int end)
{
  double x,y;
  int i;
				/* Find next fit tail to head */
  x = surf(begin, 3);
  y = surf(begin, 4);

  for (i=begin+1; i<=end; i++) {
    if (fabs(surf(i, 1)-x) < TOL*fabs(x) && 
	fabs(surf(i, 2)-y) < TOL*fabs(y)) return i;
    if (fabs(surf(i, 3)-x) < TOL*fabs(x) && 
	fabs(surf(i, 4)-y) < TOL*fabs(y)) return -i;
  }

  return 0;
}



/* Swap segment i1 with segment i2 */

void level_surface_sort_swap(Eigen::MatrixXd& surf,int i1, int i2)
{
  double temp[4];

  for (int i=0; i<4; i++)
    temp[i] = surf(abs(i1), i+1);

  for (int i=0; i<4; i++)
    surf(abs(i1), i+1) = surf(abs(i2), i+1);

  for (int i=0; i<4; i++)
    surf(abs(i2), i+1) = temp[i];

}



/* Invert direction of line segment at index i */
    
void level_surface_exch(Eigen::MatrixXd& surf, int i)
{
  double xtmp,ytmp;

  if (i<0) {
    xtmp = surf(abs(i), 1);
    ytmp = surf(abs(i), 2);
    surf(abs(i), 1) = surf(abs(i), 3);
    surf(abs(i), 2) = surf(abs(i), 4);
    surf(abs(i), 3) = xtmp;
    surf(abs(i), 4) = ytmp;
  }
}


/* Debug dump routines */
#ifdef DEBUG

void level_surface_dump_raw(Eigen::MatrixXd& surf, int k)
{
  int i;

  ofstream out("raw.dat");

  for (i=1; i<=k; i++)
    out << surf(i, 1] << " " << surf(i, 2] << " " 
	<< surf(i, 3] << " " << surf(i, 4] << endl;

}


void level_surface_dump_cooked(Eigen::MatrixXd& surf, int k)
{
  int i;

  ofstream out("cooked.dat");

  for (i=1; i<=k; i++)
    out << surf(i, 1] << " " << surf(i, 2] << " " 
	<< surf(i, 3] << " " << surf(i, 4] << endl;

}

#endif
