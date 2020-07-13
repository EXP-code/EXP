#ifndef Grid2D_H
#define Grid2D_H

#include <vector>

/**
   Create, store, and restore a two-dimensional binned grid
 */
class Grid2D
{
private:

  //@{
  //! Parameters
  double R, T, dr;
  int N;
  //@}
  
  //! Data store
  std::vector< std::vector<double> > dens;

public:
  //! Constructor
  Grid2D(double Rmax, int num, double t) : R(Rmax), N(num), T(t)
  {
    dr = 2.0*Rmax/N;
    dens.resize(N);
    for (auto & a : dens) a.resize(N, 0.0);
  }

  //! Constructor from istream (restore from cache file)
  void Grid2d(std::istream& in)
  {
    read(in);
  }

  //! Add
  void addPoint(double mass, double x, double y)
  {
    if (x<-R or x>=R or y<-R or y>=R) return;
    int iX = (x + R)/dr;
    int iY = (y + R)/dr;
    dens[iX][iY] += mass;
  }

  //! Get density
  double operator()(double x, double y)
  {
    if (x<-R or x>=R or y<-R or y>=R) return 0.0;

				// Grid indices
    int iX1 = (x + R)/dr, iX2;
    int iY1 = (y + R)/dr, iY2;
    
    // Sanity check
    //
    iX1 = std::min<int>(iX1, N-1);
    iY1 = std::min<int>(iY1, N-1);

    // Closed bin center
    //
    double x1 = -R + (iX1 + 0.5)*dr;
    double y1 = -R + (iY1 + 0.5)*dr;

    if (x<x1) {
      if (iX1>0)   iX2 = iX1 - 1;
      else         iX2 = iX1 + 1;
    } else {
      if (iX1<N-1) iX2 = iX1 + 1;
      else         iX2 = iX1 - 1;
    }

    if (y<y1) {
      if (iY1>0)   iY2 = iY1 - 1;
      else         iY2 = iY1 + 1;
    } else {
      if (iY1<N-1) iY2 = iY1 + 1;
      else         iY2 = iY1 - 1;
    }

    // Next bin center
    //
    double x2 = -R + (iX2 + 0.5)*dr;
    double y2 = -R + (iY2 + 0.5)*dr;

    // Perform the bilinear interpolation
    //
    double a = (x - x1)/(x2 - x1), b = (x - x2)/(x1 - x2);
    double c = (y - x1)/(y2 - y1), d = (y - y2)/(y1 - y2);

    double mass = 
      a*c*dens[iX1][iY1] +
      a*d*dens[iX1][iY2] +
      b*c*dens[iX2][iY1] +
      b*d*dens[iX2][iY2] ;
    
    return mass/(dr*dr);
  }

  //! Get time
  double Time() { return T; }

  //! Synchonize arrays
  void sync()
  {
    for (auto & a : dens) {
      MPI_Allreduce(MPI_IN_PLACE, a.data(), N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
  }

  //! Write to file
  void write(std::ostream& out)
  {
    if (myid==0) {
      out.write((const char *)&N,  sizeof(int));
      out.write((const char *)&R,  sizeof(double));
      out.write((const char *)&T,  sizeof(double));
      out.write((const char *)&dr, sizeof(double));
      for (auto & a : dens)
	out.write((const char *)a.data(), N*sizeof(double));
    }
  }

  //! Read from file
  void read(std::istream& in)
  {
    in.read((char *)&N,  sizeof(int));
    in.read((char *)&R,  sizeof(double));
    in.read((char *)&T,  sizeof(double));
    in.read((char *)&dr, sizeof(double));

    dens.resize(N);
    for (auto & a : dens) {
      a.resize(N);
      in.read((char *)a.data(), N*sizeof(double));
    }
  }

};

#endif
