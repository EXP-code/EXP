#ifndef Grid2D_H
#define Grid2D_H

#include <fstream>
#include <vector>

// Define NO_MPI to compile without the MPI libraries
//
#ifdef NO_MPI
int myid=0;
#endif

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
  
  //! Data store (flattened into a single vector)
  std::vector<double> dens;

public:
  //! Constructor
  Grid2D(double Rmax, int num, double t) : R(Rmax), N(num), T(t)
  {
    dr = 2.0*Rmax/N;
    dens.resize(N*N, 0.0);
  }

  //! Constructor from istream (restore from cache file)
  Grid2D(std::istream& in) { read(in); }

  //! Add
  void addPoint(double mass, double x, double y)
  {
    if (x<-R or x>=R or y<-R or y>=R) return;
    int iX = (x + R)/dr;
    int iY = (y + R)/dr;
    dens[iX*N+iY] += mass;
  }

  //! Get density
  double operator()(double x, double y)
  {
    // Bounds check
    //
    if (x<-R or x>=R or y<-R or y>=R) return 0.0;

    // Grid indices
    //
    int iX1 = (x + R)/dr, iX2;
    int iY1 = (y + R)/dr, iY2;
    
    // Sanity check for enclosing bin
    //
    iX1 = std::min<int>(iX1, N-1);
    iY1 = std::min<int>(iY1, N-1);

    // Closest bin center
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

    // The adjacent bin center
    //
    double x2 = -R + (iX2 + 0.5)*dr;
    double y2 = -R + (iY2 + 0.5)*dr;

    // Perform the bilinear interpolation
    //
    double a = (x2 - x)/(x2 - x1), b = (x - x1)/(x2 - x1);
    double c = (y2 - y)/(y2 - y1), d = (y - y1)/(y2 - y1);

    double mass = 
      a*(c*dens[iX1*N+iY1] + d*dens[iX1*N+iY2]) +
      b*(c*dens[iX2*N+iY1] + d*dens[iX2*N+iY2]) ;
    
    // Compute and return the density
    //
    return mass/(dr*dr);
  }

  //! Get time
  double Time() { return T; }

  //! Synchonize arrays
  void sync()
  {
#ifndef NO_MPI
    MPI_Allreduce(MPI_IN_PLACE, dens.data(), N*N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  }

  //! Write to file
  void write(std::ostream& out)
  {
    if (myid==0) {
      out.write((const char *)&N,  sizeof(int));
      out.write((const char *)&R,  sizeof(double));
      out.write((const char *)&T,  sizeof(double));
      out.write((const char *)&dr, sizeof(double));
      out.write((const char *)dens.data(), N*N*sizeof(double));
    }
  }

  //! Read from file
  void read(std::istream& in)
  {
    in.read((char *)&N,  sizeof(int));
    in.read((char *)&R,  sizeof(double));
    in.read((char *)&T,  sizeof(double));
    in.read((char *)&dr, sizeof(double));

    dens.resize(N*N);
    in.read((char *)dens.data(), N*N*sizeof(double));
  }

};

#endif
