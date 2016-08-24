#include <exception>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>

#include "Quantile.H"

using namespace NTC;

// Count instances for debugging
unsigned Quantile::instance = 0;

#include "Quantile.H"

Quantile::Quantile()
{
  // Type not yet identified
  t = Box::none;

  // No data to start
  M = 0;

  // No markers to start
  N = 0;

  // For debugging
  instance++;
}

Quantile::Quantile(double p)
{
  quantile(p);

  // For debugging
  instance++;
}

Quantile::Quantile(unsigned n)
{
  histogram(n);

  // For debugging
  instance++;
}

Quantile::Quantile(const Quantile & p)
{
  t  = p.t;
  q  = p.q;
  dn = p.dn;
  np = p.np;
  n  = p.n;
  M  = p.M;
  N  = p.N;

  // For debugging
  instance++;
}

void Quantile::init(void)
{
  // Initialize storage
  q .resize(N, 0);
  dn.resize(N, 0);
  np.resize(N, 0);
  n .resize(N);

  // Set markers
  for (unsigned i=0; i<N; i++) n[i] = i;
}

void Quantile::add(double x)
{
  if (t==Box::quantile)  addBox1(x);
  if (t==Box::histogram) addBox2(x);
  if (t==Box::none) {
    throw std::runtime_error("Quantile: trying to operate with BoxType=none");
  }
}

void Quantile::quantile(double p)
{
  // We are a quantile
  t = Box::quantile;

  // Number of markers
  N = 5;

  // No data to start
  M = 0;

  // Initialize storage and set marker positions
  // (Box 1.A)
  init();

  // Set desired marker positions
  // (Box 1.A)
  np[0] = 0.0;
  np[1] = 2.0*p;
  np[2] = 4.0*p;
  np[3] = 2.0 + 2.0*p;
  np[4] = 4.0;

  // Set marker increments
  // (Box 1.A)
  dn[0] = 0.0;
  dn[1] = 0.5*p;
  dn[2] = p;
  dn[3] = 0.5*(1.0 + p);
  dn[4] = 1.0;
}

void Quantile::histogram(unsigned n)
{
  // We are a histogram
  t = Box::histogram;

  // No data to start
  M = 0;

  // Number of markers
  N = n;

  // Initialize storage
  init();

  // Set marker increments
  for (unsigned i=0; i<N; i++) dn[i] = static_cast<double>(i)/(N-1);
}

double Quantile::parabolic(int i, int d)
{
  return q[i] + 
    d / (n[i+1] - n[i-1]) * 
    ( (n[i] - n[i-1] + d) * (q[i+1] - q[i] ) / (n[i+1] - n[i]) + 
      (n[i+1] - n[i] - d) * (q[i] - q[i-1] ) / (n[i] - n[i-1]) );
}

double Quantile::linear(int i, int d)
{
  return q[i] + d * (q[i+d] - q[i]) / (n[i+d] - n[i]);
}

void Quantile::addBox1(double x)
{
  int k = 0;

  // Algorithm from Box 1 on page 1079 of Jain & Chlamtac
  
  if (M >= N) {
    
    M++;			// Update data count

    // Box 1.B.1
    //
    if(x < q[0]) {
      q[0] = x;			// Left of left end
      k = 0;
    } else if (x >= q[N - 1]) {
      q[N - 1] = x;		// Right of right end
      k = N - 2;
    } else {			// Look for marker
      for(unsigned i=0; i<N-1; i++) {
	if (x < q[i+1]) {
	  k = i;
	  break;
	}
      }
    }

    // Box 1.B.2
    //
    for (unsigned i=k+1; i<N; i++) n[i]++;
    for (unsigned i=0; i<N; i++) np[i] += dn[i];
    
    // Box 1.B.3
    //
    for (unsigned i=1; i<N-1; i++) {
      double d = np[i] - n[i];
      if ( (d >=  1.0 && n[i+1] - n[i] >  1) or
	   (d <= -1.0 && n[i-1] - n[i] < -1)) {

	double tq = parabolic(i, std::copysign(1.0, d));
	if (q[i-1] < tq && tq < q[i+1]) {
	  q[i] = tq;
	} else {
	  q[i] = linear(i, std::copysign(1.0, d));
	}
	n[i] += std::copysign(1.0, d);
      }
    }
    
  } else {

    // Box 1.A
    //

    q[M] = x;

    M++;			// Increment data count

    if (M == N) {
      // We have enough to start the algorithm, sort and initialize
      // position vector
      std::sort(q.begin(), q.end());
      for (unsigned i=0; i<N; i++) n[i] = i + 1;
    }
  }
}


void Quantile::addBox2(double x)
{
  int k = 0;

  // Algorithm from Box 2 on page 1084 of Jain & Chlamtac

  if (M >= N) {
    
    M++;			// Update data count
    
    // Box 2: B.1
    //
    if (x < q[0]) {		
      q[0] = x;			// Left of left end
      k = 0;
    } else if (x > q[N - 1]) {
      q[N - 1] = x;	// Right of right end
      k = N - 2;
    } else if (x >= q[N - 1] and x <= q[N -2]) {
      k = N - 2;
    } else {			// Look for marker
      for (unsigned i=0; i<N-1; i++) {
	if (x < q[i+1]) {
	  k = i;
	  break;
	}
      }
    }

    // Box 2: B.2
    //
    for (unsigned i=k+1; i<N; i++) n[i]++;
    for (unsigned i=0; i<N; i++) np[i] += dn[i];
    
    // Box 2: B.3
    //
    for (unsigned i=1; i<N-1; i++) {
      double d = np[i] - n[i];
      if ( (d >=  1.0 && n[i+1] - n[i] >  1) ||
	   (d <= -1.0 && n[i-1] - n[i] < -1)) {

	double tq = parabolic(i, std::copysign(1.0, d));
	if (q[i-1] < tq && tq < q[i+1]) {
	  q[i] = tq;
	} else {
	  q[i] = linear(i, std::copysign(1.0, d));
	}
	n[i] += std::copysign(1.0, d);
      }
    }

  } else {

    // Box 2.A
    //

    q[M] = x;

    M++;			// Increment data count

    if (M == N) {
      // We have enough to start the algorithm, sort and initialize
      // position vector
      std::sort(q.begin(), q.end());
      for (unsigned i=0; i<N; i++) n[i] = i + 1;
    }
  }
}

double Quantile::operator()()
{
  if (N != 5) {
    throw std::runtime_error("You can only use the () operator for a single quantile");
  }
  return (*this)(dn[(N-1)/2]);
}

double Quantile::operator()(double p)
{
  if (M < N) {
    std::sort(q.begin(), q.end());

    // Use simplest empirical CDF
    //
    unsigned best = 1;
    for (unsigned i=2; i<M; i++) {
      if (fabs(static_cast<double>(i)/M - p) < fabs(static_cast<double>(best)/N - p)) {
	best = i;
      }
    }
    return q[best];

  } else {

    // Find the closest quantile value to p
    //
    int best = 1;
    for (unsigned i=2; i<N-1; i++) {
      if (fabs(dn[i] - p) < fabs(dn[best] - p)) {
	best = i;
      }
    }
    return q[best];
  }
}

// Node sends its internal data to root
void Quantile::send()
{
  unsigned sz = q.size();

  // Send vector size
  MPI_Send(&sz,    1, MPI_UNSIGNED,  0, 1100, MPI_COMM_WORLD);

  // Send data
  MPI_Send(&q[0],  sz, MPI_DOUBLE,        0, 1101, MPI_COMM_WORLD);

  // Send marker values
  MPI_Send(&dn[0], sz, MPI_DOUBLE,        0, 1102, MPI_COMM_WORLD);

  // Send position values
  MPI_Send(&np[0], sz, MPI_DOUBLE,        0, 1103, MPI_COMM_WORLD);

  // Send current indices
  MPI_Send(&n[0],  sz, MPI_DOUBLE,        0, 1104, MPI_COMM_WORLD);

  // Send number of data processed so far
  MPI_Send(&M,      1, MPI_UNSIGNED_LONG, 0, 1105, MPI_COMM_WORLD);

  // Send number of markers
  MPI_Send(&N,      1, MPI_UNSIGNED,      0, 1106, MPI_COMM_WORLD);
}
    
// Root intializes itself from node's data
void Quantile::recv(int id)
{
  unsigned sz;			// Data size

  MPI_Status status;		// MPI return values
  int count;

  // Send vector size
  MPI_Recv(&sz,     1, MPI_UNSIGNED,      id, 1100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  q .resize(sz);
  dn.resize(sz);
  np.resize(sz);
  n .resize(sz);

  // Receive data values
  MPI_Recv(&q[0],  sz, MPI_DOUBLE,        id, 1101, MPI_COMM_WORLD, &status);

  // Sanity
  MPI_Get_count(&status, MPI_DOUBLE, &count);
  if (DBG_VERBOSE && static_cast<unsigned>(count) != sz) {
    int myid; MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    std::cout << "Bad count q [" << myid << "] count=" << count
	      << ", expected=" << sz << std::endl;
  }
  // Receive marker values
  MPI_Recv(&dn[0], sz, MPI_DOUBLE,        id, 1102, MPI_COMM_WORLD, &status);

  // Sanity
  MPI_Get_count(&status, MPI_DOUBLE, &count);
  if (DBG_VERBOSE && static_cast<unsigned>(count) != sz) {
    int myid; MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    std::cout << "Bad count dn [" << myid << "] count=" << count
	      << ", expected=" << sz << std::endl;
  }

  // Receive position values
  MPI_Recv(&np[0], sz, MPI_DOUBLE,        id, 1103, MPI_COMM_WORLD, &status);

  // Sanity
  MPI_Get_count(&status, MPI_DOUBLE, &count);
  if (DBG_VERBOSE && static_cast<unsigned>(count) != sz) {
    int myid; MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    std::cout << "Bad count np [" << myid << "] count=" << count
	      << ", expected=" << sz << std::endl;
  }

  // Receive the rank position of data
  MPI_Recv(&n[0],  sz, MPI_DOUBLE,        id, 1104, MPI_COMM_WORLD, &status);

  // Sanity
  MPI_Get_count(&status, MPI_INT, &count);
  if (DBG_VERBOSE && static_cast<unsigned>(count) != sz) {
    int myid; MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    std::cout << "Bad count n [" << myid << "] count=" << count
	      << ", expected=" << sz << std::endl;
  }

  // Receice the probability value
  MPI_Recv(&M,      1, MPI_UNSIGNED_LONG, id, 1105, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  // Number of markers
  MPI_Recv(&N,      1, MPI_UNSIGNED,      id, 1106, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  if (DBG_VERBOSE && instance % 1000 == 0)
    std::cout << "I=" << instance << std::endl;
}

Quantile::~Quantile()
{
  instance--;			// Count instances for debugging only
}

double Quantile::inverse(double x)
{
  if (M < N) std::sort(q.begin(), q.end());

  // Deal with end points as special cases
  //
  if (x <= q.front()) return 0.0;
  if (x >= q.back ()) return 1.0;

  unsigned K = N; if (M < N) K = M;
  double ret = q[(q.size()-1)/2];

  // Linear interpolation to get x
  //
  for (unsigned i=1; i<K; i++) {
    if (x > q[i-1] and x <= q[i]) {
      if (M<N) 
	ret = ( (q[i] - x)*i + (x - q[i-1])*(i-1) ) / ( (q[i] - q[i-1])*M );
      else     
	ret = ( (q[i] - x)*dn[i] + (x - q[i-1])*dn[i-1] ) / (q[i] - q[i-1]);
    }
  }
  
  return ret;
}

double Quantile::xmin()
{
  if (M < N) std::sort(q.begin(), q.end());
  return q.front();
}

double Quantile::xmax()
{
  if (M < N) std::sort(q.begin(), q.end());
  return q.back();
}

void Quantile::dump(std::ostream& out)
{
  size_t prc = out.precision(3);
  out << std::left
      << std::setw(12) << "q"
      << std::setw(12) << "dn"
      << std::setw(12) << "np"
      << std::setw( 6) << "n"
      << std::endl
      << std::setw(12) << "---"
      << std::setw(12) << "---"
      << std::setw(12) << "---"
      << std::setw( 6) << "---"
      << std::endl;
  for (size_t i=0; i<q.size(); i++)
    out << std::setw(12) << q [i]
	<< std::setw(12) << dn[i]
	<< std::setw(12) << np[i]
	<< std::setw( 6) << n [i]
	<< std::endl;
  out << std::string(42, '-') << std::endl;
  out.precision(prc);
}

