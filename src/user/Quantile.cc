#include <stdexcept>
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
  // No data to start
  M = 0;

  // No markers to start
  N = 0;

  // Initialize storage
  init();
}

Quantile::Quantile(double p)
{
  // No data to start
  M = 0;

  // Initialize storage, will set markers
  init();

  // Add a quantile
  newQ(p);
}

Quantile::Quantile(const Quantile & p)
{
  q  = p.q;
  dn = p.dn;
  np = p.np;
  n  = p.n;
  M  = p.M;
  N  = p.N;
}

void Quantile::init(void)
{
  // Set end points
  N = 2;

  // Initialize storage for end points
  q .resize(N);
  dn.resize(N);
  np.resize(N);
  n .resize(N);

  // Add end points
  dn[0] = 0.0;
  dn[1] = 1.0;

  update();
}

double * Quantile::extend(int count)
{
  q .resize(N + count);
  dn.resize(N + count);
  np.resize(N + count);
  n .resize(N + count);

  int last = N;
  N += count;

  return &dn[last];
}

void Quantile::update()
{
  std::sort(dn.begin(), dn.end());

  // Reset markers to uniform distribution
  for (int i=0; i<N; i++)  
    np[i] = dn[i]*(N - 1) + 1;
}

void Quantile::newQ(double p)
{
  double *markers = extend(3);

  // Add new dn markers for this quantile
  markers[0] = p;
  markers[1] = 0.5*p;
  markers[2] = 0.5*(1.0 + p);
  
  update();
}

void Quantile::histogram(int n)
{
  double *markers = extend(n - 1);

  // Add in the new markers for a uniform distribution
  for (int i=1; i<n; i++) 
    markers[i-1] = static_cast<double>(i) / n;
  
  update();
}

double Quantile::parabolic(int i, int d)
{
  return q[i] + 
    d / static_cast<double>(n[i+1] - n[i-1]) * 
    ( (n[i] - n[i-1] + d) * (q[i+1] - q[i] ) / (n[i+1] - n[i]) + 
      (n[i+1] - n[i] - d) * (q[i] - q[i-1] ) / (n[i] - n[i-1]) );
}

double Quantile::linear(int i, int d)
{
  return q[i] + d * (q[i+d] - q[i]) / (n[i+d] - n[i]);
}

void Quantile::add(double x)
{
  int k = 0;

  // Algorithm from Box 2 on page 1084 of Jain & Chlamtac

  if (M >= N) {
    
    M++;			// Update data count

    // Box 2: B.1
    //
    if(x < q[0]) {		
      q[0] = x;			// Left of left end
      k = 1;
    } else if (x >= q[N - 1]) {
      q[N - 1] = x;	// Right of right end
      k = N - 1;
    } else {			// Look for marker
      for(int i=1; i<N; i++) {
	if (x < q[i]) {
	  k = i;
	  break;
	}
      }
    }

    // Box 2: B.2
    //
    for (int i=k; i<N; i++) {
      n[i]++;
      np[i] += dn[i];
    }
    for (int i=0; i<k; i++) np[i] += dn[i];
    
    // Box 2: B.3
    //
    for (int i=1; i<N-1; i++) {
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

    // Box 2: A
    //

    q[M] = x;

    M++;			// Increment data count

    if (M == N) {
      // We have enough to start the algorithm, sort and initialize
      // position vector
      std::sort(q.begin(), q.end());
      for (int i=0; i<N; i++) n[i] = i + 1;
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
    int best = 1;
    for (int i=2; i<M; i++) {
      if (fabs(static_cast<double>(i)/M - p) < fabs(static_cast<double>(best)/N - p)) {
	best = i;
      }
    }
    return q[best];

  } else {

    // Find the closest quantile value to p
    //
    int best = 1;
    for (int i=2; i<N-1; i++) {
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
  int sz = q.size();

  // Send vector size
  MPI_Send(&sz,    1, MPI_INT,     0, 1100, MPI_COMM_WORLD);

  // Send data
  MPI_Send(&q[0],  sz, MPI_DOUBLE, 0, 1101, MPI_COMM_WORLD);

  // Send marker values
  MPI_Send(&dn[0], sz, MPI_DOUBLE, 0, 1102, MPI_COMM_WORLD);

  // Send position values
  MPI_Send(&np[0], sz, MPI_DOUBLE, 0, 1103, MPI_COMM_WORLD);

  // Send current indices
  MPI_Send(&n[0],  sz, MPI_INT,    0, 1104, MPI_COMM_WORLD);

  // Send number of data processed so far
  MPI_Send(&M,      1, MPI_INT,    0, 1105, MPI_COMM_WORLD);

  // Send number of markers
  MPI_Send(&N,      1, MPI_INT,    0, 1106, MPI_COMM_WORLD);
}
    
// Root intializes itself from node's data
void Quantile::recv(int id)
{
  int sz;			// Data size

  MPI_Status status;		// MPI return values
  int count;

  // Send vector size
  MPI_Recv(&sz,     1, MPI_INT,    id, 1100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  q .resize(sz);
  dn.resize(sz);
  np.resize(sz);
  n .resize(sz);

  // Receive data values
  MPI_Recv(&q[0],  sz, MPI_DOUBLE, id, 1101, MPI_COMM_WORLD, &status);

  // Sanity
  MPI_Get_count(&status, MPI_DOUBLE, &count);
  if (DBG_VERBOSE && count != sz) {
    int myid; MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    std::cout << "Bad count q [" << myid << "] count=" << count
	      << ", expected=" << sz << std::endl;
  }
  // Receive marker values
  MPI_Recv(&dn[0], sz, MPI_DOUBLE, id, 1102, MPI_COMM_WORLD, &status);

  // Sanity
  MPI_Get_count(&status, MPI_DOUBLE, &count);
  if (DBG_VERBOSE && count != sz) {
    int myid; MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    std::cout << "Bad count dn [" << myid << "] count=" << count
	      << ", expected=" << sz << std::endl;
  }

  // Receive position values
  MPI_Recv(&np[0], sz, MPI_DOUBLE, id, 1103, MPI_COMM_WORLD, &status);

  // Sanity
  MPI_Get_count(&status, MPI_DOUBLE, &count);
  if (DBG_VERBOSE && count != sz) {
    int myid; MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    std::cout << "Bad count np [" << myid << "] count=" << count
	      << ", expected=" << sz << std::endl;
  }

  // Receive the rank position of data
  MPI_Recv(&n[0],  sz, MPI_INT,    id, 1104, MPI_COMM_WORLD, &status);

  // Sanity
  MPI_Get_count(&status, MPI_INT, &count);
  if (DBG_VERBOSE && count != sz) {
    int myid; MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    std::cout << "Bad count n [" << myid << "] count=" << count
	      << ", expected=" << sz << std::endl;
  }

  // Receice the probability value
  MPI_Recv(&M,      1, MPI_INT,    id, 1105, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  // Number of markers
  MPI_Recv(&N,      1, MPI_INT,     id, 1106, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

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

  int K = std::min(M, N);
  double ret = q[(q.size()-1)/2];

  // Linear interpolation to get x
  //
  for (int i=1; i<K; i++) {
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

