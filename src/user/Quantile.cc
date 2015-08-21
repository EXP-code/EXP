#include <iostream>
#include <iomanip>
#include <sstream>

#include "Quantile.H"

using namespace NTC;

unsigned Quantile::instance = 0;

Quantile::Quantile(const Quantile& q)
{
  dn   = q.dn;
  npos = q.npos;
  hgt  = q.hgt;
  pos  = q.pos;
  p    = q.p;
  num  = q.num;
  full = q.full;

  instance++;
}


double Quantile::P2
(double d,
 double qp1, double q, double qm1, 
 double np1, double n, double nm1)
{
  double outer       = d / (np1 - nm1);
  double inner_left  = (n - nm1 + d) * (qp1 - q ) / (np1 - n);
  double inner_right = (np1 - n - d) * (q - qm1 ) / (n - nm1);
    
  return q + outer * (inner_left + inner_right);
}

void Quantile::update()
{
  for (size_t i=1; i<ssize-1; i++) {
    int    N = pos[i];
    double n = N;
    double q = hgt[i];
	
    double d = npos[i] - n;
	
    if ( (d >=  1.0 and pos[i+1] - n >  1) or 
	 (d <= -1.0 and pos[i-1] - n < -1)  )
      {
	d = floor(copysign(1.0, d));
	int D = static_cast<int>(d);
	
	double qp1 = hgt[i+1];
	double qm1 = hgt[i-1];
	double np1 = pos[i+1];
	double nm1 = pos[i-1];
	double qn  = P2(d, qp1, q, qm1, np1, n, nm1);
	    
	// Accept the parabolic value
	if (qm1 < qn and qn < qp1)
	  hgt[i] = qn;
	else			// Linear rather than parabolic
	  hgt[i] = q + d * (hgt[i+D] - q) / (pos[i+D] - N);
	
	pos[i] = N + D;
      }
  }
}

void Quantile::reset(double P)
{
  // Register p value
  p = P;

  // Initial marker values
  double DN[] = {0.0, 0.5*p, p, 0.5*(1.0 + p), 1.0};
  dn = std::vector<double>(DN, DN + sizeof(DN)/sizeof(double) );
      
  // Initial position values
  double NPOS[] = {1.0, 1.0 + 2.0*p, 1.0 + 4.0*p, 3.0 + 2.0*p, 5.0};
  npos = std::vector<double>(NPOS, NPOS + sizeof(NPOS)/sizeof(double) );

  pos.clear();
  for (size_t i=0; i<ssize; i++) pos.push_back(i+1);
      
  num = 0;

  full = false;
  
  hgt.clear();
}

void Quantile::operator()(double item)
{
  num++;

  // Get the first ssize values
  if (hgt.size() != ssize) {

    hgt.push_back(item);

  } else {

    if (!full) {
      std::sort(hgt.begin(), hgt.end());
      full = true;
    }
    
    // Which interval?
    //
    size_t k = 0;		// Force upper value check
    
    if (item < hgt[0]) {	// Replace lower value
      hgt[0] = item;
      k = 1;
    } else {
      for (size_t i=1; i<ssize; i++) {
	if (hgt[i-1] <= item and item < hgt[i]) {
	  k = i;
	  break;
	}
      }
      
      if (k==0) {
	k = 4;
	if (hgt[k] < item)	// Replace upper value
	  hgt[k] = item;
      }
    }
    
    // Increment all positions greater than k
    //
    for (size_t i=0; i<ssize; i++) {
      if (i >= k) pos[i]++;
      npos[i] += dn[i];
    }
    
    update();
  }
}


double Quantile::operator()()
{
  if (!full) {
    std::sort(hgt.begin(), hgt.end());
    int l = hgt.size();
    // make sure we don't overflow on p == 1 
    // or underflow on p == 0
    return hgt[int(std::min<double>(p * l, std::max<int>(l - 1, 0)))];
  } else {
    return hgt[2];
  }
}


// Node sends its internal data to root
void Quantile::send()
{
  MPI_Send(&dn[0],   ssize, MPI_DOUBLE,        0, 1101, MPI_COMM_WORLD);
  MPI_Send(&npos[0], ssize, MPI_DOUBLE,        0, 1102, MPI_COMM_WORLD);
  int hsz = hgt.size();
  MPI_Send(&hsz,         1, MPI_INT,           0, 1103, MPI_COMM_WORLD);
  MPI_Send(&hgt[0],    hsz, MPI_DOUBLE,        0, 1104, MPI_COMM_WORLD);
  MPI_Send(&pos[0],  ssize, MPI_INT,           0, 1105, MPI_COMM_WORLD);
  MPI_Send(&p,           1, MPI_DOUBLE,        0, 1106, MPI_COMM_WORLD);
  MPI_Send(&num,         1, MPI_UNSIGNED_LONG, 0, 1107, MPI_COMM_WORLD);
  int fl = full ? 1 : 0;
  MPI_Send(&fl,          1, MPI_INT,           0, 1108, MPI_COMM_WORLD);
}
    
// Root intializes itself from node's data
void Quantile::recv(int id)
{
  int ii;
  MPI_Recv(&dn[0],   ssize, MPI_DOUBLE,        id, 1101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(&npos[0], ssize, MPI_DOUBLE,        id, 1102, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(&ii,          1, MPI_INT,           id, 1103, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  hgt.resize(ii);
  MPI_Recv(&hgt[0],     ii, MPI_DOUBLE,        id, 1104, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(&pos[0],  ssize, MPI_DOUBLE,        id, 1105, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(&p,           1, MPI_DOUBLE,        id, 1106, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(&num,         1, MPI_UNSIGNED_LONG, id, 1107, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(&ii,          1, MPI_INT,           id, 1108, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  full = ii ? true : false;
}

Quantile::~Quantile()
{
  instance--;
}
