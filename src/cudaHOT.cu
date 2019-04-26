#include <cuda.h>
#include <cuda_runtime.h>

#include <thrust/binary_search.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>

#include <cudaUtil.cuH>

#include <pHOT.H>

// Debug keys across all nodes (uses MPI calls)
static bool DEBUG_KEYS = false;


__device__
inline uint64_t split3( unsigned int a )
{
  // we only use the first 21 bits
  uint64_t x = a & 0x1fffff;
  x = (x | x << 32) & 0x001f00000000ffff; // shift left 32 bits, OR with self, and 00011111000000000000000000000000000000001111111111111111
  x = (x | x << 16) & 0x001f0000ff0000ff; // shift left 32 bits, OR with self, and 00011111000000000000000011111111000000000000000011111111
  x = (x | x << 8 ) & 0x100f00f00f00f00f; // shift left 32 bits, OR with self, and 0001000000001111000000001111000000001111000000001111000000000000
  x = (x | x << 4 ) & 0x10c30c30c30c30c3; // shift left 32 bits, OR with self, and 0001000011000011000011000011000011000011000011000011000100000000
  x = (x | x << 2 ) & 0x1249249249249249;
  return x;
}

__device__
inline uint64_t mortonEncode_mask( unsigned int* u )
{
  uint64_t answer = 0;
  answer |= split3(u[0]) | split3(u[1]) << 1 | split3(u[2]) << 2;
  return answer;
}

__device__
inline uint64_t genKey( double x, double y, double z, int nbits )
{
  const unsigned int maxI = 0x1fffff;
  unsigned u[3];
  u[2] = x * maxI;
  u[1] = y * maxI;
  u[0] = z * maxI;
  
  uint64_t nkey = mortonEncode_mask(&u[0]);

  uint64_t place(1u);
  place <<= 3*nbits;
  unsigned shift = 63 - 3*nbits;
  nkey >>= shift;
  nkey += place;

  return nkey;
}


__global__ void computeKeysKernel(dArray<double>   x,
				  dArray<double>   y,
				  dArray<double>   z,
				  dArray<uint64_t> k,
				  int nbits, int stride)
{
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;
  const int N   = k._s;

  for (int s=0; s<stride; s++) {

    int n = tid*stride + s;	// Current interaction pair

    if (n < N) {
      k._v[n] = genKey(x._v[n], y._v[n], z._v[n], nbits);
    }
  }
}


void pHOT::keyProcessCuda(std::vector<key_wght> keywght)
{
				// Cuda propertices
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, cc->cudaDevice);

				// For diagnostics
  Timer *timer_debug;
  if (DEBUG_KEYS && myid==0) {
    timer_debug = new Timer();
    timer_debug->start();
  }

  timer_keygenr.start();

  size_t pN = cc->Particles().size();

  keywght.resize(pN);		// Return key-weight list

  std::vector<double> x_h, y_h, z_h;
  std::vector<bool> mask(pN, true);

  int I=0;
  for (auto & v : cc->Particles()) {
    
    double X = (v.second->pos[0] + offset[0])/sides[0];
    double Y = (v.second->pos[1] + offset[1])/sides[1];
    double Z = (v.second->pos[2] + offset[2])/sides[2];

    // Compute keys for in-bounds particles
    if (X>0.0 and Y>0.0 and Z>0.0 and X<1.0 and Y<1.0 and Z<1.0) {
      x_h.push_back(X);
      y_h.push_back(Y);
      z_h.push_back(Z);
    } else {
      mask[I] = false;
    }
    I++;
  }

  size_t N = x_h.size();
  thrust::device_vector<double>   x_d = x_h;
  thrust::device_vector<double>   y_d = y_h;
  thrust::device_vector<double>   z_d = z_h;
  thrust::device_vector<uint64_t> k_d(N);

  int stride   = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
  int gridSize = (N+BLOCK_SIZE*stride-1)/(BLOCK_SIZE*stride);

  computeKeysKernel<<<gridSize, BLOCK_SIZE>>>
    (toKernel(x_d), toKernel(y_d), toKernel(z_d), toKernel(k_d), 
     nbits, stride);

  thrust::host_vector<uint64_t> k_h(k_d);
  
  I = 0;
  oob.clear();
  for (auto & v : cc->Particles()) {
    if (mask[I]) {
      v.second->key = k_h[I];
    } else {
      v.second->key = 0ull;
      oob.insert(v.first);
    }
    keywght[I] = key_wght(k_h[I], 1.0);
    I++;
  }

  timer_keygenr.stop();
  timer_keysort.start();

  spreadOOB();

				// Sort the keys on the device
  thrust::sort(k_d.begin(), k_d.end());
  k_h = k_d;			// Now, get the sorted keys

  std::vector<uint64_t> keys1(k_h.size()), keys;
  std::copy(k_h.begin(), k_h.end(), keys1.begin());

  parallelMergeOne(keys1, keys);

  vector<key_wght> keylist(keys.size());

  for (unsigned i=0; i<keylist.size(); i++) {
    keylist[i].first  = keys[i];
    keylist[i].second = 1.0;
  }
				// Accumulate weight list
  for (unsigned i=1; i<keylist.size(); i++) 
    keylist[i].second += keylist[i-1].second;
    
				// Normalize weight list
  double ktop = keylist.back().second;
  for (unsigned i=0; i<keylist.size(); i++)  keylist[i].second /= ktop;

  vector<double> frate(numprocs);

				// Use an even rate
  frate[0] = 1.0;
  for (unsigned i=1; i<numprocs; i++) 
    frate[i] = frate[i-1] + 1.0;
  //                        ^
  //                        |
  //                        |
  // Replace with the computation rate for the node <i>
  // for load balancing
  //
    
  struct wghtDBL
  {
    bool operator()(const std::pair<key_type, double>& a, 
		    const std::pair<key_type, double>& b)
    { 
      return (a.second < b.second); 
    }
  };

  
  if (myid==0) {

    // The overhead for computing these is small so
    // no matter if they are not used below

    vector<double>   wbeg(numprocs), wfin(numprocs); // Weights for debugging
    vector<unsigned> pbeg(numprocs), pfin(numprocs); // Counts  for debugging
    
				// Compute the key boundaries in the partition
				//
    for (unsigned i=0; i<numprocs-1; i++) {
      if (keylist.size()) {
	key_wght k(0u, frate[i]/frate[numprocs-1]);
	vector<key_wght>::iterator
	  ret = lower_bound(keylist.begin(), keylist.end(), k, wghtDBL());
	kfin[i] = ret->first;
	wfin[i] = ret->second;
	pfin[i] = ret - keylist.begin();
      }
      else {
	kfin[i] = key_min;
	wfin[i] = 0.0;
	pfin[i] = 0;
      }
    }

    kfin[numprocs-1] = key_max;
    wfin[numprocs-1] = 1.0;
    pfin[numprocs-1] = keylist.size();

    kbeg[0] = key_min;
    wbeg[0] = 0.0;
    pbeg[0] = 0;
    for (unsigned i=1; i<numprocs; i++) {
      kbeg[i] = kfin[i-1];
      wbeg[i] = wfin[i-1];
      pbeg[i] = pfin[i-1];
    }
      
    if (DEBUG_KEYS) {	 // If true, print key ranges for each process
      ofstream out(debugf.c_str(), ios::app);
      unsigned nhead = 5 + 3*15 + 3*10 + 3*klen;
      out << setw(nhead) << setfill('-') << '-' << endl
	  << "---- partitionKeys: keys in list="<< keylist.size() << endl
	  << setw(nhead) << setfill('-') << '-' << endl << setfill(' ')
	  << left << setw(5) << "proc" << right << setw(klen) << "kbeg"
	  << setw(klen) << "kfin" << setw(klen) << "# keys" 
	  << setw(15) << "wbeg" << setw(15) << "wend" << setw(15) << "wdif" 
	  << setw(10) << "pbeg" << setw(10) << "pend" << setw(10) << "pdif" 
	  << endl
	  << left << setw(5) << "----" << right << setw(klen) << "----"
	  << setw(klen) << "----" << setw(klen) << "------" 
	  << setw(15) << "----" << setw(15) << "----" << setw(15) << "----"  
	  << setw(10) << "----" << setw(10) << "----" << setw(10) << "----"  
	  << endl;
      double mdif=0.0, mdif2=0.0;
      for (int i=0; i<numprocs; i++) {
	out << left << setw(5) << i << right
	    << hex << setw(klen) << kbeg[i]
	    << setw(klen) << kfin[i] << dec
	    << setw(klen) << (kfin[i] - kbeg[i])
	    << setw(15) << wbeg[i]
	    << setw(15) << wfin[i]
	    << setw(15) << wfin[i] - wbeg[i]
	    << setw(10) << pbeg[i]
	    << setw(10) << pfin[i]
	    << setw(10) << pfin[i] - pbeg[i]
	    << endl;
	mdif  += wfin[i] - wbeg[i];
	mdif2 += (wfin[i] - wbeg[i])*(wfin[i] - wbeg[i]);
      }
      out << setw(nhead) << setfill('-') << '-' << endl << setfill(' ') 
	  << endl;
      if (numprocs>1) {
	mdif /= numprocs;
	mdif2 = sqrt((mdif2 - mdif*mdif*numprocs)/(numprocs-1));
	out << "----  mean wght = " << setw(10) << keylist.size() << endl
	    << "----  std. dev. = " << setw(10) << keylist.size() << endl;
      }
    }
  }

  MPI_Bcast(&kbeg[0], numprocs, MPI_EXP_KEYTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&kfin[0], numprocs, MPI_EXP_KEYTYPE, 0, MPI_COMM_WORLD);


  if (DEBUG_KEYS) {		// If true, print key totals
    unsigned oobn = oobNumber();
    unsigned tkey1 = keys.size(), tkey0 = 0;
    MPI_Reduce(&tkey1, &tkey0, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
    if (myid==0) {
      ofstream out(debugf.c_str(), ios::app);
      out << endl
	  << setfill('-') << setw(60) << '-' << setfill(' ') << endl
	  << "---- partitionKeys" << endl
	  << setfill('-') << setw(60) << '-' << setfill(' ') << endl
	  << "----  list size = " << setw(10) << keylist.size() << endl
	  << "---- total keys = " << setw(10) << tkey0 << endl
	  << "----  total oob = " << setw(10) << oobn << endl
	  << "----      TOTAL = " << setw(10) << tkey0 + oobn << endl
	  << setfill('-') << setw(60) << '-' << setfill(' ') << endl
	  << "----   Time (s) = " << setw(10) << timer_debug->stop()
	  << endl
	  << setfill('-') << setw(60) << '-' << setfill(' ') << endl
	  << endl;
      delete timer_debug;
    }
  }

  timer_keysort.stop();

  if (DEBUG_KEYS) {		// Deep debug output
    static unsigned count = 0;
    std::ostringstream sout;
    sout << debugf << "." << myid << "." << count++;
    ofstream out(sout.str());

    for (unsigned i=0; i<keylist.size(); i++) {
      out << std::setw( 5) << i
	  << std::setw(18) << hex << keylist[i].first
	  << std::setw(18) << dec << keylist[i].second
	  << std::endl;
      if (i>0 and keylist[i].first < keylist[i-1].first)
	out << "####" << std::endl;
    }
  } // END: DEBUG_KEYS

}

//
// This routine combines two sorted vectors into one
// larger sorted vector
//
void pHOT::sortCombineOne(vector<key_type>& one, vector<key_type>& two,
			  vector<key_type>& comb)
{
  int i=0, j=0;
  int n = one.size()-1;
  int m = two.size()-1;
  
  comb = vector<key_type>(one.size()+two.size());

  for (int k=0; k<n+m+2; k++) {
    if (i > n)
      comb[k] = two[j++];
    else if(j > m)
      comb[k] = one[i++];
    else {
      if(one[i] < two[j])
	comb[k] = one[i++];
      else
	comb[k] = two[j++];
    }
  }
}

void pHOT::parallelMergeOne(std::vector<key_type>& initl,
			    std::vector<key_type>& final)
{
  MPI_Status status;
  std::vector<key_type> work;
  unsigned n;

  // Find the largest power of two smaller than
  // the number of processors
  // 
  int M2 = 1;
  while (M2*2 < numprocs) M2 = M2*2;

  // Combine the particles of the high nodes
  // with those of the lower nodes so that
  // all particles are within M2 nodes
  //
  // NB: if M2 == numprocs, no particles
  // will be sent or received
  //
  if (myid >= M2) {
    n = initl.size();
    MPI_Send(&n, 1, MPI_UNSIGNED, myid-M2, 11, MPI_COMM_WORLD);
    if (n) {
      MPI_Send(&initl[0], n, MPI_EXP_KEYTYPE, myid-M2, 12,
	       MPI_COMM_WORLD);
    }
  }

  std::vector<key_type> data = initl;

  //
  // Retrieve the excess particles
  //
  if (myid + M2 < numprocs) {
    MPI_Recv(&n, 1, MPI_UNSIGNED, myid+M2, 11, MPI_COMM_WORLD, &status);
    if (n) {
      std::vector<key_type> recv(n);

      MPI_Recv(&recv[0], n, MPI_EXP_KEYTYPE, status.MPI_SOURCE, 12,
	       MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				// data=data+new_data
      sortCombineOne(initl, recv, data);
    }
  }

  if (myid < M2) {

    //
    // Now do the iterative binary merge
    //
    while (M2 > 1) {

      M2 = M2/2;

      // When M2 = 1, we are on the the last iteration.
      // The final node left will be the root with the entire sorted array.
      
      //
      // The upper half of the nodes send to the lower half and is done
      //
      if (myid >= M2) {
	n = data.size();
	MPI_Send(&n, 1, MPI_UNSIGNED, myid-M2, 11, MPI_COMM_WORLD);
	if (n) {
	  MPI_Send(&data[0], n, MPI_EXP_KEYTYPE, myid-M2, 12, 
		   MPI_COMM_WORLD);
	}
	break;

      } else {
	MPI_Recv(&n, 1, MPI_UNSIGNED, MPI_ANY_SOURCE, 11, MPI_COMM_WORLD, 
		 &status);
	if (n) {
	  
	  vector<key_type> recv(n);

	  
	  MPI_Recv(&recv[0], n, MPI_EXP_KEYTYPE, status.MPI_SOURCE, 12, 
		   MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	  //
	  // The lower half sorts and loop again
	  //
	  sortCombineOne(data, recv, work);
	  data = work;
	}
      }
    }
  }


  //
  // We are done, return the result
  //

  if (myid==0) final = data;

  unsigned sz = final.size();
  MPI_Bcast(&sz, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  if (myid) final.resize(sz);
  MPI_Bcast(&final[0], sz, MPI_EXP_KEYTYPE, 0, MPI_COMM_WORLD);

  std::cout << "[" << myid
	    << "] pHOT::parallelMerge: data size=" << final.size() << std::endl;

  return;
}

