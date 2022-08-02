#include <random>
#include <chrono>

#include <KMeans.H>

void KMeans::kMeansClustering::iterate(KMeans::KMeansDistance& distance,
				       int niter, int k, int s, bool verbose)
{
  // Compute initial cen
  //
  cen.clear();

  make_sane();			// Remove null points from classes
				// vector

  if (s>0) {			// Seed centers by stride
    for (int i=0; i<classes.size(); i+=s) {
      if (cen.size()>=k) break;
      cen.push_back(classes.at(i)->x);
    }
    k = cen.size();
  } else {		// obtain a seed from the system clock:
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 gen(seed);
			// Seed centers randomly
    for (int i=0; i<k; ++i) {
      cen.push_back(classes.at(gen() % classes.size())->x);
    }
  }

  // Save last center
  //
  auto last(cen);

  // Iterations
  //
  for (int n=0; n<niter; n++) {

    for (int id=0; id<k; id++) {
      for (auto p : classes) {
	double dist = distance(p->x, cen[id]);
	if (dist < p->minDist) {
	  p->minDist = dist;
	  p->cid = id;
	}
      }
    }
    
    // Initialize for update
    //
    std::vector<int> nPoints(k, 0);
    std::vector< std::vector<double> > sumX(k);
    for (int id=0; id<k; id++) sumX[id].resize(ndim, 0);
    
    // Iterate over points to append data to new centroid
    //
    for (auto p : classes) {
      if (p->cid>=0) {
	nPoints[p->cid] += 1;
	for (int i=0; i<ndim; i++) sumX[p->cid][i] += p->x[i];
	p->minDist = std::numeric_limits<double>::max();
      }
    }
    
    // Compute the new centroid
    //
    last = cen;
    for (int id=0; id<k; id++) {
      if (nPoints[id]) {
	for (int i=0; i<ndim; i++)
	  cen[id][i] = sumX[id][i]/nPoints[id];
      }
    }

    std::vector<double> cdiff(k, 0.0);
    total = 0.0;
    for (int id=0; id<k; id++) {
      for (int i=0; i<ndim; i++)
	cdiff[id] += (cen[id][i] - last[id][i])*(cen[id][i] - last[id][i]);
      total += cdiff[id];
    }

    if (verbose) {
      std::cout << "Iteration " << n << ", total=" << total << std::endl;
      for (int id=0; id<k; id++) {
	std::cout << std::setw(12) << cdiff[id];
	for (int i=0; i<ndim; i++)
	  std::cout << std::setw(12) << cen[id][i];
	std::cout << std::endl;
      }
    }

    if (total<=0.0) break;
  }

}


double KMeans::EuclideanDistance::operator()
  (const std::vector<double>& x, const std::vector<double>& y)
{
  double dist = 0.0;
  int ndim = x.size();
  for (int i=0; i<ndim; i++) dist += (x[i] - y[i]) * (x[i] - y[i]);
  return dist;
}


double KMeans::WcorrDistance::operator()
  (const std::vector<double>& x, const std::vector<double>& y)
{
  // A Lambda for the weight function
  auto w = [&](int i) {
	     if      (i < Lstar) return i;
	     else if (i < Kstar) return Lstar;
	     else                return numT - i + 1;
	   };
  
  double corr = 0.0, nrmx = 0.0, nrmy = 0.0;
  for (int i=0; i<numT; i++) {
    corr += w(i) * x[i] * y[i];
    nrmx += w(i) * x[i] * x[i];
    nrmy += w(i) * y[i] * y[i];
  }
  
  double ret = 1.0;
  if (nrmx*nrmy > 0.0) ret -= sqrt(corr/sqrt(nrmx*nrmy));
  
  return ret;
}

double KMeans::WcorrDistMulti::operator()
  (const std::vector<double>& x, const std::vector<double>& y)
{
  // A Lambda for the weight function
  auto w = [&](int i) {
	     if      (i < Lstar) return i;
	     else if (i < Kstar) return Lstar;
	     else                return numT - i + 1;
	   };
  
  double corr = 0.0, nrmx = 0.0, nrmy = 0.0;
  for (int n=0; n<nchn; n++) {
    for (int i=0; i<numT; i++) {
      corr += w(i) * x[i+n*numT] * y[i+n*numT];
      nrmx += w(i) * x[i+n*numT] * x[i+n*numT];
      nrmy += w(i) * y[i+n*numT] * y[i+n*numT];
    }
  }
  
  double ret = 1.0;
  if (nrmx*nrmy > 0.0) ret -= sqrt(corr/sqrt(nrmx*nrmy));
  
  return ret;
}

