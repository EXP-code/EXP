
#include "interactSelect.H"

using namespace std;
using namespace boost;

void InteractSelect::setSeed() 
{
  // generator.seed(static_cast<unsigned int>(time(0)));
  
  srand(unsigned(time(NULL)));
  // ^
  // |
  // + ---- is most probably not a good idea to use the built-in generator
  //
}

InteractSelect::InteractSelect() 
{
  setSeed();
  // uni_dist = dist_ptr(new dist_type(0, 1));
  // uni = gen_ptr(new gen_type(generator, *uni_dist));
}


double InteractSelect::selectCEInteract
(Ion& a, std::vector< std::pair< double, double > > cumCross) 
{
  typedef std::vector< std::pair< double, double > >::iterator cIter;
				// Location in the cumulative distribution
  double rn = (double)rand()/(double)RAND_MAX * cumCross.back().first;

				// Find the energy
  for (cIter i=cumCross.begin(); i!=cumCross.end(); i++) {
    if (rn <= i->first) return i->second;
  }

  return 0.0;
}

double InteractSelect::selectFFInteract(Ion& a, double E) 
{
  std::vector< double > dum;
  std::vector< double > normed;
  std::vector< double > kgrid_log;
  
  double hbc  = 197.327; // Value of h-bar * c in eV nm
  
  double Emin = hbc*pow(10, a.kgrid[0]);

  if (E < Emin) return 0.0;

  for (int i = 0; i < a.kffsteps; i++) {
    double de = a.egrid[1]-a.egrid[0];
    int e_1 = int(floor(E/de));
    int e_2 = e_1 + 1;

    // Extrapolate? [Need asymptotic expression here?]
    if (e_2 >= static_cast<int>(a.egrid.size())) {
      e_2 = a.egrid.size()-1;
      e_1 = e_2 - 1;
    }

    double y1 = a.ffCumCross[e_1][i];
    double y2 = a.ffCumCross[e_2][i];
    double x1 = a.egrid[e_1];
    double x2 = a.egrid[e_2];
    double y_tmp = y1 + ((y2-y1)/(x2-x1))*(E-x1);

    if (E > hbc*pow(10, a.kgrid[i]))
      dum.push_back(y_tmp);
  }

  std::vector <double>::iterator max = max_element(dum.begin(), dum.end());
  for (int i = 0; i < a.kffsteps; i++) {
    if (E > hbc*pow(10, a.kgrid[i])) {
      normed.push_back(dum[i]/(*max));
      kgrid_log.push_back(a.kgrid[i]);
    }
  }
  
  double k = a.kgrid[0];	// Minimum grid value (default)
  size_t n = normed.size();
				// Linear
  if (n==2) {
    double rn = (double)rand()/(double)RAND_MAX;
    k = ( (normed[1] - rn)*kgrid_log[0] + (rn - normed[0])*kgrid_log[1] ) /
      (normed[1] - normed[0]);
  }
				// Spline
  if (n>2) {
    double rn = (double)rand()/(double)RAND_MAX;
    Cspline<double, double> interp(normed, kgrid_log);
    k = interp(rn);
  }

  k = pow(10, k);

  return k*hbc;
}

double InteractSelect::DIInterLoss(chdata& ch, Ion& a) 
{
  int Z1 = a.getZ(); 
  int C1 = a.getC();

  if (C1 > Z1) {
    std::cout << "ERROR: IONIZING PAST BARE NUCLEUS" << std::endl;
    return 0;
  }

  double ip1 = ch.ipdata[Z1-1][C1-1];
  double E   = ip1;
  
  return E;
}

double InteractSelect::selectRRInteract
(Ion& a, std::vector<double> cumCross, double Ee) 
{
  std::vector<double> normed;
  std::vector<double> kgrid_log;
  const double hbc = 197.327; // value of h-bar * c in eV nm
  
  int N = a.kffsteps;
  for (int i = 0; i < a.kffsteps; i++) {
    if (cumCross[a.kffsteps-1] > 0.0) {
      if (Ee > hbc*pow(10, a.kgrid[i])) {
	normed.push_back(cumCross[i]/cumCross[a.kffsteps-1]);
	kgrid_log.push_back(a.kgrid[i]);
      } else {
	N = i-1;
	break;
      }
    }
    else return 0.0;
  }

  // Sanity check
  //
  for (int i = 0; i < N; i++) {
    if (isnan(normed[i])) 
      std::cout << "Cross nan at i= " << i 
		<< " cross = " << cumCross[i] << std::endl;
    if (i >= 1) assert(normed[i] >= normed[i-1]);
  }
  
  Cspline<double, double> interp(normed, kgrid_log);
  
  double rn = (double)rand()/(double)RAND_MAX;
  double k = interp(rn);
  k = pow(10.0, k);
  double E = k*hbc;
  
  return E;
}

