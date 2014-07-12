
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
(const Ion& a, const Ion::collType& cumCross) 
{
  typedef std::vector< std::pair< double, double > >::const_iterator cIter;
				// Location in the cumulative distribution
  double rn = (double)rand()/(double)RAND_MAX * cumCross.back().first;

				// Find the energy
  for (cIter i=cumCross.begin(); i!=cumCross.end(); i++) {
    if (rn <= i->first) return i->second;
  }

  return 0.0;
}

double InteractSelect::selectFFInteract(const Ion& a, int id)
{
  std::map<int, double>::const_iterator it = a.ffWaveCrossN.find(id);
  if (it != a.ffWaveCrossN.end()) return it->second;
  else return 0.0;
}


double InteractSelect::DIInterLoss(const Ion& a) 
{
  const unsigned short Z1 = a.getZ(); 
  const unsigned short C1 = a.getC();
  const lQ Q(Z1, C1);

  if (C1 > Z1) {
    std::cout << "ERROR: IONIZING PAST BARE NUCLEUS" << std::endl;
    return 0;
  }

  double ip1 = a.getIP(Q);
  double E   = ip1;
  
  return E;
}

double InteractSelect::selectRRInteract
(const Ion& a, const std::vector<double>& cumCross, double Ee) 
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
  
  Cspline interp(normed, kgrid_log);
  
  double rn = (double)rand()/(double)RAND_MAX;
  double k = interp(rn);
  k = pow(10.0, k);
  double E = k*hbc;
  
  return E;
}

