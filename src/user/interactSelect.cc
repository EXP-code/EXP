
#include "interactSelect.H"

using namespace std;
using namespace boost;

void InteractSelect::setSeed() {
	//generator.seed(static_cast<unsigned int>(time(0)));
	srand(unsigned(time(NULL)));
}

InteractSelect::InteractSelect() {
	setSeed();
	//uni_dist = dist_ptr(new dist_type(0, 1));
	//uni = gen_ptr(new gen_type(generator, *uni_dist));
}


double InteractSelect::selectCEInteract(Ion& a, std::vector< std::pair< double, double > > cumCross) {
	double eVtoRyd = 1.0/13.60569253;
	double RydtoeV = 1.0/eVtoRyd;

	int N = cumCross.size();


	std::vector< double > normed;
	normed.push_back(0);
	for(int i = 0; i < N; i++) {
		double tmp = cumCross[i].first/cumCross[N-1].first;
		normed.push_back(tmp);
	}

	double rn = (double)rand()/(double)RAND_MAX;
	double E;
	for(int i = 0; i < N; i++) {
		if (rn <= normed[i]) {
			E = cumCross[i].second;
			break;
		}
	}
	//cout << "Excite cool - dE = " << E << endl;

	return E;
}
double InteractSelect::selectFFInteract(Ion& a, double E) {
	std::vector< double > dum;
	std::vector< double > normed;
	std::vector< double > kgrid_log;


	double hbc = 197.327; //value of h-bar * c in eV nm
	if (E < a.egrid[0]) return 0;
	for(int i = 0; i < a.kffsteps; i++) {
		double de = a.egrid[1]-a.egrid[0];
		int e_1 = int(floor(E/de));
		int e_2 = e_1+1;
		//cout << e_1 << "\t" << e_2 << endl;
		double y1 = a.ffCumCross[e_1][i];
		double y2 = a.ffCumCross[e_2][i];
		double x1 = a.egrid[e_1];
		double x2 = a.egrid[e_2];
		double y_tmp = y1 + ((y2-y1)/(x2-x1))*(E-x1);
		if (E > hbc*pow(10, a.kgrid[i]))
			dum.push_back(y_tmp);
	}
	std::vector <double>::iterator max = max_element(dum.begin(), dum.end());
	for(int i = 0; i < a.kffsteps; i++) {
		if (E > hbc*pow(10, a.kgrid[i])) {
			normed.push_back(dum[i]/(*max));
			kgrid_log.push_back(a.kgrid[i]);
		}
	}
	
	//cout << "FF: " << normed.size() << "\t" << kgrid_log.size() << endl;
	Cspline<double, double> interp(normed, kgrid_log);
	double rn = (double)rand()/(double)RAND_MAX;
	double k = interp(rn);
	k = pow(10, k);
	double EF = k*hbc;

	//cout << "Ff cool: k = " << k << " EF = " << EF << endl;

	return EF;
}
double InteractSelect::DIInterLoss(chdata& ch, Ion& a) {
	int Z1, Z2, C1, C2;
	Z1 = a.getZ(); C1 = a.getC();
	Z2 = Z1; C2 = C1+1;
	if(C1 < 0) {
		cout << "ERROR: IONIZING PAST BARE NUCLEUS" << endl;
		return 0;
	}
	double ip1 = ch.ipdata[Z1-1][C1-1];
	double ip2 = ch.ipdata[Z2-1][C2+1];
	double E = ip1;

	return E;
}
double InteractSelect::selectRRInteract(Ion& a, std::vector<double> cumCross, double Ee) {

	std::vector< double > normed;
	std::vector< double > kgrid_log;
	double hbc = 197.327; //value of h-bar * c in eV nm


	int N=a.kffsteps;
	for(int i = 0; i < a.kffsteps; i++) {
		if (cumCross[a.kffsteps-1] != 0) {
			if (Ee > hbc*pow(10, a.kgrid[i])) {
				normed.push_back(cumCross[i]/cumCross[a.kffsteps-1]);
				kgrid_log.push_back(a.kgrid[i]);
			}
			else {
				N = i-1;
				break;
			}
		}
		else return 0.;
	}
	for(int i = 0; i < N; i++) {
		if(isnan(normed[i])) cout << "Cross nan at i= " << i << " cross = " << cumCross[i] << endl;
		if(i >= 1) assert(normed[i] >= normed[i-1]);
	}
	
	//cout << "RR: " << normed.size() << "\t" << kgrid_log.size() << endl;	
	Cspline<double, double> interp(normed, kgrid_log);

	double rn = (double)rand()/(double)RAND_MAX;
	double k = interp(rn);
	k = pow(10, k);
	double E = k*hbc;


	return E;
}

