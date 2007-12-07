#include "UserTreeDSMC.H"
#include "CollideLTE.H"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

using namespace std;

				// Proton mass (g)
const double mp = 1.67262158e-24;
				// Electron mass (g)
const double me =  9.1094e-28;
				// Speed of light (cm/s)
const double c =  2.99792458e+10;
				// electron volt in (cgs)
const double eV =  1.6020e-12;
				// Boltzmann constant (cgs)
const double boltz = 1.3810e-16;
				// Planck's constant (cgs)
const double planck  = 6.6262e-27;
				// Electron charge (cgs)
const double esu = 4.8032e-10;
				// Bohr radius
const double a0 = planck*planck/(4.0*M_PI*M_PI*me*esu*esu);

				// Neutral Hydrogen binding energy
const double E1 = 2.0*M_PI*M_PI*me*pow(esu, 4.0)/(planck*planck);

				// Hydrogen fraction
const double f_H = 0.76;


unsigned CollideLTE::trhocnt = 0;

CollideLTE::CollideLTE(double diameter, int Nth) : Collide(diameter, Nth)
{

  HeatCool::initialize();

  if (myid==0) {
    const double levels [13] = 
      {0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99};

    // Attempt open for read     
    fstream fs("collideL.debug", ios_base::in);
    if (!fs)
      {
	// File doesn't exist; create a new one
	fs.open("collideL.debug", ios_base::out); 
	fs << "#" << endl
	   << "# " << setw(13) << "Time"
	   << setw(15) << "Min T"
	   << setw(15) << "Avg T"
	   << setw(15) << "Max T"
	   << setw(15) << "Disp T";
	for (int n=0; n<13; n++) fs << setw(15) << levels[n];
	fs << endl << "# " << setw(13) << 1;
	for (int n=1; n<18; n++) fs << "| " << setw(13) << n+1;
	fs << endl << "#" << endl;
      }
    
  }

  cellcnt = vector<unsigned>(nthrds, 0);
  minT = vector<double>(nthrds, 1e30);
  maxT = vector<double>(nthrds, 0.0);
  avgT = vector<double>(nthrds, 0.0);
  dispT = vector<double>(nthrds, 0.0);
  tlist = vector< vector<double> >(nthrds);

  numt = 20;
  tmin = log(5000.0);
  tmax = log(200000.0);
  dtmp = (tmax - tmin)/(numt - 1);
  for (unsigned n=0; n<numt; n++) {
    thisto1.push_back(tmin + dtmp*n);
    thisto2.push_back(0.0);
  }

  numn = 20;
  nmin = log(75.0);
  nmax = log(150.0);
  ntmp = (nmax - nmin)/(numn - 1);
  for (unsigned n=0; n<numn; n++) {
    nhisto1.push_back(nmin + ntmp*n);
    nhisto2.push_back(0.0);
  }

  trho = vector< vector<double> >(numt);
  for (unsigned n=0; n<numt; n++) trho[n] = vector<double>(numn, 0);

  lostSoFar = vector<double>(nthrds, 0.0);
  deltaE = vector<double>(nthrds);
}


void CollideLTE::initialize_cell(pCell* cell, 
				 double rvmax, double tau, double number, 
				 int id)
{
				// Cell temperature and mass (cgs)
				// 
  double KEtot, KEdsp;
  cell->KE(KEtot, KEdsp);	// These are already specific in mass
  double mass = cell->Mass();
  
  double rvrel = sqrt(4.0*KEdsp);
  KEdsp *= UserTreeDSMC::Eunit;

  double Mass = mass * UserTreeDSMC::Munit;
  double mm = f_H*mp + (1.0-f_H)*4.0*mp;
  double T = 2.0*KEdsp/3.0 * mm/UserTreeDSMC::Munit/boltz;
  double volume = cell->Volume();
  double Volume = volume * pow(UserTreeDSMC::Lunit, 3);
  double Density = Mass/Volume;
  double n0 = Density/mp;

  if (log(T)>tmin && log(T)<tmax) {
    int indx = (int)float( (log(T) - tmin)/dtmp );
    thisto2[indx] += cell->Mass();
  }

  if (log(n0)>nmin && log(n0)<nmax) {
    int indx = (int)float( (log(n0) - nmin)/ntmp );
    nhisto2[indx] += cell->Mass();
  }

  if (log(T) >tmin && log(T) <tmax &&
      log(n0)>nmin && log(n0)<nmax) {
    int indx1 = (int)float( (log(T)  - tmin)/dtmp );
    int indx2 = (int)float( (log(n0) - nmin)/ntmp );
    trho[indx1][indx2] += cell->Mass();
  }

  cellcnt[id]++;
  minT[id] = min<double>(minT[id], T);
  maxT[id] = max<double>(maxT[id], T);
  avgT[id] += T;
  dispT[id] += T*T;
  tlist[id].push_back(T);

  number *= min<double>(1.0, rvrel/rvmax);

  double n_h = n0*f_H;

  HeatCool heatcool(n0, T);

  deltaE[id] = heatcool.CoolRate()/UserTreeDSMC::Eunit * n_h*n_h * Volume * tau * UserTreeDSMC::Tunit / mass / number;

				// Assign temp and/or density to particles
  if (use_temp>=0 || use_dens>=0) {
    
    unsigned nbods = cell->bods.size();
    for (unsigned j=0; j<nbods; j++) {
      if (cell->bods[j] == 0) {
	cout << "proc=" << myid << " id=" << id 
	     << " ptr=" << hex << cell << dec
	     << " j=" << j << " indx=" << cell->bods[j] 
	     << "/" << cell->bods.size() << endl;
	continue;
      }

      int sz = cell->Body(j)->dattrib.size();
      if (use_temp>=0 && use_temp<sz) 
	cell->Body(j)->dattrib[use_temp] = T;
      if (use_dens>=0 && use_dens<sz) 
	cell->Body(j)->dattrib[use_dens] = mass/volume;
    }
  }

}


int CollideLTE::inelastic(pHOT *tree, Partstruct* p1, Partstruct* p2,
			  double *cr, int id)
{
  int ret = 0;			// No error (flag)

				// Reduced mass in system units
  double Mu = p1->mass*p2->mass/(p1->mass + p2->mass);

				// Consistent: KE in coll. frame is
  if ((*cr)*(*cr) > 2.0*deltaE[id]/Mu) {
    lostSoFar[id] += deltaE[id];	// larger than the energy radiated
    (*cr) = sqrt((*cr)*(*cr) - 2.0*deltaE[id]/Mu);
    ret = 0;			// No error
  }
  else {			// Inconsistent: too much energy lost!
    lostSoFar[id] += 0.5*Mu*(*cr)*(*cr);
    (*cr) = 0.0;
    ret = 1;			// Set error flag
  }

  return ret;
}


double CollideLTE::Elost()
{ 
  double ret=0.0;
  for (int n=0; n<nthrds; n++) {
    ret += lostSoFar[n];
    lostSoFar[n] = 0.0; 
  }
  return ret; 
}


void CollideLTE::Debug(double t)
{
  unsigned cellcnt0;
  double minT0, maxT0, avgT0, dispT0;
  const double levels [13] = {0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99};
  double values[13], values0[13];
  
  // Combine from threads . . .

  unsigned cellcnt1 = 0;
  double minT1 = 1e20;
  double maxT1 = 0.0;
  double avgT1 = 0.0;
  double dispT1 = 0.0;
  vector<double> tlist1;

  for (int n=0; n<nthrds; n++) {
    cellcnt1 += cellcnt[n];
    minT1 = min<double>(minT[n], minT1);
    maxT1 = min<double>(maxT[n], maxT1);
    avgT1 += avgT[n];
    dispT1 += dispT[n];
    tlist1.insert(tlist1.end(), tlist[n].begin(), tlist[n].end());
  }

  sort(tlist1.begin(), tlist1.end());
  for (int n=0; n<13; n++) {
    if (tlist1.size())
      values[n] = tlist1[(int)floor(levels[n]*tlist1.size()+0.5)];
    else
      values[n] = 0.0;
  }

  vector<double> thist20(numt, 0);
  MPI_Reduce(&thisto2[0], &thist20[0], numt, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  vector<double> nhist20(numt, 0);
  MPI_Reduce(&nhisto2[0], &nhist20[0], numn, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  MPI_Reduce(&values,   &values0, 13, MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&minT1,    &minT0, 1,    MPI_DOUBLE,   MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&maxT1,    &maxT0, 1,    MPI_DOUBLE,   MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&avgT1,    &avgT0, 1,    MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&dispT1,   &dispT0, 1,   MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&cellcnt1, &cellcnt0, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

  if (myid==0) {
    ofstream out("collideL.debug", ios::app);
    double avg=0.0, disp=0.0;
    if (cellcnt0>0) avg  = avgT0/cellcnt0;
    if (cellcnt0>1) disp = sqrt((dispT0 - avgT0*avgT0/cellcnt0)/(cellcnt0-1));
    out << setw(15) << t
	<< setw(15) << minT0       
	<< setw(15) << avg
	<< setw(15) << maxT0
	<< setw(15) << disp;
    for (int n=0; n<13; n++) 
      out << setw(15) << values0[n]/numprocs;
    out << endl;

    double mast = 0.0;
    for (unsigned n=0; n<numt; n++) mast += thist20[n];


    double cum = 0.0;

    ofstream out2("collideH.thisto", ios::app);
    for (unsigned n=0; n<numt; n++) {
      cum += thist20[n];
      out2 << setw(18) << t
	   << setw(18) << exp(thisto1[n])
	   << setw(18) << thist20[n]/mast
	   << setw(18) << cum/mast
	   << endl;
    }
    out2 << endl;

    ofstream out3("collideH.nhisto", ios::app);
    cum = 0.0;
    for (unsigned n=0; n<numn; n++) {
      cum += nhist20[n];
      out3 << setw(18) << t
	   << setw(18) << exp(nhisto1[n])
	   << setw(18) << nhist20[n]/mast
	   << setw(18) << cum/mast
	   << endl;
    }
    out3 << endl;

    if ( (trhocnt % 10) == 0) {
      ostringstream sout;
      sout << "collideH.trho." << trhocnt;
      ofstream out4(sout.str().c_str());
      out4 << "# T=" << t << endl;
      for (unsigned n=0; n<numt; n++) {
	for (unsigned m=0; m<numn; m++) {
	  out4 << setw(18) << exp(thisto1[n])
	       << setw(18) << exp(nhisto1[m])
	       << setw(18) << trho[n][m]
	       << endl;
	}
	out4 << endl;
      }
    }
    trhocnt++;
  }

  for (int n=0; n<nthrds; n++) {
    cellcnt[n] = 0;
    minT[n] = 1e33;
    maxT[n] = 0.0;
    avgT[n] = 0.0;
    dispT[n] = 0.0;
    tlist[n].clear();
  }

  for (unsigned n=0; n<numt; n++) thisto2[n] = 0.0;
  for (unsigned m=0; m<numn; m++) nhisto2[m] = 0.0;
  for (unsigned n=0; n<numt; n++) 
    for (unsigned m=0; m<numn; m++) trho[n][m] = 0.0;
  
}
