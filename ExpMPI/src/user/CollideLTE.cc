#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <cfloat>

#include "global.H"
#include "UserTreeDSMC.H"
#include "CollideLTE.H"

using namespace std;

				// Proton mass (g)
const double mp      = 1.67262158e-24;
				// Electron mass (g)
const double me      =  9.1094e-28;
				// Speed of light (cm/s)
const double c       =  2.99792458e+10;
				// electron volt in (cgs)
const double eV      =  1.6020e-12;
				// Boltzmann constant (cgs)
const double boltz   = 1.3810e-16;
				// Planck's constant (cgs)
const double planck  = 6.6262e-27;
				// Electron charge (cgs)
const double esu     = 4.8032e-10;
				// Bohr radius
const double a0      = planck*planck/(4.0*M_PI*M_PI*me*esu*esu);

				// Neutral Hydrogen binding energy
const double E1      = 2.0*M_PI*M_PI*me*pow(esu, 4.0)/(planck*planck);

				// Hydrogen fraction
const double f_H     = 0.76;
				// Time step control for energy solution
const double rfac    = 0.01;
				// Minimum kinetic temperature for cooling
const double Tlow    = 100.0;

				// Default global class parameters
double   CollideLTE::Nmin    = 1.0e-08;
double   CollideLTE::Nmax    = 1.0e+25;
double   CollideLTE::Tmin    = 1.0e+03;
double   CollideLTE::Tmax    = 1.0e+08;
double   CollideLTE::TolV    = 1.0e-03;
unsigned CollideLTE::Nnum    = 400;
unsigned CollideLTE::Tnum    = 200;
string   CollideLTE::cache   = ".HeatCool";
unsigned CollideLTE::trhocnt = 0;
bool     CollideLTE::frost_warning = false;

CollideLTE::CollideLTE(ExternalForce *force, double diameter, int Nth) : 
  Collide(force, diameter, Nth)
{

  HeatCool::initialize();
  hc = new HeatCool(Nmin, Nmax, Tmin, Tmax, Nnum, Tnum, cache);

  cellcnt = vector<unsigned>(nthrds, 0);
  minT    = vector<double>(nthrds, 1e30);
  maxT    = vector<double>(nthrds, 0.0);
  avgT    = vector<double>(nthrds, 0.0);
  dispT   = vector<double>(nthrds, 0.0);
  tlist   = vector< vector<double> >(nthrds);

  debug_enabled = true;

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

  deltaE = vector<double>(nthrds);

				// Energy diagnostics
  totalSoFar = 0.0;
  massSoFar = 0.0;
  lostSoFar  = vector<double>(nthrds, 0.0);

  prec = vector<Precord>(nthrds);
  for (int n=0; n<nthrds; n++)
    prec[n].second = vector<double>(Nphase);
}


CollideLTE::~CollideLTE()
{
  delete hc;
}

void CollideLTE::initialize_cell(pCell* cell, 
				 double rvmax, double tau, double number, 
				 int id)
{
  sCell *samp = cell->sample;
				// Cell temperature and mass (cgs)
				// 
  double KEtot, KEdspS, KEdspC;
  samp->KE(KEtot, KEdspS);	// These are already in specific mass
  cell->KE(KEtot, KEdspC);

  double massC = cell->Mass();	// Mass in real cell
  double massS = samp->Mass();	// Mass in sample cell

  totalSoFar += massC * KEdspC;
  massSoFar  += massC;

				// Mean mass per particle
  double Mass = massC * UserTreeDSMC::Munit;
  double mm   = f_H*mp + (1.0-f_H)*4.0*mp;
  double T    = 2.0*KEdspS*UserTreeDSMC::Eunit/3.0 * mm/UserTreeDSMC::Munit/boltz;

				// Volume in cells
  double volumeC = cell->Volume();
  double volumeS = samp->Volume();

  double Volume  = volumeC * pow(UserTreeDSMC::Lunit, 3);
  double Density = Mass/Volume;
  double n0      = Density/mp;
  double T0      = T;
  double h0      = 0.0;

				// Volume in real cell
  double CellVolume = volumeC * pow(UserTreeDSMC::Lunit, 3);

  if (T>0.0) {
    if (log(T)>tmin && log(T)<tmax) {
      int indx = (int)float( (log(T) - tmin)/dtmp );
      thisto2[indx] += massC;
    }
  }

  if (n0>0.0) {
    if (log(n0)>nmin && log(n0)<nmax) {
      int indx = (int)float( (log(n0) - nmin)/ntmp );
      nhisto2[indx] += massC;
    }
  }

  if (T>0.0 && n0>0.0) {
    if (log(T) >tmin && log(T) <tmax &&
	log(n0)>nmin && log(n0)<nmax) {
      int indx1 = (int)float( (log(T)  - tmin)/dtmp );
      int indx2 = (int)float( (log(n0) - nmin)/ntmp );
      trho[indx1][indx2] += massC;
    }
  }
    
  if (debug_enabled) {
    cellcnt[id]++;
    minT[id] = min<double>(minT[id], T);
    maxT[id] = max<double>(maxT[id], T);
    avgT[id] += T;
    dispT[id] += T*T;
    tlist[id].push_back(T);
  }

  // Hydrogen number density
  //
  // double n_h = n0*f_H;

  coolTime[id].start();

  // CoolRate has units erg*cm^3/t
  // So CoolRate*n_h*n_h*Volume*Time = ergs

  double KEdspF = 0.0;		// Final KE, used for ESOL only

  // Total energy lost (for both collisions and EPSM)
  //
  if (NOCOOL || n0<=0.0 || T<=Tlow)
    coolheat[id] = 0.0;
  else {
				// Convert to energy rate in system units
    double Cfac = CellVolume * n0 * n0 *
      UserTreeDSMC::Tunit / UserTreeDSMC::Eunit;

				// Energy to lose/gain
    coolheat[id] = hc->CoolRate(n0, T) * Cfac * tau * ENHANCE;

    if (ESOL) {
				// Instantaneous loss rate
      double h = tau*cell->Mass()*KEdspS/(coolheat[id]+DBL_MIN);
      h0 = h = rfac * min<double>(tau, h);

				// Energy to temperature
      double Tfac = 3.0*UserTreeDSMC::Munit/UserTreeDSMC::Eunit/2.0 *
	cell->Mass()*boltz/mm;

				// Initial energy
      double E = cell->Mass()*KEdspS;
      double E0 = E;


#ifdef DEBUG
      int icnt=0;
#endif
      double k1, k2, k3, k4, t=0.0, TT;
      while (t<tau*(1.0-2.0*DBL_MIN)) {
	h = min<double>(h, tau-t);
				// Safety
	h = max<double>(h, 1.0e-4*tau);
                                // RK4
	TT = max<double>(1.0e-8, E/Tfac);
	k1 = -h * Cfac * hc->CoolRate(n0, TT);
	k2 = -h * Cfac * hc->CoolRate(n0, TT+k1*0.5/Tfac);
	k3 = -h * Cfac * hc->CoolRate(n0, TT+k2*0.5/Tfac);
	k4 = -h * Cfac * hc->CoolRate(n0, TT+k3/Tfac);

	E += (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
	t += h;
#ifdef DEBUG
	icnt++;
#endif
      }
				// Sanity: failure of explicit solution
      if (isnan(E)) E = E0;
      E = max<double>(E, 3.0*Tfac);
				// Final energy per unit mass in the cell
      KEdspF = E/cell->Mass();
				// Effective final temperature
      T      = E/Tfac;

#ifdef DEBUG
      if (icnt > 10000) {
	cout << "Process " << setw(4) << myid 
	     << " [" << setw(2) << id << "]:" << left
	     << " #="      << setw(8) << icnt
	     << " T_0="    << setw(15) << T0
	     << " T_F="    << setw(15) << T 
	     << endl << right;
      }
      if (T > T0*1.001) {
	cout << "Process " << setw(4) << myid 
	     << " [" << setw(2) << id << "]:" << left
	     << " #="      << setw(8) << icnt
	     << " E_0="    << setw(15) << E0
	     << " E_F="    << setw(15) << E
	     << " T_0="    << setw(15) << T0
	     << " T_F="    << setw(15) << T 
	     << " k_1="    << setw(15) << k1/Tfac 
	     << " k_4="    << setw(15) << k4/Tfac 
	     << endl << right;
	T = T0;
      }
#endif
    }
				// Sanity: failure of implicit solution
    if (isnan(coolheat[id])) coolheat[id] = 0.0;
  }

  coolSoFar[id] = coolTime[id].stop();

  // Energy per encounter
  //
  if (number<=0.0 || NOCOOL)
    deltaE[id] = 0.0;
  else {
    if (ESOL)
      deltaE[id] = cell->Mass()*(KEdspS - KEdspF) / number;
    else
      deltaE[id] = coolheat[id] / number;
  }

  if (fabs((KEdspS - KEdspF)/KEdspS)>1.0) {
    cout << "|deltaE/E|=" << (KEdspS - KEdspF)/KEdspS << " > 1"
	 << ", T=" << T << ", T0=" << T0 << ", h=" << h0 << ", tau=" << tau
	 << ", KEdspS=" << KEdspS << ", KEdspF=" << KEdspF << endl;
  } else if ( (KEdspS - KEdspF)/KEdspS < -0.05 ) {
    cout << " deltaE/E =" << (KEdspS - KEdspF)/KEdspS << " < -0.05"
	 << ", T=" << T << ", T0=" << T0 << ", h=" << h0 << ", tau=" << tau
	 << ", KEdspS=" << KEdspS << ", KEdspF=" << KEdspF << endl;
  }

  if (frost_warning && T<1000.0) {
    vector<double> pos(3);
    cell->MeanPos(pos);

    ostringstream sout;
    sout << outdir << "CollideLTE_diag." << runtag << "." << id << "." << myid;
    ofstream out(sout.str().c_str(), ios::app);
    out << setw(16) << tnow
	<< setw(16) << T
	<< setw(16) << coolheat[id];
    if (ESOL && !NOCOOL) out << setw(16) << cell->Mass()*(KEdspS - KEdspF);
    out << setw(16) << Density
	<< setw(16) << massC
	<< setw(16) << massC/volumeC
	<< setw(16) << massS/volumeS;
    for (int i=0; i<3; i++) out << setw(16) << pos[i];
    out << endl;
  }

  if (MFPDIAG) {
    ttempT[id].push_back(T);
    tdeltT[id].push_back(deltaE[id]);

    prec[id].first = massC/volumeC;
    prec[id].second[0] = T;
    prec[id].second[1] = cell->bods.size();
    prec[id].second[2] = cell->Mass();
    prec[id].second[3] = cell->Volume();
    tphaseT[id].push_back(prec[id]);
  }

  // Assign temp and/or density to particles
  //
  if (use_temp>=0 || use_dens>=0) {
    
    double dens = massC/volumeC;
    // set<unsigned long>::iterator j = cell->bods.begin();
    vector<unsigned long>::iterator j = cell->bods.begin();
    while (j != cell->bods.end()) {
      if (*j == 0) {
	cout << "proc=" << myid << " id=" << id 
	     << " ptr=" << hex << cell << dec
	     << " indx=" << *j << "/" << cell->bods.size() << endl;
	j++;
	continue;
      }

      int sz = cell->Body(j)->dattrib.size();
      if (use_temp>=0 && use_temp<sz) 
	cell->Body(j)->dattrib[use_temp] = T;
      if (use_dens>=0 && use_dens<sz) 
	cell->Body(j)->dattrib[use_dens] = dens;
      j++;
    }
  }

				// Energy ratio for time step estimation
  if (use_delt>=0) {

    double Ctime =  tau*cell->Mass()*KEdspC/(coolheat[id]+DBL_MIN);
    //                                                     ^
    // to prevent inf values ------------------------------|
    //
				// Diagnose cooling time step in this cell
    int indx = 0;
    if (Ctime>0.0) indx = (int)floor(log(Ctime/tau)/log(4.0) + 5);
    if (indx<0 )   indx = 0;
    if (indx>10)   indx = 10;
    tcoolT[id][indx]++;

				// Assign per body time step requests

    vector<unsigned long>::iterator j;
    for (j=cell->bods.begin(); j!=cell->bods.end(); j++) {
      if (*j == 0) {
	cout << "proc=" << myid << " id=" << id 
	     << " ptr=" << hex << cell << dec
	     << " indx=" << *j
	     << "/" << cell->bods.size() << endl;
	continue;
      }

      int sz = cell->Body(j)->dattrib.size();
      if (use_delt>=0 && use_delt<sz) 
	cell->Body(j)->dattrib[use_delt] = Ctime;
    }

  }

}


int CollideLTE::inelastic(pHOT *tree, Particle* p1, Particle* p2, 
			  double *cr, int id)
{
  int ret = 0;			// No error (flag)

				// Reduced mass in system units
  double Mt = p1->mass + p2->mass;
  if (Mt<=0.0) return ret;

  double Mu = p1->mass*p2->mass/Mt;
  if (Mu<=0.0) return ret;

				// Energy floor
  double kE = 0.5*Mu*(*cr)*(*cr);
  double dE = kE*TolV*TolV;
  double remE = kE - dE;
  double delE = deltaE[id];

  if (remE<=0.0 || delE<=0.0) return ret;

				// Cooling rate diagnostic
  if (TSDIAG) {
    if (delE>0.0) {
      int indx = (int)floor(log(remE/delE)/(log(2.0)*TSPOW) + 5);
      if (indx<0 ) indx = 0;
      if (indx>10) indx = 10;
      
      EoverT[id][indx] += Mt;
    }
  }

  if (use_exes>=0) {
				// (-/+) value means under/overcooled: 
				// positive/negative increment to delE
				// NB: delE may be < 0 if too much energy 
				// was radiated previously . . .
				//
    delE -= p1->dattrib[use_exes] + p2->dattrib[use_exes];
  }

				// Consistent: KE in coll. frame is
  if (remE >= delE) {
    lostSoFar[id] += delE;	// larger than the energy radiated
    decelT[id]    += delE;
    (*cr) = sqrt( 2.0*(kE - delE)/Mu );
    ret = 0;			// No error

				// Zero out internal energy excess
    if (use_exes>=0)		// since excess is now used up
      p1->dattrib[use_exes] = p2->dattrib[use_exes] = 0.0;
  }
  else {			// Inconsistent: too much energy lost!
    lostSoFar[id] += remE;
    decolT[id]    += remE - delE;
    (*cr) *= TolV;
    ret = 1;			// Set error flag

				// Reset internal energy excess
    if (use_exes>=0) {
      if (ENSEXES) 		// Energy will be spread later
	p1->dattrib[use_exes] = p2->dattrib[use_exes] = 0.0;
      else {			// Energy excess incorporated now
	p1->dattrib[use_exes] =  p1->mass*(remE - delE)/Mt;
	p2->dattrib[use_exes] =  p2->mass*(remE - delE)/Mt;
      }
    }
  }

  return ret;
}


double CollideLTE::Etotal()
{ 
  double ret = totalSoFar;
  totalSoFar = 0.0;
  return ret; 
}

double CollideLTE::Mtotal()
{ 
  double ret = massSoFar;
  massSoFar = 0.0;
  return ret; 
}

void CollideLTE::Elost(double* collide, double* epsm)
{ 
  double ret1=0.0, ret2=0.0;
  for (int n=0; n<nthrds; n++) {
    ret1 += lostSoFar[n];
    ret2 += lostSoFar_EPSM[n];
    lostSoFar[n] = 0.0; 
    lostSoFar_EPSM[n] = 0.0; 
  }
  *collide = ret1;
  *epsm = ret2;
}


void CollideLTE::Debug(double t)
{
  unsigned cellcnt0;
  double minT0, maxT0, avgT0, dispT0;
  const double levels [13] = 
    {0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99};
  double values[13], values0[13];
  
				// Will begin accumulating now . . .
  if (!debug_enabled) {
    debug_enabled = true;
    return;
  }

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

    // Attempt open for read     
    {
      fstream fs("collideL.debug", ios_base::in);
      if (!fs) {
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
    
    ofstream out(string(outdir + "collideL.debug").c_str(), ios::app);
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


    const char marker[] = "XXXXXX";
    double cum = 0.0;

    ofstream out2(string(outdir + runtag + ".collideH.thisto").c_str(), 
		  ios::app);
    for (unsigned n=0; n<numt; n++) {
      cum += thist20[n];
      out2 << setw(18) << t
	   << setw(18) << exp(thisto1[n]);
      if (mast>0.0)
	out2 << setw(18) << thist20[n]/mast
	     << setw(18) << cum/mast
	     << endl;
      else
	out2 << setw(18) << marker
	     << setw(18) << marker
	     << endl;
    }
    out2 << endl;

    ofstream out3(string(outdir + runtag + ".collideH.nhisto").c_str(), 
		  ios::app);
    cum = 0.0;
    for (unsigned n=0; n<numn; n++) {
      cum += nhist20[n];
      out3 << setw(18) << t
	   << setw(18) << exp(nhisto1[n]);
      if (mast>0.0)
	out3 << setw(18) << nhist20[n]/mast
	     << setw(18) << cum/mast
	     << endl;
      else
	out3 << setw(18) << marker
	     << setw(18) << marker
	     << endl;
    }
    out3 << endl;

    if ( (trhocnt % 10) == 0) {
      ostringstream sout;
      sout << outdir << runtag << ".collideH.trho." << trhocnt;
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


void CollideLTE::list_sizes()
{
  string sname = outdir + runtag + ".collide_storage";
  for (int n=0; n<numprocs; n++) {
    if (myid==n) {
      ofstream out(sname.c_str(), ios::app);
      if (out) {
	out << setw(18) << tnow
	    << setw(6)  << myid;
	list_sizes_proc(&out);
	Collide::list_sizes_proc(&out);
	out << endl;
	if (myid==numprocs-1) out << endl;
	out.close();
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

void CollideLTE::list_sizes_proc(ostream* out)
{
  *out << setw(12) << prec.size()
       << setw(12) << cellcnt.size()
       << setw(12) << minT.size()
       << setw(12) << maxT.size()
       << setw(12) << avgT.size()
       << setw(12) << dispT.size()
       << setw(12) << thisto1.size()
       << setw(12) << thisto2.size()
       << setw(12) << nhisto1.size()
       << setw(12) << nhisto2.size()
       << setw(12) << deltaE.size()
       << setw(12) << lostSoFar.size()
       << setw(12) << (trho.size()  ? trho[0].size()  : (size_t)0)
       << setw(12) << (tlist.size() ? tlist[0].size() : (size_t)0);
}
