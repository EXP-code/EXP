#ifndef _CollideLTE_H
#define _CollideLTE_H

#include "HeatCool.H"
#include "Collide.H"

class CollideLTE : public Collide
{
protected:
  // For debugging
  vector<unsigned> cellcnt;
  vector<double> minT, maxT, avgT, dispT;
  vector< vector<double> > tlist;

  unsigned numt;
  double tmin, tmax, dtmp;
  vector<double> thisto1, thisto2;

  unsigned numn;
  double nmin, nmax, ntmp;
  vector<double> nhisto1, nhisto2;

  static unsigned trhocnt;
  vector< vector<double> > trho;

protected:
  vector<double> deltaE;

  void initialize_cell(pCell* cell, double vrelmax, 
		       double dt, double number, int id=0);

  int inelastic(pHOT *tree, Partstruct* p1, Partstruct* p2, double *crm,
		int id=0);

  // Diagnostic output
  vector<double> lostSoFar;

public:

  CollideLTE(double diameter, int nth=1);
  double Elost();

  void Debug(double t);

};

#endif
