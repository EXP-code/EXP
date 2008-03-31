#ifndef _CollideLTE_H
#define _CollideLTE_H

#include "HeatCool.H"
#include "Collide.H"

class CollideLTE : public Collide
{
private:
  vector<Precord> prec;

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

  bool debug_enabled;

protected:
  vector<double> deltaE;

  void initialize_cell(pCell* cell, double vrelmax, 
		       double dt, double number, int id=0);

  int inelastic(pHOT *tree, Particle* p1, Particle* p2, double *crm, int id=0);

  // Diagnostic output
  //
  double totalSoFar;
  vector<double> lostSoFar;

  virtual void list_sizes();
  void list_sizes_proc(ostream*);

public:

  CollideLTE(double diameter, int nth=1);

  virtual void set_timestep(int DTpos) { use_delt=DTpos; }

  virtual double Etotal();
  virtual void Elost(double* collide, double* epsm);

  void Debug(double t);

};

#endif
