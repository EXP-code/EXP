 // This may look like C code, but it is really -*- C++ -*-

#ifndef _orbit_h
#define _orbit_h 1

const char rcsid_orbit[] = "$Id$";

#include <logic.h>
#include <biorth.h>
#include <gaussQ.h>

class AxiSymModel;

class RegularOrbit
{
public:

  int dof;
  string OrbitID;
  
  bool model_defined;
  bool action_defined;
  bool freq_defined;
  bool angle_defined;
  bool biorth_defined;
  
  //  MassModel *model;
  
  Vector action;
  Vector freq;
  
  virtual ~RegularOrbit() {}

  virtual double get_action(const int) = 0;
  virtual double get_freq(const int) = 0;
  virtual double get_angle(const int, const double) = 0;
  
  // Error function
  
  void bomb(string s) {
    cerr << "ERROR from " << OrbitID << ": " << s << '\n';
#ifdef DEBUG
    abort();
#endif
    exit(-1);
  }
  
  void bomb(char *s) {
    cerr << "ERROR from " << OrbitID << ": " << s << '\n';
#ifdef DEBUG
    abort();
#endif
    exit(-1);
  }
  
};

struct ANGLE_GRID {
  Matrix t;  // RECS
  Matrix w1;
  Matrix dw1dt;
  Matrix f;
  Matrix r;
  Matrix fr; // NMAX x RECS
  int num;
};



class SphericalOrbit : public RegularOrbit
{
  
private:
  double energy;
  double kappa;
  double beta;
  
  double jmax;
  double r_peri;
  double r_apo;
  double r_circ;
  
  AxiSymModel *model;
  AxiSymBiorth *biorth;
  int RECUR;
  
  void compute_action(void) { compute_freq(); }
  void compute_freq(void);
  void compute_angles(void);
  void compute_freq_epi(void);
  void compute_angles_epi(void);
  void compute_biorth(void);
  
  struct ANGLE_GRID angle_grid;
  int recs;
  int nmax;
  int nbsct;
  
  int l;
  int l1s, l2s;
  Vector cosvec;
  
  LegeQuad *Gkn;
  
public:
  
  // Constructors
  
  SphericalOrbit(void);
  SphericalOrbit(AxiSymModel *model);
  SphericalOrbit(AxiSymModel *model, double Energy, double kappa,
		 double beta = 0.0);
  virtual ~SphericalOrbit();
  
  //  SphericalOrbit(SphericalOrbit &orb);
  SphericalOrbit &operator=(const SphericalOrbit &);
  
  // Required functions
  
  double get_action(const int i) {
    if (i>dof) bomb("actions not defined!");
    else if (!action_defined) compute_action();
    return action[i]; }
  
  double get_freq(const int i) {
    if (i>dof) bomb("frequencies not defined!");
    if (!freq_defined) compute_freq();
    return freq[i]; }
  
  double get_angle(const int i, const double time);
  
  // Specific functions
  
  void new_orbit(double Energy, double kappa, double beta = 0.0);
  
  void set_biorth(AxiSymBiorth &type, int ll, int max, int recur = 0) 
  { nmax = max; l = ll; biorth = &type; RECUR = recur; 
  biorth_defined = false; }
  
  void set_numerical_params(int RECS=64, int NMAX=40, int NBSCT=3)
  { recs=RECS; nmax=NMAX; nbsct=NBSCT;}
  
  double pot_trans(int l1, int l2, double (*func)(double) );
  double pot_trans(int l1, int l2, int n );
  void pot_trans(int l1, int l2, Vector& t);
  
  // Access to underlying grid for pot_trans
  struct ANGLE_GRID * get_angle_grid(void) { return &angle_grid; }
  
  // Safe access
  
  double Energy(void) { return energy; }
  double AngMom(void) { 
    if (!freq_defined) compute_freq();
    return jmax*kappa; }
  double Jmax(void) {
    if (!freq_defined) compute_freq();
    return jmax; }
  double Kappa(void) { return kappa; }
  double Beta(void) { return beta; }
  double peri(void) { 
    if (!freq_defined) compute_freq();
    return r_peri; }
  double apo(void) {
    if (!freq_defined) compute_freq();
    return r_apo; }
  double circ(void) {
    if (!freq_defined) compute_freq();
    return r_circ; }
  
  AxiSymModel& modl(void) { return *model; }
  
  AxiSymBiorth& orth(void) { return *biorth; }
  
};

#include <massmodel.h>

#endif
