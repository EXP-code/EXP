// This may look like C code, but it is really -*- C++ -*-

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <math.h>

#include <CVector.h>
#include <clinalg.h>

#include <biorth.h>
#include <model3d.h>
#include <isothermal.h>
#include <hernquist.h>
#include <gaussQ.h>

#include <Perturbation.H>


int Perturbation::L1MAX = 6;
int Perturbation::NINT = 500;
bool Perturbation::selfgrav = false;
bool Perturbation::verbose = false;


Perturbation::Perturbation(int NMAX)
{
  ID = "Perturbation";
  computed = false;
  omega_computed = false;

  OMPI = default_OMPI;
  l = 2;

  nmax = NMAX;

  Epts = default_Epts;
  Kpts = default_Kpts;
  Recs = default_RECS;
}

Perturbation::~Perturbation()
{
  // Default stuff
}

void Perturbation::set_respmat_parameters(int epts, int kpts, int recs)
{
  Epts = epts;
  Kpts = kpts;
  Recs = recs;
}

void Perturbation::compute_perturbation(AxiSymModel *halo_model, 
				      AxiSymBiorth *halo_ortho,
				      CVector& total, CVector& pp )
{
  model = halo_model;
  biorth = halo_ortho;

  CVector Response(1, nmax);
  CVector Response_total(1, nmax);

  compute_coefficients();

  CMatrix ident(1, nmax, 1, nmax);
  ident.zero();
  for (int i=1; i<=ident.getnrows(); i++) ident[i][i] = 1.0;

  make_response();

  if (selfgrav) {
    total = current * asymp * bcoef;
    pp = current * Re(asymp) * bcoef;
  }
  else {
    total = asymp * bcoef;
    pp = Re(asymp) * bcoef;
  }
}

//
// Linearly interpolate on grid of response matrices
//

void Perturbation::make_response()
{
  if (computed) return;
				// Get frequency
  compute_omega();
  
  asymp = CMatrix(1, nmax, 1, nmax);
  asymp.zero();
  for (int i=1; i<=asymp.getnrows(); i++) asymp[i][i] = 1.0;
  working = asymp;

  computed = true;
}

