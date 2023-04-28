// This may look like C code, but it is really -*- C++ -*-

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>

#include <biorth.H>
#include <model3d.H>
#include <isothermal.H>
#include <hernquist_model.H>
#include <gaussQ.H>

#include <Perturbation.H>


int Perturbation::L1MAX = 6;
int Perturbation::NINT = 500;
bool Perturbation::selfgrav = false;
bool Perturbation::verbose = false;
double Perturbation::default_OMPI=0.03;


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

void Perturbation::compute_perturbation(AxiSymModPtr halo_model, 
					AxiSymBioPtr halo_ortho,
					Eigen::VectorXcd& total,
					Eigen::VectorXcd& pp )
{
  model  = halo_model;
  biorth = halo_ortho;

  Eigen::MatrixXcd Response(nmax, nmax);
  Eigen::MatrixXcd Response_total(nmax, nmax);

  compute_coefficients();

  Eigen::MatrixXcd ident = Eigen::MatrixXcd::Identity(nmax, nmax);

  make_response();

  if (selfgrav) {
    total = current * asymp * bcoef;
    pp = current * asymp.real() * bcoef;
  }
  else {
    total = asymp * bcoef;
    pp = asymp.real() * bcoef;
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
  
  asymp = Eigen::MatrixXcd::Identity(nmax, nmax);
  working = asymp;

  computed = true;
}

