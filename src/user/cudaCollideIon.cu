// -*- C++ -*-

#include <Component.H>
#include <CollideIon.H>
#include <UserTreeDSMC.H>
#include <cudaIon.cuH>
#include <cudaElastic.cuH>
#include <EXPException.H>

#include <curand.h>
#include <curand_kernel.h>

enum cudaInterTypes { 
  nothing    = 0,
  neut_neut  = 1,
  neut_elec  = 2,
  neut_prot  = 3,
  ion_ion    = 4,
  ion_elec   = 5,
  free_free  = 6,
  col_excite = 7,
  col_ionize = 8,
  recombine  = 9
};

const int maxAtomicNumber = 15;
__constant__ cuFP_t cuda_atomic_weights[maxAtomicNumber], cuFloorEV, cuEsu;
__constant__ cuFP_t cuVunit, cuMunit, cuTunit, cuLunit, cuEunit;
__constant__ cuFP_t cuAmu, cuEV, cuLogL, cuCrossfac;
__constant__ bool   cuMeanKE, cuMeanMass, cuNewRecombAlg;

const int coulSelNumT = 2000;
__constant__ cuFP_t coulSelA[coulSelNumT];
__constant__ cuFP_t coulSelTau_i, coulSelTau_f, coulSelTau_z, coulSelDel;

// Coulombic algorithm initialization
//
double cuCA_f(double x)
{
  return 1.0/tanh(x) - 1.0/x;
}

double cuCA_deriv(double x)
{
  if (x>100.0) return 1.0/(x*x);
  double z = sinh(x);
  return 1.0/(x*x) - 1.0/(z*z);
}

double cuCA_func(cuFP_t tau, cuFP_t x)
{
  const double tol = 1.0e-12;
  const int itmax  = 1000;
  double fac       = exp(-2.0*tau), corr;
  
  for (int i=0; i<itmax; i++) {
    corr  = (cuCA_f(x) - fac)/cuCA_deriv(x);
    x    -= corr;
    if (fabs(corr)<tol) break;
  }
  
  return x;
}

// Link static parameters from CollideIon.cc (decide how to link these
// later)
//
static double FloorEv      = 0.05;
static bool   newRecombAlg = false;

void CollideIon::cuda_initialize()
{
  static bool done = false;
  if (done) return;
  done = true;

  thrust::host_vector<cuIonElement> elems;

  int minSp = std::numeric_limits<int>::max();

  for (auto s : SpList) {
    speciesKey k = s.first;
    int Z = k.first;
    int C = k.second;
    // Scan
    bool found = false;
    for (auto & E : cuIonElem) {
      if (E.Z == Z and E.C == C) {
	E.I = s.second;
	minSp = std::min<int>(minSp, s.second);
	elems.push_back(E);
	break;
      }
    }
    if (not found) {
      std::cout << "CollideIon::cuda_initialize: [Z, C] = ["
		<< Z << ", " << C << "] not found" << std::endl;
    }
  }

  for (auto & E : elems) E.I -= minSp;

  cuElems = elems;

  // Coulombic velocity selection
  //
  cuFP_t tau_i = 0.0001, tau_f = 4.0, tau_z = 40.0;
  std::vector<cuFP_t> hA(coulSelNumT);

  cuFP_t del = (log(tau_f) - log(tau_i))/(coulSelNumT-1);
  cuFP_t A   = 0.5/tau_i, T;
  for (int i=0; i<coulSelNumT; i++) {
    T = tau_i*exp(del*i);
    hA[i] = A = cuCA_func(T, A);
  }

  cuda_safe_call(cudaMemcpyToSymbol(coulSelTau_i, &tau_i, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying coulSelTau_i");

  cuda_safe_call(cudaMemcpyToSymbol(coulSelTau_f, &tau_f, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying coulSelTau_f");

  cuda_safe_call(cudaMemcpyToSymbol(coulSelTau_z, &tau_z, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying coulSelTau_z");

  cuda_safe_call(cudaMemcpyToSymbol(coulSelDel, &del, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying coulSelDel");

  cuda_safe_call(cudaMemcpyToSymbol(coulSelA, &hA[0], coulSelNumT*sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying coulSelA");
}

void CollideIon::cuda_atomic_weights_init()
{
  cudaElasticInit();

  std::vector<cuFP_t> weights(maxAtomicNumber);

  weights[0]  = 0.000548579909; // Mass of electron
  weights[1]  = 1.0079;	       // Hydrogen
  weights[2]  = 4.0026;	       // Helium
  weights[3]  = 6.941;	       // Lithum
  weights[4]  = 9.0122;	       // Beryllium
  weights[5]  = 10.811;	       // Boron
  weights[6]  = 12.011;	       // Carbon
  weights[7]  = 14.007;	       // Nitrogen
  weights[8]  = 15.999;	       // Oxygen
  weights[9]  = 18.998;	       // Florine
  weights[10] = 20.180;	       // Neon
  weights[11] = 22.990;	       // Sodium
  weights[12] = 24.305;	       // Magnesium
  weights[13] = 26.982;	       // Aluminium
  weights[14] = 28.085;	       // Silicon

  cuda_safe_call(cudaMemcpyToSymbol(cuda_atomic_weights, &weights[0], sizeof(cuFP_t)*maxAtomicNumber), 
		 __FILE__, __LINE__, "Error copying cuda_atomic_weights");

  cuFP_t v = UserTreeDSMC::Vunit;
  cuda_safe_call(cudaMemcpyToSymbol(cuVunit, &v, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying cuVunit");
  v = UserTreeDSMC::Lunit;
  cuda_safe_call(cudaMemcpyToSymbol(cuLunit, &v, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying cuLunit");
  v = UserTreeDSMC::Munit;
  cuda_safe_call(cudaMemcpyToSymbol(cuMunit, &v, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying cuMunit");
  v = UserTreeDSMC::Tunit;
  cuda_safe_call(cudaMemcpyToSymbol(cuTunit, &v, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying cuTunit");
  v = UserTreeDSMC::Eunit;
  cuda_safe_call(cudaMemcpyToSymbol(cuEunit, &v, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying cuEunit");
  v = amu;
  cuda_safe_call(cudaMemcpyToSymbol(cuAmu, &v, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying cuAmu");
  v = eV;
  cuda_safe_call(cudaMemcpyToSymbol(cuEV, &v, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying cuEV");
  v = logL;
  cuda_safe_call(cudaMemcpyToSymbol(cuLogL, &v, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying cuLogL");
  v = crossfac;
  cuda_safe_call(cudaMemcpyToSymbol(cuCrossfac, &v, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying cuCrossfac");
  v = FloorEv;
  cuda_safe_call(cudaMemcpyToSymbol(cuFloorEV, &v, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying cuFloorEV");
  v = esu;
  cuda_safe_call(cudaMemcpyToSymbol(cuEsu, &v, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying cuEsu");

  cuda_safe_call(cudaMemcpyToSymbol(cuMeanMass, &MeanMass, sizeof(bool)), 
		 __FILE__, __LINE__, "Error copying cuMeanMass");

  cuda_safe_call(cudaMemcpyToSymbol(cuNewRecombAlg, &newRecombAlg, sizeof(bool)), 
		 __FILE__, __LINE__, "Error copying cuNewRecombAlg");
}  



// CURAND initialization
//
__global__ void initCurand(dArray<curandState> state, unsigned long seed)
{
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;
  if (tid < state._s) {
    curand_init(seed, tid, 0, &state._v[tid]);
  }
}

__global__ void cellInitKernel(dArray<cudaParticle> in,
			       dArray<cuFP_t> Ivel2,
			       dArray<cuFP_t> Evel2,
			       dArray<cuFP_t> PiProb,
			       dArray<cuFP_t> ABrate,
			       dArray<cuFP_t> volC,
			       dArray<cuFP_t> tauC,
			       dArray<cuFP_t> selC,
			       dArray<int>    cellI,
			       dArray<int>    cellN,
			       dArray<cuIonElement> elems,
			       int epos, unsigned int stride)
{
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  const cuFP_t dfac = cuMunit/cuAmu / (cuLunit*cuLunit*cuLunit);

  for (int s=0; s<stride; s++) {

    int c = tid*stride + s;

    if (c < cellI._s) {

      int nbods = cellN._v[c];

      cuFP_t massP = 0.0, numbP = 0.0, massE = 0.0;
      cuFP_t evel2 = 0.0, ivel2 = 0.0, numQ2 = 0.0;
      cuFP_t densI = 0.0, densQ = 0.0, densE = 0.0;
      cuFP_t molPP = 0.0, meanM = 0.0;

      for (size_t i=0; i<nbods; i++) {
	
	// The particle
	cudaParticle & p = in._v[i + cellI._v[c]];

	// Mass per cell
	massP += p.mass;
	
	// Mass-weighted trace fraction
	cuFP_t ee = 0.0;
	for (int k=0; k<elems._s; k++) {
	  cuIonElement& E = elems._v[k];
	  cuFP_t ff = p.datr[E.I];
	  cuFP_t ww = ff/cuda_atomic_weights[E.Z];
	  cuFP_t qq = E.C - 1;
	  // Mean number
	  numbP += p.mass * ww;
	  // Electron fraction
	  ee += ww * qq;
	  // Ion dens
	  densI += p.mass * ww;
	  // Charge
	  densE += p.mass * ww * qq;
	  // Charge squared
	  numQ2 += p.mass * ww * qq*qq;
	  if (E.C) densQ += p.mass * ww;
	}
	
	cuFP_t eVel2 = 0.0, iVel2 = 0.0;
	for (int l=0; l<3; l++) {
	  cuFP_t ve  = p.datr[epos+l];
	  eVel2 += ve*ve;
	  cuFP_t vi  = p.vel[l];
	  iVel2 += vi*vi;
	}
	
	evel2 += p.mass * ee * eVel2;
	ivel2 += p.mass * iVel2;
	massE += p.mass * ee;
      }
  
      if (numbP>0.0) meanM       = massP/numbP;
      if (massP>0.0) Ivel2._v[c] = ivel2/massP;
      if (massE>0.0) Evel2._v[c] = evel2/massE;
      if (densQ>0.0) numQ2      /= densQ;
      if (densI>0.0) molPP       = massP/densI;
      
      cuFP_t ddfac = dfac/volC._v[c];

      densI *= ddfac;
      densQ *= ddfac;
      densE *= ddfac;

      // Compute per channel Coulombic probabilities
      //
      // Ion probabilities
      //
      cuFP_t muii = meanM/2.0;
      cuFP_t muee = cuda_atomic_weights[0]/2.0;
      cuFP_t muie = cuda_atomic_weights[0] * meanM/(cuda_atomic_weights[0] + meanM);
      
      // Ion-Ion
      PiProb._v[c*4 + 0] =
	densQ +
	densE * pow(2.0*Ivel2._v[c], 1.5) * muii*muii /
	(pow(Ivel2._v[c] + Evel2._v[c], 1.5) * muie*muie * numQ2);
      //                                               ^
      //                                               |
      // The density is weighted by q^2 for each species
      
      // Ion-Electron
      PiProb._v[c*4 + 1] =
	densQ * pow(Ivel2._v[c] + Evel2._v[c], 1.5) * muie*muie * numQ2 /
	(pow(2.0*Ivel2._v[c], 1.5) * muii*muii) +
	densE ;
      
      // Electron-Ion
      PiProb._v[c*4 + 2] =
	densQ +
	densE * pow(Ivel2._v[c] + Evel2._v[c], 1.5) * muie*muie /
	(pow(2.0*Evel2._v[c], 1.5) * muee*muee * numQ2);
      
      // Electron-Electron
      PiProb._v[c*4 + 3] =
	densQ * pow(2.0*Evel2._v[c], 1.5) * muee*muee * numQ2 /
	(pow(Ivel2._v[c] + Evel2._v[c], 1.5) * muie*muie) +
	densE;
      
      // Rate coefficients
      ABrate._v[c*4 + 0] = 2.0*M_PI * PiProb._v[c*4 + 0] * cuLogL * pow(numQ2*numQ2, 2.0);
      
      ABrate._v[c*4 + 1] = 2.0*M_PI * PiProb._v[c*4 + 1] * cuLogL * pow(numQ2, 2.0);
      
      ABrate._v[c*4 + 2] = 2.0*M_PI * PiProb._v[c*4 + 2] * cuLogL * pow(numQ2, 2.0);
      
      ABrate._v[c*4 + 3] = 2.0*M_PI * PiProb._v[c*4 + 3] * cuLogL ;

      // Selection scaling based on plasma rate
      cuFP_t dens = massP * ddfac;
      cuFP_t tot_pairs =  dens * nbods / molPP * 
	(1.0/PiProb._v[c*4+0] + 1.0/PiProb._v[c*4+1] + 1.0/PiProb._v[c*4+2]);
      cuFP_t Prob = dens * cuMunit/cuAmu * sqrt(Ivel2._v[c]) * tauC._v[c];

      selC._v[c] = Prob/tot_pairs;

    } // END: cell

  } // END: stride

}

__global__ void crossSectionKernel(dArray<cudaParticle> in,
				   dArray<curandState> randS,
				   dArray<cuFP_t> cross,
				   dArray<cuFP_t> delph,
				   dArray<uchar3> xspc1,
				   dArray<uchar3> xspc2,
				   dArray<cudaInterTypes> xtype,
				   dArray<cuFP_t> xctot,
				   dArray<int>    i1,
				   dArray<int>    i2,
				   dArray<int>    cc,
				   dArray<cuFP_t> volC,
				   dArray<cuFP_t> Ivel2,
				   dArray<cuFP_t> Evel2,
				   dArray<cuFP_t> PiProb,
				   dArray<cuFP_t> ABrate,
				   dArray<int>    flagI, 
				   dArray<cuFP_t> xsc_H,
				   dArray<cuFP_t> xsc_He,
				   dArray<cuFP_t> xsc_pH,
				   dArray<cuFP_t> xsc_pHe,
				   dArray<cuIonElement> elems,
				   int numxc,
				   int epos,
				   unsigned int stride)
{
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;
  const int Nsp = elems._s;

  for (int s=0; s<stride; s++) {

    int n = tid*stride + s;

    if (n < i1._s) {

      cudaParticle& p1 = in._v[i1._v[n]];
      cudaParticle& p2 = in._v[i2._v[n]];
      int          cid = cc._v[n];
      
      // Zero the xtype
      for (int j=0; j<numxc; j++) xtype._v[n*numxc+j] = nothing;

      cuFP_t Eta1=0.0, Eta2=0.0, Mu1=0.0, Mu2=0.0, Sum1=0.0, Sum2=0.0;

      xctot._v[n] = 0.0;	// Total cross section accumulator

      int J = 0;		// Cross section position counter

      for (int k=0; k<Nsp; k++) {
	cuIonElement& E = elems._v[k];

				// Number fraction of ions
	cuFP_t one = p1.datr[E.I] / cuda_atomic_weights[E.Z];
	cuFP_t two = p2.datr[E.I] / cuda_atomic_weights[E.Z];

				// Electron number fraction
	Eta1 += one * (E.C - 1);
	Eta2 += two * (E.C - 1);

	Sum1 += one;
	Sum2 += two;
      }
				// The number of electrons per particle
      Eta1 /= Sum1;
      Eta2 /= Sum2;
				// The molecular weight
      Mu1 = 1.0/Sum1;
      Mu2 = 1.0/Sum2;


      // Number of atoms in each super particle
      //
      // cuFP_t N1 = p1.mass*cuMunit/(Mu1*cuAmu);
      // cuFP_t N2 = p2.mass*cuMunit/(Mu2*cuAmu);

      cuFP_t vel   = 0.0;
      cuFP_t eVel0 = 0.0, eVel1 = 0.0, eVel2 = 0.0;
      cuFP_t gVel0 = 0.0, gVel1 = 0.0, gVel2 = 0.0;
      for (int k=0; k<3; k++) {
	cuFP_t del = p1.vel[k] - p2.vel[k];
	vel += del*del;
	del = p1.vel[k] - p2.vel[k];

	cuFP_t rvel0 = p1.datr[epos+k] - p2.datr[epos+k];
	cuFP_t rvel1 = p1.datr[epos+k] - p2.vel[k];
	cuFP_t rvel2 = p2.datr[epos+k] - p1.vel[k];

	eVel0 += rvel0*rvel0;
	eVel1 += rvel1*rvel1;
	eVel2 += rvel2*rvel2;

	// Scaled electron relative velocity
	if (cuMeanMass) {
	  rvel0 = p1.datr[epos+k]*sqrt(Eta1) - p2.datr[epos+k]*sqrt(Eta2);
	  rvel1 = p1.datr[epos+k]*sqrt(Eta1) - p2.vel[k];
	  rvel2 = p2.datr[epos+k]*sqrt(Eta2) - p1.vel[k];
	  
	  gVel0 += rvel0*rvel0;
	  gVel1 += rvel1*rvel1;
	  gVel2 += rvel2*rvel2;
	}
      }

      // Energy available in the center of mass of the atomic collision
      //
      vel   = sqrt(vel)   * cuVunit;
      eVel0 = sqrt(eVel0) * cuVunit;
      eVel1 = sqrt(eVel1) * cuVunit;
      eVel2 = sqrt(eVel2) * cuVunit;

      // Pick scaled relative velocities for mean-mass algorithm
      if (cuMeanMass) {
	gVel0 = sqrt(gVel0) * cuVunit / vel;
	gVel1 = sqrt(gVel1) * cuVunit / vel;
	gVel2 = sqrt(gVel2) * cuVunit / vel;
      }
      // Pick true relative velocity for all other algorithms
      else {
	gVel0 = eVel0 / vel;
	gVel1 = eVel1 / vel;
	gVel2 = eVel2 / vel;
      }

      cuFP_t  m1 = Mu1 * cuAmu;
      cuFP_t  m2 = Mu2 * cuAmu;
      cuFP_t  me = cuda_atomic_weights[0] * cuAmu;
      cuFP_t mu0 = m1 * m2 / (m1 + m2);
      cuFP_t mu1 = m1 * me / (m1 + me);
      cuFP_t mu2 = m2 * me / (m2 + me);

      cuFP_t kEi = 0.5 * mu0 * vel * vel;

      cuFP_t facE = 0.5  * cuAmu * cuVunit * cuVunit / cuEV;
      cuFP_t Eion = facE * Ivel2._v[n] * 0.5*(Mu1 + Mu2);
      cuFP_t Eelc = facE * Evel2._v[n] * cuda_atomic_weights[0];
      cuFP_t kEe1 = 0.5  * mu1 * eVel2*eVel2;
      cuFP_t kEe2 = 0.5  * mu2 * eVel1*eVel1;
      // cuFP_t kEee = 0.25 * me  * eVel0*eVel0;

      for (int k1=0; k1<Nsp; k1++) {

	cuIonElement& elem = elems._v[k1];

	int Z = elem.Z;
	int C = elem.C;
	int P = elem.C - 1;
	int I = elem.I;

	cuFP_t fac1 = p1.datr[I] / cuda_atomic_weights[Z] / Sum1;
	cuFP_t fac2 = p2.datr[I] / cuda_atomic_weights[Z] / Sum2;

	for (int kk=0; kk<Nsp; kk++) {
	  
	  cuIonElement& elem2 = elems._v[kk];
	  
	  int ZZ = elem2.Z;
	  int CC = elem2.C;
	  int PP = elem2.C - 1;
	  int II = elem2.I;

	  //--------------------------------------------------
	  // Particle 1 interacts with Particle 2
	  //--------------------------------------------------
	  
	  cuFP_t cfac = p1.datr[I] * p2.datr[II];

	  //-------------------------------
	  // *** Both particles neutral
	  //-------------------------------
    
	  if (P==0 and PP==0) {
	    
	    // Geometric cross sections based on
	    // atomic radius
	    cuFP_t crs = (cudaGeometric(Z) + cudaGeometric(ZZ)) * cfac;
	
	    // Double counting
	    if (Z == ZZ) crs *= 0.5;

	    cross._v[n*numxc+J] = crs*cuCrossfac;
	    xspc1._v[n*numxc+J] = make_uchar3(Z,  C,  I );
	    xspc2._v[n*numxc+J] = make_uchar3(ZZ, CC, II);
	    xtype._v[n*numxc+J] = neut_neut;
	    xctot._v[n]        += crs*cuCrossfac;
	    
	    J++;		// Increment interaction counter
	  }

	  // --------------------------------------
	  // *** Neutral atom-proton scattering
	  // --------------------------------------

	  cuFP_t crs1 = 0;

	  if (ZZ==1 and CC==2) {
	    if (Z==1 and P==0) crs1 = cudaElasticInterp(kEi, cuPH_Emin,  cuH_H,   xsc_pH );
	    if (Z==2 and P==0) crs1 = cudaElasticInterp(kEi, cuPHe_Emin, cuPHe_H, xsc_pHe);
	    crs1 *= cuCrossfac * cfac;
	  }
	  
	  if (Z==1 and C==2) {
	    if (ZZ==1 and PP==0) crs1 = cudaElasticInterp(kEi, cuPH_Emin,  cuPH_H,  xsc_pH );
	    if (ZZ==2 and PP==0) crs1 = cudaElasticInterp(kEi, cuPHe_Emin, cuPHe_H, xsc_pHe);
	    crs1 *= cuCrossfac * cfac;
	  }

	  if (crs1>0.0) {
	    cross._v[n*numxc+J] = crs1;
	    xspc1._v[n*numxc+J] = make_uchar3(Z,  C,  I );
	    xspc2._v[n*numxc+J] = make_uchar3(ZZ, CC, II);
	    xtype._v[n*numxc+J] = neut_prot;
	    xctot._v[n]        += crs1;
	    
	    J++;
	  }

	  // --------------------------------------
	  // *** Ion-ion scattering
	  // --------------------------------------
	  //
	  if (P>0 and PP>0) {
	    cuFP_t kEc  = cuMeanKE ? Eion : kEi;
	    kEc = 2.0*kEc > cuFloorEV ? 2.0*kEc*cuEV : cuFloorEV*cuEV;
	    cuFP_t afac = cuEsu*cuEsu/kEc * 1.0e7;
	    cuFP_t crs  = 2.0 * ABrate._v[cid*4+0] * afac*afac / PiProb._v[cid*4+0];

	    cross._v[n*numxc+J] = crs;
	    xspc1._v[n*numxc+J] = make_uchar3(Z,  C,  I );
	    xspc2._v[n*numxc+J] = make_uchar3(ZZ, CC, II);
	    xtype._v[n*numxc+J] = ion_ion;
	    xctot._v[n]        += crs;
	    
	    J++;
	  }

	} // End of inner species loop

	// --------------------------------------
	// *** Neutral atom-electron scattering
	// --------------------------------------
    
	// Particle 1 ION, Particle 2 ELECTRON
	//
	if (Z<=2 and P==0 and Eta2>0.0) {
	  cuFP_t crs = 0.0;
	  if (Z==1) crs = cudaElasticInterp(kEe1, cuH_Emin, cuH_H, xsc_H) * gVel2 * Eta2 *
		      cuCrossfac * fac1;
	  if (Z==2) crs = cudaElasticInterp(kEe1, cuHe_Emin, cuHe_H, xsc_He) * gVel2 * Eta2 *
		      cuCrossfac * fac1;

	  if (crs>0.0) {
	    cross._v[n*numxc+J] = crs;
	    xspc1._v[n*numxc+J] = make_uchar3(Z, C,   I);
	    xspc2._v[n*numxc+J] = make_uchar3(0, 0, 255);
	    xtype._v[n*numxc+J] = neut_elec;
	    xctot._v[n]        += crs;
	    
	    J++;
	  }
	}

	// Particle 2 ION, Particle 1 ELECTRON
	//
	if (P==0 and Eta1>0.0) {
	  cuFP_t crs = 0.0;
	  if (Z==1) crs = cudaElasticInterp(kEe2, cuH_Emin, cuH_H, xsc_H) * gVel2 * Eta2 *
		      cuCrossfac * fac1;
	  if (Z==2) crs = cudaElasticInterp(kEe2, cuHe_Emin, cuHe_H, xsc_He) * gVel2 * Eta2 *
		      cuCrossfac * fac1;

	  if (crs>0.0) {
	    cross._v[n*numxc+J] = crs;
	    xspc1._v[n*numxc+J] = make_uchar3(0, 0, 255);
	    xspc2._v[n*numxc+J] = make_uchar3(Z, C,   I);
	    xtype._v[n*numxc+J] = neut_elec;
	    xctot._v[n]        += crs;
	    
	    J++;
	  }
	}

	// --------------------------------------
	// *** Ion-electron scattering
	// --------------------------------------
	//
	if (P>0 and Eta2>0) {
	  // Particle 1 ION, Particle 2 ELECTRON
	  cuFP_t crs  = 0.0;
	  cuFP_t kEc  = cuMeanKE ? Eelc : kEe1;
	  kEc = 2.0*kEc > cuFloorEV ? 2.0*kEc*cuEV : cuFloorEV*cuEV;
	  cuFP_t afac = cuEsu*cuEsu/kEc * 1.0e7;
	    
	  crs = 2.0 * ABrate._v[cid*4+1] * afac*afac * gVel2 / PiProb._v[cid*4+1];
	  
	  cross._v[n*numxc+J] = crs;
	  xspc1._v[n*numxc+J] = make_uchar3(Z, C,   I);
	  xspc2._v[n*numxc+J] = make_uchar3(0, 0, 255);
	  xtype._v[n*numxc+J] = ion_elec;
	  xctot._v[n]        += crs;
	  
	  J++;
	  
	  // Particle 2 ION, Particle 1 ELECTRON
	  crs  = 0.0;
	  kEc  = cuMeanKE ? Eelc : kEe2;
	  afac = cuEsu*cuEsu/(2.0*kEc > cuFloorEV ? 2.0*kEc : cuFloorEV) * cuEV * 1.0e7;

	  crs = 2.0 * ABrate._v[cid*4+2] * afac*afac * gVel1 / PiProb._v[cid*4+2];

	  cross._v[n*numxc+J] = crs;
	  xspc1._v[n*numxc+J] = make_uchar3(0, 0, 255);
	  xspc2._v[n*numxc+J] = make_uchar3(Z, C,   I);
	  xtype._v[n*numxc+J] = ion_elec;
	  xctot._v[n]        += crs;
	  
	  J++;
	} // end: ion-electron scattering


	//-------------------------------
	// *** Free-free
	//-------------------------------
      
	if (Eta1>0.0 and Eta2>0.0) {

	  // Particle 1 ION, Particle 2 ELECTRON
	  cuFP_t ke = kEe1 > cuFloorEV ? kEe1 : cuFloorEV, ff, ph;
	  cuFP_t rn;
#if cuREAL == 4
	  rn = curand_uniform(&randS._v[n]);
#else
	  rn = curand_uniform_double(&randS._v[n]);
#endif
	  computeFreeFree(ke, rn, ph, ff, elem);

	  cuFP_t crs  = gVel2 * Eta2 * ff * fac1;
	
	  if (std::isinf(crs)) crs = 0.0; // Sanity check
	
	  if (crs>0.0) {

	    cross._v[n*numxc+J] = crs;
	    delph._v[n*numxc+J] = ph;
	    xspc1._v[n*numxc+J] = make_uchar3(Z, C,   I);
	    xspc2._v[n*numxc+J] = make_uchar3(0, 0, 255);
	    xtype._v[n*numxc+J] = free_free;
	    xctot._v[n]        += crs;
	  
	    J++;
	  }

	  // Particle 2 ION, Particle 1 ELECTRON
	  ke = kEe2 > cuFloorEV ? kEe2 : cuFloorEV;

#if cuREAL == 4
	  rn = curand_uniform(&randS._v[n]);
#else
	  rn = curand_uniform_double(&randS._v[n]);
#endif
	  computeFreeFree(ke, rn, ph, ff, elem);

	  crs = gVel1 * Eta1 * ff * fac2;
	  
	  if (crs>0.0) {

	    cross._v[n*numxc+J] = crs;
	    delph._v[n*numxc+J] = ph;
	    xspc1._v[n*numxc+J] = make_uchar3(0, 0, 255);
	    xspc2._v[n*numxc+J] = make_uchar3(Z, C,   I);
	    xtype._v[n*numxc+J] = free_free;
	    xctot._v[n]        += crs;
	  
	    J++;
	  }
	}
	// end: free-free 

	//-------------------------------
	// *** Collisional excitation
	//-------------------------------
    
	// Particle 1 nucleus has BOUND ELECTRON, Particle 2 has FREE ELECTRON
	//
	//  +--- Charge of the current subspecies
	//  |
	//  |       +--- Electron fraction of partner
	//  |       |
	//  V       V
	if (P<Z and Eta2>0.0) {
	  cuFP_t ke = kEe1 > cuFloorEV ? kEe1 : cuFloorEV, ph, xc;
	  computeColExcite(ke, ph, xc, elem);

	  cuFP_t crs = gVel2 * Eta2 * xc * fac1;
      
	  if (crs > 0.0) {
	    cross._v[n*numxc+J] = crs;
	    delph._v[n*numxc+J] = ph;
	    xspc1._v[n*numxc+J] = make_uchar3(Z, C,   I);
	    xspc2._v[n*numxc+J] = make_uchar3(0, 0, 255);
	    xtype._v[n*numxc+J] = col_excite;
	    xctot._v[n]        += crs;
	  
	    J++;
	  }
	}

    
	// Particle 2 nucleus has BOUND ELECTRON, Particle 1 has FREE ELECTRON
	//
	//  +--- Charge of the current subspecies
	//  |
	//  |       +--- Electron fraction of partner
	//  |       |
	//  V       V
	if (P<Z and Eta1>0) {

	  cuFP_t ke = kEe2 > cuFloorEV ? kEe2 : cuFloorEV, ph, xc;
	  computeColExcite(ke, ph, xc, elem);

	  cuFP_t crs = gVel1 * Eta1 * xc * fac2;
      
	  if (crs > 0.0) {
	    cross._v[n*numxc+J] = crs;
	    delph._v[n*numxc+J] = ph;
	    xspc1._v[n*numxc+J] = make_uchar3(0, 0, 255);
	    xspc2._v[n*numxc+J] = make_uchar3(Z, C,   I);
	    xtype._v[n*numxc+J] = col_excite;
	    xctot._v[n]        += crs;
	  
	    J++;
	  }
	}	  
	// end: colexcite

	//-------------------------------
	// *** Ionization cross section
	//-------------------------------
	
	// Particle 1 nucleus has BOUND ELECTRON, Particle 2 has FREE ELECTRON
	//
	//  +--- Charge of the current subspecies
	//  |
	//  |       +--- Electron fraction of partner
	//  |       |
	//  V       V
	if (P<Z and Eta2>0) {
      
	  cuFP_t ke = kEe1 > cuFloorEV ? kEe1 : cuFloorEV, xc;
	  computeColIonize(ke, xc, elem);

	  cuFP_t crs = gVel2 * Eta2 * xc * fac1;
      
	  if (crs > 0.0) {
	    cross._v[n*numxc+J] = crs;
	    xspc1._v[n*numxc+J] = make_uchar3(Z, C,   I);
	    xspc2._v[n*numxc+J] = make_uchar3(0, 0, 255);
	    xtype._v[n*numxc+J] = col_ionize;
	    xctot._v[n]        += crs;
	  
	    J++;
	  }
	}
    
	// Particle 2 nucleus has BOUND ELECTRON, Particle 1 has FREE ELECTRON
	//
	//  +--- Charge of the current subspecies
	//  |
	//  |       +--- Electron fraction of partner
	//  |       |
	//  V       V
	if (P<Z and Eta1) {
	  
	  cuFP_t ke = kEe2 > cuFloorEV ? kEe2 : cuFloorEV, xc;
	  computeColIonize(ke, xc, elem);

	  cuFP_t crs = gVel1 * Eta1 * xc * fac2;
      
	  if (crs > 0.0) {
	    cross._v[n*numxc+J] = crs;
	    xspc1._v[n*numxc+J] = make_uchar3(0, 0, 255);
	    xspc2._v[n*numxc+J] = make_uchar3(Z, C,   I);
	    xtype._v[n*numxc+J] = col_ionize;
	    xctot._v[n]        += crs;
	  
	    J++;
	  }
	}
	// end: ionize

	//-------------------------------
	// *** Radiative recombination
	//-------------------------------
	
	// The "new" algorithm uses the electron energy of the ion's
	// electron rather than the standard particle partner.
	//
	
	if (cuNewRecombAlg) {
	  
	  // Particle 1 is ION, Particle 2 has ELECTRON
	  //
	  //  +--- Ion charge
	  //  |
	  //  v
	  if (P>0) {
	    cuFP_t ke = kEe1 > cuFloorEV ? kEe1 : cuFloorEV, xc;
	    computeRadRecomb(ke, xc, elem);
	    
	    cuFP_t crs = gVel1 * Eta1 * xc * fac1;
	
	    if (crs > 0.0) {
	      cross._v[n*numxc+J] = crs;
	      xspc1._v[n*numxc+J] = make_uchar3(Z, C,   I);
	      xspc2._v[n*numxc+J] = make_uchar3(0, 0, 255);
	      xtype._v[n*numxc+J] = recombine;
	      xctot._v[n]        += crs;
	      
	      J++;
	    }
	  }
      
	  // Particle 2 is ION, Particle 1 has ELECTRON
	  //
	  //  +--- Ion charge
	  //  |
	  //  v
	  if (P>0) {
	    cuFP_t ke = kEe2 > cuFloorEV ? kEe2 : cuFloorEV, xc;
	    computeRadRecomb(ke, xc, elem);
	    
	    cuFP_t crs = gVel2 * Eta2 * xc * fac2;
	    
	    if (crs > 0.0) {
	      cross._v[n*numxc+J] = crs;
	      xspc1._v[n*numxc+J] = make_uchar3(0, 0, 255);
	      xspc2._v[n*numxc+J] = make_uchar3(Z, C,   I);
	      xtype._v[n*numxc+J] = recombine;
	      xctot._v[n] += crs;
	      
	      J++;
	    }
	  }
	} // end: new recomb algorithm
	else {
	  // Particle 1 is ION, Particle 2 has ELECTRON
	  //
	  //  +--- Charge of the current subspecies
	  //  |
	  //  |       +--- Electron fraction of partner
	  //  |       |
	  //  V       V
	  if (P>0 and Eta2>0.0) {
	    cuFP_t ke = kEe1 > cuFloorEV ? kEe1 : cuFloorEV, xc;
	    computeRadRecomb(ke, xc, elem);
	    
	    cuFP_t crs = gVel2 * Eta2 * xc * fac1;
	
	    if (crs > 0.0) {
	      cross._v[n*numxc+J] = crs;
	      xspc1._v[n*numxc+J] = make_uchar3(Z, C,   I);
	      xspc2._v[n*numxc+J] = make_uchar3(0, 0, 255);
	      xtype._v[n*numxc+J] = recombine;
	      xctot._v[n]        += crs;
	      
	      J++;
	    }
	  }

	  // Particle 2 is ION, Particle 1 has ELECTRON
	  //
	  //  +--- Charge of the current subspecies
	  //  |
	  //  |       +--- Electron fraction of partner
	  //  |       |
	  //  V       V
	  if (P>0 and Eta1>0.0) {

	    cuFP_t ke = kEe2 > cuFloorEV ? kEe2 : cuFloorEV, xc;
	    computeRadRecomb(ke, xc, elem);
	    
	    cuFP_t crs = gVel1 * Eta1 * xc * fac2;
	    
	    if (crs > 0.0) {
	      cross._v[n*numxc+J] = crs;
	      xspc1._v[n*numxc+J] = make_uchar3(0, 0, 255);
	      xspc2._v[n*numxc+J] = make_uchar3(Z, C,   I);
	      xtype._v[n*numxc+J] = recombine;
	      xctot._v[n]        += crs;
	      
	      J++;
	    }
	  }
	} // end: original recomb algorithm
      
      } // end: outer ionization state loop

    } // end: particle loop

  } // end: stride

}

// Return random 3d unit vector
__device__
void cudaUnitVector(cuFP_t *ret, curandState* state)
{
  enum UV {trig, gaus};		// Method choice
  static const UV uv(gaus);	// Pick the method
  
  if (uv == trig) {
#if cuREAL == 4
    cuFP_t cos_th = 1.0 - 2.0*curand_uniform(state);
    cuFP_t phi    = 2.0*M_PI*curand_uniform(state);
#else
    cuFP_t cos_th = 1.0 - 2.0*curand_uniform_double(state);
    cuFP_t phi    = 2.0*M_PI*curand_uniform_double(state);
#endif

    cuFP_t sin_th = sqrt(fabs(1.0 - cos_th*cos_th));

    ret[0] = sin_th*cos(phi);
    ret[1] = sin_th*sin(phi);
    ret[2] = cos_th;
  } else {
    cuFP_t nrm = 0.0;
    for (int i=0; i<3; i++) {
#if cuREAL == 4
      cuFP_t rn = curand_normal(state);
#else
      cuFP_t rn = curand_normal_double(state);
#endif
      ret[i] = rn;
      nrm += ret[i]*ret[i];
    }
    nrm = sqrt(nrm);
    for (int i=0; i<3; i++) ret[i] /= nrm;
  } 
}



__device__
cuFP_t cuCA_eval(cuFP_t tau)
{
  if      (tau <= coulSelTau_i) return 1.0/(2.0*tau);
  else if (tau >= coulSelTau_z) return 0.0;
  else if (tau >= coulSelTau_f) return 3.0*exp(-2.0*tau);
  else {
    cuFP_t lo  = log(coulSelTau_i);
    cuFP_t hi  = log(coulSelTau_f);
    cuFP_t lv  = log(tau), lo_t, hi_t;
    int indx   = floor( (lv - lo)/coulSelDel );

    if (indx<0) indx = 0;
    if (indx>coulSelNumT-2) indx = coulSelNumT - 2;

    lo_t = log(coulSelTau_i) + coulSelDel*indx;
    hi_t = log(coulSelTau_i) + coulSelDel*(indx+1);

    cuFP_t A = (hi_t - lv)/coulSelDel;
    cuFP_t B = (lv - lo_t)/coulSelDel;
    
    return A*coulSelA[indx-1] + B*coulSelA[indx];
  }
}

//! Select tau given random number U in [0,1)
__device__
cuFP_t cuCA_get(cuFP_t tau, cuFP_t U)
{
  cuFP_t A = cuCA_eval(tau);
  if (U<1.0e-14)
    return -1.0;
  else if (A<1.0e-10)
    return 2.0*U - 1.0;
  else if (A>40.0)
    return 1.0 + log(U)/A;
  else
    return log(exp(-A) + 2.0*U*sinh(A))/A;
}

// Return 3d Colombic scattering vector
//
__device__
void cudaCoulombVector(cuFP_t *ret, cuFP_t *rel, cuFP_t W1, cuFP_t W2, cuFP_t Tau,
		       curandState *state)
{
				// Normalize
  cuFP_t rel2 = 0.0;
  for (int i=0; i<3; i++) rel2 += rel[i]*rel[i];
  cuFP_t vfac = sqrt(rel2);
  if (vfac>0.0) for (int i=0; i<3; i++) rel[i] /= vfac;

				// Random generation
#if cuREAL == 4
  cuFP_t rn  = curand_uniform_double(state);
  cuFP_t phi = 2.0*M_PI*curand_uniform(state);
#else
  cuFP_t rn  = curand_uniform_double(state);
  cuFP_t phi = 2.0*M_PI*curand_uniform_double(state);
#endif


  cuFP_t fac    = (W1>W2 ? W1 : W2)/(W1<W2 ? W1 : W2);
  cuFP_t tau    = fac * Tau;
  cuFP_t cosx   = cuCA_get(tau, rn);
  cuFP_t sinx   = sqrt(fabs(1.0 - cosx*cosx));
  cuFP_t cosp   = cos(phi);
  cuFP_t sinp   = sin(phi);
  cuFP_t g_perp = sqrt(rel[1]*rel[1] + rel[2]*rel[2]);

				// Compute randomly-oriented
				// perpendicular vector
  cuFP_t h[3];
  if (g_perp>0.0) {
    h[0] = g_perp * cosp;
    h[1] = -(rel[1]*rel[0]*cosp + rel[2]*sinp)/g_perp;
    h[2] = -(rel[2]*rel[0]*cosp - rel[1]*sinp)/g_perp;
  } else {
    h[0] = 0.0;
    h[1] = cosp;
    h[2] = sinp;
  }
  
  for (int i=0; i<3; i++) rel[i] = rel[i]*cosx - h[i]*sinx;
}

__device__
void cudaDeferredEnergy
(
 cuFP_t E,
 cuFP_t m1,   cuFP_t m2,
 cuFP_t a,    cuFP_t b,
 cuFP_t *E1,  cuFP_t *E2
)
{
  if (m1 < 1.0) {
    E1[1] += a*E/(a + b);
    E2[0] += b*E/(a + b);
  }
  else if (m2 < 1.0) {
    E1[0] += a*E/(a + b);
    E2[1] += b*E/(a + b);
  }
  else {
    E1[0]  += a*E/(a + b);
    E2[0]  += b*E/(a + b);
  }
}

__device__
void cudaScatterTrace
(cuFP_t m1,    cuFP_t m2,
 cuFP_t eta1,  cuFP_t eta2,
 cuFP_t W1,    cuFP_t W2,
 cuFP_t *E1,   cuFP_t *E2,  cuFP_t &totE,
 cuFP_t *v1,   cuFP_t *v2,
 cuFP_t delE,  cuFP_t Tau,
 curandState *state, bool Coulombic
 )
{
  if (m1<1.0) m1 *= eta1;
  if (m2<1.0) m2 *= eta2;

  // Total effective mass in the collision
  //
  double mt = m1 + m2;

  // Reduced mass (atomic mass units)
  //
  double mu = m1 * m2 / mt;

  // Set COM frame
  //
  cuFP_t vcom[3], vrel[3];
  cuFP_t vi = 0.0;

  for (size_t k=0; k<3; k++) {
    vcom[k] = (m1*v1[k] + m2*v2[k])/mt;
    vrel[k] = v1[k] - v2[k];
    vi += vrel[k] * vrel[k];
  }
				// Energy in COM
  cuFP_t kE = 0.5*W2*mu*vi;
				// Energy reduced by loss
  cuFP_t vfac = 1.0;
  totE = kE - delE;

  // KE is positive
  //
  if (kE>0.0) {
    // More loss energy requested than available?
    //
    if (totE < 0.0) {
      // Add to energy bucket for these particles
      //
      cudaDeferredEnergy(-totE, m1, m2, W1, W2, E1, E2);
      totE = 0.0;
    }
    // Update the outgoing energy in COM
    vfac = sqrt(totE/kE);
  }
  // KE is zero (limiting case)
  //
  else {
    if (delE>0.0) {
      // Defer all energy loss
      //
      cudaDeferredEnergy(delE, m1, m2, W1, W2, E1, E2);
    } else {
      // Apply delE to COM
      //
      vi = -2.0*delE/(W1*mu);
    }
  }

  // Assign interaction energy variables
  //
  if (Coulombic)
    cudaCoulombVector(vrel, vrel, W1, W2, Tau, state);
  else
    cudaUnitVector(vrel, state);
  
  vi   = sqrt(vi);
  for (auto & v : vrel) v *= vi;
  //                         ^
  //                         |
  // Velocity in center of mass, computed from v1, v2 and adjusted
  // according to the inelastic energy loss
  //

  for (size_t k=0; k<3; k++) {
    v1[k] = vcom[k] + m2/mt*vrel[k] * vfac;
    v2[k] = vcom[k] - m1/mt*vrel[k] * vfac;
  }
      
} // END: cudaScatterTrace

__global__ void partInteractions(dArray<cudaParticle> in,
				 dArray<curandState> randS,
				 dArray<cuFP_t> cross,
				 dArray<cuFP_t> delph,
				 dArray<uchar3> xspc1,
				 dArray<uchar3> xspc2,
				 dArray<cudaInterTypes> xtype,
				 dArray<cuFP_t> xctot,
				 dArray<int>    i1,
				 dArray<int>    i2,
				 dArray<int>    cc,
				 dArray<cuFP_t> selC,
				 dArray<cuFP_t> Ivel2,
				 dArray<cuFP_t> Evel2,
				 dArray<cuFP_t> PiProb,
				 dArray<cuFP_t> ABrate,
				 dArray<cuFP_t> deltT,
				 dArray<int>    flagI, 
				 dArray<cuIonElement> elems,
				 cuFP_t spTau,
				 int numxc,
				 int epos,
				 unsigned int stride)
{
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;
  const int Nsp = elems._s;

  cuFP_t *FF1 = new cuFP_t [Nsp];
  cuFP_t *FF2 = new cuFP_t [Nsp];

  for (int s=0; s<stride; s++) {

    int n = tid*stride + s;

    if (n < i1._s) {

      cudaParticle& p1    = in._v[i1._v[n]];
      cudaParticle& p2    = in._v[i2._v[n]];
      curandState*  state = &randS._v[n];

      // Cell position in arrays
      //
      int cid = cc._v[n];
      
      cuFP_t PE[3] = {0, 0, 0}, EE[3] = {0, 0, 0};

      // Cross section in system units
      //
      cuFP_t totalXS = xctot._v[n];

      // Upscale rate
      //
      cuFP_t colCf = selC._v[cid] * totalXS;

      // Energy computation
      //
      cuFP_t Eta1=0.0, Eta2=0.0, Mu1=0.0, Mu2=0.0, Sum1=0.0, Sum2=0.0;
      cuFP_t kEe1 = 0.0, kEe2 = 0.0, kEi = 0.0, Mue = cuda_atomic_weights[0];

      // Electron KE
      //
      cuFP_t iE1 = 0.0, iE2 = 0.0;
      for (size_t k=0; k<3; k++) {
	iE1 += p1.datr[epos+k] * p1.datr[epos+k];
	iE2 += p2.datr[epos+k] * p2.datr[epos+k];

	cuFP_t rvel0 = p1.datr[epos+k] - p2.datr[epos+k];
	cuFP_t rvel1 = p1.datr[epos+k] - p2.vel[k];
	cuFP_t rvel2 = p2.datr[epos+k] - p1.vel[k];

	kEi  += rvel0*rvel0;
	kEe1 += rvel1*rvel1;
	kEe2 += rvel2*rvel2;
      }

      // Electron and molecular weight
      //
      for (int k=0; k<Nsp; k++) {
	cuIonElement& E = elems._v[k];

	FF1[k] = p1.datr[E.I];
	FF2[k] = p2.datr[E.I];

				// Number fraction of ions
	cuFP_t one = p1.datr[E.I] / cuda_atomic_weights[E.Z];
	cuFP_t two = p2.datr[E.I] / cuda_atomic_weights[E.Z];

				// Electron number fraction
	Eta1 += one * (E.C - 1);
	Eta2 += two * (E.C - 1);

	Sum1 += one;
	Sum2 += two;
      }
				// The number of electrons per particle
      Eta1 /= Sum1;
      Eta2 /= Sum2;
				// The molecular weight
      Mu1 = 1.0/Sum1;
      Mu2 = 1.0/Sum2;

      cuFP_t m1  = Mu1 * cuAmu;
      cuFP_t m2  = Mu2 * cuAmu;
      cuFP_t me  = Mue * cuAmu;
  
      cuFP_t mu0 = m1 * m2 / (m1 + m2);
      cuFP_t mu1 = m1 * me / (m1 + me);
      cuFP_t mu2 = me * m2 / (me + m2);

				// Finalize the KEs
      iE1 *= 0.5*p1.mass * cuda_atomic_weights[0];
      iE2 *= 0.5*p2.mass * cuda_atomic_weights[0];

      kEi  *= 0.5 * mu0;
      kEe1 *= 0.5 * mu1;
      kEe2 *= 0.5 * mu2;

      cuFP_t maxP = 0.0;
      cudaInterTypes maxT = nothing;

      struct Itype
      {
	unsigned char Z1, Z2, C1, C2, P1, P2, I1, I2;
	cuFP_t W1, W2, N1, N2;
      } IT;

      for (int J=0; J<numxc; J++) {
	cudaInterTypes T = xtype._v[n*numxc+J];
	if (T == nothing) break;
    
	cuFP_t XS          = cross._v[n*numxc+J];
	cuFP_t XE          = delph._v[n*numxc+J];
	double Prob        = XS/totalXS;

	if (Prob < 1.0e-14) continue;

	// Atomic number
	//
	IT.Z1 = xspc1._v[n*numxc+J].x;
	IT.C1 = xspc1._v[n*numxc+J].y;
	IT.I1 = xspc1._v[n*numxc+J].z;
	IT.Z2 = xspc2._v[n*numxc+J].x;
	IT.C2 = xspc2._v[n*numxc+J].y;
	IT.I2 = xspc2._v[n*numxc+J].z;

	// Traditional ionization state (e.g. C1=1 is neutral)
	//
	IT.P1 = IT.C1 - 1;
	IT.P2 = IT.C2 - 1;

	// Energy loss
	//
	cuFP_t dE = 0.0;

	// Number interacting atoms (255 is electron)
	//
	if (IT.I1<255) IT.W1 = p1.datr[IT.I1] / cuda_atomic_weights[IT.Z1];
	else           IT.W1 = Eta1;

	if (IT.I2<255) IT.W2 = p1.datr[IT.I1] / cuda_atomic_weights[IT.Z2];
	else           IT.W2 = Eta1;

	IT.N1 = IT.W1 * cuMunit / cuAmu;
	IT.N2 = IT.W2 * cuMunit / cuAmu;

	cuFP_t N0 = IT.N1 > IT.N2 ? IT.N2 : IT.N1;

	// Select the maximum probability channel
	//
	if (Prob > maxP) {
	  maxT  = T;
	  maxP  = Prob;
	}

	//-----------------------------
	// Parse each interaction type
	//-----------------------------

	if (T == neut_neut) {
	  PE[0] += Prob;
	}

	if (T == neut_elec) {
	  if (IT.I1) PE[1] += Prob;
	  else       PE[2] += Prob;
	}
	
	if (T == neut_prot) {
	  PE[0] += Prob;
	}

	if (T == ion_elec) {
	  if (IT.I1) PE[1] += Prob;
	  else       PE[2] += Prob;
	}
      
	// Upscale inelastic cross sections by ratio of cross-section
	// and plasma rate scattering
	//
	Prob *= colCf;
	
	if (T == free_free) {
	  
	  if (IT.I1) {		// Ion is p1

	    dE = XE * Prob;
	    
	    PE[1] += Prob;
	    EE[1] += dE;
	  } else {		// Ion is p2
	    dE = XE * Prob;
	    
	    PE[2] += Prob;
	    EE[2] += dE;
	  }
	}

	if (T == col_excite) {
	  
	  if (IT.I1) {		// Ion is p1
	    dE = XE * Prob;

	    PE[1] += Prob;
	    EE[1] += dE;
	  } else {		// Ion is p2
	    dE = XE * Prob;
	  
	    PE[2] += Prob;
	    EE[2] += dE;
	  }
	} // END: col_excite

	if (T == col_ionize) {

	  if (IT.I1) {		// Ion is p1
	    cuFP_t WW = Prob;	    

	    if (WW < FF1[J]) {
	      FF1[J]   -= WW;
	      FF1[J+1] += WW;
	    } else {
	      WW = FF1[J];
	      FF1[J]    = 0.0;
	      FF1[J+1] += WW;
	    }

	    Prob = WW;
	  
	    dE = XE * Prob;
	  
	    // The kinetic energy of the ionized electron is lost
	    // from the COM KE
	    //
	    cuFP_t Echg = iE1 * Prob / cuda_atomic_weights[IT.Z1];

	    // Energy for ionized electron comes from COM
	    dE += Echg * cuEunit / (N0*eV);
	    
	    PE[1] += Prob;
	    EE[1] += dE;
	    
	  } // END: ion-electron
	  else {		// Ion is p2
	    cuFP_t WW = Prob;

	    if (WW < FF2[J]) {
	      FF2[J]   -= WW;
	      FF2[J+1] += WW;
	    } else {
	      WW = FF2[J];
	      FF2[J]    = 0.0;
	      FF2[J+1] += WW;
	    }

	    Prob = WW;
	  
	    dE = XE * Prob;

	    // The kinetic energy of the ionized electron is lost
	    // from the COM KE
	    //
	    cuFP_t Echg = iE2 * Prob / cuda_atomic_weights[IT.Z2];

	    // Energy for ionized electron comes from COM
	    dE += Echg * cuEunit / (N0*eV);

	    PE[2] += Prob;
	    EE[2] += dE;
	  } // END: electron-ion
	  
	} // END: ionize
	
	if (T == recombine) {

	  if (IT.I1) {		// Ion is p1
	    cuFP_t WW = Prob;

	    if (WW < FF1[J]) {
	      FF1[J]   -= WW;
	      FF1[J-1] += WW;
	    } else {
	      WW = FF1[J];
	      FF1[J]    = 0.0;
	      FF1[J-1] += WW;
	    }

	    Prob = WW;		// Update to truncated value

	    // Electron KE lost in recombination is radiated by does not
	    // change COM energy
	    //
	    // cuFP_t Echg = iE1 * Prob / cuda_atomic_weights[IT.Z1];

	    // Electron KE fraction in recombination
	    //
	    // cuFP_t eE = Echg * cuEunit / (N0*cuEV);
	    
	    PE[1] += Prob;
	    EE[1] += dE;
	  } // END: ion-electron
	  else {		// Ion is p2
	    cuFP_t WW = Prob;

	    if (WW < FF2[J]) {
	      FF2[J]   -= WW;
	      FF2[J-1] += WW;
	    } else {
	      WW = FF2[J];
	      FF2[J]    = 0.0;
	      FF2[J-1] += WW;
	    }
	  
	    Prob = WW;		// Update to truncated value

	    // Electron KE lost in recombination is radiated by does not
	    // change COM energy
	    //
	    // cuFP_t Echg = iE2 * Prob / cuda_atomic_weights[IT.Z2];
	    
	    // Electron KE radiated in recombination
	    // cuFP_t eE = Echg * cuEunit / (N0*eV);
	    
	    PE[2] += Prob;
	    EE[2] += dE;
	  } // END: electro-ion

	} // END: recomb

      } // END: interaction loop

      // Convert energy loss from eV to system units and
      // total energy change for all interation
      //
      cuFP_t totalDE = 0.0;
      for (int i=0; i<3; i++) {
	EE[i] *= cuEV / cuEunit;
	totalDE += EE[i];
      }

      // Normalize probabilities and sum inelastic energy changes
      //
      cuFP_t probTot = 0.0;
      for (int i=0; i<3; i++) probTot += PE[i];

      if (probTot > 0.0) {
	for (int i=0; i<3; i++) PE[i] /= probTot;
      }

      //
      // Select interaction
      //
      cuFP_t Pr;
#if cuREAL == 4
      Pr = curand_uniform(state);
#else
      Pr = curand_uniform_double(state);
#endif
      
      unsigned int J = 2;
      if      (Pr < PE[0])         J = 0;
      else if (Pr < PE[0] + PE[1]) J = 1;

      // Deferred energy
      //
      cuFP_t E1[2] = {0.0, 0.0};
      cuFP_t E2[2] = {0.0, 0.0};
      cuFP_t totE  = 0.0;

      // Parse for Coulombic interaction
      //
      bool Coulombic = false;
      cuFP_t Tau = 0.0;
      if (maxT==ion_elec or maxT==ion_ion) {
    	Coulombic = true;

	cuFP_t dT  = spTau * cuTunit;

	if (maxT == ion_elec) {
	  if (IT.I1) {
	    cuFP_t pVel = sqrt(2.0 * kEe1 / mu1);
	    cuFP_t afac = esu*esu/(2.0*kEe1 > cuFloorEV*cuEV ? 2.0*kEe1 : cuFloorEV*cuEV);
	    Tau = ABrate._v[cid*4+1]*afac*afac*pVel * dT;
	  } else {
	    cuFP_t pVel = sqrt(2.0 * kEe2 / mu2);
	    cuFP_t afac = esu*esu/(2.0*kEe2 > cuFloorEV*cuEV ? 2.0*kEe2 : cuFloorEV*cuEV);
	    Tau = ABrate._v[cid*4+2]*afac*afac*pVel * dT;
	  }
	} else {
	  cuFP_t pVel = sqrt(2.0 * kEi / mu0);
	  cuFP_t afac = esu*esu/(2.0*kEi > cuFloorEV*cuEV ? 2.0*kEi : cuFloorEV*cuEV);
	  Tau = ABrate._v[cid*4+0]*afac*afac*pVel * dT;
	}
      }

      
      //
      // Apply neutral-neutral scattering and energy loss
      //
      if (J==0) {
	cuFP_t v1[3], v2[3];

	for (int k=0; k<3; k++) {
	  // Both particles are neutrals or ions
	  v1[k]  = p1.vel[k];
	  v2[k]  = p2.vel[k];
	}

	if (IT.W1 >= IT.W2)
	  cudaScatterTrace(Mu1, Mu2, Eta1, Eta2, IT.W1, IT.W2, E1, E2, totE,
			   v1, v2, totalDE, Tau, state, Coulombic);
	else
	  cudaScatterTrace(Mu2, Mu1, Eta2, Eta1, IT.W2, IT.W1, E2, E1, totE,
			   v2, v1, totalDE, Tau, state, Coulombic);

	// Time-step computation
	//
	{
	  cuFP_t dt = totE/totalDE;
	  deltT._v[i1._v[n]] = dt < deltT._v[i1._v[n]] ? dt : deltT._v[i1._v[n]];
	  deltT._v[i2._v[n]] = dt < deltT._v[i2._v[n]] ? dt : deltT._v[i2._v[n]];
	}
      
	for (int k=0; k<3; k++) {
	  // Both particles are ions
	  p1.vel[k] = v1[k];
	  p2.vel[k] = v2[k];
	}
      
      } // END: Atom-atom interaction


      // Apply ion/neutral-electron scattering and energy loss
      // Ion is Particle 1, Electron is Particle 2
      //
      if (J==1) {
	cuFP_t v1[3], v2[3];
      
	for (int k=0; k<3; k++) {
	  // Particle 1 is the ion
	  v1[k]  = p1.vel[k];
	  // Particle 2 is the electron
	  v2[k]  = p2.datr[epos+k];
	}

	PE[1] = totalDE;
      
	if (IT.W1 >= IT.W2)
	  cudaScatterTrace(Mu1, Mue*Eta2, Eta1, Eta2, IT.W1, IT.W2, E1, E2, totE,
			   v1, v2, totalDE, Tau, state, Coulombic);
	else
	  cudaScatterTrace(Mue*Eta2, Mu1, Eta2, Eta1, IT.W2, IT.W1, E2, E1, totE,
			   v2, v1, totalDE, Tau, state, Coulombic);
			   

	// Time-step computation
	//
	{
	  cuFP_t dt = totE/totalDE;
	  deltT._v[i1._v[n]] = dt < deltT._v[i1._v[n]] ? dt : deltT._v[i1._v[n]];
	  deltT._v[i2._v[n]] = dt < deltT._v[i2._v[n]] ? dt : deltT._v[i2._v[n]];
	}

      
	for (int k=0; k<3; k++) {
	  // Particle 1 is the ion
	  p1.vel[k] = v1[k];

	  // Particle 2 is the elctron
	  p2.datr[epos+k] = v2[k];
	}
	
      } // END: PE[1] (Ion-electron interaction)

      //
      // Apply ion/neutral-electron scattering and energy loss
      // Ion is Particle 2, Electron is Particle 1
      //
      if (J==2) {
	cuFP_t v1[3], v2[3];

	for (int k=0; k<3; k++) {
	  // Particle 1 is the elctron
	  v1[k]  = p1.datr[epos+k];
	  // Particle 2 is the ion
	  v2[k]  = p2.vel[k];
	}

	PE[2]  = totalDE;
      
	if (IT.W1 >= IT.W2)
	  cudaScatterTrace(Mue*Eta1, Mu2, Eta1, Eta2, IT.W1, IT.W2, E1, E2, totE,
			   v1, v2, totalDE, Tau, state, Coulombic);
	else
	  cudaScatterTrace(Mu2, Mue*Eta1, Eta2, Eta1, IT.W2, IT.W1, E2, E1, totE,
			   v2, v1, totalDE, Tau, state, Coulombic);

	// Time-step computation
	//
	{
	  cuFP_t dt = totE/totalDE;
	  deltT._v[i1._v[n]] = dt < deltT._v[i1._v[n]] ? dt : deltT._v[i1._v[n]];
	  deltT._v[i2._v[n]] = dt < deltT._v[i2._v[n]] ? dt : deltT._v[i2._v[n]];
	}
	
	for (int k=0; k<3; k++) {
	  // Particle 1 is the electron
	  p1.datr[epos+k] = v1[k];

	  // Particle 2 is the ion
	  p2.vel[k] = v2[k];
	}
      
      } // END: Electron-Ion interaction

    } // END: interactions with atoms AND electrons
    
  } // END: stride

  delete [] FF1;
  delete [] FF2;
}

void * CollideIon::collide_thread_cuda(void * arg)
{
  cuda_initialize();

  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, c0->cudaDevice);

  int id = static_cast<int>(((thrd_pass_arguments*)arg)->id);
  
  thread_timing_beg(id);
  
  if (id==0) {
    std::ostringstream sout;
    sout << "Collide::collide: ENTERING cuda collide_thread, T=" << tnow;
    (*barrier)(sout.str(), __FILE__, __LINE__);
  }

  // Initialize cell loop diagnostics
  //
  pre_cell_loop(id);

  // Start execution timer
  //
  cellTime[id].start();
  
  // Number of cells to process
  //
  size_t Ncells = cellist[id].size();

  // Structures for cell boundaries and counts
  //
  thrust::host_vector<int>    cellI, cellN, pairs, i1, i2, cc;
  thrust::host_vector<char>   flagI;
  thrust::host_vector<cuFP_t> h_volC(Ncells), h_tauC(Ncells);

  size_t Pcount = 0, Npairs = 0, Count = 0;

  // Loop over cells to get count of bodies and cells to process
  //
  for (unsigned j=0; j<cellist[id].size(); j++ ) {
    
    // The current cell
    //
    pCell *c = cellist[id][j];

    // Skip cell if this time has already been computed
    if (c->time >= tnow) {
      continue;
    }

    h_volC[j] = c->Volume();
    h_tauC[j] = dtime / (1<<c->maxplev);

    auto number = c->bods.size();
    if (number>1) {

      cellI.push_back(Pcount);
      cellN.push_back(number);
      pairs.push_back(Npairs);
      flagI.push_back(1);
      Pcount += number;

      for (size_t b=0; b<number/2; b++) {
	i1.push_back(c->bods[2*b+0]);
	i2.push_back(c->bods[2*b+1]);
	cc.push_back(Count);
	Npairs++;
      }
      if ((number/2)*2 != number) {
	flagI.back() = 2;
	i1.push_back(c->bods.front());
	i2.push_back(c->bods.back());
	flagI.push_back(0);
	Npairs++;
      }
      Count++;
    }
  }

  // Prepare for cudaParticle staging
  //
  if (c0->host_particles.capacity()<Pcount)
    c0->host_particles.reserve(Pcount);
  c0->host_particles.resize(Pcount);

  Component::hostPartItr hit = c0->host_particles.begin();

  // Species map info
  int minSp = std::numeric_limits<int>::max(), maxSp = 0;
  int numNeut = 0, numProt = 0, numIon = 0;
  for (auto s : SpList) {
    speciesKey k = s.first;
    minSp = std::min<int>(minSp, s.second);
    maxSp = std::max<int>(maxSp, s.second);
    if (k.second==1) numNeut++;
    else             numIon++;
    if (k.first==1)  numProt++;
  }
  numProt *= numNeut;
  maxSp++;

  if (use_elec>=0)
    maxSp = std::max<int>(maxSp, use_elec+3);
  else
    throw GenericError("use_elec must be set to use CUDA Trace implementation",  __FILE__, __LINE__);

  // Copy particles to DEVICE
  //
  thrust::device_vector<cudaParticle> d_part = c0->host_particles;

  // Copy cell boundaries and counts to DEVICE
  //
  thrust::device_vector<int>    d_cellI = cellI;
  thrust::device_vector<int>    d_cellN = cellN;
  thrust::device_vector<int>    d_i1    = i1;
  thrust::device_vector<int>    d_i2    = i2;
  thrust::device_vector<int>    d_cc    = cc;
  thrust::device_vector<int>    d_pairs = pairs;

  // Grid size computation
  //
  int N        = cellI.size();	// Number of cells
  int stride   = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
  int gridSize = (N+BLOCK_SIZE*stride-1)/(BLOCK_SIZE*stride);

				// These do not need copying back
  thrust::device_vector<cuFP_t> d_Ivel2(N), d_Evel2(N);
  thrust::device_vector<cuFP_t> d_PiProb(N*4), d_ABrate(N*4);
  thrust::device_vector<cuFP_t> d_volC(h_volC), d_tauC(h_tauC), d_selC(N);


  // Initialize per cell info
  //
  cellInitKernel<<<gridSize, BLOCK_SIZE>>>
    (toKernel(d_part),		// Particle array (input)
     toKernel(d_Ivel2),		// Mean squared ion velocity (output)
     toKernel(d_Evel2),		// Mean squared electron velocity (output)
     toKernel(d_PiProb),	// For BN algorithm (output)
     toKernel(d_ABrate),	// For BN algorithm (output)
     toKernel(d_volC),		// Cell volume (input)
     toKernel(d_tauC),		// Cell time step (input)
     toKernel(d_selC),		// Prefactor for inelastic normalization(output)
     toKernel(d_cellI),		// Cell index (input)
     toKernel(d_cellN),		// Cell body count (output)
     toKernel(cuElems),		// Ionization state info (input)
     minSp,			// Location of state data in particle attribute array (input)
     stride);			// Stride for processing (input)

  // Compute the cross sections for each interaction en masse
  //
  N        = i1.size();	// Number of pairs
  stride   = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
  gridSize = (N+BLOCK_SIZE*stride-1)/(BLOCK_SIZE*stride);

				// These do not need copying back
  unsigned int totalXCsize =
    numNeut*numNeut + numNeut*2 + numNeut*numProt*2 + numIon*numIon + 
    numIon*5*2;

  thrust::device_vector<cuFP_t>         d_cross(N*totalXCsize);
  thrust::device_vector<cuFP_t>         d_delph(N*totalXCsize);
  thrust::device_vector<uchar3>         d_xspc1(N*totalXCsize);
  thrust::device_vector<uchar3>         d_xspc2(N*totalXCsize);
  thrust::device_vector<cudaInterTypes> d_xtype(N*totalXCsize);
  thrust::device_vector<int>            d_flagI(flagI);
  thrust::device_vector<curandState>    d_randS(N);

  initCurand<<<gridSize, BLOCK_SIZE>>>(toKernel(d_randS), seed);

  crossSectionKernel<<<gridSize, BLOCK_SIZE>>>
    (toKernel(d_part),   toKernel(d_randS),
     toKernel(d_cross),  toKernel(d_delph), toKernel(d_xspc1),  toKernel(d_xspc2),
     toKernel(d_xtype),  toKernel(d_cross),
     toKernel(d_i1),     toKernel(d_i2),     toKernel(d_cc),    toKernel(d_volC),
     toKernel(d_Ivel2),  toKernel(d_Evel2),
     toKernel(d_PiProb), toKernel(d_ABrate), toKernel(d_flagI), 
     toKernel(xsc_H),    toKernel(xsc_He),   toKernel(xsc_pH),  toKernel(xsc_pHe),
     toKernel(cuElems),
     totalXCsize, use_elec, stride);

  // Do the interactions
  //
  thrust::device_vector<cuFP_t> d_deltT(Pcount, std::numeric_limits<cuFP_t>::max());

  partInteractions<<<gridSize, BLOCK_SIZE>>>
    (toKernel(d_part),   toKernel(d_randS),
     toKernel(d_cross),  toKernel(d_delph),  toKernel(d_xspc1),  toKernel(d_xspc2),
     toKernel(d_xtype),  toKernel(d_cross),
     toKernel(d_i1),     toKernel(d_i2),     toKernel(d_cc),    toKernel(d_selC),
     toKernel(d_Ivel2),  toKernel(d_Evel2),
     toKernel(d_PiProb), toKernel(d_ABrate), toKernel(d_deltT), toKernel(d_flagI), 
     toKernel(cuElems),
     spTau[id], totalXCsize, use_elec, stride);

  // Finally, copy back particles to host
  // 

  /*

  // Loop over cells, processing collisions in each cell
  //
  for (unsigned j=0; j<cellist[id].size(); j++ ) {

    // Cache acceptance fraction for scaling MFP for time step
    // selection
    //
    if (acceptCount > 0) {
      cuFP_t scale = static_cast<cuFP_t>(totalCount)/acceptCount; 
      meanLambda *= scale;
      meanCollP  /= scale;
      if (MFPTS) {
	pthread_mutex_lock(&tlock);
	selMFP[c] = meanLambda;
	pthread_mutex_unlock(&tlock);
      }
    }

    collSoFar[id] = collTime[id].stop();
  
    // Update collision counter in the cell
    //
    c->iattrib["collCount"] = acceptCount;

    // Compute dispersion diagnostics
    //
    stat3Time[id].start();
  
    cuFP_t tmass = 0.0;
    vector<cuFP_t> velm(3, 0.0), velm2(3, 0.0);
    for (unsigned j=0; j<number; j++) {
      Particle* p = tree->Body(c->bods[j]);
      for (unsigned k=0; k<3; k++) {
	velm[k]  += p.mass*p.vel[k];
	velm2[k] += p.mass*p.vel[k]*p.vel[k];
      }
      tmass += p.mass;
    }
  
    if (tmass>0.0) {
      for (unsigned k=0; k<3; k++) {
	
	velm[k] /= tmass;
	velm2[k] = velm2[k] - velm[k]*velm[k]*tmass;
	if (velm2[k]>0.0) {
	  tdispT[id][k] += velm2[k];
	  tmassT[id]    += tmass;
	}
      }
    }
    
    //
    // General hook for the derived classes for specific computations
    // and diagnostics
    //
  
    finalize_cell(c, Fn, kedsp, tau, id);
  
    // Update cell time
    //
    c->time = tnow;

    stat3SoFar[id] = stat3Time[id].stop();
  
    //
    // Compute Knudsen and/or Strouhal number
    //
    if (use_Kn>=0 || use_St>=0) {
      cuFP_t cL = pow(volc, 0.33333333);
      cuFP_t Kn = meanLambda/cL;
      cuFP_t St = cL/fabs(tau*sqrt(fabs(kedsp))+1.0e-18);
      for (unsigned j=0; j<number; j++) {
	Particle* p = tree->Body(c->bods[j]);
	if (use_Kn>=0) p.datr[use_Kn] = Kn;
	if (use_St>=0) p.datr[use_St] = St;
      }
    }
    
    // Record effort per particle in microseconds
    //
    curcSoFar[id] = curcTime[id].stop();
    cuFP_t tt = curcSoFar[id];
    if (EFFORT) {
      if (effortAccum) 
	effortNumber[id].push_back(pair<long, unsigned>(tt, number));
      cuFP_t effort = static_cast<cuFP_t>(tt)/number;
      for (unsigned k=0; k<number; k++) 
	tree->Body(c->bods[k])->effort += effort;
    }
  
    // Usage debuging
    //
    if (minUsage[id*2+EPSMused] > tt) {
      minUsage[id*2+EPSMused] = tt;
      minPart [id*2+EPSMused] = number;
      minCollP[id*2+EPSMused] = meanCollP;
    }
    if (maxUsage[id*2+EPSMused] < tt) {
      maxUsage[id*2+EPSMused] = tt;
      maxPart [id*2+EPSMused] = number;
      maxCollP[id*2+EPSMused] = meanCollP;
    }
    
  } // Loop over cells

  */

  if (id==0) {
    std::ostringstream sout;
    sout << "Collide::collide: AFTER cell loop, T=" << tnow;
    (*barrier)(sout.str(), __FILE__, __LINE__);
  }

  cellSoFar[id] = cellTime[id].stop();

  // Diagnostics at end of cell loop
  //
  post_cell_loop(id);

  thread_timing_end(id);
  
  return (NULL);
}
