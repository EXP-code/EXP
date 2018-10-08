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
  neut_neut  = 0,
  neut_elec  = 1,
  neut_prot  = 2,
  ion_ion    = 3,
  ion_elec   = 4,
  free_free  = 5,
  col_excite = 6,
  col_ionize = 7,
  recombine  = 8
};

const int maxAtomicNumber = 15;
__constant__ cuFP_t cuda_atomic_weights[maxAtomicNumber], cuFloorEV, cuEsu;
__constant__ cuFP_t cuVunit, cuMunit, cuTunit, cuLunit, cuAmu, cuEV, cuLogL, cuCrossfac;
__constant__ bool   cuMeanKE, cuMeanMass, cuNewRecombAlg;

// Link static parameters from CollideIon.cc
//
extern double FloorEv;
extern bool   MEAN_KE;
extern bool   newRecombAlg;;

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

  cuda_safe_call(cudaMemcpyToSymbol(cuMeanKE, &MEAN_KE, sizeof(bool)), 
		 __FILE__, __LINE__, "Error copying cuMeanKE");

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
			       dArray<cuFP_t> meanM,
			       dArray<cuFP_t> Ivel2,
			       dArray<cuFP_t> Evel2,
			       dArray<cuFP_t> PiProb,
			       dArray<cuFP_t> ABrate,
			       dArray<cuFP_t> volC,
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

      meanM._v[c] = 0.0;

      int nbods = cellN._v[c];

      cuFP_t massP = 0.0, numbP = 0.0, massE = 0.0;
      cuFP_t evel2 = 0.0, ivel2 = 0.0, numQ2 = 0.0;
      cuFP_t densI = 0.0, densQ = 0.0, densE = 0.0;

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
  
      if (numbP>0.0) meanM._v[c] = massP/numbP;
      if (massP>0.0) Ivel2._v[c] = ivel2/massP;
      if (massE>0.0) Evel2._v[c] = evel2/massE;
      if (densQ>0.0) numQ2      /= densQ;
      
      cuFP_t ddfac = dfac/volC._v[c];

      densI *= ddfac;
      densQ *= ddfac;
      densE *= ddfac;

      // Compute per channel Coulombic probabilities
      //
      // Ion probabilities
      //
      cuFP_t muii = meanM._v[c]/2.0;
      cuFP_t muee = cuda_atomic_weights[0]/2.0;
      cuFP_t muie = cuda_atomic_weights[0] * meanM._v[c]/(cuda_atomic_weights[0] + meanM._v[c]);
      
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

    } // END: cell

  } // END: stride

}

__global__ void crossSectionKernel(dArray<cudaParticle> in,
				   dArray<curandState> randS,
				   dArray<cuFP_t> cross,
				   dArray<uchar4> xspec,
				   dArray<cudaInterTypes> xtype,
				   dArray<cuFP_t> xctot,
				   dArray<int>    i1,
				   dArray<int>    i2,
				   dArray<int>    cc,
				   dArray<cuFP_t> volC,
				   dArray<cuFP_t> meanM,
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

    if (n < meanM._s) {

      cudaParticle& p1 = in._v[i1._v[n]];
      cudaParticle& p2 = in._v[i2._v[n]];
      int          cid = cc._v[n];
      
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
	    xspec._v[n*numxc+J] = make_uchar4(Z, C, ZZ, CC);
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
	    xspec._v[n*numxc+J] = make_uchar4(Z, C, ZZ, CC);
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
	    xspec._v[n*numxc+J] = make_uchar4(Z, C, ZZ, CC);
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
	    xspec._v[n*numxc+J] = make_uchar4(Z, C, 0, 0);
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
	    xspec._v[n*numxc+J] = make_uchar4(0, 0, Z, C);
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
	  xspec._v[n*numxc+J] = make_uchar4(Z, C, 0, 0);
	  xtype._v[n*numxc+J] = ion_elec;
	  xctot._v[n]        += crs;
	  
	  J++;
	  
	  // Particle 2 ION, Particle 1 ELECTRON
	  crs  = 0.0;
	  kEc  = cuMeanKE ? Eelc : kEe2;
	  afac = cuEsu*cuEsu/(2.0*kEc > cuFloorEV ? 2.0*kEc : cuFloorEV) * cuEV * 1.0e7;

	  crs = 2.0 * ABrate._v[cid*4+2] * afac*afac * gVel1 / PiProb._v[cid*4+2];

	  cross._v[n*numxc+J] = crs;
	  xspec._v[n*numxc+J] = make_uchar4(0, 0, Z, C);
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
	    xspec._v[n*numxc+J] = make_uchar4(Z, C, 0, 0);
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
	    xspec._v[n*numxc+J] = make_uchar4(0, 0, Z, C);
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
	    xspec._v[n*numxc+J] = make_uchar4(Z, C, 0, 0);
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
	    xspec._v[n*numxc+J] = make_uchar4(0, 0, Z, C);
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
	    xspec._v[n*numxc+J] = make_uchar4(Z, C, 0, 0);
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
	    xspec._v[n*numxc+J] = make_uchar4(0, 0, Z, C);
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
	      xspec._v[n*numxc+J] = make_uchar4(Z, C, 0, 0);
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
	      xspec._v[n*numxc+J] = make_uchar4(0, 0, Z, C);
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
	      xspec._v[n*numxc+J] = make_uchar4(Z, C, 0, 0);
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
	      xspec._v[n*numxc+J] = make_uchar4(0, 0, Z, C);
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

void * CollideIon::collide_thread_cuda(void * arg)
{
  cuda_initialize();

  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, c0->cudaDevice);

  sKeyDmap *Fn  = static_cast<sKeyDmap*>(((thrd_pass_arguments*)arg)->fn  );
  int id        = static_cast<int>      (((thrd_pass_arguments*)arg)->id  );
  
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
  
  // Structures for cell boundaries and counts
  //
  thrust::host_vector<int>    cellI, cellN, pairs, i1, i2, cc;
  thrust::host_vector<char>   flagI;
  thrust::host_vector<cuFP_t> h_volC;

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
  thrust::device_vector<cuFP_t> d_meanM(N), d_Ivel2(N), d_Evel2(N);
  thrust::device_vector<cuFP_t> d_PiProb(N*4), d_ABrate(N*4);
  thrust::device_vector<cuFP_t> d_volC(h_volC);


  // Initialize per cell info
  //
  cellInitKernel<<<gridSize, BLOCK_SIZE>>>
    (toKernel(d_part),
     toKernel(d_meanM),  toKernel(d_Ivel2),  toKernel(d_Evel2),
     toKernel(d_PiProb), toKernel(d_ABrate), toKernel(d_volC),
     toKernel(d_cellI),  toKernel(d_cellN),  toKernel(cuElems),
     minSp, stride);



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
  thrust::device_vector<uchar4>         d_xspec(N*totalXCsize);
  thrust::device_vector<cudaInterTypes> d_xtype(N*totalXCsize);
  thrust::device_vector<int>            d_flagI(flagI);
  thrust::device_vector<curandState>    d_randS(N);

  initCurand<<<gridSize, BLOCK_SIZE>>>(toKernel(d_randS), seed);

  crossSectionKernel<<<gridSize, BLOCK_SIZE>>>
    (toKernel(d_part),   toKernel(d_randS),
     toKernel(d_cross),  toKernel(d_xspec),  toKernel(d_xtype), toKernel(d_cross),
     toKernel(d_i1),     toKernel(d_i2),     toKernel(d_cc),    toKernel(d_volC),
     toKernel(d_meanM),  toKernel(d_Ivel2),  toKernel(d_Evel2),
     toKernel(d_PiProb), toKernel(d_ABrate), toKernel(d_flagI), 
     toKernel(xsc_H),    toKernel(xsc_He),   toKernel(xsc_pH),  toKernel(xsc_pHe),
     toKernel(cuElems),
     totalXCsize, use_elec, stride);

  // Do the interactions
  //

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
