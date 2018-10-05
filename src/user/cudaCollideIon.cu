// -*- C++ -*-

#include <Component.H>
#include <CollideIon.H>
#include <cudaElastic.cuH>

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
__constant__ cuFP_t cuda_atomic_weights[MaxAtomicNumber], cuFloorEV;
__constant__ cuFP_t cuVunit, cuMunit, cuTunit, cuLunit, cuAmu, cuEV, cuLogL;
__constant__ bool   cuMeanKE, cuMeanMass;

void CollideIon::cuda_atomic_weights_init()
{
  cudaElasticInit();

  std::vector<cuFP_t> weights(MaxAtomicNumber);

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

  cuda_safe_call(cudaMemcpyToSymbol(cuda_atomic_weights, &weights[0], sizeof(cuFP_t)*MaxAtomicNumber), 
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
  v = FloorEv;
  cuda_safe_call(cudaMemcpyToSymbol(cuFloorEV, &v, sizeof(cuFP_t)), 
		 __FILE__, __LINE__, "Error copying cuFloorEV");

  cuda_safe_call(cudaMemcpyToSymbol(cuMeanKE, &MEAN_KE, sizeof(bool)), 
		 __FILE__, __LINE__, "Error copying cuMeanKE");

  cuda_safe_call(cudaMemcpyToSymbol(cuMeanMass, &MeanMass, sizeof(bool)), 
		 __FILE__, __LINE__, "Error copying cuMeanMass");
}  

__global__ void cellInitKernel(dArray<cudaParticle> in,
			       dArray<cuFP_t> meanM,
			       dArray<cuFP_t> Ivel2,
			       dArray<cuFP_t> Evel2,
			       dArray<cuFP_t> PiProb,
			       dArray<cuFP_t> ABrate,
			       dArray<cuFP_t> cellI,
			       dArray<cuFP_t> cellN,
			       dArray<char> Z,
			       dArray<char> C,
			       dArray<char> I,
			       int minSp,
			       unsigned int stride)
{
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;

  for (int s=0; s<stride; s++) {

    int c = tid*stride + s;

    if (c < cellI._s) {

      meanM._v[c] = 0.0;

      int nbods = cellN._v[c];

      cuFP_t massP = 0.0, numbP = 0.0, massE = 0.0;
      cuFP_t evel2 = 0.0, ivel2 = 0.0, numQ2 = 0.0, densQ = 0.0;

      for (size_t i=0; i<nbods; i++) {
	
	// The particle
	cudaParticle *p = in._v[i + cellI._v[c]];

	// Mass per cell
	massP += p->mass;
	
	// Mass-weighted trace fraction
	cuFP_t ee = 0.0;
	for (int k=0; k<Z._s; k++) {
	  cuFP_t ff = p->datr[I._v[k] - minSp];
	  cuFP_t ww = ff/cuda_atomic_weights[Z._v[k]];
	  cuFP_t qq = C._v[k] - 1;
	  // Mean number
	  numbP += p->mass * ww;
	  // Electron fraction
	  ee += ww * qq;
	  // Charge squared
	  numQ2 += p->mass * ww * qq*qq;
	  if (C._v[k]) densQ += p->mass * ww;
	}
	
	cuFP_t eVel2 = 0.0, iVel2 = 0.0;
	for (int l=0; l<3; l++) {
	  cuFP_t ve  = p->dattrib[use_elec+l];
	  eVel2 += ve*ve;
	  cuFP_t vi  = p->vel[l];
	  iVel2 += vi*vi;
	}
	
	evel2 += p->mass * ee * eVel2;
	ivel2 += p->mass * iVel2;
	massE += p->mass * ee;
      }
    }
  
    if (numbP>0.0) meanM._v[c] = massP/numbP;
    if (massP>0.0) Ivel2._v[c] = ivel2/massP;
    if (massE>0.0) Evel2._v[c] = evel2/massE;
    if (densQ>0.0) numQ2      /= densQ;
    
    // Compute per channel Coulombic probabilities
    //
    // Ion probabilities
    //
    cuFP_t muii = meanM._v[c]/2.0;
    cuFP_t muee = cuda_atomic_weights[0]/2.0;
    cuFP_t muie = cuda_atomic_weights[0] * meanM._v[c]/(cuda_atomic_weights[0] + meanM._v[c]);
    
    // Ion-Ion
    PiProb._v[c*4 + 0] =
      densQtot +
      densEtot * pow(2.0*Ivel2._v[c], 1.5) * muii*muii /
      (pow(Ivel2._v[c] + Evel2._v[c], 1.5) * muie*muie * numQ2);
    //                                               ^
    //                                               |
    // The density is weighted by q^2 for each species
    
    // Ion-Electron
    PiProb._v[c*4 + 1] =
      densQtot * pow(Ivel2._v[c] + Evel2._v[c], 1.5) * muie*muie * numQ2 /
      (pow(2.0*Ivel2._v[c], 1.5) * muii*muii) +
      densEtot ;
    
    // Electron-Ion
    PiProb._v[c*4 + 2] =
      densQtot +
      densEtot * pow(Ivel2._v[c] + Evel2._v[c], 1.5) * muie*muie /
      (pow(2.0*Evel2._v[c], 1.5) * muee*muee * numQ2);
    
    // Electron-Electron
    PiProb._v[c*4 + 3] =
      densQtot * pow(2.0*Evel2[id], 1.5) * muee*muee * numQ2[id] /
      (pow(Ivel2[id] + Evel2[id], 1.5) * muie*muie) +
      densEtot;
    
    // Rate coefficients
    ABrate._v[c*4 + 0] = 2.0*M_PI * PiProb[c*4 + 0] * logL * pow(numQ2[id]*numQ2[id], 2.0);
    
    ABrate._v[c*4 + 1] = 2.0*M_PI * PiProb[c*4 + 1] * logL * pow(numQ2[id], 2.0);
  
    ABrate._v[c*4 + 2] = 2.0*M_PI * PiProb[c*4 + 2] * logL * pow(numQ2[id], 2.0);
    
    ABrate._v[c*4 + 3] = 2.0*M_PI * PiProb[c*4 + 3] * logL ;
    
   }
}

__global__ void crossSectionKernel(dArray<cudaParticle> in,
				   dArray<cuFP_t> cross,
				   dArray<uchar4> xspec,
				   dArray<uchar>  xtype,
				   dArray<cuFP_t> xctot,
				   dArray<int>    i1,
				   dArray<int>    i2,
				   dArray<int>    cc,
				   dArray<cuFP_t> meanM,
				   dArray<cuFP_t> Ivel2,
				   dArray<cuFP_t> Evel2,
				   dArray<cuFP_t> PiProb,
				   dArray<cuFP_t> ABrate,
				   dArray<int>    flagI, 
				   dArray<char>   dZ,
				   dArray<char>   dC,
				   dArray<char>   dI,
				   dArray<cuFP_t> xsc_H,
				   dArray<cuFP_t> xsc_He,
				   dArray<cuFP_t> xsc_pH,
				   dArray<cuFP_t> xsc_pHe,
				   int numxc,
				   int minSp,
				   int eePos,
				   unsigned int stride)
{
  const int tid = blockDim.x * blockIdx.x + threadIdx.x;
  const int Nsp = dZ._s;

  for (int s=0; s<stride; s++) {

    int n = tid*stride + s;

    if (n < meanM._s) {

      cudaParticle *p1 = i1._v[n];
      cudaParticle *p2 = i2._v[n];
      int          cid = cc._v[n];

      cuFP_t Eta1=0.0, Eta2=0.0, Mu1=0.0, Mu2=0.0, Sum1=0.0, Sum2=0.0;

      xctot._v[n] = 0.0;	// Total cross section accumulator

      int J = 0;		// Cross section position counter

      for (int k=0; k<Nsp; k++) {

				// Number fraction of ions
	cuFP_t one = p1->datr[dI._v[k]-minSp] / cuda_atomic_weights[dZ._v[k]];
	cuFP_t two = p2->datr[dI._v[k]-minSp] / cuda_atomic_weights[dZ._v[k]];

				// Electron number fraction
	Eta1 += one * (dC._v[k] - 1);
	Eta2 += two * (dC._v[k] - 1);

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
      cuFP_t N1 = p1->mass*cuMunit/(Mu1*cuAmu);
      cuFP_t N2 = p2->mass*cuMunit/(Mu2*cuAmu);

      cuFP_t vel = 0.0;
      cuFP_t eVel0 = 0.0, eVel1 = 0.0, eVel2 = 0.0;
      cuFP_t gVel0 = 0.0, gVel1 = 0.0, gVel2 = 0.0;
      for (int k=0; k<3; k++) {
	cuFP_t del = p1->vel[k] - p2->vel[k];
	vel += del*del;
	del = p1->vel[k] - p2->vel[k];

	double rvel0 = p1->datr[eePos+i] - p2->datr[eePos+i];
	double rvel1 = p1->datr[eePos+i] - p2->vel[i];
	double rvel2 = p2->datr[eePos+i] - p1->vel[i];

	eVel0 += rvel0*rvel0;
	eVel1 += rvel1*rvel1;
	eVel2 += rvel2*rvel2;

	// Scaled electron relative velocity
	if (cuMeanMass) {
	  rvel0 = p1->datr[eePos+i]*sqrt(Eta1) - p2->datr[eePos+i]*sqrt(Eta2);
	  rvel1 = p1->datr[eePos+i]*sqrt(Eta1) - p2->vel[i];
	  rvel2 = p2->datr[eePos+i]*sqrt(Eta2) - p1->vel[i];
	  
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
	gVel0 = eVel0;
	gVel1 = eVel1;
	gVel2 = eVel2;
      }

      cuFP_t  m1 = Mu1 * cuAmu;
      cuFP_t  m2 = Mu2 * cuAmu;
      cuFP_t  me = cuda_atomic_weights[0] * cuAmu;
      cuFP_t mu0 = m1 * m2 / (m1 + m2);
      cuFP_t mu1 = m1 * me / (m1 + me);
      cuFP_t mu2 = m2 * me / (m2 + me);

      cuFP_t kEi = 0.5 * mu0 * vel * vel;

      cuFP_t facE = 0.5 * cuAmu * cuVunit * cuVunit / cuEV;
      cuFP_t Eion = fac * Ivel2._v[c] * 2.0/(Mu1 + Mu2);
      cuFP_t kEe1 = 0.5  * me * eVel2*eVel2;
      cuFP_t kEe2 = 0.5  * me * eVel1*eVel1;
      cuFP_t kEee = 0.25 * me * eVel0*eVel0;

      int nZ1 = 0;
      for (int k1=0; k1<Nsp; k1++) {

	unsigned char Z = dZ._v[k];
	unsigned char C = dC._v[k];
	unsigned char P = C1 - 1;
	int I = dI._v[k] - minSp;

	cuFP_t fac1 = p1->datr[I] / cuda_atomic_weights[Z] / Sum1;
	cuFP_t fac2 = p2->datr[I] / atomic_weights[Z] / Sum2;


	int NZ2 = 0;
	for (int kk=0; kk<Nsp; kk++) {
	  
	  unsigned char ZZ = dZ._v[kk];
	  unsigned char CC = dC._v[kk];
	  unsigned char PP = CC - 1;
	  int II = dI._v[kk] - minSp;

	  //--------------------------------------------------
	  // Particle 1 interacts with Particle 2
	  //--------------------------------------------------
	  
	  cuFP_t cfac = p1->datr[I] * p2->datr[II];

	  //-------------------------------
	  // *** Both particles neutral
	  //-------------------------------
    
	  if (P==0 and PP==0) {
	    
	    double cross = 0.0;
	    // Geometric cross sections based on
	    // atomic radius
	    double crs = (cudaGeometric(Z) + cudaGeometric(ZZ)) * cfac;
	
	    // Double counting
	    if (Z == ZZ) crs *= 0.5;

	    cross += crs*crossfac;

	    cross._v[n*numxc+J] = cross;
	    xspec._v[n*numxc+J] = make_char4(Z, C, ZZ, CC);
	    xtype._v[n*numxc+J] = neut_neut;
	    xctot._v[n] += cross;
	    
	    J++;		// Increment interaction counter
	  }

	  // --------------------------------------
	  // *** Neutral atom-proton scattering
	  // --------------------------------------

	  double crs1 = 0;

	  if (ZZ==1 and CC==2) {
	    if (Z==1 and P==0) crs1 = cudaElasticInterp(kEi, cuPH_Emin,  cuH_PH,  xsc_PH );
	    if (Z==2 and P==0) crs1 = cudaElasticInterp(kEi, cuPHe_Emin, cuH_PHe, xsc_PHe);
	    crs1 *= crossfac * cfac;
	  }
	  
	  if (Z==1 and C==2) {
	    if (ZZ==1 and PP==0) crs1 = cudaElasticInterp(kEi, cuPH_Emin,  cuH_PH,  xsc_PH );
	    if (ZZ==2 and PP==0) crs1 = cudaElasticInterp(kEi, cuPHe_Emin, cuH_PHe, xsc_PHe);
	    crs1 *= crossfac * cfac;
	  }

	  if (crs1>0.0) {
	    cross._v[n*numxc+J] = crs1;
	    xspec._v[n*numxc+J] = make_char4(Z, C, ZZ, CC);
	    xtype._v[n*numxc+J] = neut_prot;
	    xctot._v[n] += crs1;
	    
	    J++;
	  }

	  // --------------------------------------
	  // *** Ion-ion scattering
	  // --------------------------------------
	  //
	  if (P>0 and PP>0) {
	    double kEc  = cuMeanKE ? Eion : kEi;
	    double afac = cuEsu*cuEsu/std::max<double>(2.0*kEc*eV, cuFloorEV*cuEV) * 1.0e7;
	    double crs  = 2.0 * ABrate._v[cid*4+0] * afac*afac / PiProb[cid*4+0];

	    cross._v[n*numxc+J] = crs;
	    xspec._v[n*numxc+J] = make_char4(Z, C, ZZ, CC);
	    xtype._v[n*numxc+J] = ion_ion;
	    xctot._v[n] += crs;
	    
	    J++;
	  }

	} // End of inner species loop

	// --------------------------------------
	// *** Neutral atom-electron scattering
	// --------------------------------------
    
	// Particle 1 ION, Particle 2 ELECTRON
	//
	if (Z<=2 and P==0 and Eta2>0.0) {
	  double crs = 0.0;
	  if (Z==1) crs = cudaElasticInterp(kEe1, cuH_Emin, cuH_H, xsc_H) * gVel2 * Eta2 *
		      crossfac * cscl_[Z] * fac1;
	  if (Z==2) crs = cudaElasticInterp(kEe1, cuHe_Emin, cuH_He, xsc_He) * gVel2 * Eta2 *
		      crossfac * cscl_[Z] * fac1;

	  if (crs>0.0) {
	    cross._v[n*numxc+J] = crs;
	    xspec._v[n*numxc+J] = make_char4(Z, C, 0, 0);
	    xtype._v[n*numxc+J] = neut_elec;
	    xctot._v[n] += crs;
	    
	    J++;
	  }
	}

	// Particle 2 ION, Particle 1 ELECTRON
	//
	if (P==0 and Eta1>0.0) {
	  double crs = 0.0;
	  if (Z==1) crs = cudaElasticInterp(kEe2, cuH_Emin, cuH_H, xsc_H) * gVel2 * Eta2 *
		      crossfac * cscl_[Z] * fac1;
	  if (Z==2) crs = cudaElasticInterp(kEe2, cuHe_Emin, cuH_He, xsc_He) * gVel2 * Eta2 *
		      crossfac * cscl_[Z] * fac1;

	  if (crs>0.0) {
	    cross._v[n*numxc+J] = crs;
	    xspec._v[n*numxc+J] = make_char4(0, 0, Z, C);
	    xtype._v[n*numxc+J] = neut_elec;
	    xctot._v[n] += crs;
	    
	    J++;
	  }
	}

	// --------------------------------------
	// *** Ion-electron scattering
	// --------------------------------------
	//
	if (P>0 and Eta2>0) {
	  // Particle 1 ION, Particle 2 ELECTRON
	  double crs  = 0.0;
	  double kEc  = cuMeanKE ? Eelc : kEe1;
	  double afac = cuEsu*cuEsu/std::max<double>(2.0*kEc*cuEV, cuFloorEV*eV) * 1.0e7;
	    
	  crs = 2.0 * ABrate[cid*4+1] * afac*afac * gVel2 / PiProb[cid*4+1];
	  
	  cross._v[n*numxc+J] = crs;
	  xspec._v[n*numxc+J] = make_char4(Z, C, 0, 0);
	  xtype._v[n*numxc+J] = ion_elec;
	  xctot._v[n] += crs;
	  
	  J++;
	  
	  // Particle 2 ION, Particle 1 ELECTRON
	  crs  = 0.0;
	  kEc  = MEAN_KE ? Eelc  : kEe2;
	  afac = esu*esu/std::max<double>(2.0*kEc*eV, FloorEv*eV) * 1.0e7;

	  crs = 2.0 * ABrate[id*4+2] * afac*afac * gVel1 / PiProb[cid*4+2];

	  cross._v[n*numxc+J] = crs;
	  xspec._v[n*numxc+J] = make_char4(0, 0, Z, C);
	  xtype._v[n*numxc+J] = ion_elec;
	  xctot._v[n] += crs;
	  
	  J++;
	} // end: ion-electron scattering


	//-------------------------------
	// *** Free-free
	//-------------------------------
      
	if (Eta1>0.0 and Eta2>0.0) {

	  // Particle 1 ION, Particle 2 ELECTRON
	  double   ke  = std::max<double>(kEe1, cuFloorEV);
	  CFreturn ff  = ch.IonList[Q]->freeFreeCross(ke, id);

	double crs  = gVel2 * Eta2 * ff.first * fac1;
	
	if (std::isinf(crs)) crs = 0.0; // Sanity check
	
	if (crs>0.0) {

	  Interact::T t { free_free, Ion, Interact::edef };
	    
	  hCross[id].push_back(XStup(t));
	  hCross[id].back().crs = crs;
	  hCross[id].back().CF  = ff;

	  CProb[id][1] += crs;
	}
      }

      // Particle 2 ION, Particle 1 ELECTRON
      {
	double    ke = std::max<double>(kEe2[id], FloorEv);
	CFreturn  ff = ch.IonList[Q]->freeFreeCross(ke, id);

	double crs  = gVel1 * Eta1 * ff.first * fac2;
	
	if (std::isinf(crs)) crs = 0.0; // Sanity check
	  
	if (crs>0.0) {
	  
	  Interact::T t { free_free, Interact::edef, Ion };
	    
	  hCross[id].push_back(XStup(t));
	  hCross[id].back().crs = crs;
	  hCross[id].back().CF  = ff;

	  CProb[id][2] += crs;
	}
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
      double    ke = std::max<double>(kEe1[id], FloorEv);
      CEvector  CE = ch.IonList[Q]->collExciteCross(ke, id);

      double   crs = gVel2 * Eta2 * CE.back().first * fac1;
      
      if (DEBUG_CRS) trap_crs(crs);
      
      if (crs > 0.0) {
	Interact::T t { colexcite, Ion, Interact::edef };
	
	hCross[id].push_back(XStup(t));
	hCross[id].back().crs = crs;
	hCross[id].back().CE  = CE;

	CProb[id][1] += crs;
      }
    }
    // end: colexcite
    
    // Particle 2 nucleus has BOUND ELECTRON, Particle 1 has FREE ELECTRON
    //
    //  +--- Charge of the current subspecies
    //  |
    //  |       +--- Electron fraction of partner
    //  |       |
    //  V       V
    if (P<Z and Eta1>0) {
      double    ke = std::max<double>(kEe2[id], FloorEv);
      CEvector  CE = ch.IonList[Q]->collExciteCross(ke, id);

      double   crs = gVel1 * Eta1 * CE.back().first * fac2;
      
      if (DEBUG_CRS) trap_crs(crs);
      
      if (crs > 0.0) {
	Interact::T t { colexcite, Interact::edef, Ion };
	
	hCross[id].push_back(XStup(t));
	hCross[id].back().crs = crs;
	hCross[id].back().CE  = CE;
	
	CProb[id][2] += crs;
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
      
      double ke    = std::max<double>(kEe1[id], FloorEv);
      double DI    = ch.IonList[Q]->directIonCross(ke, id);

      double crs   = gVel2 * Eta2 * DI * fac1;
      
      if (DEBUG_CRS) trap_crs(crs);
      
      if (crs > 0.0) {
	Interact::T t { ionize, Ion, Interact::edef };
	
	hCross[id].push_back(XStup(t));
	hCross[id].back().crs = crs;
	
	CProb[id][1] += crs;
      }
    }
    // end: ionize
    
    // Particle 2 nucleus has BOUND ELECTRON, Particle 1 has FREE ELECTRON
    //
    //  +--- Charge of the current subspecies
    //  |
    //  |       +--- Electron fraction of partner
    //  |       |
    //  V       V
    if (P<Z and Eta1) {
      
      double ke    = std::max<double>(kEe2[id], FloorEv);
      double DI    = ch.IonList[Q]->directIonCross(ke, id);

      double crs   = gVel1 * Eta1 * DI * fac2;
      
      if (DEBUG_CRS) trap_crs(crs);
      
      if (crs > 0.0) {
	Interact::T t { ionize, Interact::edef, Ion };
	
	hCross[id].push_back(XStup(t));
	hCross[id].back().crs = crs;
	
	CProb[id][2] += crs;
      }
    }
    // end: ionize

    //-------------------------------
    // *** Radiative recombination
    //-------------------------------

    // The "new" algorithm uses the electron energy of the ion's
    // electron rather than the standard particle partner.
    //

    if (newRecombAlg) {

      // Particle 1 is ION, Particle 2 has ELECTRON
      //
      //  +--- Ion charge
      //  |
      //  v
      if (P>0) {
	double ke              = std::max<double>(kE1s[id], FloorEv);
	std::vector<double> RE = ch.IonList[Q]->radRecombCross(ke, id);

	double crs = sVel1 * Eta1 * RE.back() * fac1;
	
	if (scatter_check and recomb_check) {
	  double val = sVel1 * vel * 1.0e-14 * RE.back();
	  recombA[id].add(k, Eta1, val);
	}

	if (DEBUG_CRS) trap_crs(crs);
	
	if (crs > 0.0) {
	  Interact::T t { recomb, Ion, Interact::edef };
	  
	  hCross[id].push_back(XStup(t));
	  hCross[id].back().crs = crs;
	  
	  CProb[id][1] += crs;
	}
      }
      
      // Particle 2 is ION, Particle 1 has ELECTRON
      //
      //  +--- Ion charge
      //  |
      //  v
      if (P>0) {
	double ke              = std::max<double>(kE2s[id], FloorEv);
	std::vector<double> RE = ch.IonList[Q]->radRecombCross(ke, id);

	double crs = sVel2 * Eta2 * RE.back() * fac2;
	
	if (scatter_check and recomb_check) {
	  double val = sVel2 * vel * 1.0e-14 * RE.back();
	  recombA[id].add(k, Eta2, val);
	}

	if (DEBUG_CRS) trap_crs(crs);
	  
	if (crs > 0.0) {
	  Interact::T t { recomb, Interact::edef, Ion };
	  
	  hCross[id].push_back(XStup(t));
	  hCross[id].back().crs = crs;
	  
	  CProb[id][2] += crs;
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
	  double ke              = std::max<double>(kEe1[id], FloorEv);
	  std::vector<double> RE = ch.IonList[Q]->radRecombCross(ke, id);

	  double crs = gVel2 * Eta2 * RE.back() * fac1;
	  
	  if (scatter_check and recomb_check) {
	    double val = sVel2 * vel * 1.0e-14 * RE.back();
	    recombA[id].add(k, Eta2, val);
	  }

	  if (DEBUG_CRS) trap_crs(crs);
	  
	  if (crs > 0.0) {
	    Interact::T t { recomb, Ion, Interact::edef };
	    
	    hCross[id].push_back(XStup(t));
	    hCross[id].back().crs = crs;
	    
	    CProb[id][1] += crs;
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
	double ke = std::max<double>(kEe2[id], FloorEv);
	  std::vector<double> RE = ch.IonList[Q]->radRecombCross(ke, id);
	  
	  double crs = gVel1 * Eta1 * RE.back() * fac2;
	  
	  if (scatter_check and recomb_check) {
	    double val = sVel1 * vel * 1.0e-14 * RE.back();
	    recombA[id].add(k, Eta2, val);
	  }
	  
	  if (DEBUG_CRS) trap_crs(crs);
	  
	  if (crs > 0.0) {
	    Interact::T t { recomb, Interact::edef, Ion };
	    
	    hCross[id].push_back(XStup(t));
	    hCross[id].back().crs = crs;
	    
	    CProb[id][2] += crs;
	  }
	  
      } // end: original recomb algorithm
      
    } // end: recombination

  } // end: outer ionization state loop

  double totalXS = 0.0;
  for (auto & v : CProb[id]) {
    v *= crs_units;
    totalXS += v;
  }

  return totalXS;

  }


  __global__ void energyKernel
(dArray<cudaParticle> in, 
 dArray<cuFP_t> ionion,
 dArray<cuFP_t> ionelc,
 dArray<cuFP_t> elcion,
 dArray<cuFP_t> elcelc,
 dArray<cuFP_t> ionelS,
 dArray<cuFP_t> elionS,
 dArray<cuFP_t> eta1,
 dArray<cuFP_t> eta2,
 dArray<cuFP_t> mol1,
 dArray<cuFP_t> mol2,
 dArray<int>    i1,
 dArray<int>    i2,
 dArray<char>   dZ,
 dArray<char>   dC,
 dArray<char>   dI,
 int pos, int minSP,
 unsigned int stride)
{
  const int tid   = blockDim.x * blockIdx.x + threadIdx.x;

  for (int n=0; n<stride; n++) {
    int i = tid*stride + n;

    if (i < i_1._s) {

      cudaParticle p1 = in._v[i1._v[i]];
      cudaParticle p2 = in._v[i2._v[i]];
    
      // Electron fraction and mean molecular weight for each particle
      //
      cuFP_t Eta1=0.0, Eta2=0.0, Mu1=0.0, Mu2=0.0, Sum1=0.0, Sum2=0.0;

      for (int k=0; k<dZ._s; k++) {
				// Number fraction of ions
	cuFP_t one = p1->datr[dI._v[k]-minSp] / cuda_atomic_weights[dZ._v[k]];
	cuFP_t two = p2->datr[dI._v[k]-minSp] / cuda_atomic_weights[dZ._v[k]];
	
	Eta1 += one * (dC._v[k] - 1);
	Eta2 += two * (dC._v[k] - 1);

	Sum1 += one;
	Sum2 += two;
      }
				// The number of electrons per particle
      Eta1 /= Sum1;
      Eta2 /= Sum2;
				// The molecular weight
      Mu1 = 1.0/Sum1;
      Mu2 = 1.0/Sum2;

      // Ion-ion, ion-electron, and electron-electron relative velocities
      //
      
      cuFP_t eVel0 = 0.0;		// Electron relative velocities
      cuFP_t eVel1 = 0.0;
      cuFP_t eVel2 = 0.0;
      cuFP_t gVel0 = 0.0;		// Scaled velocities for mean-mass algorithm
      cuFP_t gVel1 = 0.0;
      cuFP_t gVel2 = 0.0;
      cuFP_t eVelI = 0.0;		// Ion relative velocity
  
      cuFP_t rel;
      for (int k=0; k<3; k++) {
				// Ion-ion
	rel = p1->vel[i] - p2->vel[i];
	eVelI += rel * rel;
    
				// Electron-electron
	rel = p1->datr[pos+k] - p2->datr[pos+k];
	eVel0 += rel * rel;

				// Electron (p2) and Ion (p1)
	rel = p2->datr[pos+k] - p1->vel[k];
	eVel1 += rel * rel;
				// Electron (p1) and Ion (p2)
	rel = p1->datr[pos+k] - p2->vel[k];
	eVel2 += rel * rel;

				// Scaled electron relative velocity
	if (MeanMass) {
	  rel = p1->datr[pos+k]*sqrt(Eta1) - p2->datr[pos+k]*sqrt(Eta2);
	  gVel0 += rel * rel;
	  rel = p1->datr[pos+k]*sqrt(Eta1) - p2->vel[k];
	  gVel1 += rel * rel;
	  rel = p2->datr[pos+k]*sqrt(Eta2) - p1->vel[k];
	  gVel2 += rel * rel;
	}
      }
    }

    eVel0 = sqrt(eVel0) * cuVunit;
    eVel1 = sqrt(eVel1) * cuVunit;
    eVel2 = sqrt(eVel2) * cuVunit;
    
    // Pick scaled relative velocities for mean-mass algorithm
    if (MeanMass) {
      gVel0 = sqrt(gVel0) * cuVunit;
      gVel1 = sqrt(gVel1) * cuVunit;
      gVel2 = sqrt(gVel2) * cuVunit;
    }

    cuFP_t m1  = Mu1 * cuAmu;
    cuFP_t m2  = Mu2 * cuAmu;
    cuFP_t me  = cuda atomic_weights[0] * cuAmu;
    
    cuFP_t mu0 = m1 * m2 / (m1 + m2);
    cuFP_t mu1 = m1 * me / (m1 + me);
    cuFP_t mu2 = me * m2 / (me + m2);

    ionion[i] = 0.5 * mu0 * eVelI / cuEV;

    ionelc[i] = 0.5 * mu1 * eVel1 / cuEV;
    elcion[i] = 0.5 * mu2 * eVel2 / cuEV;
    elcelc[i] = 0.25 * me * eVel0 / cuEV;

    ionelS[i] = 0.5 * mu1 * gVel1 / cuEV;
    elionS[i] = 0.5 * mu2 * gVel2 / cuEV;

    eta1[i]   = Eta1;
    eta2[i]   = Eta2;

    mol1[i]   = Mu1;
    mol2[i]   = Mu2;
  }
}

void * Collide::collide_thread_cuda(void * arg)
{
  sKeyDmap *Fn  = static_cast<sKeyDmap*>(((thrd_pass_arguments*)arg)->fn  );
  int id        = static_cast<int>      (((thrd_pass_arguments*)arg)->id  );
  
  thread_timing_beg(id);
  
  if (id==0) {
    std::ostringstream sout;
    sout << "Collide::collide: ENTERING collide_thread, T=" << tnow;
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
  thrust::host_vector<int>  cellI, cellN, pairs, i1, i2, cc;
  thrust::host_vector<char> flagI;
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

  hostPartItr hit = host_particles.begin();

  // Species map info
  thrust::host_vector<int> Z, C, I;
  int minSp = std::numeric_limits<int>::max(), maxSp = 0;
  int numNeut = 0, numProt = 0, numIon = 0;
  for (auto s : SpList) {
    speciesKey k = s.first;
    Z.push_back[k.first ];
    C.push_back[k.second];
    I.push_back[s.second];
    minSp = std::min<int>(minSp, s.second);
    maxSp = std::max<int>(maxSp, s.second);
    if (k.second==1) numNeut++;
    else             numIon++;
    if (k.first==1)  numProt++;
  }
  numProt *= numNeut;
  maxSp++;

  // Loop over cells, processing collisions in each cell
  //
  for (unsigned j=0; j<cellist[id].size(); j++ ) {
    
    // The current cell
    //
    pCell *c = cellist[id][j];

    // Skip cell if this time has already been computed
    if (c->time >= tnow) {
      continue;
    }

    int EPSMused = 0;
    
    // Start the effort time
    //
    curcTime[id].reset();
    curcTime[id].start();
    listTime[id].start();
    
    // Number of particles in this cell
    //
    unsigned number = c->bods.size();
    numcntT[id].push_back(number);
    
    // Skip cells with only one particle
    //
    if ( number < 2 ) {
      colcntT[id].push_back(0);
      // Stop timers
      curcTime[id].stop();
      listSoFar[id] = listTime[id].stop();
      // Skip to the next cell
      continue;
    }
    
    listSoFar[id] = listTime[id].stop();
    stat1Time[id].start();
    
    // Energy lost in this cell
    //
    decolT[id] = 0.0;
    decelT[id] = 0.0;
    
    // Compute 1.5 times the mean relative velocity in each MACRO cell
    //
    pCell *samp = c->sample;

    sKeyUmap::iterator it1, it2;

    //
    // Sanity check
    //
    if (samp == 0x0) {
      cout << "Process "  << myid << " in collide: no sample cell"
	   << ", owner="   << c->owner << hex
	   << ", mykey="   << c->mykey
	   << ", mask="    << c->mask  << dec
	   << ", level="   << c->level    
	   << ", Count="   << c->ctotal
	   << ", maxplev=" << c->maxplev;
      if (tree->onFrontier(c->mykey)) cout << ", ON frontier" << endl;
      else cout << ", NOT on frontier" << endl;
      
    }

    cuFP_t crm = 0.0;
    
    if (samp) {
      if (samp->stotal[0]>0.0) {
	for (unsigned k=0; k<3; k++) {
	  crm += (samp->stotal[1+k] - 
		  samp->stotal[4+k]*samp->stotal[4+k]/samp->stotal[0])
	    /samp->stotal[0];}
      }
      crm  = crm>0.0 ? sqrt(2.0*crm) : 0.0;
    }
    
    stat1SoFar[id] = stat1Time[id].stop();
    stat2Time [id].start();
    
    // KE in the cell
    //
    cuFP_t kedsp=0.0;
    if (MFPDIAG) {
      if (c->stotal[0]>0.0) {
	for (unsigned k=0; k<3; k++) 
	  kedsp += 
	    0.5*(c->stotal[1+k] - c->stotal[4+k]*c->stotal[4+k]/c->stotal[0]);
      }
    }
    
    // Timestep for this cell
    //
    cuFP_t tau = dtime / (1<<c->maxplev);

    // Volume in the cell
    //
    cuFP_t volc = c->Volume();
    
    // Mass in the cell
    //
    cuFP_t mass = c->Mass();
    
    // Mass density in the cell
    cuFP_t dens = mass/volc;
    
    if (mass <= 0.0) continue;
    
    // Cell length
    //
    cuFP_t cL = pow(volc, 0.33333333);
    
    
    // Cell initialization (generate cross sections)
    //
    initialize_cell(c, crm, id);

    // Per species quantities
    //
    cuFP_t            meanLambda, meanCollP, totalNsel;
    sKey2Amap         nselM   = generateSelection(c, Fn, crm, tau, id, 
						  meanLambda, meanCollP, 
						  totalNsel);
    
    // Move computation of selection inline here from generateSelectionTrace


    // Ratio of cell size to mean free path
    //
    double mfpCL = cL/meanLambda;

    if (MFPCL) mfpCLdata[id].push_back(mfpCL);
    
    if (TSDIAG) {		// Diagnose time step in this cell
      double vmass;
      vector<double> V1, V2;
      c->Vel(vmass, V1, V2);
      double scale = c->Scale();
      double taudiag = 1.0e40;
      for (int k=0; k<3; k++) {	// Time of flight
	taudiag = min<double>
	  (pHOT::sides[k]*scale/(fabs(V1[k]/vmass)+sqrt(V2[k]/vmass)+1.0e-40), 
	   taudiag);
      }
      
      int indx = (int)floor(log(taudiag/tau)/log(4.0) + 5);
      if (indx<0 ) indx = 0;
      if (indx>10) indx = 10;
      tdiagT[id][indx]++;
    }
    
    if (VOLDIAG) {
      if (c->level<nbits) {
	VcntT[id][c->level]++;
	VdblT[id][c->level*nvold+0] += dens;
	VdblT[id][c->level*nvold+1] += 1.0 / mfpCL;
	VdblT[id][c->level*nvold+2] += meanCollP;
	VdblT[id][c->level*nvold+3] += crm*tau / cL;
	VdblT[id][c->level*nvold+4] += number;
	VdblT[id][c->level*nvold+5] += number*number;
      }
    }
    // Selection number per particle
    if (MFPDIAG)
      tselnT[id].push_back(totalNsel/number);
    
    stat2SoFar[id] = stat2Time[id].stop();
    
    // Species map for collisions

    unsigned colc = 0;
    std::map<speciesKey, std::vector<unsigned long> > bmap;

    // No collisions, primarily for testing . . .
    //
    if (DRYRUN) continue;
    
    collTime[id].start();
    collCnt[id]++;
    
    initTime[id].start();
    initSoFar[id] = initTime[id].stop();
      
    // Need to load bodies into particle host_vector
    //
    for (size_t k=0; k<c->bods.size(); k++) {
      ParticlesHtoD(Particles()[c->bods[k]], (*hit++), minSp, maxSp);
    }
  }

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

  // Create interaction energy files
  //
  thrust::device_vector<cuFP_t> ionion(pairs.size());
  thrust::device_vector<cuFP_t> ionelc(pairs.size());
  thrust::device_vector<cuFP_t> elcion(pairs.size());
  thrust::device_vector<cuFP_t> elcelc(pairs.size());

  // Okay, now compute all of the pieces
  //
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, c0->cudaDevice);

  unsigned int N        = i1.size(); // Number of pairs
  unsigned int stride   = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
  unsigned int gridSize = N/BLOCK_SIZE/stride;

  if (N > gridSize*BLOCK_SIZE*stride) gridSize++;

  // Compute energies
  //
  energyKernel<<<gridSize, BLOCK_SIZE>>>
    (toKernel(d_part),
     toKernel(ionion), toKernel(ionelc),
     toKernel(elcion), toKernel(elcelc),
     toKernel(d_i1), toKernel(d_i2), eePos, stride);


  // Create ion-electron cross-section lists
  //
  thrust::device_vector<cuFP_t> ieXC(pairs.size()*4);
  thrust::device_vector<cuFP_t> eiXC(pairs.size()*4);

  // Compute ion-electron cross sections
  //
  ionElecKernel<<<gridSize, BLOCK_SIZE>>>
    (toKernel(ionion), toKernel(ionelc),
     toKernel(elcion), toKernel(elcelc),
     toKernel(ieXC),   toKernel(eiXC),
     eePos, stride);

  // Compute molecular weight and eta
  //
  N        = d_part.size();	// Number of particles
  stride   = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
  gridSize = N/BLOCK_SIZE/stride;

				// These do not need copying back
  thrust::device_vector<int> d_ZZ), d_C(C), d_I(I);
  thrust::device_vector<cuFP_t> d_mol(N), d_eta(N);

  etaKernel<<<gridSize, BLOCK_SIZE>>>
    (toKernel(d_part), toKernel(d_mol), toKernel(d_eta),
     toKernel(d_Z), toKernel(d_C), toKernel(d_I), minSp, maxSp, stride);


  // Initialize cells
  //
  N        = d_cellI.size();	// Number of cells
  stride   = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
  gridSize = N/BLOCK_SIZE/stride;

				// These do not need copying back
  thrust::device_vector<cuFP_t> d_meanM(N), d_Ivel2(N), d_Evel2(N);
  thrust::device_vector<cuFP_t> d_PiProb(N*4), d_ABrate(N*4);


  cellInitKernel<<<gridSize, BLOCK_SIZE>>>
    (toKernel(d_part),
     toKernel(d_meanM), toKernel(d_Ivel2), toKernel(d_Evel2),
     toKernel(d_PiProb), toKernel(d_ABrate),
     toKernel(d_cellI), toKernel(d_cellN),
     toKernel(d_Z), toKernel(d_C), toKernel(d_I),
     minSp, stride);

  // Compute the cross sections for each interaction en masse
  //
  N        = i1.size();	// Number of pairs
  stride   = N/BLOCK_SIZE/deviceProp.maxGridSize[0] + 1;
  gridSize = N/BLOCK_SIZE/stride;

				// These do not need copying back
  unsigned int totalXCsize =
    numCudaInterTypesNeut*numNeut +
    numCudaInterTypesProt*numProt +
    numCudaInterTypesIon*numIon   ;
  thrust::device_vector<cuFP_t> d_cross(N*totalXCsize);
  thrust::device_vector<uchar4> d_xspec(N*totalXCsize);
  thrust::device_vector<uchar>  d_xtype(N*totalXCsize);
  thrust::device_vector<int>    d_flagI(flagI);

  crossSectionKernel<<<gridSize, BLOCK_SIZE>>>
    (toKernel(d_part),
     toKernel(d_cross),  toKernel(d_xspec),  toKernel(d_xtype), toKernel(d_crossT),
     toKernel(d_i1),     toKernel(d_i2),     toKernel(d_cc),
     toKernel(d_meanM),  toKernel(d_Ivel2),  toKernel(d_Evel2),
     toKernel(d_PiProb), toKernel(d_ABrate), toKernel(d_flagI), 
     toKernel(d_Z), toKernel(d_C), toKernel(d_I),
     toKernel(xsc_H), toKernel(xsc_He), toKernel(xsc_pH), toKernel(xsc_pHe),
     totalXCsize, minSp, eePos, stride);

  // Compute the total cross section per cell
  //

  // Do the interactions
  //

  // Finally, copy back particles to host
  // 

  // Loop over cells, processing collisions in each cell
  //
  for (unsigned j=0; j<cellist[id].size(); j++ ) {

    // Cache acceptance fraction for scaling MFP for time step
    // selection
    //
    if (acceptCount > 0) {
      double scale = static_cast<double>(totalCount)/acceptCount; 
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
  
    double tmass = 0.0;
    vector<double> velm(3, 0.0), velm2(3, 0.0);
    for (unsigned j=0; j<number; j++) {
      Particle* p = tree->Body(c->bods[j]);
      for (unsigned k=0; k<3; k++) {
	velm[k]  += p->mass*p->vel[k];
	velm2[k] += p->mass*p->vel[k]*p->vel[k];
      }
      tmass += p->mass;
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
      double cL = pow(volc, 0.33333333);
      double Kn = meanLambda/cL;
      double St = cL/fabs(tau*sqrt(fabs(kedsp))+1.0e-18);
      for (unsigned j=0; j<number; j++) {
	Particle* p = tree->Body(c->bods[j]);
	if (use_Kn>=0) p->datr[use_Kn] = Kn;
	if (use_St>=0) p->datr[use_St] = St;
      }
    }
    
    // Record effort per particle in microseconds
    //
    curcSoFar[id] = curcTime[id].stop();
    double tt = curcSoFar[id];
    if (EFFORT) {
      if (effortAccum) 
	effortNumber[id].push_back(pair<long, unsigned>(tt, number));
      double effort = static_cast<double>(tt)/number;
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
