#include <filesystem>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <memory>
#include <limits>
#include <string>

#include <CylEXP.H>

CylEXP::CylEXP(int numr, int lmax, int mmax, int nord,
	       double ascale, double hscale, int Nodd,
	       std::string cachename) :
  EmpCylSL(numr, lmax, mmax, nord, ascale, hscale, Nodd, cachename)
{
  // Initialize storage for the multistep arrays
  //
  differC1.resize(nthrds);
  differS1.resize(nthrds);

  for (auto & v : differC1) v.resize(multistep+1);
  for (auto & v : differS1) v.resize(multistep+1);
  
  for (int nth=0; nth<nthrds; nth++) {
    for (unsigned M=0; M<=multistep; M++) {
      differC1[nth][M].resize(mmax+1, nord);
      differS1[nth][M].resize(mmax+1, nord);
      differC1[nth][M].setZero();
      differS1[nth][M].setZero();
    }
  }
    
  unsigned sz = (multistep+1)*(mmax+1)*nord;
  workC1.resize(sz);
  workC .resize(sz);
  workS1.resize(sz);
  workS .resize(sz);
}

void CylEXP::multistep_reset()
{
}

void CylEXP::multistep_update_begin()
{
				// Clear the update matricies
  for (int nth=0; nth<nthrds; nth++) {
    for (unsigned M=mfirst[mdrft]; M<=multistep; M++) {
      differC1[nth][M].setZero();
      differS1[nth][M].setZero();
    }
  }
}

void CylEXP::multistep_update_finish()
{
  unsigned offset0, offset1;
  unsigned sz = (multistep - mfirst[mdrft] + 1)*(MMAX+1)*rank3;
  for (unsigned j=0; j<sz; j++) 
    workC1[j] = workC[j] = workS1[j] = workS[j] = 0.0;

				// Combine the update matricies
  for (unsigned M=mfirst[mdrft]; M<=multistep; M++) {

    offset0 = (M - mfirst[mdrft])*(MMAX+1)*rank3;

    for (int mm=0; mm<=MMAX; mm++) {
      
      offset1 = mm*rank3;

      for (int nth=0; nth<nthrds; nth++)
	for (int k=0; k<rank3; k++) 
	  workC1[offset0+offset1+k] += differC1[nth][M](mm, k);

      if (mm) {
	for (int nth=0; nth<nthrds; nth++)
	  for (int k=0; k<rank3; k++) 
	    workS1[offset0+offset1+k] += differS1[nth][M](mm, k);
      }
    }
  }

  MPI_Allreduce (&workC1[0], &workC[0], sz,
		 MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce (&workS1[0], &workS[0], sz, 
		 MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  //  +--- Deep debugging
  //  |
  //  v
  if (false and myid==0) {
    std::string filename = runtag + ".differ_cyl";
    std::ofstream out(filename, ios::app);
    if (out) {
      out << std::string(13+16*rank3, '-') << std::endl;
      out << "# T=" << tnow << " mstep=" << mstep << std::endl;
      for (unsigned M=mfirst[mdrft]; M<=multistep; M++) {
	offset0 = (M - mfirst[mdrft])*(MMAX+1)*rank3;
	for (int mm=0; mm<=MMAX; mm++) {
	  offset1 = mm*rank3;
	  out << std::setw(5) << M << " C " << std::setw(5) << mm;
	  for (int nn=0; nn<rank3; nn++) 
	    out << std::setw(16) << workC[offset0+offset1+nn];
	  out << std::endl;
	  if (mm) {
	    out << std::setw(5) << M << " S " << std::setw(5) << mm;
	    for (int nn=0; nn<rank3; nn++) 
	      out << std::setw(16) << workS[offset0+offset1+nn];
	    out << std::endl;
	  }
	  break;
	}
      }
      out << std::string(13+16*rank3, '-') << std::endl;
      for (int mm=0; mm<=MMAX; mm++) {
	offset1 = mm*rank3;
	out << std::setw(5) << " *** " << " C " << std::setw(5) << mm;
	for (int nn=0; nn<rank3; nn++) 
	  out << std::setw(16) << accum_cos[mm][nn];
	out << std::endl;
	if (mm) {
	  out << std::setw(5) << " *** " << " S " << std::setw(5) << mm;
	  for (int nn=0; nn<rank3; nn++) 
	    out << std::setw(16) << accum_sin[mm][nn];
	  out << std::endl;
	}
	break;
      }
      out << std::string(13+16*rank3, '-') << std::endl;
      out << std::string(13+16*rank3, '-') << std::endl;
    } else {
      std::cout << "Error opening test file <test_differ.sph>"
		<< std::endl;
    }
  }


  for (unsigned M=mfirst[mdrft]; M<=multistep; M++) {

    offset0 = (M - mfirst[mdrft])*(MMAX+1)*rank3;

    for (int mm=0; mm<=MMAX; mm++) {
      
      offset1 = mm*rank3;

      for (int nn=0; nn<rank3; nn++) 
	cosN(M)[0][mm][nn] += workC[offset0+offset1+nn];

      if (mm) {
	for (int nn=0; nn<rank3; nn++) 
	  sinN(M)[0][mm][nn] += workS[offset0+offset1+nn];
      }
    }
  }
}

void CylEXP::multistep_update(int from, int to, double r, double z, double phi, double mass, int id)
{
  double rr = sqrt(r*r+z*z);
  if (rr/ASCALE>Rtable) return;

  double msin, mcos;

  double norm = -4.0*M_PI;
  
  get_pot(vc[id], vs[id], r, z);

  for (int mm=0; mm<=MMAX; mm++) {

    mcos = cos(phi*mm);
    msin = sin(phi*mm);

    for (int nn=0; nn<rank3; nn++) {
      double hold = norm * mass * mcos * vc[id](mm, nn);
      differC1[id][from](mm, nn) -= hold;
      differC1[id][to  ](mm, nn) += hold;

      if (mm>0) {
	hold = norm * mass * msin * vs[id](mm, nn);
	differS1[id][from](mm, nn) -= hold;
	differS1[id][to  ](mm, nn) += hold;
      }
    }
  }

}



void CylEXP::compute_multistep_coefficients(unsigned mlevel)
{
				// Clean coefficient matrix
				// 
  for (int mm=0; mm<=MMAX; mm++) {
    for (int nn=0; nn<rank3; nn++) {
      accum_cos[mm][nn] = 0.0;
      if (mm) accum_sin[mm][nn] = 0.0;
    }
  }
  
				// Interpolate to get coefficients above the
				// current active level
  for (unsigned M=0; M<mfirst[mdrft]; M++) {
    
    double numer = static_cast<double>(mdrft            - dstepL[M][mdrft]);
    double denom = static_cast<double>(dstepN[M][mdrft] - dstepL[M][mdrft]);

    double b = numer/denom;	// Interpolation weights
    double a = 1.0 - b;

    //  +--- Deep debugging
    //  |
    //  v
    if (false and myid==0) {
      std::cout << std::left << std::fixed
		<< "CYL INTERP M=" << std::setw(2) << M
		<< " mstep=" << std::setw(3) << mstep
		<< " mdrft=" << std::setw(3) << mdrft
		<< " a=" << std::setw(16) << a
		<< " b=" << std::setw(16) << b
		<< std::endl << std::right;
    }

    for (int mm=0; mm<=MMAX; mm++) {
      for (int nn=0; nn<rank3; nn++) {
	accum_cos[mm][nn] += a*cosL(M)[0][mm][nn] + b*cosN(M)[0][mm][nn];
	if (mm)
	  accum_sin[mm][nn] += a*sinL(M)[0][mm][nn] + b*sinN(M)[0][mm][nn];
      }
    }

    if (false and myid==0) {
      std::cout << "CYL interpolate:"
		<< " M="     << std::setw( 4) << M
		<< " mstep=" << std::setw( 4) << mstep 
		<< " mdrft=" << std::setw( 4) << mdrft 
		<< " minS="  << std::setw( 4) << dstepL[M][mdrft]
		<< " maxS="  << std::setw( 4) << dstepN[M][mdrft]
		<< " a="     << std::setw(14) << a 
		<< " b="     << std::setw(14) << b 
		<< " L00="   << std::setw(14) << cosL(M)[0][0][0]
		<< " N00="   << std::setw(14) << cosN(M)[0][0][0]
		<< " c00="   << std::setw(14) << accum_cos[0][0]
		<< std::endl;
    }

    // Sanity debug check
    if (a<0.0 && a>1.0) {
      cout << "Process " << myid << ": interpolation error in multistep [a]" << endl;
    }
    if (b<0.0 && b>1.0) {
      cout << "Process " << myid << ": interpolation error in multistep [b]" << endl;
    }
  }
				// Add coefficients at or above this level
				// 
  for (unsigned M=mfirst[mdrft]; M<=multistep; M++) {

    //  +--- Deep debugging
    //  |
    //  v
    if (false and myid==0) {
      std::cout << std::left << std::fixed
		<< "CYL FULVAL M=" << std::setw(2) << M
		<< " mstep=" << std::setw(3) << mstep
		<< " mdrft=" << std::setw(3) << mdrft
		<< std::endl << std::right;
    }

    for (int mm=0; mm<=MMAX; mm++) {
      for (int nn=0; nn<rank3; nn++) {
	accum_cos[mm][nn] += cosN(M)[0][mm][nn];
	if (mm)
	  accum_sin[mm][nn] += sinN(M)[0][mm][nn];
      }
    }
  }

  coefs_made = vector<short>(multistep+1, true);
}


void CylEXP::multistep_debug()
{
  for (int n=0; n<numprocs; n++) {
    if (n==myid) {
      cout << "Process " << myid << endl
	   << "   accum_cos[0][0]="
	   << accum_cos[0][0] << endl
	   << "   accum_sin[1][0]="
	   << accum_sin[1][0] << endl;

      cout.setf(ios::scientific);
      int c = cout.precision(2);
      for (unsigned M=0; M<=multistep; M++) {
	cout << "   M=" << M << ": #part=" << howmany[M] << endl;

	cout << "   M=" << M << ", m=0, C_N: ";
	for (int j=0; j<NORDER; j++) 
	  cout << setprecision(2) << setw(10) << cosN(M)[0][0][j];
	cout << endl;

	cout << "   M=" << M << ", m=0, C_L: ";
	for (int j=0; j<NORDER; j++) 
	  cout << setprecision(2) << setw(10) << cosL(M)[0][0][j];
	cout << endl;

	cout << "   M=" << M << ", m=1, C_N: ";
	for (int j=0; j<NORDER; j++) 
	  cout << setprecision(2) << setw(10) << cosN(M)[0][1][j];
	cout << endl;

	cout << "   M=" << M << ", m=1, C_L: ";
	for (int j=0; j<NORDER; j++) 
	  cout << setprecision(2) << setw(10) << cosL(M)[0][1][j];
	cout << endl;

	cout << "   M=" << M << ", m=1, S_N: ";
	for (int j=0; j<NORDER; j++) 
	  cout << setprecision(2) << setw(10) << sinN(M)[0][1][j];
	cout << endl;

	cout << "   M=" << M << ", m=1, S_L: ";
	for (int j=0; j<NORDER; j++) 
	  cout << setprecision(2) << setw(10) << sinL(M)[0][1][j];
	cout << endl;

	cout.precision(c);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

