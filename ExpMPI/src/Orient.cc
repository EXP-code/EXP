#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>

#include "expand.h"
#include <localmpi.h>

#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

#include <Orient.H>


Matrix return_euler_slater(double PHI, double THETA, double PSI, int BODY);

Orient::Orient(int n, int nwant, double Einit, unsigned Oflg, unsigned Cflg,
	       string Logfile) : gen(11, 20), gauss(0.0, 1.0, &gen) 
{
  keep = n;
  current = 0;
  many = nwant;
  Egrad = 0.0;
  Ecurr = Einit;
  oflags = Oflg;
  cflags = Cflg;
  logfile = Logfile;
  Nlast = 0;
  linear = false;

				// Work vectors
  axis1.setsize(1, 3);
  center1.setsize(1, 3);
  sumY.setsize(1, 3);
  sumXY.setsize(1, 3);
  sumY2.setsize(1, 3);
  slope.setsize(1, 3);

				// Center and axis
  axis.setsize(1, 3);
  center.setsize(1, 3);
  center0.setsize(1, 3);
  cenvel0.setsize(1, 3);
  axis.zero();
  center.zero();
  center0.zero();
  cenvel0.zero();

				// Set up identity
  body.setsize(1, 3, 1, 3);
  body.zero();
  body[1][1] = body[2][2] = body[3][3] = 1.0;
  orig = body;

				// Check for previous state
  int in_ok;
  double *in1 = new double [3];
  double *in2 = new double [3];

  if (myid==0) {		// Master does the reading

    ifstream in(logfile.c_str());
    

				// If the logfile is there, read it
    if (in) {
      in.close();

				// Backup old file
      string backupfile = logfile + ".bak";
      if (rename(logfile.c_str(), backupfile.c_str())) {
	cerr << "Orient: error making backup file <" 
	     << backupfile << ">\n";
	MPI_Abort(MPI_COMM_WORLD, 42);
      }

				// Open new output stream for writing
      ofstream out(logfile.c_str());
      if (!out) {
	cerr << "Orient: error opening new log file <" 
	     << logfile << "> for writing\n";
	MPI_Abort(MPI_COMM_WORLD, 43);
      }
	  
      // Open old file for reading
      in.open(backupfile.c_str());
      if (!in) {
	ostringstream message;
	cerr << "Orient: error opening original log file <" 
	     << backupfile << "> for reading\n";
	MPI_Abort(MPI_COMM_WORLD, 44);
      }


      double time;
      int nused;

      in_ok = 1;		// Signal slave: OK

      MPI_Bcast(&in_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);

      const int cbufsiz = 16384;
      char *cbuffer = new char [cbufsiz];
	
				// Look for data and write it while
				// accumlating data for averaging
      while (in) {

	in.getline(cbuffer, cbufsiz);
	if (in.rdstate() & (ios::failbit | ios::eofbit)) break;

	istringstream line(cbuffer);

	line >> time;
	if (tnow < time || fabs(tnow-time)<1.0e-8) break;

	out << cbuffer << "\n";

	line >> Ecurr;
	line >> nused;
	line >> axis[1];
	line >> axis[2];
	line >> axis[3];
	line >> center0[1];
	line >> center0[2];
	line >> center0[3];

	if (oflags & AXIS) {
	  sumsA.push_back(axis);
	  if (sumsA.size() > keep) sumsA.pop_front();
	}

	if (oflags & CENTER) {
	  sumsC.push_back(center0);
	  if (sumsC.size() > keep) sumsC.pop_front();
	}

      }

      cout << " Orient: current log=" << logfile << "  backup=" << backupfile << endl;

      cout << " Orient: cached time=" << time << "  Ecurr= " << Ecurr << endl;

      cout << " Orient: axis master (cache size=" << sumsA.size() << "): " 
	   << axis[1] << ", "
	   << axis[2] << ", "
	   << axis[3] << endl;

      cout << " Orient: center master (cache size=" << sumsC.size() << "): " 
	   << center0[1] << ", "
	   << center0[2] << ", "
	   << center0[3] << endl;

      MPI_Bcast(&Ecurr, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      int howmany = max<int>(sumsA.size(), sumsC.size());
      MPI_Bcast(&howmany, 1, MPI_INT, 0, MPI_COMM_WORLD);

      for (int k=0; k<howmany; k++) {

	if (oflags & AXIS) {
	  for (int j=0; j<3; j++) in1[j] = sumsA[k][j+1];
	  MPI_Bcast(in1, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}

	if (oflags & CENTER) {
	  for (int j=0; j<3; j++) in2[j] = sumsC[k][j+1];
	  MPI_Bcast(in2, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
      }


    } else {

      in_ok = 0;		// Signal slave: NO VALUES

      MPI_Bcast(&in_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }


  } else {

				// Get state from Master

    MPI_Bcast(&in_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (in_ok) {

      MPI_Bcast(&Ecurr, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      int howmany;
      MPI_Bcast(&howmany, 1, MPI_INT, 0, MPI_COMM_WORLD);

      for (int k=0; k<howmany; k++) {

	if (oflags & AXIS) {
	  MPI_Bcast(in1, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	  for (int j=0; j<3; j++) axis[j+1] = in1[j];
	  sumsA.push_back(axis);
	}

	if (oflags & CENTER) {
	  MPI_Bcast(in2, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	  for (int j=0; j<3; j++) center0[j+1] = in2[j];
	  sumsC.push_back(center0);
	}

      }

    }

  }

  delete [] in1;
  delete [] in2;

  if (in_ok) {

    if (oflags & CENTER && sumsC.size()>0)
      center = sumsC[sumsC.size()-1];
	
    if (oflags & AXIS && sumsA.size()>0) {

      axis = sumsA[sumsA.size()-1];

      double phi = atan2(axis[2], axis[1]);
      double theta = -acos(axis[3]/sqrt(axis*axis));
      double psi = 0.0;

      body = return_euler_slater(phi, theta, psi, 0);
      orig = return_euler_slater(phi, theta, psi, 1);
    }

  }
      
}


void Orient::accumulate(double time, vector<Particle> *p, 
			double *cen, double *com, double *cov)
{
  if (linear) {
      center = center0;
      center0 += cenvel0*dtime;
      if (myid==0) write_log(time, 0.0, 0.0, com);
      return;
  }

  double energy, mass;
  double Emin1= 1.0e20, Emin0= 1.0e20;
  double Emax1=-1.0e20, Emax0=-1.0e20;

  angm.clear();

  int nbodies = p->size();
  for (int i=0; i<nbodies; i++) {

    energy = (*p)[i].pot;
    
    if (cflags & KE) energy += 0.5*((*p)[i].vel[0]*(*p)[i].vel[0] + 
				    (*p)[i].vel[1]*(*p)[i].vel[1] + 
				    (*p)[i].vel[2]*(*p)[i].vel[2]);

    if (cflags & EXTERNAL) energy += (*p)[i].potext;

    Emin1 = min<double>(energy, Emin1);
    Emax1 = max<double>(energy, Emax1);
    
    if (energy < Ecurr) {

      mass = (*p)[i].mass;

      t.E = energy;

      t.M = mass;

      t.L[1] = mass*(((*p)[i].pos[1]-com[1]-cen[1])*((*p)[i].vel[2] - cov[2]) -
		     ((*p)[i].pos[2]-com[2]-cen[2])*((*p)[i].vel[1] - cov[1]));

      t.L[2] = mass*(((*p)[i].pos[2]-com[2]-cen[2])*((*p)[i].vel[0] - cov[0]) -
		     ((*p)[i].pos[0]-com[0]-cen[0])*((*p)[i].vel[2] - cov[2]));

      t.L[3] = mass*(((*p)[i].pos[0]-com[0]-cen[0])*((*p)[i].vel[1] - cov[1]) -
		     ((*p)[i].pos[1]-com[1]-cen[1])*((*p)[i].vel[0] - cov[0]));

      /*
	t.R[1] = (*p)[i].pot*((*p)[i].pos[0] - com[0]);
	t.R[2] = (*p)[i].pot*((*p)[i].pos[1] - com[1]);
	t.R[3] = (*p)[i].pot*((*p)[i].pos[2] - com[2]);
      */

      t.R[1] = mass*((*p)[i].pos[0] - com[0]);
      t.R[2] = mass*((*p)[i].pos[1] - com[1]);
      t.R[3] = mass*((*p)[i].pos[2] - com[2]);
      angm.push_back(t);
    }
  }

  // Sanity check

  int size0=0, size1=(int)angm.size();

  // Propagate minimum energy and current cached low energy particles
  // with nodes
  MPI_Allreduce(&Emin1, &Emin0, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&Emax1, &Emax0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&size1, &size0, 1, MPI_INT,    MPI_SUM, MPI_COMM_WORLD);

  if (size0 == 0) {
    
    if (myid==0) {
      cout << " Orient: Ecurr error, current value=" << Ecurr 
	   << " but Min(Phi)=" << Emin0 << endl;
    }
    
    Ecurr = Emin0*0.98;
    
    for (int i=0; i<nbodies; i++) {
      
      energy = (*p)[i].pot;
      
      if (cflags & KE) 
	energy += 0.5*(
		       ((*p)[i].vel[0] - cov[0])*((*p)[i].vel[0] - cov[0]) +
		       ((*p)[i].vel[1] - cov[1])*((*p)[i].vel[1] - cov[1]) + 
		       ((*p)[i].vel[2] - cov[2])*((*p)[i].vel[2] - cov[2])
		       );
      
      if (cflags & EXTERNAL) energy += (*p)[i].potext;

      if (energy < Ecurr) {
	
	mass = (*p)[i].mass;
	
	t.E = energy;
	
	t.M = mass;

	t.L[1] = mass*(((*p)[i].pos[1]-com[1]-cen[1])*((*p)[i].vel[2]-cov[2]) -
		       ((*p)[i].pos[2]-com[2]-cen[2])*((*p)[i].vel[1]-cov[1]));

	t.L[2] = mass*(((*p)[i].pos[2]-com[2]-cen[2])*((*p)[i].vel[0]-cov[0]) -
		       ((*p)[i].pos[0]-com[0]-cen[0])*((*p)[i].vel[2]-cov[2]));

	t.L[3] = mass*(((*p)[i].pos[0]-com[0]-cen[0])*((*p)[i].vel[1]-cov[1]) -
		       ((*p)[i].pos[1]-com[1]-cen[1])*((*p)[i].vel[0]-cov[0]));
	
	/*
	  t.R[1] = (*p)[i].pot*((*p)[i].pos[0] - com[0]);
	  t.R[2] = (*p)[i].pot*((*p)[i].pos[1] - com[1]);
	  t.R[3] = (*p)[i].pot*((*p)[i].pos[2] - com[2]);
	*/
	  
	t.R[1] = mass*((*p)[i].pos[0] - com[0]);
	t.R[2] = mass*((*p)[i].pos[1] - com[1]);
	t.R[3] = mass*((*p)[i].pos[2] - com[2]);
	angm.push_back(t);
      }
    }
  }
				// Compute values for this step
  axis1.zero();
  center1.zero();
  double mtot=0.0, mtot1=0.0, dE=1e20;
  vector<EL3>::iterator i;
  for (i = angm.begin(); i != angm.end(); i++) {
    axis1 += i->L;
    center1 += i->R;
    mtot1 += i->M;
  }

  int cnum = angm.size();
  Vector inA = axis1;
  Vector inC = center1;

  used = 0;
  axis1.zero();
  center1.zero();

				// Share stuff between nodes
  MPI_Allreduce(&cnum, &used, 1, 
		MPI_INT,    MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(inA.array(0, 2), axis1.array(0, 2), 3, 
		MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce(&mtot1, &mtot, 1, 
		MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(inC.array(0, 2), center1.array(0, 2), 3, 
		MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (used && Nlast) {
				// Don't divide by zero
    if (used != Nlast) 
      Egrad = (Ecurr - Elast)/(used - Nlast);

    dE = (double)(many - used) * Egrad;

    if (fabs(dE) <= 1.0e-10)
      dE = (Ecurr - Emin0)*0.01*gauss();

  }

				// Push current value onto stack
  if (used) {

    Nlast = used;
    Elast = Ecurr;

    axis1   /= mtot;
    center1 /= mtot;
    if (oflags & AXIS)   sumsA.push_back(axis1);
    if (oflags & CENTER) sumsC.push_back(center1);


    if ((cflags & DIAG) && myid==0) {
      cout << " Orient info: " << used << " particles used, Ecurr=" << Ecurr 
	   << " Center=" 
	   << center1[1] << ", "
	   << center1[2] << ", "
	   << center1[3] << endl;
    }

				// Compute delta energy for next step


    if (sumsA.size() > keep) {

      sumsA.pop_front();

      double x;
      int i=0;

      sumX = 0.0;
      sumX2 = 0.0;
      sumY.zero();
      sumXY.zero();
      sumY2.zero();

      int N = sumsC.size();

      deque<Vector>::iterator j;
      for (j = sumsA.begin(); j != sumsA.end(); j++) {
	x = time - dtime*(N-1-i);
	sumX += x;
	sumX2 += x*x;
	sumY += *j;
	sumXY += *j * x;
	sumY2 += *j & *j;

	i++;
      }

				// Linear least squares estimate for axis

      slope = (sumXY*N - sumX*sumY)/(sumX2*N - sumX*sumX);
      axis = (sumX2*sumY - sumX*sumXY)/(sumX2*N - sumX*sumX) + slope*time;

      i = 0;
      sigA = 0.0;
      for (j = sumsA.begin(); j != sumsA.end(); j++) {
	x = -dtime*(N-1-i);
	sigA += (*j - axis - slope*x)*(*j - axis - slope*x);

	i++;
      }
      sigA /= i;

      double phi = atan2(axis[2], axis[1]);
      double theta = -acos(axis[3]/sqrt(axis*axis));
      double psi = 0.0;

      body = return_euler_slater(phi, theta, psi, 0);
      orig = return_euler_slater(phi, theta, psi, 1);
    }


    if (sumsC.size() > 1) {

      if (sumsC.size() > keep) sumsC.pop_front();

      double x;
      int i=0;
      sumX = 0.0;
      sumX2 = 0.0;
      sumY.zero();
      sumXY.zero();
      sumY2.zero();

      int N = sumsC.size();

      deque<Vector>::iterator j;
      for (j = sumsC.begin(); j != sumsC.end(); j++) {
	x = time - dtime*(N-1-i);
	sumX += x;
	sumX2 += x*x;
	sumY += *j;
	sumXY += *j * x;
	sumY2 += *j & *j;

	i++;

	if ((cflags & DIAG) && myid==0)
	  cout << " Orient debug i=" << i << ":"
	       << "      t=" << setw(15) << x
	       << "      x=" << setw(15) << (*j)[1]
	       << "      y=" << setw(15) << (*j)[2]
	       << "      z=" << setw(15) << (*j)[3]
	       << "   SumX=" << setw(15) << sumX 
	       << "  SumX2=" << setw(15) << sumX2 
	       << "  Delta=" << setw(15) << sumX2*i - sumX*sumX 
	       << endl;
      }
				// Linear least squares estimate for center

      slope = (sumXY*N - sumX*sumY)/(sumX2*N - sumX*sumX);
      center = (sumX2*sumY - sumX*sumXY)/(sumX2*N - sumX*sumX) + slope*time;

      i = 0;
      sigC = 0.0;
      sigCz = 0.0;
      for (j = sumsC.begin(); j != sumsC.end(); j++) {
	x = -dtime*(N-1-i);
	sigC += (*j - center - slope*x)*(*j - center - slope*x);
	sigCz += 
	  ( (*j)[3] - center[3] - slope[3]*x ) *
	  ( (*j)[3] - center[3] - slope[3]*x );
	i++;
      }
      sigC  /= i;
      sigCz /= i;

    }

    if (keep>1) {
      if (sumsC.size()>1) {
	double factor = (double)((int)sumsC.size() - keep)/keep;
	factor = factor*factor;
	center = center0*factor + center*(1.0 - factor);
      } else
	center = center0;
    } else
      center = center1;

				// Increment initial center according 
				// to user specified velocity
    center0 += cenvel0*dtime;

    if ((cflags & DIAG) && myid==0) {
      cout << "===================================================" << endl
	   << " Orient info: size=" << sumsC.size()
	   << "  SumX=" << sumX << " SumX2=" << sumX2 << endl
	   << "  SumY="
	   << sumY[1] << " "
	   << sumY[2] << " "
	   << sumY[3] << endl
	   << "  SumXY="
	   << sumXY[1] << " "
	   << sumXY[2] << " "
	   << sumXY[3] << endl
	   << "  SumY2="
	   << sumY2[1] << " "
	   << sumY2[2] << " "
	   << sumY2[3] << endl
	   << "  slope="
	   << slope[1] << " "
	   << slope[2] << " "
	   << slope[3] << endl
	   << "  center=" 
	   << center[1] << " "
	   << center[2] << " "
	   << center[3] << endl
	   << "===================================================" << endl;
    }
    
    // Will force center to be last minimum energy average
    // rather than pooled average of the last "keep" states
    if (keep == 0) {
      for (int i=1; i<=3; i++) center[i] = (*sumsC.begin())[i];
      sumsC.pop_front();
    }
    
  }

  // Energy for next iteration
  // =======================================================
  // Will use the secant method computation from above to
  // update.  This makes sure that estimate stays in bounds
  // =======================================================
  //
  if (Ecurr+dE < Emin0 || Ecurr+dE > Emax0) {

    int dtype = 0;
				// Ecurr error, algorithm failure!!
    if (Ecurr < Emin0 || Ecurr > Emax0) {
      if (myid==0)
	cerr << "\nOrient: Ecurr pre-step out of bounds!!! dE=" << dE 
	     << "  Emin0=" << Emin0 
	     << "  Emax0=" << Emax0 
	     << "  Ecurr=" << Ecurr << "\n";
      Ecurr = 0.99*Emin0;
    }
				// Try to take a smaller step
				// using secant gradient
    else if (Ecurr+0.2*dE > Emin0 && Ecurr+0.2*dE < Emax0) {
      Ecurr += 0.2*dE;
      dtype = 1;
    }
    
    else {
				// Final strategy: take small steps
				// based on distance to minimum
      if (used > many)
	Ecurr += 0.06*(Emin0-Ecurr);
      else
	Ecurr -= 0.10*(Emin0-Ecurr);

      dtype = 2;
    } 
    
				// Ecurr error, algorithm failure!!
    if (Ecurr < Emin0 && myid==0) {
      cerr << "Orient: Ecurr post-step out of bounds: "
	   << "  dE=" << dE 
	   << "  Emin0=" << Emin0 
	   << "  Emax0=" << Emax0 
	   << "  Ecurr=" << Ecurr;

      switch (dtype) {
      case 0:
	cout << "  Case: out of bounds to start\n";
	break;
      case 1: 
	cout << "  Case: small gradient step\n";
	break;
      case 2:
	cout << "  Case: small distance from center step\n";
	break;
      }
    }

  }
  else				// Secant method: use approx. gradient
    Ecurr += dE;
  

    if (myid==0) write_log(time, Egrad, dE, com);
}

void Orient::write_log(double time, double Egrad, double dE, double *com)
{
  ofstream outl(logfile.c_str(), ios::app);
  if (outl) {
    outl << setw(15) << time << setw(15) << Ecurr << setw(15) << used;
    for (int k=0; k<3; k++) outl << setw(15) << axis[k+1];
    for (int k=0; k<3; k++) outl << setw(15) << center[k+1] + com[k];
    outl << setw(15) << Egrad << setw(15) << dE;
    outl << endl;
  }
}

