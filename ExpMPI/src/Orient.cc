#include <stdlib.h>
#include <math.h>

#include <iostream>
#include <fstream>
#include <string>

#include <localmpi.h>

#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

#include <Orient.H>


const string axis_cache = "axis.cache";

Matrix return_euler_slater(double PHI, double THETA, double PSI, int BODY);

Orient::Orient(int n, int nwant, double Einit)
{
  keep = n;
  current = 0;
  many = nwant;
  Egrad = 0.0;
  Ecurr = Einit;
  Nlast = 0;

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
  axis.zero();
  center.zero();

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

    ifstream in(axis_cache.c_str());
    
    if (in) {

      in_ok = 1;		// Signal slave: OK

      MPI_Bcast(&in_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);

      in >> axis[1];
      in >> axis[2];
      in >> axis[3];
      in >> center[1];
      in >> center[2];
      in >> center[3];

      cout << "Orient: cached axis master: " 
	   << axis[1] << ", "
	   << axis[2] << ", "
	   << axis[3] << "\n";
      cout << "Orient: cached center master: " 
	   << center[1] << ", "
	   << center[2] << ", "
	   << center[3] << "\n";

      for (int j=0; j<3; j++) {
	in1[j] = axis[j+1];
	in2[j] = center[j+1];
      }
      MPI_Bcast(in1, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(in2, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    } else {

      in_ok = 0;		// Signal slave: NO VALUES

      MPI_Bcast(&in_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

  } else {

				// Get state from Master

    MPI_Bcast(&in_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (in_ok) {

      MPI_Bcast(in1, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(in2, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      for (int j=0; j<3; j++) {
	axis[j+1] = in1[j];
	center[j+1] = in2[j];
      }

      if (myid==1) {
	cerr << "Orient: cached axis slave: " 
	     << axis[1] << ", "
	     << axis[2] << ", "
	     << axis[3] << "\n";
	cerr << "Orient: cached center slave: " 
	     << center[1] << ", "
	     << center[2] << ", "
	     << center[3] << "\n";
      }

    }
  }
  delete [] in1;
  delete [] in2;

  if (in_ok) {

    double phi = atan2(axis[2], axis[1]);
    double theta = -acos(axis[3]/sqrt(axis*axis));
    double psi = 0.0;

    body = return_euler_slater(phi, theta, psi, 0);
    orig = return_euler_slater(phi, theta, psi, 1);
  }
      
}


void Orient::accumulate(vector<Particle> *p, double *com)
{
  double mass;

  angm.clear();

  int nbodies = p->size();
  for (int i=0; i<nbodies; i++) {
    if ((*p)[i].pot < Ecurr) {
      mass = (*p)[i].mass;

      t.E = (*p)[i].pot + 0.5*((*p)[i].vel[0]*(*p)[i].vel[0] + (*p)[i].vel[1]*(*p)[i].vel[1] + (*p)[i].vel[2]*(*p)[i].vel[2]);
      /*
	t.E = (*p)[i].pot;
      */
      t.M = mass;
      t.L[1] = mass*(((*p)[i].pos[1]-com[1])*(*p)[i].vel[2]-((*p)[i].pos[2]-com[2])*(*p)[i].vel[1]);
      t.L[2] = mass*(((*p)[i].pos[2]-com[2])*(*p)[i].vel[0]-((*p)[i].pos[0]-com[0])*(*p)[i].vel[2]);
      t.L[3] = mass*(((*p)[i].pos[0]-com[0])*(*p)[i].vel[1]-((*p)[i].pos[1]-com[1])*(*p)[i].vel[0]);
      /*
	t.R[1] = (*p)[i].pot*(*p)[i].pos[0];
	t.R[2] = (*p)[i].pot*(*p)[i].pos[1];
	t.R[3] = (*p)[i].pot*(*p)[i].pos[2];
      */
      t.R[1] = mass*(*p)[i].pos[0];
      t.R[2] = mass*(*p)[i].pos[1];
      t.R[3] = mass*(*p)[i].pos[2];
      angm.push_back(t);
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


				// Compute delta energy for next step
  if (used && Nlast) {
				// Don't divide by zero
    if (used != Nlast) 
      Egrad = (Ecurr - Elast)/(used - Nlast);

				// Don't believe a zero gradient
    if (Egrad != 0.0)
      dE = (double)(many - used) * Egrad;
  }

				// Push current value onto stack
  if (used) {

    Nlast = used;
    Elast = Ecurr;

    axis1   /= mtot;
    center1 /= mtot;
    sumsA.push_back(axis1);
    sumsC.push_back(center1);

    if (sumsA.size() > keep) {

      sumsA.pop_front();

      double x;
      int i=0;

      sumX = 0.0;
      sumX2 = 0.0;
      sumY.zero();
      sumXY.zero();
      sumY2.zero();

      deque<Vector>::iterator j;
      for (j = sumsA.begin(); j != sumsA.end(); j++) {
	x = sumsA.size() - 1 - (double)i;
	sumX += x;
	sumX2 += x*x;
	sumY += *j;
	sumXY += *j * x;
	sumY2 += *j & *j;

	i++;
      }

				// Linear least squares estimate for axis

      axis = (sumX2*sumY - sumX*sumXY)/(sumX2*i - sumX*sumX);
      slope = (sumXY*i - sumX*sumY)/(sumX2*i - sumX*sumX);

      i = 0;
      sigA = 0.0;
      for (j = sumsA.begin(); j != sumsA.end(); j++) {
	x = sumsA.size() - 1 - (double)i;
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

    /*
    cout << "Process " << myid << ": size=" << sumsC.size()
	 << " com=" 
	 << center[1] << " "
	 << center[2] << " "
	 << center[3] << "\n";
    */

    if (sumsC.size() > 2) {

      if (sumsC.size() > keep) sumsC.pop_front();

      double x;
      int i=0;
      sumX = 0.0;
      sumX2 = 0.0;
      sumY.zero();
      sumXY.zero();
      sumY2.zero();
      deque<Vector>::iterator j;
      for (j = sumsC.begin(); j != sumsC.end(); j++) {
	x = sumsC.size() - 1 - i;
	sumX += x;
	sumX2 += x*x;
	sumY += *j;
	sumXY += *j * x;
	sumY2 += *j & *j;

	i++;
      }
				// Linear least squares estimate for center

      center = (sumX2*sumY - sumX*sumXY)/(sumX2*i - sumX*sumX);	
      slope = (sumXY*i - sumX*sumY)/(sumX2*i - sumX*sumX);

      i = 0;
      sigC = 0.0;
      sigCz = 0.0;
      for (j = sumsC.begin(); j != sumsC.end(); j++) {
	x = sumsC.size() - 1 - (double)i;
	sigC += (*j - center - slope*x)*(*j - center - slope*x);
	sigCz += 
	  ( (*j)[3] - center[3] - slope[3]*x ) *
	  ( (*j)[3] - center[3] - slope[3]*x );
	i++;
      }
      sigC  /= i;
      sigCz /= i;

    }

    double factor = (double)(sumsC.size() - keep)/keep;
    factor = 1.0 - factor*factor;

    center *= factor*factor;

  }

				// Energy for next iteration
  if (fabs(dE)>fabs(0.1*Ecurr)) {
				// Suppress large step or gradient computation
				// failure
    if (used > many)
      Ecurr *= 1.01;
    else
      Ecurr *= 0.99;
  } 
  else				// Use gradient computed shift
    Ecurr += dE;
  
}
