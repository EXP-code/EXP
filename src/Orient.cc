#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>

#include "expand.H"

#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

#include <Orient.H>


void EL3::debug() const 
{
  if (myid==0) {
    cerr << left << setfill('-');
    ostringstream ostr;
    ostr << "--- EL3 [" << myid << "] ";
    cerr << setw(60) << ostr.str() << endl << setfill(' ');
    cerr.precision(4);
    cerr << setw(12) << T
	 << setw(12) << M
	 << setw(12) << E;
    cerr << " [";
    for (int i=1; i<=3; i++) cerr << setw(12) << L[i];
    cerr << "] [";
    for (int i=1; i<=3; i++) cerr << setw(12) << R[i];
    cerr << "]" << endl;
    cerr << left << setfill('-') 
	 << setw(60) << '-' << endl << setfill(' ');
  }
}

Eigen::MatrixXd
return_euler_slater(double PHI, double THETA, double PSI, int BODY);

Orient::Orient(int n, int nwant, unsigned Oflg, unsigned Cflg,
	       string Logfile, double dt, double damping)
{
  keep    = n;
  current = 0;
  many    = nwant;
  oflags  = Oflg;
  cflags  = Cflg;
  logfile = Logfile;
  deltaT  = dt;
  Nlast   = 0;
  damp    = damping;
  linear  = false;

  pos = vector<double>(3);
  psa = vector<double>(3);
  vel = vector<double>(3);

  center .setZero();
  center0.setZero();
  cenvel0.setZero();
  axis   .setZero();

				// Initialize last time to something
				// in the near infinite past
  lasttime = -std::numeric_limits<double>::max();

  axis[2] = 1;			// This sets the axis to the z-axis

  used = 0;			// No particles used to start

				// Set up identity
  body.setIdentity();
  orig = body;

				// Check for previous state on
				// a restart
  int in_ok;
  std::vector<double> in1(4), in2(4);

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
      int tused;

      in_ok = 1;		// Signal slave: OK
	
      MPI_Bcast(&in_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);

      const int cbufsiz = 16384;
      char *cbuffer = new char [cbufsiz];
	
				// Look for data and write it while
				// accumlating data for averaging
      while (in && restart) {

	in.getline(cbuffer, cbufsiz);
	if (in.rdstate() & (ios::failbit | ios::eofbit)) break;

				// Skip comment lines
				//
	if (cbuffer[0] == '#') continue;

	istringstream line(cbuffer);

				// Read until current time is reached
	line >> time;		// 
	if (tnow+0.1*dtime/Mstep < time) break;

	out << cbuffer << "\n";

	line >> Ecurr;
	line >> tused;
	line >> axis[0];	// Last computed axis from regression
	line >> axis[1];
	line >> axis[2];
	line >> axis1[0];	// Last axis from particle algorithm
	line >> axis1[1];
	line >> axis1[2];
	line >> center[0];	// Last computed center from regression
	line >> center[1];
	line >> center[2];
	line >> center0[0];	// Analytic center (from velocity intgration)
	line >> center0[1];
	line >> center0[2];
	line >> center1[0];	// Last center from particle algorithm
	line >> center1[1];
	line >> center1[2];
	  
	if (oflags & AXIS) {
	  sumsA.push_back(DV(time, axis1));
	  if (static_cast<int>(sumsA.size()) > keep) sumsA.pop_front();
	}

	if (oflags & CENTER) {
	  sumsC.push_back(DV(time, center1));
	  if (static_cast<int>(sumsC.size()) > keep) sumsC.pop_front();
	}
	
      }

      cout << "  -- Orient: current log=" << logfile << "  backup=" << backupfile << endl;

      cout << "  -- Orient: cached time=" << time << "  Ecurr= " << Ecurr << endl;

      cout << "  -- Orient: axis master (cache size=" << sumsA.size() << "): " 
	   << axis[0] << ", "
	   << axis[1] << ", "
	   << axis[2] << endl;
      
      cout << "  -- Orient: center master (cache size=" << sumsC.size() << "): " 
	   << center[0] << ", "
	   << center[1] << ", "
	   << center[2] << endl;
      
      MPI_Bcast(&Ecurr,         1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(axis.data(),    3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(center.data(),  3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(center0.data(), 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      
      int howmany = max<int>(sumsA.size(), sumsC.size());
      MPI_Bcast(&howmany, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
      for (int k=0; k<howmany; k++) {

	if (oflags & AXIS) {
	  in1[0] = sumsA[k].first;
	  for (int j=1; j<=3; j++) in1[j] = sumsA[k].second[j-1];
	  MPI_Bcast(in1.data(), 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}

	if (oflags & CENTER) {
	  in2[0] = sumsC[k].first;
	  for (int j=1; j<=3; j++) in2[j] = sumsC[k].second[j-1];
	  MPI_Bcast(in2.data(), 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
      }


    } else {
	
      in_ok = 0;		// Signal slave: NO VALUES
	
      MPI_Bcast(&in_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);

				// Write header
      ofstream out(logfile.c_str());
      if (out) {
	out.setf(ios::left);
	out << setw(15) << "# Time"		// 1
	    << setw(15) << "| E_curr"		// 2
	    << setw(15) << "| Used"		// 3
	    << setw(15) << "| X-axis(reg)"	// 4
	    << setw(15) << "| Y-axis(reg)"	// 5
	    << setw(15) << "| Z-axis(reg)"	// 6
	    << setw(15) << "| X-axis(cur)"	// 7
	    << setw(15) << "| Y-axis(cur)"	// 8
	    << setw(15) << "| Z-axis(cur)"	// 9
	    << setw(15) << "| X-center(anl)"	// 10
	    << setw(15) << "| Y-center(anl)"	// 11
	    << setw(15) << "| Z-center(anl)"	// 12
	    << setw(15) << "| X-center(reg)"	// 13
	    << setw(15) << "| Y-center(reg)"	// 14
	    << setw(15) << "| Z-center(reg)"	// 15
	    << setw(15) << "| X-center(cur)"	// 16
	    << setw(15) << "| Y-center(cur)"	// 17
	    << setw(15) << "| Z-center(cur)"	// 18
	    << setw(15) << "| X-com(cur)"	// 19
	    << setw(15) << "| Y-com(cur)"	// 20
	    << setw(15) << "| Z-com(cur)"	// 21
	    << setw(15) << "| X-com(dif)"	// 22
	    << setw(15) << "| Y-com(dif)"	// 23
	    << setw(15) << "| Z-com(dif)"	// 24
	    << endl;
	out.fill('-');

	int icnt = 1;
	out << "# " << setw(13) << icnt++;
	for (int i=0; i<23; i++) out << "| " << setw(13) << icnt++;
	out << endl;

	out.close();

      } else {
	cerr << "Orient: error opening log file <" << logfile << ">\n";
      }
    }
    
  } else {

				// Get state from Master

    MPI_Bcast(&in_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (in_ok) {

      MPI_Bcast(&Ecurr,         1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(axis.data(),    3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(center.data(),  3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(center0.data(), 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      int howmany;
      MPI_Bcast(&howmany, 1, MPI_INT, 0, MPI_COMM_WORLD);

      for (int k=0; k<howmany; k++) {

	if (oflags & AXIS) {
	  MPI_Bcast(in1.data(), 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	  for (int j=0; j<3; j++) axis1[j] = in1[j+1];
	  sumsA.push_back(DV(in1[0], axis1));
	}

	if (oflags & CENTER) {
	  MPI_Bcast(in2.data(), 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	  for (int j=0; j<3; j++) center1[j] = in2[j+1];
	  sumsC.push_back(DV(in2[0], center1));
	}

      }

    }

  }

  if (in_ok) {

    if (oflags & AXIS) {

      double phi = atan2(axis[0], axis[0]);
      double theta = -acos(axis[2]/sqrt(axis.dot(axis)));
      double psi = 0.0;

      body = return_euler_slater(phi, theta, psi, 0);
      orig = return_euler_slater(phi, theta, psi, 1);
    }

  }
      
}

void Orient::accumulate_cpu(double time, Component *c)
{
  double energy, mass, v2;
  unsigned nbodies = c->Number();
  PartMapItr it = c->Particles().begin();
  unsigned tkeep = many/numprocs;
  set<EL3, ltEL3>::reverse_iterator el3last = angm.rbegin();

  for (unsigned q=0; q<nbodies; q++) {

    unsigned long i = (it++)->first;
    Particle     *p = c->Part(i);

    v2 = 0.0;
    for (int k=0; k<3; k++) {
      pos[k] = c->Pos(i, k, Component::Local);
      if (std::isnan(pos[k])) {
	cerr << "Orient: process " << myid << " index=" << i
	     << " has NaN on component ";
	for (int s=0; s<3; s++)
	  cerr << setw(16) << p->pos[s];
	for (int s=0; s<3; s++)
	  cerr << setw(16) << p->vel[s];
	for (int s=0; s<3; s++)
	  cerr << setw(16) << p->acc[s];
	cerr << endl;
      }
      vel[k] = c->Vel(i, k, Component::Local);
      psa[k] = pos[k] - center[k];
      v2 += vel[k]*vel[k];
    }

    energy = p->pot;
    
    if (cflags & KE) energy += 0.5*v2;

    if (cflags & EXTERNAL) energy += p->potext;

    unsigned size0 = angm.size();
    bool test1 = (size0 <= tkeep);
    bool test2 = true;

    if (size0) test2 = (energy < el3last->E);
    
    if (test1 || test2) {

      mass = p->mass;

      t.E = energy;
      t.T = time;
      t.M = mass;

      t.L[0] = mass*(psa[1]*vel[2] - psa[2]*vel[1]);
      t.L[1] = mass*(psa[2]*vel[0] - psa[0]*vel[2]);
      t.L[2] = mass*(psa[0]*vel[1] - psa[1]*vel[0]);

      t.R[0] = mass*pos[0];
      t.R[1] = mass*pos[1];
      t.R[2] = mass*pos[2];

      // Trim the list by removing the element with the largest energy
      //
      if (test2 && !test1) angm.erase(*el3last);

      // Insert the new element
      //
      angm.insert(t);

      // Reset the iterator at the top of the heap
      //
      el3last = angm.rbegin();

#ifdef DEBUG      
      t.debug();
#endif
    }
  }
}


void Orient::accumulate(double time, Component *c)
{
				// Do not register a duplicate entry
  if (fabs(lasttime - time) < 1.0e-12) return;
				// Space entries by at least deltaT
  if (time - deltaT - lasttime < 0.0 ) return;

  lasttime = time;

  if (linear) {
      center = center0;
      center0 += cenvel0*dtime;
      return;
  }

  comp->timer_orient.start();	// Time the accumulate step
  
  angm.clear();

#if HAVE_LIBCUDA
  if (use_cuda)
    accumulate_gpu(time, c);
  else
#endif
    accumulate_cpu(time, c);

  comp->timer_orient.stop();

  unsigned tkeep = many/numprocs;

  std::vector<double> ee;
  for (auto it = angm.begin(); it != angm.end(); it++) ee.push_back(it->E);
  
  for (int n=1; n<numprocs; n++) {
    if (n==myid) {
      unsigned nsiz = ee.size();
      // Size trim approximation; 3x the target size
      nsiz = std::min<unsigned>(nsiz, 3.0*tkeep);
      // Send to root node for sorting
      MPI_Send(&nsiz,     1, MPI_UNSIGNED, 0, 331, MPI_COMM_WORLD);
      MPI_Send(&ee[0], nsiz, MPI_DOUBLE,   0, 332, MPI_COMM_WORLD);
    }
    if (0==myid) {
      unsigned nsiz;
      MPI_Recv(&nsiz, 1, MPI_UNSIGNED, n, 331, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      std::vector<double> tt(nsiz);
      MPI_Recv(&tt[0], nsiz, MPI_DOUBLE, n, 332, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      ee.insert(ee.end(), tt.begin(), tt.end());
    }
  }

  if (myid==0) {
    std::sort(ee.begin(), ee.end());
    if (ee.size()<=many) Ecurr = ee.back();
    else                 Ecurr = *(ee.begin()+many);
  }

  // Propagate minimum energy and current cached low energy particles
  // within nodes

  MPI_Bcast(&Ecurr, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

				// Compute values for this step
  axis1  .setZero();
  center1.setZero();

  double mtot=0.0, mtot1=0.0;
  int cnum = 0;
  for (set<EL3, ltEL3>::iterator i=angm.begin(); i!=angm.end() && i->E<Ecurr; i++) {
    axis1   += i->L;
    center1 += i->R;
    mtot1   += i->M;
    cnum++;
  }

  if (0) {
    for (int n=0; n<numprocs; n++) {
      if (n==myid) {
	if (myid==0) cout << "------------------------" << endl
			  << "Center check in Orient  " << endl
			  << setw(15) << left << "Component: "
			  << c->name << endl
			  << setw(15) << left << "Time: "
			  << time << endl
			  << "------------------------" << endl 
			  << right;
	cout << setw(4) << myid << setw(4);

	for (int k=1; k<=3; k++) {
	  if (mtot1>0.0)
	    cout << setw(18) << center1[k]/mtot1;
	  else
	    cout << setw(18) << 0.0;
	}

	for (int k=1; k<=3; k++) 
	  cout << setw(18) << angm.begin()->R[k]/angm.begin()->M;

	cout << setw(18) << angm.begin()->E << setw(18) << Ecurr
	     << setw(18) << angm.end()  ->E << setw(18) << endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==0) cout << endl;
  }

  Eigen::Vector3d inA = axis1;
  Eigen::Vector3d inC = center1;

  axis1.setZero();
  center1.setZero();

				// Share stuff between nodes
  MPI_Allreduce(&cnum, &used, 1, 
		MPI_INT,    MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(inA.data(), axis1.data(), 3, 
		MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce(&mtot1, &mtot, 1, 
		MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(inC.data(), center1.data(), 3, 
		MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  // Push current value onto stack

  Nlast = used;
  Elast = Ecurr;

  if (mtot>0.0) {
    axis1   /= mtot;
    center1 /= mtot;
    if (oflags & AXIS)   sumsA.push_back(DV(time, axis1  ));
    if (oflags & CENTER) sumsC.push_back(DV(time, center1));
  }

  if ((cflags & DIAG) && myid==0) {
    cout << " Orient info [" << time << ", " << c->name << "]: " 
	 << used << " particles used, Ecurr=" << Ecurr 
	 << " Center=" 
	 << center1[0] << ", "
	 << center1[1] << ", "
	 << center1[2] << endl;
  }

  if (static_cast<int>(sumsA.size()) > keep + 1) {

    sumsA.pop_front();

    double x;
    int i=0;
    
    sumX = 0.0;
    sumX2 = 0.0;
    sumY.setZero();
    sumXY.setZero();
    sumY2.setZero();
    
    int N = sumsC.size();
    
    deque<DV>::iterator j;
    for (j = sumsA.begin(); j != sumsA.end(); j++) {
      x      = j->first;
      sumX  += x;
      sumX2 += x * x;
      sumY  += j->second;
      sumXY += j->second * x;
      for (int k=0; k<3; k++) 
	sumY2[k] += j->second[k] * j->second[k];
      
      i++;
    }
    
    // Linear least squares estimate for axis

    slope     = (sumXY*N - sumX*sumY)/(sumX2*N - sumX*sumX);
    intercept = (sumX2*sumY - sumX*sumXY)/(sumX2*N - sumX*sumX);
    axis      = intercept + slope*(damp*time + (1.0 - damp)*sumsA.front().first);
    
    i = 0;
    sigA = 0.0;
    for (j = sumsA.begin(); j != sumsA.end(); j++) {
      sigA += 
	(j->second - intercept - slope*j->first).adjoint() *
	(j->second - intercept - slope*j->first) ;
      i++;
    }
    sigA /= i;
    
    double phi = atan2(axis[2], axis[0]);
    double theta = -acos(axis[2]/sqrt(axis.dot(axis)));
    double psi = 0.0;
    
    body = return_euler_slater(phi, theta, psi, 0);
    orig = return_euler_slater(phi, theta, psi, 1);
  }


  if (sumsC.size() > 1) {
    
    if (static_cast<int>(sumsC.size()) > keep+1) sumsC.pop_front();

    double x;
    int i = 0;
    sumX  = 0.0;
    sumX2 = 0.0;
    sumY. setZero();
    sumXY.setZero();
    sumY2.setZero();
    
    int N = sumsC.size();
      
    deque<DV>::iterator j;
    for (j = sumsC.begin(); j != sumsC.end(); j++) {
      x      = j->first;
      sumX  += x;
      sumX2 += x * x;
      sumY  += j->second;
      sumXY += j->second * x;
      for (int k=0; k<3; k++)
	sumY2[k] += j->second[k] * j->second[k];

      i++;

      if ((cflags & DIAG) && myid==0)
	cout << " Orient debug [" << time << ", " << c->name << "] i=" 
	     << i << ":"
	     << "      t=" << setw(15) << x
	     << "      x=" << setw(15) << j->second[0]
	     << "      y=" << setw(15) << j->second[1]
	     << "      z=" << setw(15) << j->second[2]
	     << "   SumX=" << setw(15) << sumX 
	     << "  SumX2=" << setw(15) << sumX2 
	     << "  Delta=" << setw(15) << sumX2*i - sumX*sumX 
	     << endl;
    }
    // Linear least squares estimate for center
    
    slope     = (sumXY*N - sumX*sumY)/(sumX2*N - sumX*sumX);
    intercept = (sumX2*sumY - sumX*sumXY)/(sumX2*N - sumX*sumX);
    center    = intercept + slope*(damp*time + (1.0 - damp)*sumsC.front().first);
    
    i = 0;
    sigC = 0.0;
    sigCz = 0.0;
    for (j = sumsC.begin(); j != sumsC.end(); j++) {
      sigC += 
	(j->second - intercept - slope*j->first).adjoint() *
	(j->second - intercept - slope*j->first) ;
      sigCz += 
	( j->second[2] - intercept[2] - slope[2]*j->first ) *
	( j->second[2] - intercept[2] - slope[2]*j->first ) ;
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
	 << " Orient info [" << time << ", " << c->name << "]:"   << endl
	 << " size=" << sumsC.size() << " sigC=" << sigC  
	 << " sigCz=" << sigCz << endl
	 << "  SumX=" << sumX << " SumX2=" << sumX2 << endl
	 << "  SumY="
	 << sumY[0] << " "
	 << sumY[1] << " "
	 << sumY[2] << endl
	 << "  SumXY="
	 << sumXY[0] << " "
	 << sumXY[1] << " "
	 << sumXY[2] << endl
	 << "  SumY2="
	 << sumY2[0] << " "
	 << sumY2[1] << " "
	 << sumY2[2] << endl
	 << "  slope="
	 << slope[0] << " "
	 << slope[1] << " "
	 << slope[2] << endl
	 << "  center=" 
	 << center[0] << " "
	 << center[1] << " "
	 << center[2] << endl
	 << "===================================================" << endl;
  }
  
  // Will force center to be last minimum energy average
  // rather than pooled average of the last "keep" states
  if (keep == 0) {
    for (int i=1; i<=3; i++) center[i] = sumsC.begin()->second[i];
    sumsC.pop_front();
  }
    
}

void Orient::logEntry(double time, Component *c)
{
  if (myid) return;

  ofstream outl(logfile.c_str(), ios::app);
  if (outl) {
    // Columns 1 - 3
    outl << setw(15) << time << setw(15) << Ecurr << setw(15) << used;

    // Columns 4 - 6
    for (int k=0; k<3; k++) outl << setw(15) << axis[k];

    // Columns 7 - 9
    for (int k=0; k<3; k++) outl << setw(15) << axis1[k];

    // Columns 10 - 12
    for (int k=0; k<3; k++) outl << setw(15) << center[k];

    // Columns 13 - 15
    for (int k=0; k<3; k++) outl << setw(15) << center0[k];

    // Columns 16 - 18
    for (int k=0; k<3; k++) outl << setw(15) << center1[k];

    // Columns 19 - 21
    for (int k=0; k<3; k++) outl << setw(15) << c->com[k];

    // Columns 22 - 24
    for (int k=0; k<3; k++) outl << setw(15) << c->com0[k];

    outl << endl;
  }
}


Orient::~Orient()
{
}
