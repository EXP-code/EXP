#include "expand.h"
#include "localmpi.h"

void set_global_com(void);

/*
  Tags:

   1 = particle count
   2 = buffer count
   3 = particle buffer
  52 = buffer count on gather
  53 = particle gather buffer

*/


#include <stdlib.h>
#include <unistd.h>
#include <math.h>

				// STL stuff
#include <iostream>
#include <iomanip>
#include <fstream>
#include <strstream>
#include <string>
#include <vector>
#include <algorithm>

#include <Vector.h>

void parallel_gather_coefficients(Matrix& expcoef, Matrix& expcoef1,
				  Matrix*& cc, Matrix*& cc1,
				  int lmax)
{
  int Ldim, L0, loffset, moffset, l, m, n, nn;

#ifdef MPE_PROFILE
  MPE_Log_event(5, myid, "b_gather_c");
#endif

  if (dof==3) {
    Ldim = lmax*(lmax + 2) + 1;
    L0 = 0;
  }
  else {
    Ldim = 2*lmax + 1;
    L0 = lmax;
  }    

  if (myid == 0) {

    for (l=L0, loffset=0; l<=lmax; loffset+=(2*l+1), l++) {

      for (m=0, moffset=0; m<=l; m++) {

	if (m==0) {
	  for (n=1; n<=nmax; ++n) {
	    expcoef[loffset+moffset][n] = 0.0;

	    for (nn=n; nn<=nmax; nn++)
	      cc[loffset+moffset][n][nn] = 0.0;
	  }
	  moffset++;
	}
	else {
	  for (n=1; n<=nmax; ++n) {
	    expcoef[loffset+moffset][n] = 0.0;
	    expcoef[loffset+moffset+1][n] = 0.0;

	    for (nn=n; nn<=nmax; nn++) {
	      cc[loffset+moffset][n][nn] = 0.0;
	      cc[loffset+moffset+1][n][nn] = 0.0;
	    }
	  }
	  moffset+=2;
	}
      }
    }
  }


  for (l=L0, loffset=0; l<=lmax; loffset+=(2*l+1), l++) {

    for (m=0, moffset=0; m<=l; m++) {

      if (m==0) {
	MPI_Reduce(&expcoef1[loffset+moffset][1], 
		   &expcoef[loffset+moffset][1], nmax, 
		   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	for (n=1; n<=nmax; n++)
	  MPI_Reduce(&cc1[loffset+moffset][n][n],
		     &cc[loffset+moffset][n][n], nmax-n+1, 
		     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	moffset++;
      }
      else {
	MPI_Reduce(&expcoef1[loffset+moffset][1], 
		   &expcoef[loffset+moffset][1], nmax, 
		   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	MPI_Reduce(&expcoef1[loffset+moffset+1][1],
		   &expcoef[loffset+moffset+1][1], nmax, 
		   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	
	for (n=1; n<=nmax; n++) {
	  MPI_Reduce(&cc1[loffset+moffset][n][n],
		     &cc[loffset+moffset][n][n], nmax-n+1, 
		     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&cc1[loffset+moffset+1][n][n],
		     &cc[loffset+moffset+1][n][n], nmax-n+1, 
		     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}

	moffset+=2;
      }
    }
  }

#ifdef MPE_PROFILE
  MPE_Log_event(6, myid, "e_gather_c");
#endif

}

void parallel_distribute_coefficients(Matrix& expcoef, int lmax)
{
  int Ldim, L0, loffset, moffset, l, m;

#ifdef MPE_PROFILE
  MPE_Log_event(7, myid, "b_distrib_c");
#endif

  if (dof==3) {
    Ldim = lmax*(lmax + 2) + 1;
    L0 = 0;
  }
  else {
    Ldim = 2*lmax + 1;
    L0 = lmax;
  }    

  for (l=L0, loffset=0; l<=lmax; loffset+=(2*l+1), l++) {

      for (m=0, moffset=0; m<=l; m++) {

	if (m==0) {
	  MPI_Bcast(&expcoef[loffset+moffset][1], nmax, MPI_DOUBLE,
		    0, MPI_COMM_WORLD);
	  moffset++;
	}
	else {
	  MPI_Bcast(&expcoef[loffset+moffset][1], nmax, MPI_DOUBLE,
		     0, MPI_COMM_WORLD);
	  MPI_Bcast(&expcoef[loffset+moffset+1][1], nmax, MPI_DOUBLE,
		    0, MPI_COMM_WORLD);
	  moffset+=2;
	}
      }
  }

#ifdef MPE_PROFILE
  MPE_Log_event(8, myid, "e_distrib_c");
#endif

}

static vector<int> plist, plist2;
static vector<int> ncount;

/*
  The debug Partstruct carries a sequence number for strong checking
  of the redistribution algorithm and implementation
*/
  
#ifdef SEQCHECK
struct Partstruct
{				// Offsets:
  int component;		// 0
  unsigned int sequence;	// 1
  double mass;			// 2
  double pos[3];		// 3,4,5
  double vel[3];		// 6,7,8
  double pot;			// 9
  double potext;		// 10
  double esave;			// 11
  double mfp;			// 12
};
#else
struct Partstruct
{				// Offsets:
  int component;		// 0
  double mass;			// 1
  double pos[3];		// 2,3,4
  double vel[3];		// 5,6,7
  double pot;			// 8
  double potext;		// 9
  double esave;			// 10
  double mfp;			// 11
};
#endif

MPI_Datatype Particletype;

static const int nbuf = 2000;
static int seq_beg, seq_end;

extern "C" void read_bodies_and_distribute(void)
{
				// MPI buffers
  MPI_Status status0;
#ifdef SEQCHECK
  struct Partstruct buf[nbuf];

				// Make MPI datatype
  MPI_Datatype	type[9] = {MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
			   MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};

				// Get displacements
  MPI_Aint	disp[9];
  MPI_Get_address(&buf[0].component,	&disp[0]);
  MPI_Get_address(&buf[0].sequence,	&disp[1]);
  MPI_Get_address(&buf[0].mass,		&disp[2]);
  MPI_Get_address(&buf[0].pos,		&disp[3]);
  MPI_Get_address(&buf[0].vel,		&disp[4]);
  MPI_Get_address(&buf[0].pot,		&disp[5]);
  MPI_Get_address(&buf[0].potext,	&disp[6]);
  MPI_Get_address(&buf[0].esave,	&disp[7]);
  MPI_Get_address(&buf[0].mfp,		&disp[8]);

  for (int i=8; i>=0; i--) disp[i] -= disp[0];
  
				// Block offsets
  int		blocklen[9] = {1, 1, 1, 3, 3, 1, 1, 1, 1};
  
  MPI_Type_create_struct(9, blocklen, disp, type, &Particletype);
  MPI_Type_commit(&Particletype);

#else
  struct Partstruct buf[nbuf];

				// Make MPI datatype
  MPI_Datatype	type[8] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
			   MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};

				// Get displacements
  MPI_Aint	disp[8];
  MPI_Get_address(&buf[0].component,	&disp[0]);
  MPI_Get_address(&buf[0].mass,		&disp[1]);
  MPI_Get_address(&buf[0].pos,		&disp[2]);
  MPI_Get_address(&buf[0].vel,		&disp[3]);
  MPI_Get_address(&buf[0].pot,		&disp[4]);
  MPI_Get_address(&buf[0].potext,	&disp[5]);
  MPI_Get_address(&buf[0].esave,	&disp[6]);
  MPI_Get_address(&buf[0].mfp,		&disp[7]);

  for (int i=7; i>=0; i--) disp[i] -= disp[0];
  
				// Block offsets
  int		blocklen[8] = {1, 1, 3, 3, 1, 1, 1, 1};
  
  MPI_Type_create_struct(8, blocklen, disp, type, &Particletype);
  MPI_Type_commit(&Particletype);
#endif

  double tmp, rmax1, r2;

  ifstream *fin;
  const int nline = 256;
  char line[nline];
  
  if (myid == 0) {
				// Read in phase space

    fin = new ifstream(infile);

    if (!*fin) {
      cerr << "Couldn't open " << infile << " . . . quitting\n";
      exit(-1);
    }

    fin->getline(line, nline);
    istrstream ins(line);
    
    ins >> nbodies_tot;
    ins >> tnow;
  }

  MPI_Bcast(&nbodies_tot, 1, MPI_INT, 0, MPI_COMM_WORLD);


  if (nbodies_tot > nbodmax*numprocs) {
    if (myid==0) {
      cerr << "Not enough space on all processors to hold phase space\n";
      cerr << "nbodmax should be at least "
	   << (int)( (double)nbodies_tot/numprocs + 1) << endl;
    }
    MPI_Finalize();
    exit(-1);
  }

  is_init = 1;
  setup_distribution();
  is_init = 0;

				// Form cumulative and differential bodies list
  
  plist = vector<int>(numprocs);
  for (int i=0; i<numprocs; i++) plist[i] = nbodies_index[i];
  ncount = plist;
  for (int i=1; i<numprocs; i++) ncount[i] -= plist[i-1];

  if (myid==0) seq_beg = 1;
  else seq_beg = plist[myid-1]+1;
  seq_end = plist[myid];

  if (myid==0) {

				// Read root node particles

				// First line
    {
      fin->getline(line, nline);
      istrstream ins(line);

      ins >> component[1];
#ifdef SEQCHECK
      sequence[1] = 1;
#endif
      ins >> mass[1];
      ins >> x[1];
      ins >> y[1];
      ins >> z[1];
      ins >> vx[1];
      ins >> vy[1];
      ins >> vz[1];
      ins >> tmp;
      if (ins.rdstate() & ios::eofbit) restart = 0;
      else restart = 1;
      pot[1] = potext[1] = 0.0;

      rmax1 = x[1]*x[1] + y[1]*y[1] + z[1]*z[1];
    }

    if (restart)
      cout << "Restart . . .\n";
    else
      cout << "New run . . .\n";
    cout << flush;

				// Remainder of Node 0's particles

    for (int i=2; i<=ncount[0]; i++) {
      fin->getline(line, nline);
      istrstream ins(line);
      ins >> component[i];
#ifdef SEQCHECK
      sequence[i] = i;
#endif
      ins >> mass[i];
      ins >> x[i];
      ins >> y[i];
      ins >> z[i];
      ins >> vx[i];
      ins >> vy[i];
      ins >> vz[i];

      pot[i] = potext[i] = 0.0;

      r2 = x[i]*x[i] + y[i]*y[i] + z[i]*z[i];
      rmax1 = max<double>(r2, rmax1);
    }

    nbodies = ncount[0];

    int icount, ibufcount;
    for (int n=1; n<numprocs; n++) {

      MPI_Send(&ncount[n], 1, MPI_INT, n, 1, MPI_COMM_WORLD);

      icount = 0;
      ibufcount = 0;
      while (icount < ncount[n]) {

	fin->getline(line, nline);
	istrstream ins(line);
	ins >> buf[ibufcount].component;
#ifdef SEQCHECK
	buf[ibufcount].sequence = icount + nbodies_index[n-1] + 1;
#endif
	ins >> buf[ibufcount].mass;
	for (int k=0; k<3; k++) ins >> buf[ibufcount].pos[k];
	for (int k=0; k<3; k++) ins >> buf[ibufcount].vel[k];

	r2 = 0.0;
	for (int k=0; k<3; k++) 
	  r2 += buf[ibufcount].pos[k]*buf[ibufcount].pos[k];
	rmax1 = max<double>(r2, rmax1);

	ibufcount++;
	icount++;

	if (ibufcount == nbuf || icount == ncount[n]) {
	  MPI_Send(&ibufcount, 1, MPI_INT, n, 2, MPI_COMM_WORLD);
	  MPI_Send(buf, ibufcount, Particletype, n, 3, MPI_COMM_WORLD);
	  ibufcount = 0;
	}

      }

    }

  } else {

      MPI_Recv(&nbodies, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status0);
      
      int icount = 0, isent;
      while (icount < nbodies) {
	MPI_Recv(&isent, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &status0);
	MPI_Recv(buf, isent, Particletype, 0, 3, MPI_COMM_WORLD, &status0);
	for (int i=0; i<isent; i++) {
	  component[icount+1] = buf[i].component;
#ifdef SEQCHECK
	  sequence[icount+1] = buf[i].sequence;
#endif
	  mass[icount+1] = buf[i].mass;
	  x[icount+1] = buf[i].pos[0];
	  y[icount+1] = buf[i].pos[1];
	  z[icount+1] = buf[i].pos[2];
	  vx[icount+1] = buf[i].vel[0];
	  vy[icount+1] = buf[i].vel[1];
	  vz[icount+1] = buf[i].vel[2];
	  pot[icount+1] = potext[icount+1] = 0.0;
	  
	  icount++;
	}
      }
  }

#ifdef SEQCHECK			// Sanity check
  if (seq_beg != sequence[1] || seq_end != sequence[nbodies]) {
    cout << "Process " << myid << ": sequence error on init,"
	 << " seq_beg=" << seq_beg
	 << " seq_end=" << seq_end
	 << " seq[1]=" << sequence[1]
	 << " seq[N]=" << sequence[nbodies]
	 << endl << flush;
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
#endif
				// Default: set to max radius

  if (rmax <= 0.0) rmax = sqrt(rmax1);

  MPI_Bcast(&restart, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&rmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&tnow, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  

				// Initialize time
  tpos = tvel = tnow;
  
				// Set center of mass and com velocity
				// to zero for each component

  if (!restart && (zerocom || zerovel)) set_global_com();

				// Look for point mass particles
  make_pointmass_list();

  if (myid==0) delete fin;

}

extern "C" Partstruct* get_particles(int* number)
{
  MPI_Status status;

  static bool firsttime = true;
  static int counter = 0;
  static int node = 0;

  static struct Partstruct *buf;
  if (firsttime) {
    buf = new Partstruct [nbuf];
    firsttime = false;
  }

  if (*number < 0) {
    counter = 0;
    node = 0;
  }

  if (node==numprocs) {
    *number = 0;
    return 0;
  }

  int icount, indx, curnode = node;
  if (node == 0) {
    
    if (myid == 0) {
      
      icount = 0;
      while (icount < nbuf && counter < nbodies_index[node]) {

				// Pack structure

	buf[icount].component = component[counter+1];
#ifdef SEQCHECK
	buf[icount].sequence = sequence[counter+1];
#endif
	buf[icount].mass = mass[counter+1];
	buf[icount].pos[0] = x[counter+1];
	buf[icount].pos[1] = y[counter+1];
	buf[icount].pos[2] = z[counter+1];
	buf[icount].vel[0] = vx[counter+1];
	buf[icount].vel[1] = vy[counter+1];
	buf[icount].vel[2] = vz[counter+1];
	buf[icount].pot = pot[counter+1];
	buf[icount].potext = potext[counter+1];
	if (relx) buf[icount].esave = esave[counter+1];
	if (scatter) buf[icount].mfp = mfp[counter+1];

	icount++;
	counter++;
      }

      if (counter == nbodies_index[node]) node++;

#ifdef DEBUG
    cout << "Process " << myid 
	 << ": buffered " << icount << " particles from master"
	 << endl << flush;
#endif    

    }

  } else {

    if (myid == 0) {

      MPI_Recv(&icount, 1, MPI_INT, node, 52, MPI_COMM_WORLD, &status);
      MPI_Recv(buf, icount, Particletype, node, 53, MPI_COMM_WORLD, &status);

#ifdef DEBUG
    cout << "Process " << myid 
	 << ": received " << icount << " particles from slave"
	 << endl << flush;
#endif    

    } else if (myid == node) {
      
      icount = 0;
      while (icount < nbuf && counter < nbodies_index[node]) {

	indx = counter - nbodies_index[node-1] + 1;
	
				// Pack structure
	buf[icount].component = component[indx];
#ifdef SEQCHECK
	buf[icount].sequence = sequence[indx];
#endif
	buf[icount].mass = mass[indx];
	buf[icount].pos[0] = x[indx];
	buf[icount].pos[1] = y[indx];
	buf[icount].pos[2] = z[indx];
	buf[icount].vel[0] = vx[indx];
	buf[icount].vel[1] = vy[indx];
	buf[icount].vel[2] = vz[indx];
	buf[icount].pot = pot[indx];
	buf[icount].potext = potext[indx];
	if (relx) buf[icount].esave = esave[indx];
	if (scatter) buf[icount].mfp = mfp[indx];

	icount++;
	counter++;
      }

      MPI_Send(&icount, 1, MPI_INT, 0, 52, MPI_COMM_WORLD);
      MPI_Send(buf, icount, Particletype, 0, 53, MPI_COMM_WORLD);
#ifdef DEBUG
    cout << "Process " << myid 
	 << ": sent " << icount << " particles to master"
	 << ", counter value=" << counter
	 << ", nbodies_index=" << nbodies_index[node]
	 << endl << flush;
#endif    

      if (counter == nbodies_index[node]) node++;
      
    }

  }
  
  MPI_Bcast(&counter, 1, MPI_INT, curnode, MPI_COMM_WORLD);
  MPI_Bcast(&node, 1, MPI_INT, curnode, MPI_COMM_WORLD);
  MPI_Bcast(&icount, 1, MPI_INT, curnode, MPI_COMM_WORLD);

				// Return values
  *number = icount;
  return buf;

#ifdef DEBUG
  cout << "Process " << myid 
       << ": received counter=" << counter << " node=" << node
       << " number=" << *number
       << " address=" << buf
       << endl << flush;
#endif    

}


static double tdump = -1.0e20;

#include "tipsydefs.h"

extern "C" void create_tipsy(void)
{
  static bool firsttime = true;
  const int ncomps = 2;
  static ofstream **ocmp;
  if (firsttime) {
    ocmp = new ofstream* [ncomps+1];
    firsttime = false;
  }

  int number=-1;
  Partstruct *p;
  
  int ngas = 0;
  int ndark = 0;
  int nstar = 0;

  gas_particle gas;
  dark_particle dark;
  star_particle star;

				// Check for temp files for this step . . .
  if (fabs(tdump-tnow)<1.0e-10) return;

  bool ok = true;
  if (myid==0) {
    
    char fnames[3][50] = { "particle.header\0", 
			 "dark.particles\0", 
			 "star.particles\0" };

    for (int i=0; i<3; i++) {
      ocmp[i] = new ofstream(fnames[i], ios::out | ios::trunc);
      if (!*ocmp[i]) {
	cerr << "Error opening <" << fnames[i] << ">\n";
	ok = false;
      }
    }

    cout << "Process " << myid << ": files opened\n";

  }

  if (!ok) return;

#ifdef SEQCHECK
  int seq_check = 1;
  int seq_error = 0;
#endif

  while (1) {

    p = get_particles(&number);

    if (number == 0) break;

    if (myid==0) {
  
#ifdef DEBUG
    cout << "Process " << myid 
	 << ": received " << number << " from get_particles"
	 << endl << flush;
#endif    

      for (int i=0; i<number; i++) {

#ifdef SEQCHECK
	if (seq_check != p[i].sequence) {
	  if (!seq_error)
	    cout << "Process " << myid << ": tipsy dump out of order, "
		 << "seq_check=" << seq_check << " found=" << p[i].sequence
		 << endl << flush;
	  seq_error++;
	}
	seq_check++;
#endif
	switch(p[i].component) {

	case 0:			/* Point particles */
	case 1:			/* Dark particles */
	  dark.mass = p[i].mass;
	  for (int k=0; k<3; k++) dark.pos[k] = p[i].pos[k];
	  for (int k=0; k<3; k++) dark.vel[k] = p[i].vel[k];
	  dark.eps = 0;
	  dark.phi = p[i].pot;
	
	  ocmp[1]->write((char *)&dark, sizeof(dark));
	  ndark++;
	  break;
	  
	case -1:		/* Slab particles */
	case 2:			/* Disk particles */

	  star.mass = p[i].mass;
	  for (int k=0; k<3; k++) star.pos[k] = p[i].pos[k];
	  for (int k=0; k<3; k++) star.vel[k] = p[i].vel[k];
	  star.metals = 0;
	  star.tform = 0;
	  star.eps = 0;
	  star.phi = p[i].pot;
	
	  ocmp[2]->write((char *)&star, sizeof(star));
	  nstar++;
	  break;
	default:
	  cerr << "create_tipsy: Uh, oh, no such component number, " 
	       << p[i].component << " for index=" << i << "\n";
	  break;
	} /* No gas, yet . . . */
      }
    }
  }

  if (myid==0) {

    header.time    = tnow;
    header.ndim    = 3;
    header.nbodies = ngas+ndark+nstar;
    
    header.nsph    = ngas;
    header.ndark   = ndark;
    header.nstar   = nstar;
    
#ifdef DEBUG
    cerr << "T=" << tnow << "   ndark, nstar, ngas: "
	 << ndark << ", "
	 << nstar << ", "
	 << ngas  << endl << flush;

    if (seq_error)
      cerr << "Master: " << seq_error << " particles out of order"
	   << endl << flush;
#endif

    ocmp[0]->write((char *)&header, sizeof(header));
    
    for (int n=0; n<=ncomps; n++) {
      ocmp[n]->close();
      delete ocmp[n];
    }
    
  }

  tdump = tnow;

}

#ifdef SEQCHECK
void sequence_check(void)
{
  unsigned int counter;

  if (myid==0) counter = 1;
  else counter = nbodies_index[myid-1]+1;

  for (int i=1; i<=nbodies; i++) {
    if (sequence[i] != counter) {
      cerr << endl << "Process " << myid << ": sequence error, "
	   << sequence[i] << " != " << counter 
	   << endl << flush;
      break;
    }
    counter++;
  }

#ifdef DEBUG
  cout << "Process " << myid 
       << ": [seq_beg, seq_end]=[" 
       << seq_beg << ", " << seq_end 
       << "] [seq(1), seq(N)]=[" 
       << sequence[1] << ", " << sequence[nbodies] 
       << "]" << endl << flush;
#endif
}
#endif // SEQCHECK

				/* Normalization of rates to optimize 
				   n-body throughput */

extern "C" void recurse_norm(double maxrate, double* r)
{
  pair<double, int> element;
  vector< pair<double, int> > array;

  for (int i=0; i<numprocs; i++) {
    element.first = r[i];
    element.second = i;
    array.push_back(element);
  }


  double excess, norm;
  int iover = numprocs - 1;
  while (1) {

    sort(array.begin(), &array[iover]);
  
    excess = 0.0;
    for (; iover>=0; iover--) {
      if (array[iover].first < maxrate) break;
      excess += array[iover].first - maxrate;
      array[iover].first = maxrate;
    }

    
    if (excess == 0.0) break;

    norm = 0.0;
    for (int j=iover; j>=0; j--) norm += array[j].first;
    for (int j=iover; j>=0; j--) array[j].first += excess*array[j].first/norm;

  }

  if (excess > 0.0) {
    cerr << "recurse_norm: error renormalizing rates\n";
    cerr << "Rates:\n";
    for (int i=0; i<numprocs; i++)
      cerr << setw(5)  << i
	   << setw(15) << r[i]
	   << setw(15) << array[i].first
	   << setw(5)  << array[i].second
	   << endl;
  }

  for (int i=0; i<numprocs; i++)
    r[array[i].second] = array[i].first;
}


void shift_particles(int beg, int offset)
{

  if (offset>0) {
    for (int i=nbodies; i>=beg; i--) {
      component[i+offset] = component[i];
#ifdef SEQCHECK
      sequence[i+offset] = sequence[i];
#endif
      mass[i+offset] = mass[i];
      x[i+offset] = x[i];
      y[i+offset] = y[i];
      z[i+offset] = z[i];
      vx[i+offset] = vx[i];
      vy[i+offset] = vy[i];
      vz[i+offset] = vz[i];
      pot[i+offset] = pot[i];
      potext[i+offset] = potext[i];
      
      if (relx) esave[i+offset] = esave[i];
      if (scatter) mfp[i+offset] = mfp[i];
    }
  }

  if (offset<0) {
    for (int i=beg; i<=nbodies; i++) {
      component[i+offset] = component[i];
#ifdef SEQCHECK
      sequence[i+offset] = sequence[i];
#endif
      mass[i+offset] = mass[i];
      x[i+offset] = x[i];
      y[i+offset] = y[i];
      z[i+offset] = z[i];
      vx[i+offset] = vx[i];
      vy[i+offset] = vy[i];
      vz[i+offset] = vz[i];
      pot[i+offset] = pot[i];
      potext[i+offset] = potext[i];
      
      if (relx) esave[i+offset] = esave[i];
      if (scatter) mfp[i+offset] = mfp[i];
    }
  }

}


void pack_and_send(int orig, int dest, int ibeg, int iend, int number)
{
  MPI_Status status;
  
  static struct Partstruct *buf = NULL;
  if (!buf) buf = new Partstruct [nbuf];
  
  int counter = ibeg - seq_beg;
  int icount, ibufcount;

  if (myid == orig) {
    
#ifdef DEBUG
    if (ibeg == seq_beg)
      cout << "Process " << myid << ": ibeg=seq_beg (" 
	   << ibeg << ")" << endl << flush;

    if (iend == seq_end)
      cout << "Process " << myid << ": iend=seq_end (" 
	   << iend << ")" << endl << flush;
				// Sanity check
    if (ibeg != seq_beg  && iend != seq_end) {
      cerr << "Process " << myid << ": out of range,"
	   << " ibeg=" << ibeg
	   << " iend=" << iend
	   << " seq_beg=" << seq_beg
	   << " seq_end=" << seq_end
	   << endl << flush;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
#endif

    icount = 0;
    while (icount < number) {
      
      ibufcount = 0;

      for (int i=0; i<nbuf; i++) {
				// Pack structure
	buf[i].component = component[counter+1];
#ifdef SEQCHECK
	buf[i].sequence = sequence[counter+1];
#endif
	buf[i].mass = mass[counter+1];
	buf[i].pos[0] = x[counter+1];
	buf[i].pos[1] = y[counter+1];
	buf[i].pos[2] = z[counter+1];
	buf[i].vel[0] = vx[counter+1];
	buf[i].vel[1] = vy[counter+1];
	buf[i].vel[2] = vz[counter+1];
	buf[i].pot = pot[counter+1];
	buf[i].potext = potext[counter+1];
	if (relx) buf[i].esave = esave[counter+1];
	if (scatter) buf[i].mfp = mfp[counter+1];

	icount++;
	counter++;
	ibufcount++;
	
	if (icount == number) break;
      }

      MPI_Send(&ibufcount, 1, MPI_INT, dest, 2, MPI_COMM_WORLD);
      MPI_Send(buf, ibufcount, Particletype, dest, 3, MPI_COMM_WORLD);
      ibufcount = 0;
    }
				// Keep track of particle sequence
    if (ibeg==seq_beg) {
      shift_particles(number+1, -number);
      seq_beg = iend+1;
#ifdef DEBUG
      cout << "Process " << myid << ": shifting particles indexed "
	   << "[" << seq_beg << ", " << seq_end << "] down by " << number 
	   << endl << flush;
#endif
    } else seq_end = ibeg-1;

				// Shave particles from count
    nbodies -= number;

  }

}

void receive_and_unpack(int orig, int dest, int ibeg, int iend, int number)
{
  MPI_Status status;

  static struct Partstruct *buf = NULL;
  if (!buf) buf = new Partstruct [nbuf];
  
  int counter;
  int icount, ibufcount;

  if (myid == dest) {

    if (ibeg == seq_end+1) {
#ifdef DEBUG
      cerr << "Process " << myid << ": adding [" 
	   << ibeg << ", " << iend 
	   << "] to end of sequence beginning " << seq_end
	   << endl << flush;
#endif
      counter = nbodies;
    }
    else if (iend == seq_beg-1) {
#ifdef DEBUG
      cerr << "Process " << myid << ": adding [" 
	   << ibeg << ", " << iend 
	   << "] to start of sequence beginning " << seq_beg
	   << endl << flush;
#endif
      shift_particles(1, number);
      counter = 0;
    }
    else if (nbodies==0) {	// Empty node
#ifdef DEBUG
      cerr << "Process " << myid << ": adding [" 
	   << ibeg << ", " << iend 
	   << "] to currently empty node"
	   << endl << flush;
#endif
      seq_beg = ibeg;
      seq_end = iend;
      counter = 0;
    }
#ifdef DEBUG
    else {			// Sanity check
      cerr << "Process " << myid << ": out of range,"
	   << " ibeg=" << ibeg
	   << " iend=" << iend
	   << " seq_beg=" << seq_beg
	   << " seq_end=" << seq_end
	   << endl << flush;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
#endif

    icount = 0;
    while (icount < number) {
     
      MPI_Recv(&ibufcount, 1, MPI_INT, orig, 2, MPI_COMM_WORLD, &status);
      MPI_Recv(buf, ibufcount, Particletype, orig, 3, MPI_COMM_WORLD, &status);

      for (int i=0; i<ibufcount; i++) {
	component[counter+1] = buf[i].component;
#ifdef SEQCHECK
	sequence[counter+1] = buf[i].sequence;
#endif
	mass[counter+1] = buf[i].mass;
	x[counter+1] = buf[i].pos[0];
	y[counter+1] = buf[i].pos[1];
	z[counter+1] = buf[i].pos[2];
	vx[counter+1] = buf[i].vel[0];
	vy[counter+1] = buf[i].vel[1];
	vz[counter+1] = buf[i].vel[2];
	pot[counter+1] = buf[i].pot;
	potext[counter+1] = buf[i].potext;
	
	if (relx) esave[counter+1] = buf[i].esave;
	if (scatter) mfp[counter+1] = buf[i].mfp;
	
	icount++;
	counter++;
      }
    }
				// Add particles to count
    nbodies += number;

				// Keep track of particle sequence
    if (ibeg == seq_end+1)
      seq_end = iend;
    else
      seq_beg = ibeg;
  }

}

extern "C" void redistribute_particles(void)
{
  // Distribute current rate and number lists to all processors

  MPI_Bcast(rates, numprocs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  plist.clear();
  plist2.clear();
  plist.push_back((int)(nbodies_table[0]));
  plist2.push_back((int)(rates[0]*nbodies_tot));
  for (int i=1; i<numprocs-1; i++) {
    plist.push_back((int)(nbodies_table[i]) + plist[i-1]);
    plist2.push_back((int)(rates[i]*nbodies_tot) + plist2[i-1]);
  }
  plist.push_back(nbodies_tot);
  plist2.push_back(nbodies_tot);

#ifdef DEBUG
  if (myid==0) {
    cout << endl;
    cout << setw(5) << "Node" 
	 << setw(10) << "Old #" 
	 << setw(10) << "New #"
	 << endl
	 << setw(5) << "----" 
	 << setw(10) << "---------"
	 << setw(10) << "---------"
	 << endl;
    for (int i=0; i<numprocs; i++)
      cout << setw(5) << i
	   << setw(10) << plist[i]
	   << setw(10) << plist2[i]
	   << endl;
  }
#endif

  //========================================================================
  // Make merged list
  //========================================================================

  pair<int, int> id;
  pair<int, pair<int, int> > item;
  vector< pair<int, pair<int, int> > > total;

  for (int i=0; i<numprocs-1; i++) {

    id.first = i;
    id.second = 1;
    item.first = plist[i];
    item.second = id;

    total.push_back(item);

    id.first = i;
    id.second = 2;
    item.first = plist2[i];
    item.second = id;

    total.push_back(item);
  }

  sort(total.begin(), total.end());

#ifdef DEBUG
  if (myid==0) {
    cout << endl;
    cout << setw(5) << "Rank"
	 << setw(10) << "Number"
	 << setw(10) << "Node"
	 << setw(10) << "Old/New"
	 << endl
	 << setw(5) << "----"
	 << setw(10) << "--------"
	 << setw(10) << "----"
	 << setw(10) << "-------"
	 << endl;
    for (int i=0; i<2*(numprocs-1); i++)
      cout << setw(5) << i
	   << setw(10) << total[i].first
	   << setw(10) << total[i].second.first
	   << setw(10) << total[i].second.second
	   << endl;
  }
#endif
    

  //========================================================================
  // Check transaction algorithm
  //========================================================================

  int icur1=0, icur2=0, lastcnt=0;
  int number, ibeg, iend;
  
  vector<int> check = plist;
  for (int i=1; i<numprocs; i++) check[i] = plist[i] - plist[i-1];

#ifdef DEBUG
  if (myid==0) cout << endl;
#endif

  for (int i=0; i<2*(numprocs-1); i++) {

    if (icur1 != icur2) {
				// Put strategy in here . . .

      number = total[i].first - lastcnt;
      ibeg = lastcnt+1;
      iend = total[i].first;

      pack_and_send(icur1, icur2, ibeg, iend, number);
      receive_and_unpack(icur1, icur2, ibeg, iend, number);

#ifdef DEBUG
      if (myid==0) {
	cout << "Node " << icur1 << " sending " << number
	     << " particles to Node " << icur2 
	     << " [" << ibeg << ", " << iend << "]" << endl << flush;
      }
#endif
      MPI_Barrier(MPI_COMM_WORLD);

      check[icur1] -= number;
      check[icur2] += number;
    }
      
    if (total[i].second.second == 2) {

      if (total[i].second.first != icur2)
	cout << "  Mismatch: " << total[i].second.first << " !=" 
	     << icur2 << endl;
      icur2++;

    } else {
	  
      if (total[i].second.first != icur1)
	cout << "  Mismatch: " << total[i].second.first << " !=" 
	     << icur2 << endl;
      icur1++;
      
    }

    lastcnt = total[i].first;

  }

#ifdef DEBUG
  if (myid==0) {
    
    cout << endl;
    for (int i=1; i<numprocs; i++) check[i] += check[i-1];
    
    bool ok = true;
    for (int i=0; i<numprocs; i++)
      if (plist2[i] != check[i]) ok = false;
      
    if (!ok) {
      cout << "Fail\n";
      for (int i=0; i<numprocs; i++)
	cout << setw(5) << i
	     << setw(10) << plist2[i]
	     << setw(10) << check[i]
	     << endl;
    }

  }
#endif

				// Reset cached list to new list
  nbodies_table[0] = nbodies_index[0] = plist2[0];
  for (int i=1; i<numprocs; i++) {
    nbodies_index[i] = plist2[i];
    nbodies_table[i] = plist2[i] - plist2[i-1];
  }

  //========================================================================

#ifdef SEQCHECK
  sequence_check();
#endif

}

