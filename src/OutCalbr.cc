#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <expand.h>
#include <Timer.h>
#include <OutCalbr.H>


OutCalbr::OutCalbr(string& line) : Output(line)
{
  nint = 10;
  filename = outdir + "OUTCALBR." + runtag;
  tcomp = NULL;
  num = 10;

  initialize();

  if (!tcomp) {
    if (myid==0) {
      cerr << "OutCalbr: no component to trace\n";
      MPI_Abort(MPI_COMM_WORLD, 112);
    }
  }
}


void OutCalbr::set_energies()
{
  if (!restart) {
				// Compute energies and angular momentum
    
    double Emin1=1e30, Emax1=-1e30, v2, E;

    PartMapItr it = tcomp->Particles().begin();
    unsigned long n;
    Particle *p = tcomp->Part(n);

    for (int q=0; q<tcomp->Number(); q++) {
      n = (it++)->first;

      v2 = 0.0;
      for (int j=0; j<3; j++) 
	v2 += p->vel[j]*p->vel[j];
				// Energy
      p->dattrib[0] = E =
	0.5*v2 + p->pot + p->potext;
				// Lx
      p->dattrib[1] = 
	p->pos[1]*p->vel[2] - 
	p->pos[2]*p->vel[1] ;
      
				// Ly
      p->dattrib[2] = 
	p->pos[2]*p->vel[0] -
	p->pos[0]*p->vel[2] ; 
      
				// Lz
      p->dattrib[3] = 
	p->pos[0]*p->vel[1] - 
	p->pos[1]*p->vel[0] ;

      Emin1 = min<double>(Emin1, E);
      Emax1 = max<double>(Emax1, E);
    }

    MPI_Allreduce(&Emin1, &Emin, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&Emax1, &Emax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    dE = (Emax - Emin)/num;
    Ec = vector<double>(num);
    for (int i=0; i<num; i++) Ec[i] = Emin + dE*(0.5+i);

    deltaE   = vector<double>(num);
    deltaLx  = vector<double>(num);
    deltaLy  = vector<double>(num);
    deltaLz  = vector<double>(num);

    deltaE1  = vector<double>(num);
    deltaLx1 = vector<double>(num);
    deltaLy1 = vector<double>(num);
    deltaLz1 = vector<double>(num);

    ncnt     = vector<unsigned>(num);
    ncnt1    = vector<unsigned>(num);
    
				// Make the data header
    if (myid==0) {
      ofstream out(filename.c_str());
      if (!out) {
	cout << "OutCalbr: can't open file <" << filename << ">\n";
      }

      out.setf(ios::left);
      out << setw(18) << "# Time";
      for (int i=0; i<num; i++) {
	ostringstream label;
	label << setprecision(3) << fixed 
	      << "| [" << Ec[i]-0.5*dE 
	      << ", "  << Ec[i] 
	      << ", "  << Ec[i]+0.5*dE << ")";
	out << setw(4*18) << label.str();
      }
      out << endl << setw(18) << "#";
      string labs[4] = {"| E", "Lx", "Ly", "Lz"};
      for (int i=0; i<num; i++) {
	for (int j=0; j<4; j++) out << setw(18) << labs[j];
      }
      out << endl;
      out.fill('-');
      out << setw(18) << "# 1 ";
      int icnt=2;
      for (int i=0; i<num*4; i++) {
	ostringstream label;
	label << "| [" << icnt++ << "] ";
	out << setw(18) << label.str();
      }
      out << endl;
    }
  }
}


void OutCalbr::initialize()
{
  string tmp;
				// Get file name
  get_value(string("filename"), filename);
  
  if (get_value(string("nint"), tmp)) 
    nint = atoi(tmp.c_str());

  if (get_value(string("N"), tmp)) 
    num = atoi(tmp.c_str());

				// Search for desired component
  if (get_value(string("name"), tmp)) {
    for (auto c : comp->components) {
      if (!(c->name.compare(tmp))) tcomp  = c;
    }
  }

}

void OutCalbr::Run(int ns, bool last)
{
  if (ns==0) set_energies();

  if (ns % nint != 0 && !last) return;

  MPI_Status status;

  Timer timer;

  if (myid==0) timer.start();

				// Open output file
  ofstream out;
  if (myid==0) {
    out.open(filename.c_str(), ios::out | ios::app);
    if (!out) {
      cout << "OutCalbr: can't open file <" << filename << ">\n";
      return;
    }
  }
  
  for (int i=0; i<num; i++) {
    deltaE1 [i] = 0.0;
    deltaLx1[i] = 0.0;
    deltaLy1[i] = 0.0;
    deltaLz1[i] = 0.0;
    ncnt[i]     = 0;
  }
				// Compute energies and angmom
  double E, Lx, Ly, Lz, v2;
  int indx;
  PartMapItr it = tcomp->Particles().begin();
  unsigned long n;

  for (int q=0; q<tcomp->Number(); q++) {
    n = (it++)->first;

    Particle *p = tcomp->Part(n);

    v2 = 0.0;
    for (int j=0; j<3; j++) 
      v2 += p->vel[j]*p->vel[j];

    E =	0.5*v2 + p->pot + p->potext;

    if (E<Emin || E>=Emax) continue;
    indx = min<int>((int)floor((E-Emin)/dE), num-1);

    Lx = 
      p->pos[1]*p->vel[2] - 
      p->pos[2]*p->vel[1] ;
      
    Ly = 
      p->pos[2]*p->vel[0] -
      p->pos[0]*p->vel[2] ; 
      
    Lz = 
      p->pos[0]*p->vel[1] - 
      p->pos[1]*p->vel[0] ;

    deltaE1[indx] += 
      (E - p->dattrib[0])*
      (E - p->dattrib[0]);

    deltaLx1[indx] += 
      (Lx - p->dattrib[1])*
      (Lx - p->dattrib[1]);

    deltaLy1[indx] += 
      (Ly - p->dattrib[2])*
      (Ly - p->dattrib[2]);

    deltaLz1[indx] += 
      (Lz - p->dattrib[3])*
      (Lz - p->dattrib[3]);

    ncnt1[indx]++;
  }

				// Send arrays to master
  MPI_Reduce(&ncnt1[0], &ncnt[0], num, MPI_UNSIGNED, MPI_SUM, 
	     0, MPI_COMM_WORLD);

  MPI_Reduce(&deltaE1[0], &deltaE[0], num, MPI_DOUBLE, MPI_SUM, 
	     0, MPI_COMM_WORLD);

  MPI_Reduce(&deltaLx1[0], &deltaLx[0], num, MPI_DOUBLE, MPI_SUM, 
	     0, MPI_COMM_WORLD);

  MPI_Reduce(&deltaLy1[0], &deltaLy[0], num, MPI_DOUBLE, MPI_SUM, 
	     0, MPI_COMM_WORLD);

  MPI_Reduce(&deltaLz1[0], &deltaLz[0], num, MPI_DOUBLE, MPI_SUM, 
	     0, MPI_COMM_WORLD);


  if (myid==0) {

    out.setf(ios::left);
    out << setw(18) << tnow;
    
    for (int i=0; i<num; i++) {
      if (ncnt[i]>0) {
	out << setw(18) << sqrt(deltaE[i] /ncnt[i])
	    << setw(18) << sqrt(deltaLx[i]/ncnt[i])
	    << setw(18) << sqrt(deltaLy[i]/ncnt[i])
	    << setw(18) << sqrt(deltaLz[i]/ncnt[i]);
      } else {
	out << setw(18) << 0.0
	    << setw(18) << 0.0
	    << setw(18) << 0.0
	    << setw(18) << 0.0;
      }
    }
    out << endl;
  }

}
