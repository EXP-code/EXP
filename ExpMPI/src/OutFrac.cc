#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <expand.h>
#include <Timer.h>
#include <OutFrac.H>


OutFrac::OutFrac(string& line) : Output(line)
{
  nint = 10;
  filename = outdir + "OUTFRAC." + runtag;
  tcomp = NULL;
  numQuant = 0;

  initialize();

  if (!tcomp) {
    if (myid==0) {
      cerr << "OutFrac: no component to trace\n";
      MPI_Abort(MPI_COMM_WORLD, 112);
    }
  }

  if (myid==0) {
    if (numQuant==0) {
      cerr << "OutFrac: no quantiles defined!\n";
      MPI_Abort(MPI_COMM_WORLD, 113);
    }
    else
      cout << "OutFrac: using " << numQuant << " quantiles\n";
  }

				// If not a restart, make quantile header
  if (!restart && myid==0) {
    ofstream out(filename.c_str());
    if (!out) {
      cout << "OutFrac: can't open file <" << filename << ">\n";
    }

    out.setf(ios::left);
    out << setw(18) << "# Time";
    for (int i=0; i<numQuant; i++) {
      ostringstream label;
      label << "| " << Quant[i];
      out << setw(18) << label.str();
    }
    out << setw(18) << "| elapsed time" << endl;
    out.fill('-');
    out << setw(18) << "# 1 ";
    for (int i=0; i<=numQuant; i++) {
      ostringstream label;
      label << "| " << i+2 << " ";
      out << setw(18) << label.str();
    }
    out << endl;
  }
}

void OutFrac::initialize()
{
  string tmp;
				// Get file name
  get_value(string("filename"), filename);
  
  if (get_value(string("nint"), tmp)) 
    nint = atoi(tmp.c_str());

				// Search for desired component
  if (get_value(string("name"), tmp)) {
    list<Component*>::iterator cc;
    Component* c;
    for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
      c = *cc;
      if (!(c->name.compare(tmp))) tcomp  = c;
    }
  }

				// Get quantiles

  string val;
  for (numQuant=0; numQuant<1000; numQuant++) {
    ostringstream count;
    count << "frac(" << numQuant+1 << ")";
    if (get_value(count.str(), val)) {
      Quant.push_back(atof(val.c_str()));
    } else break;
  }

}

void OutFrac::Run(int n, bool last)
{
  if (n % nint != 0 && !last) return;

  MPI_Status status;

  Timer timer(true);

  if (myid==0) timer.start();

				// Open output file
  ofstream out;
  if (myid==0) {
    out.open(filename.c_str(), ios::out | ios::app);
    if (!out) {
      cout << "OutFrac: can't open file <" << filename << ">\n";
      return;

    }
  }
  
				// Compute R and rank by radius
  double r, pos[3];
  vector<double> rad(tcomp->Number());
  map<unsigned long, Particle>::iterator it = tcomp->Particles().begin();
  unsigned long j;

  for (int n=0; n<tcomp->Number(); n++) {
    j = it->first;
    it++;
    tcomp->Pos(pos, j, Component::Centered);
    r = 0.0;
    for (int j=0; j<3; j++) r += pos[j]*pos[j];
    rad[n] = sqrt(r);
  }

				// Send arrays to master
  int nbodies = tcomp->Number();
  int max_nbodies, cur_bodies;
  MPI_Reduce(&nbodies, &max_nbodies,
	     1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

  vector<double> rtot;
  if (myid==0) {
    rtot = rad;
    rad = vector<double>(max_nbodies);
  }

  for (int n=1; n<numprocs; n++) {
    if (myid==n) {
      MPI_Send(&nbodies, 1, MPI_INT, 0, 153, MPI_COMM_WORLD);
      MPI_Send(&rad[0], nbodies, MPI_DOUBLE, 0, 154, MPI_COMM_WORLD);
    }
    if (myid==0) {
      MPI_Recv(&cur_bodies, 1, MPI_INT, n, 153, MPI_COMM_WORLD, &status);
      MPI_Recv(&rad[0], cur_bodies, MPI_DOUBLE, n, 154, MPI_COMM_WORLD, &status);
      rtot.insert(rtot.end(), rad.begin(), rad.begin()+cur_bodies);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  if (myid==0) {

    if (tcomp->nbodies_tot != rtot.size()) {
      cerr << "OutFrac: body count mismatch!\n";
    }

    out.setf(ios::left);
    out << setw(18) << tnow;
    
				// Sort the radii
    sort(rtot.begin(), rtot.end());
    
				// Send all radii to 
				// Put quantiles into file
    int indx;

    for (int i=0; i<numQuant; i++) {
				// get the index (nearest integer)
      indx = (int)(Quant[i]*rtot.size()+0.5);
      if (indx >= rtot.size()) indx = rtot.size()-1;

      out << setw(18) << rtot[indx];
    }
    out << setw(18) << timer.stop().getRealTime();
    out << endl;
  }

}
