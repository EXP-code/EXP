#include <iostream>
#include <iomanip>
#include <string>

#include <math.h>
#include <getopt.h>

#include <SLGridMP2.h>
#include <gaussQ.h>

int numprocs, myid, proc_namelen;
char processor_name[MPI_MAX_PROCESSOR_NAME];

MPI_Comm MPI_COMM_SLAVE;

char threading_on = 0;
pthread_mutex_t mem_lock;


int main(int argc, char** argv)
{
  bool use_mpi = false;
  int cmap = 0;
  double scale = 1.0;

  int c;
  while (1) {
    int this_option_optind = optind ? optind : 1;
    int option_index = 0;
    static struct option long_options[] = {
      {"mpi", 0, 0, 0},
      {"cmap", 0, 0, 0},
      {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "",
		     long_options, &option_index);

    if (c == -1) break;

    switch (c) {
    case 0:
      {
	string optname(long_options[option_index].name);

	if (!optname.compare("mpi")) {
	  use_mpi = true;
	} else if (!optname.compare("cmap")) {
	  cmap = 1;
	} else {
	  cout << "Option " << long_options[option_index].name;
	  if (optarg) cout << " with arg " << optarg;
	  cout << " is not defined \n";
	  exit(0);
	}
      }
      break;

    case 'c':
      cmap = 1;
      break;

    case 'm':
      use_mpi = true;
      break;
    }

  }

  //===================
  // MPI preliminaries 
  //===================
  if (use_mpi) {
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Get_processor_name(processor_name, &proc_namelen);
  }

  //===================
  // Get info
  //===================

  int Lmax, nmax;
  int numr = 1000;
  double rmin, rmax, rs;
  
  if (!use_mpi || myid==0) {

				// Get info from user
    cout << "Lmax, Nmax? ";
    cin >> Lmax;
    cin >> nmax;
    
    cout << "Rmin, Rmax, Rs? ";
    cin >> rmin;
    cin >> rmax;
    cin >> rs;
  }

  if (use_mpi) {

    MPI_Bcast(&Lmax, 1, MPI_INT,    0, MPI_COMM_WORLD);
    MPI_Bcast(&nmax, 1, MPI_INT,    0, MPI_COMM_WORLD);
    MPI_Bcast(&rmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rs,   1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast(&Lmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    SLGridSph::mpi = 1;		// Turn on MPI

  } else {

    SLGridSph::mpi = 0;		// Turn off MPI
  }

				// Generate Sturm-Liouville grid
  SLGridSph *ortho = new SLGridSph(Lmax, nmax, numr, rmin, rmax, 
				   cmap, rs);


				// Do what?
  while (1) {
    bool done=false;
    int iwhich;

    cout << "Task:\n";
    cout << "1: Print out density, potential pairs\n";
    cout << "2: Check orthogonality\n";
    cout << "3: Quit\n";
    cout << "?? ";
    cin >> iwhich;

    switch(iwhich) {
    case 1:
      {
	string filename;
	cout << "Filename? ";
	cin >> filename;
	ofstream out (filename.c_str());
	if (!out) {
	  cout << "Can't open <" << filename << "> for output\n";
	  break;
	}

	cout << "Number of points? ";
	int num;
	cin >> num;

	cout << "L, N? ";
	int L, N;
	cin >> L;
	cin >> N;

	double ximin = ortho->r_to_xi(rmin);
	double ximax = ortho->r_to_xi(rmax);

	double x, r;

	for (int i=0; i<num; i++) {
	  x = ximin + (ximax-ximin)*(0.5 + i)/num;
	  r = ortho->xi_to_r(x);
	  out << setw(15) << x
	      << setw(15) << r
	      << setw(15) << ortho->get_pot(x, L, N, 0)
	      << setw(15) << ortho->get_force(x, L, N, 0)
	      << setw(15) << ortho->get_dens(x, L, N, 0)
	      << endl;
	}
      }

      break;

    case 2:
      {
	cout << "Number of knots? ";
	int num;
	cin >> num;

	LegeQuad lw(num);

	cout << "L, N1, N2? ";
	int L, N1, N2;
	cin >> L;
	cin >> N1;
	cin >> N2;

	double ximin = ortho->r_to_xi(rmin);
	double ximax = ortho->r_to_xi(rmax);

	double x, r, ans=0.0;
	for (int i=0; i<num; i++) {

	  x = ximin + (ximax - ximin)*lw.knot(i+1);
	  r = ortho->xi_to_r(x);

	  ans += r*r*ortho->get_pot(x, L, N1, 0)*ortho->get_dens(x, L, N2, 0) /
	    ortho->d_xi_to_r(x) * (ximax - ximin)*lw.weight(i+1);

	}

	cout << "<" << N1 << "|" << N2 << "> = " << ans << endl;
      }

      break;

    default:
      done = true;
      break;
    }

    if (done) break;
  }

  delete ortho;
  return 0;
}
