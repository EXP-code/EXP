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


//===========================================================================

void usage(char *prog)
{
  cout << "Usage:" << endl << endl
       << prog << " [options]" << endl << endl
       << setw(15) << "Option" << setw(10) << "Argument" << setw(10) << " " 
       << setiosflags(ios::left)
       << setw(40) << "Description" << endl
       << resetiosflags(ios::left)
       << endl
       << setw(15) << "-m or --mpi" << setw(10) << "No" << setw(10) << " " 
       << setiosflags(ios::left)
       << setw(40) << "Turn on MPI for SL computation" << endl
       << resetiosflags(ios::left)
       << setw(15) << "-c or --cmap" << setw(10) << "Yes" << setw(10) << " " 
       << setiosflags(ios::left)
       << setw(40) << "Use mapped rather than linear coordinates (0|1|2)" << endl
       << resetiosflags(ios::left)
       << setw(15) << "-n or --numr" << setw(10) << "Yes" << setw(10) << " " 
       << setiosflags(ios::left)
       << setw(40) << "Number of points in radial table" << endl
       << resetiosflags(ios::left)
       << setw(15) << "-i or --file" << setw(10) << "SLGridSph.model" << setw(10) << " " 
       << setiosflags(ios::left)
       << setw(40) << "Model profile" << endl
       << resetiosflags(ios::left)
       << setw(15) << "-l or --logr" << setw(10) << "No" << setw(10) << " " 
       << setiosflags(ios::left)
       << setw(40) << "Turn on log spacing in radius" << endl
       << resetiosflags(ios::left)
       << endl;

  exit(0);
}

int main(int argc, char** argv)
{
  bool use_mpi = false;
  bool use_logr = false;
  int cmap = 0;
  double scale = 1.0;
  int numr = 10000;
  int diverge = 0;
  double dfac = 1.0;
  string filename = "SLGridSph.model";

  int c;
  while (1) {
    int this_option_optind = optind ? optind : 1;
    int option_index = 0;
    static struct option long_options[] = {
      {"mpi", 0, 0, 0},
      {"logr", 0, 0, 0},
      {"cmap", 1, 0, 0},
      {"numr", 1, 0, 0},
      {"dfac", 1, 0, 0},
      {"file", 1, 0, 0},
      {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "c:lmn:d:i:h",
		     long_options, &option_index);

    if (c == -1) break;

    switch (c) {
    case 0:
      {
	string optname(long_options[option_index].name);

	if (!optname.compare("mpi")) {
	  use_mpi = true;
	} else if (!optname.compare("cmap")) {
	  cmap = atoi(optarg);
	} else if (!optname.compare("logr")) {
	  use_logr = true;
	} else if (!optname.compare("numr")) {
	  numr = atoi(optarg);
	} else if (!optname.compare("dfac")) {
	  diverge = 1;
	  dfac = atof(optarg);
	} else if (!optname.compare("file")) {
	  filename = optarg;
	} else {
	  cout << "Option " << long_options[option_index].name;
	  if (optarg) cout << " with arg " << optarg;
	  cout << " is not defined " << endl;
	  exit(0);
	}
      }
      break;

    case 'c':
      cmap = atoi(optarg);
      break;

    case 'l':
      use_logr = true;
      break;

    case 'm':
      use_mpi = true;
      break;

    case 'n':
      numr = atoi(optarg);
      break;

    case 'd':
      diverge = 1;
      dfac = atof(optarg);
      break;

    case 'i':
      filename = optarg;
      break;

    case 'h':
    default:
      usage(argv[0]);
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

				// Set model file
  SLGridSph::model_file_name = filename;

  cout << "Model=" << filename << endl;
  cout << "CMAP=" << cmap << endl;

				// Generate Sturm-Liouville grid
  SLGridSph *ortho = new SLGridSph(Lmax, nmax, numr, rmin, rmax, 
				   cmap, rs, diverge, dfac);


				// Slaves exit
  if (use_mpi && myid>0) {
    MPI_Finalize();
    exit(0);
  }

				// Do what?
  while (1) {
    bool done=false;
    int iwhich;

    cout << "Task:" << endl;
    cout << "1: Print out density, potential pairs" << endl;
    cout << "2: Check orthogonality" << endl;
    cout << "3: Quit" << endl;
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
	  cout << "Can't open <" << filename << "> for output" << endl;
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

	double x, r, lrmin, lrmax;

	if (use_logr) {
	  if (rmin<1.0e-16) use_logr = false;
	  else {
	    lrmin = log(rmin);
	    lrmax = log(rmax);
	  }
	}

	for (int i=0; i<num; i++) {

	  if (use_logr)
	    x = ortho->r_to_xi(exp(lrmin + (lrmax-lrmin)*i/(num-1)));
	  else
	    x = ximin + (ximax-ximin)*i/(num-1);

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

	  ans += r*r*ortho->get_pot(x, L, N1, 0)*
	    ortho->get_dens(x, L, N2, 0) /
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

  if (use_mpi) MPI_Finalize();

  return 0;
}




