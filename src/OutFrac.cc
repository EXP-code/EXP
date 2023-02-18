#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <unistd.h>		// for unlink

#include <expand.H>
#include <Timer.H>
#include <OutFrac.H>


const double default_quant[] = {0.001, 0.003, 0.01, 0.03, 0.1, 0.2, 0.4, 0.5, 0.6, 0.8, 0.9, 0.97, 0.99, 0.993, 0.999};

const std::set<std::string>
OutFrac::valid_keys = {
  "filename",
  "nint",
  "nintsub",
  "name",
  "frac"
};

OutFrac::OutFrac(const YAML::Node& conf) : Output(conf)
{
  nint = 10;
  nintsub = std::numeric_limits<int>::max();
  filename = outdir + "OUTFRAC." + runtag;
  tcomp = NULL;
  numQuant = sizeof(default_quant)/sizeof(double);
  for (int i=0; i<numQuant; i++) Quant.push_back(default_quant[i]);

  initialize();

  if (!tcomp) {
    if (myid==0) {
      std::cerr << "OutFrac: no component to trace\n";
    }
    MPI_Finalize();
    exit(112);
  }

  if (numQuant==0) {
    if (myid==0) {
      std::cerr << "OutFrac: no quantiles defined!\n";
    }
    MPI_Finalize();
    exit(113);
  }
  else if (myid==0)
    std::cout << "OutFrac: using " << numQuant << " quantiles\n";


				// If not a restart, make quantile header
  if (myid==0) {

    if (restart) {
      string backfile = filename + ".bak";
				// Remove an old backup files
      if (unlink(backfile.c_str())) {
	perror("OutFrac::Run()");
	std::cout << "OutFrac::Run(): error unlinking old backup file <" 
		  << backfile << ">" << std::endl;
      } else {
	std::cout << "OutFrac::Run(): successfully unlinked <"
		  << backfile << ">" << std::endl;
      }
      if (rename(filename.c_str(), backfile.c_str())) {
	perror("OutFrac::Run()");
	std::ostringstream sout;
	sout << "OutFrac: error renaming the current file <"
	     << filename << "> to the backup file <" 
	     << backfile << ">";
	throw GenericError(sout.str(), __FILE__, __LINE__, 114, true);
      } else {
	std::cout << "OutFrac::Run(): successfully renamed <"
		  << filename << "> to <" << backfile << ">" << std::endl;
      }

      std::ifstream in(backfile.c_str());
      if (!in) {
	std::ostringstream sout;
	sout << "OutFrac: error opening backup file <" 
	     << backfile << "> for input";
	throw GenericError(sout.str(), __FILE__, __LINE__, 115, true);
      }

      std::ofstream out(filename.c_str());
      if (!out) {
	std::ostringstream sout;
	sout << "OutFrac: error opening new file <" 
	     << filename << "> for output";
	throw GenericError(sout.str(), __FILE__, __LINE__, 116, true);
      }

      const unsigned linesz = 4196;
      char line[linesz];
      double time;

      in.get(line, linesz);	// Copy over the header
      while (line[0] == '#') {
	out << line;
	in.get(line, linesz);
      }
				// Copy all records with times earlier
				// than the current time
      while (in.good() && !in.eof())  {
	istringstream sin(line);
	sin >> time;
	if (time>=tnow) break;
	out << line;
	in.get(line, linesz);
      }

    } else {			// Make the header

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

}

void OutFrac::initialize()
{
  // Remove matched keys
  //
  for (auto v : valid_keys) current_keys.erase(v);
  
  // Assign values from YAML
  //
  try {
				// Get file name
    if (Output::conf["filename"])
      filename = Output::conf["filename"].as<std::string>();

    if (Output::conf["nint"])
      nint     = Output::conf["nint"]    .as<int>();


    if (Output::conf["nintsub"]) {
      nintsub     = Output::conf["nintsub"]    .as<int>();
      if (nintsub <= 0) nintsub = 1;
    }

				// Search for desired component
    if (Output::conf["name"]) {
      std::string tmp = Output::conf["name"].as<std::string>();
      for (auto c : comp->components) {
	if (!(c->name.compare(tmp))) tcomp  = c;
      }
    }

    // Get quantiles
    if (Output::conf["frac"])
      Quant = Output::conf["frac"].as<std::vector<double>>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in OutFrac: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }

}

void OutFrac::Run(int n, int mstep, bool last)
{
  if (n % nint != 0 && !last) return;
  if (multistep>1 and mstep % nintsub !=0) return;

  MPI_Status status;

  Timer timer;

  if (myid==0) timer.start();

#ifdef HAVE_LIBCUDA
  if (use_cuda) {
    if (tcomp->force->cudaAware() and not comp->fetched[tcomp]) {
      comp->fetched[tcomp] = true;
      tcomp->CudaToParticles();
    }
  }
#endif
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
  PartMapItr it = tcomp->Particles().begin();
  unsigned long j;

  for (int n=0; n<tcomp->Number(); n++) {
    j = it++->first;
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

    if (tcomp->CurTotal() != rtot.size()) {
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
    out << setw(18) << timer.stop();
    out << endl;
  }

}
