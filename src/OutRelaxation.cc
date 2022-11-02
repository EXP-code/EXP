#include <iostream>
#include <iomanip>
#include <fstream>

#include <global.H>
#include <Component.H>
#include <ComponentContainer.H>
#include <OutRelaxation.H>

const std::set<std::string> OutRelaxation::valid_keys = {
  "suffix",
  "epos"
};

OutRelaxation::OutRelaxation(const YAML::Node& conf) : Output(conf)
{
  id = "OutRelaxation";

  epos = 0;			// Default

  initialize();

				// Initialize output file
  if (myid==0) {
    
    fname = "relx.";
    fname += suffix;

    ofstream out(fname.c_str(), ios::out | ios::app);
    if (!out) {
      string msg("Couldn't open <");
      msg += fname + ">";
      bomb(msg);
    }

    out << "! 1) time 2) step 3) Delta E; 4) Root variance; 5) |Delta E|\n";

    cout << "Created an OutRelaxation with dattrib index = " << epos << "\n";
  }

}

void OutRelaxation::initialize()
{
  // Remove matched keys
  for (auto v : valid_keys) current_keys.erase(v);
  
  // Assign values from YAML
  //
  try {
				// Get file name
    if (conf["suffix"])
      suffix = conf["suffix"].as<std::string>();
    else
      suffix = "out";

    if (conf["epos"]) epos = conf["epos"].as<int>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in OutRelaxation: "
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

void OutRelaxation::Run(int n, int mstep, bool final)
{
  double delta_e, e_average=0.0, e_absolute=0.0, variance=0.0;
  double e_average1=0.0, e_absolute1=0.0, variance1=0.0, esave;
  int used1 = 0, used0 = 0;
  

  int nbodies;

  for (auto c : comp->components) {

    nbodies = c->Number();

    PartMapItr it = c->Particles().begin();
    unsigned long i;

    for (int q=1; q<=nbodies; q++) {

      i = it->first; it++;

      if (c->freeze(i)) continue;


      esave = c->Part(i)->dattrib[epos];

      delta_e = 0.0;
      for (int j=0; j<3; j++) 
	delta_e += c->Vel(i, j) * c->Vel(i, j);

      delta_e = 0.5*c->Mass(i)*delta_e + 
	c->Mass(i)*(c->Part(i)->pot + c->Part(i)->potext) -
	esave;
	
      e_average1 += delta_e/esave;
      e_absolute1 += fabs(delta_e/esave);
      variance1 += (delta_e*delta_e/(esave*esave));
      used1++;
    }

  }

				/* Collect values */

  MPI_Reduce(&used1, &used0, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&e_average1, &e_average, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&e_absolute1, &e_absolute, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&variance1, &variance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (myid == 0) {

    e_average /= used0;
    e_absolute /= used0;
    variance = (variance - e_average*e_average)/(used0 - 1);

    ofstream out(fname.c_str(), ios::out | ios::app);
    if (!out) {
      cout << "OutRelaxation: Couldn't reopen <" << fname << ">\n";
      return;
    }
    out.setf(ios::scientific);
    out << setw(14) << tnow
	<< " " << setw(5) << n
	<< " " << setw(14) << e_average
	<< " " << setw(14) << sqrt(variance)
	<< " " << setw(14) << e_absolute
	<< endl;
  }

}
