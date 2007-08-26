#include <iostream>
#include <iomanip>
#include <fstream>
#include <global.H>

#include <OutRelaxation.H>

OutRelaxation::OutRelaxation(string& line) : Output(line)
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
  string tmp;

				// Get file name
  if (!get_value(string("suffix"), suffix)) {
    suffix.erase();
    suffix = "out";
  }

  if (get_value(string("epos"), tmp))  epos = atoi(tmp.c_str());
}

void OutRelaxation::Run(int n, bool final)
{
  double delta_e, e_average=0.0, e_absolute=0.0, variance=0.0;
  double e_average1=0.0, e_absolute1=0.0, variance1=0.0, esave;
  int used1 = 0, used0 = 0;
  

  int nbodies;
  list<Component*>::iterator cc;
  Component* c;

  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {

    c = *cc;

    nbodies = c->Number();

    for (int i=1; i<=nbodies; i++) {

      if (c->freeze(*(c->Part(i)))) continue;


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
