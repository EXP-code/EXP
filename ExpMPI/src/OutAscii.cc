#include <iostream>
#include <fstream>
#include <iomanip>
#include <strstream>

#include "expand.h"
#include <global.H>

#include <AxisymmetricBasis.H>
#include <OutAscii.H>

OutAscii::OutAscii(string& line) : Output(line)
{
  nint = 100;
  nbeg = 0;
  name = "";
  accel = false;
  initialize();

  if (name.size()>0) {
				// Look for the fiducial component
    bool found = false;
    list<Component*>::iterator cc;
    Component *c;
    for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
      c = *cc;
      if ( !name.compare(c->name) ) {
	c0 = c;
	found = true;
      break;
      }
    }

    if (!found) {
      cerr << "Process " << myid << ": can't find desired component <"
	   << name << ">" << endl;
      MPI_Abort(MPI_COMM_WORLD, 35);
    }

  }
  else
    c0 = NULL;
}

void OutAscii::initialize()
{
  string tmp;

  if (Output::get_value(string("nint"), tmp))  nint = atoi(tmp.c_str());
  if (Output::get_value(string("nbeg"), tmp))  nbeg = atoi(tmp.c_str());
  if (Output::get_value(string("name"), tmp))  name = tmp;
  if (!Output::get_value(string("filename"), filename)) {
    filename.erase();
    filename = "OUTASC\0";
  }
  if (Output::get_value(string("accel"), tmp)) {
    if (atoi(tmp.c_str())) accel = true;
  }
				// Determine last file

  if (restart && nbeg==0 && myid==0) {

    for (nbeg=0; nbeg<100000; nbeg++) {

				// Output name
      ostrstream fname;
      fname << filename << "." << setw(5) << setfill('0') << nbeg << '\0';

				// See if we can open file
      ofstream out(fname.str(), ios::out | ios::noreplace);

      if (out) {
	out.close();
	ostrstream command;
	command << "rm " <<  fname.str() << '\0';
	system(command.str());
	cout << "OutAscii: will begin with nbeg=" << nbeg << endl;
	break;
      }
    }
  }
}


void OutAscii::Run(int n, bool last)
{
  if (n % nint && !last) return;
  if (!c0) return;

  synchronize_velocity(1);

  ofstream *out;

  if (myid==0) {
				// Output name
    ostrstream fname;
    fname << filename << "." << setw(5) << setfill('0') << nbeg++ << '\0';

				// Open file and write master header
    out = new ofstream(fname.str(), ios::out | ios::noreplace);

    if (!*out) {
      cerr << "OutAscii: can't open file <" << fname.str() 
	   << "> . . . quitting\n";
      MPI_Abort(MPI_COMM_WORLD, 33);
    }
    
    *out << "# Time=" << tnow << "\n";
    *out << setw(10) << c0->nbodies_tot
	 << setw(10) << c0->niattrib
	 << setw(10) << c0->ndattrib << "\n";
  }

  c0->write_ascii(out);

  if (myid==0) {
    out->close();
    delete out;
  }

  synchronize_velocity(0);

}

