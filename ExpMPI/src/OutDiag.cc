#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "expand.h"

#include <AxisymmetricBasis.H>
#include <OutDiag.H>

OutDiag::OutDiag(string& line) : Output(line)
{
  if (myid) return;
				// Defaults
  RMIN = 1.0e-3;
  RMAX = 10.0;
  THETA = PHI = 1.0e-10;
  NUM = 100;

  list<Component*>::iterator cc;
  Component *c;
  
  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    c = *cc;

    if (c->force->geometry == PotAccel::sphere || 
	c->force->geometry == PotAccel::cylinder) 
      {
	lcomp.push_back(&(*c));
      }
  }
  
  names.push_back("Rho");
  names.push_back("Pot");
  names.push_back("d(Pot)/dr)");
  names.push_back("d(Pot)/d cos(theta)");
  names.push_back("d(Pot)/d phi");

  initialize();
}

void OutDiag::initialize()
{
  string tmp;
				// Get file name
  if (!Output::get_value(string("filename"), filename)) {
    filename.erase();
    filename = "ORBDIAG." + runtag;
  }

  if (Output::get_value(string("nint"), tmp))
    nint = atoi(tmp.c_str());
  else
    nint = 1;

  if (Output::get_value(string("RMIN"), tmp))
    RMIN = atof(tmp.c_str());

  if (Output::get_value(string("RMAX"), tmp))
    RMAX = atof(tmp.c_str());

  if (Output::get_value(string("THETA"), tmp))
    THETA = atof(tmp.c_str());

  if (Output::get_value(string("PHI"), tmp))
    PHI = atof(tmp.c_str());

  if (Output::get_value(string("NUM"), tmp))
    NUM = atoi(tmp.c_str());
}


void OutDiag::header(ostream& out)
{
  list<Component*>::iterator c;
  int ncur = 0;

  out << "# " << ++ncur << ": " << "Radius\n";

  for (c=lcomp.begin(); c != lcomp.end(); c++) {

    out << "# [" << (*c)->id << "]\n";
    
    for (int i=0; i<5; i++) 
      out << "# " << setw(3) << ++ncur << ": " << names[i] << endl;
  }

  out << "#" << endl;
}


void OutDiag::Run(int n, bool last)
{
  if (myid) return;
  if (n % nint && !last) return;

  double r, dr, dens;
  double potl, potr, pott, potp;
    
  ostringstream outs;
  outs << filename.c_str() << "." << n << '\0';
  ofstream out(outs.str().c_str());
  if (!out) return;

  out.setf(ios::scientific);
  out.precision(6);

  header(out);

  // Determine potential and acceleration

  list<Component*>::iterator c;

  dr = (RMAX - RMIN)/(double)NUM;
  for (int i=0; i<=NUM; i++) {
    r = RMIN + dr*i;

    out << setw(15) << r;
    
    for (c=lcomp.begin(); c != lcomp.end(); c++) {
      AxisymmetricBasis * q = (AxisymmetricBasis *)((*c)->force);
      q->determine_fields_at_point_sph(r, THETA, PHI,
				       &dens, &potl, &potr, &pott, &potp);
      out << setw(15) << dens
	  << setw(15) << potl
	  << setw(15) << potr
	  << setw(15) << pott
	  << setw(15) << potp;
    }
    out << endl;
  }

}

