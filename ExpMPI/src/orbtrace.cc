#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>

#include <expand.h>

static char rcsid[] = "$Id$";

const int norb = 5;
static int orblist[norb] = {100, 150, 200, 250, 300};

extern "C" void orb_trace(void)
{
  static bool firsttime = true;
  ofstream out(orbtracefile, ios::out | ios::app);
  if (!out) {
    if (firsttime) {
      cerr << "orb_trace: can't open file <" << orbtracefile << ">\n";
      firsttime = false;
    }
    else return;
  }

  out << setw(15) << tnow;
  for (int i=0; i<norb; i++) {
    out 
      << setw(15) << x[orblist[i]]
      << setw(15) << y[orblist[i]]
      << setw(15) << z[orblist[i]];
  }
  out << endl;
  
}
