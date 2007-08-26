/* 
  Multi-timestepping support routines
*/

#include "expand.h"
#include <sstream>

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void sync_eval_multistep(void)
{
  comp.multistep_reset();
				// This forces interpolation to evaluate to 
				// the last set of computed coefficients
  mstep = Mstep;
  for (int M=0; M<=multistep; M++) {
    stepL[M] = 0;
    stepN[M] = Mstep;
  }
}

void adjust_multistep_level(bool all)
{
  if (!multistep) return;

  // Examine all time steps at or below this level and compute timestep
  // criterion and adjust level if necessary

  list<Component*>::iterator cc;
  Component *c;
  double dt;
  int npart, lev, offgrid;
  vector<unsigned> off1, off;

  if (all) mstep=Mstep;

  //
  // Run through particles in each component
  //
  for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
    c = *cc;

    c->force->multistep_update_begin();

    npart = c->Number();
    offgrid = 0;

    for (int n=0; n<npart; n++) {

      if (all || mactive[mstep-1][c->Part(n)->level]) {

	double rtot = 0.0;
	double vtot = 0.0;
	double atot = 0.0;

	for (int i=0; i<c->dim; i++) {
	  rtot += c->Part(n)->pos[i]*c->Part(n)->pos[i];
	  vtot += c->Part(n)->vel[i]*c->Part(n)->vel[i];
	  atot += c->Part(n)->acc[i]*c->Part(n)->acc[i];
	}
	rtot = sqrt(rtot);
	vtot = sqrt(vtot) + 1.0e-18;
	atot = sqrt(atot) + 1.0e-18;

	dt = dynfrac*min<double>(rtot/vtot, sqrt(rtot/atot));
	
	if (dt>dtime) lev = 0;
	else lev = (int)floor(log(dtime/dt)/log(2.0));

	if (lev>multistep) {
	  lev = multistep;
	  offgrid++;
	}

	if ( lev != c->Part(n)->level && (all || mactive[mstep-1][lev]) ) {
	  c->force->multistep_update(c->Part(n)->level, lev, c, n);
	  c->Part(n)->level = lev;
	}
      }
    }

    c->force->multistep_update_finish();
    if (VERBOSE>0 && (this_step % 100 == 0)) {
      off1.push_back(offgrid);
      off.push_back(0);
    }

  }

  if (VERBOSE>0) {
    MPI_Reduce(&off1[0], &off[0], off1.size(), MPI_UNSIGNED, MPI_SUM, 0, 
	       MPI_COMM_WORLD);
    
    if (myid==0) {

      unsigned sum=0;
      for (unsigned i=0; i<off.size(); i++) sum += off[i];
      
      if (sum) {
	cout << setw(70) << setfill('-') << '-' << endl
	     << setw(70) << left << "--- Multistepping overrun" << endl
	     << setw(70) << setfill('-') << '-' << endl << setfill(' ') << right;

	vector<unsigned>::iterator it = off.begin();

	for (cc=comp.components.begin(); cc != comp.components.end(); cc++) {
	  if (*it > 0) {
	    ostringstream sout;
	    sout << "Component <" << (*cc)->name << ">";
	    cout << setw(20) << sout.str() << ": "
		 << *it << "/" << (*cc)->nbodies_tot << endl;
	  }
	  it++;
	}
	cout << setw(70) << setfill('-') << '-' << endl << setfill(' ');
      }
    }

  }

}

void initialize_multistep()
{
  Mstep = 1 << multistep;	// Number of substeps

				// Number of substeps at each level
  mintvl.push_back(Mstep);
  for (int n=1; n<=multistep; n++) mintvl.push_back(mintvl.back()/2);

				// Set up the step-level bool array
  for (int ms=0; ms<Mstep; ms++)
    mactive.push_back(vector<bool>(multistep+1, false));

				// Find and save the active levels at each step
  for (int ms=1; ms<=Mstep; ms++) {
    for (int mlevel=0; mlevel<=multistep; mlevel++) {
      if ( (ms % (1 << (multistep-mlevel))) == 0) 
	mactive[ms-1][mlevel] = true;
    }
  }

  stepL = vector<int>(multistep+1);
  stepN = vector<int>(multistep+1);

  mfirst = vector<int>(Mstep);
				// Lowest active level at each step
  for (int ms=1; ms<=Mstep; ms++) {
    for (int mlevel=0; mlevel<=multistep; mlevel++) {
      if (mactive[ms-1][mlevel]) {
	mfirst[ms-1] = mlevel;
	break;
      }
    }
  }

  if (VERBOSE>0 && myid==0 && multistep) {
    cout << setw(70) << setfill('-') << '-' << endl
	 << setw(70) << left << "--- Multistep level control structure " << endl
	 << setw(70) << setfill('-') << '-' << endl << setfill(' ') << right;
    cout << "Mstep=" << Mstep << "  levels=" << multistep << endl << endl;
    cout << setw(4) << "Step" << "/" << setw(3)<< "1st" << "  ";
    for (int mlevel=0; mlevel<=multistep; mlevel++) cout << setw(3) << mlevel;
    cout << endl;
    for (int ms=1; ms<=Mstep; ms++) {
      cout << setw(4) << ms << "/";
      cout << setw(3) << mfirst[ms-1] << ": ";
      for (int mlevel=0; mlevel<=multistep; mlevel++)
	if (mactive[ms-1][mlevel]) cout << setw(3) << "+";
	else                          cout << setw(3) << "-";
      cout << endl;
    }
    cout << setw(70) << setfill('-') << '-' << endl << setfill(' ');
  }

}
