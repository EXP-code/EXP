
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <memory>
#include <cmath>

#include <time.h>
#include "numerical.H"

#include "phase.H"

using namespace std;

/*
  ==============================================================================
  
  LOG I/O:
  record the vital statistics of an ensemble: conserved
  quantities, virial, system and real time.
  
  
  ==============================================================================
*/

void Ensemble::write_log(ostream& logfile)
{
  double pe, ke;
  Eigen::Vector3d xcm, vcm, J;
  time_t realtime;
  
  time(&realtime);
  logfile << "\n\t --------------------------------------- \n";
  logfile << "Log: " << ctime(&realtime);
  logfile << "time: " << t << '\n';
  xcm = CM_Position();
  vcm = CM_Velocity();
  J = total_Angular_Momentum();
  logfile << "x_CM: " 
	  << xcm[1] << " " << xcm[2] << " " << xcm[3] << '\n';
  logfile << "v_CM: "
	  << vcm[1] << " " << vcm[2] << " " << vcm[3] << '\n';
  logfile << "Jtot: "
	  << J[1] << " " << J[2] << " " << J[3] << '\n';
  pe = total_Potential();
  ke = total_Kinetic();
  logfile << "E = " << pe+ke
	  << "T = " << ke
	  << "V = " << pe
	  << "2T + V = " << 2.0*ke + pe
	  << '\n';
}








/*
  =============================================================================
  
  
  "SNAPSHOT" I/O: record system time, 
  individual masses and phase space for each star.
  Individual times not recorded.
  
  
  =============================================================================
*/



void Ensemble::write_snapshot(ostream& fp)
{
  fp << Nstars << " " << t << '\n';

  for (int i=0; i<Nstars; i++)
    {
      double p = 0.0;
      if (Phase::potential_defined) 
	p = (*Phase::Potential)(t, stars[i].x);

      fp << stars[i].m << " " 
	 << stars[i].x[1] << " " 
	 << stars[i].x[2] << " "
	 << stars[i].x[3] << " " 
	 << stars[i].v[1] << " " 
	 << stars[i].v[2] << " " 
	 << stars[i].v[3] << '\n';
    }
}



void Ensemble::read_snapshot(istream& fp)
{
  char buffer[512];
  
  fp.get(buffer, 512);
  {
    istringstream sbuf(buffer);
    sbuf >> Nstars;
    sbuf >> t;
  }
  
  stars.resize(Nstars);
  
  for (int i=0; i<Nstars; i++)
    {
      fp.get(buffer, 512);
      std::istringstream sbuf(buffer);
      sbuf >> stars[i].m;
      sbuf >> stars[i].x[1];
      sbuf >> stars[i].x[2];
      sbuf >> stars[i].x[3];
      sbuf >> stars[i].v[1];
      sbuf >> stars[i].v[2];
      sbuf >> stars[i].v[3];
      sbuf >> stars[i].t;
      stars[i].work = 0.0;
    }
}



/* Variant to read MDW's and DFC's tides phase space */

void Ensemble::read_tidesfile(istream& fp)
{
  char tidbuf[512];
  
  // Determine the number of stars in the dump
  // Admittedly, this is wasteful . . . 
  
  Nstars = 0;
  double tmp;		// Translation of c construct; may have
  // bugs . . .
  while (fp.get(tidbuf, 512)) {
    std::istringstream sbuf(tidbuf);
    int icnt = 0;
    while ((sbuf >> tmp)) icnt++;
    if (icnt == 9) Nstars++;
  }
  fp.seekg(ios::beg);
  
  stars.resize(Nstars);
  
  double m = 1.0/Nstars;
  for (int i=0; i<Nstars; i++)
    {
      fp.get(tidbuf, 512);
      istringstream sbuf(tidbuf);
      stars[i].m = m;
      sbuf >> stars[i].x[1];
      sbuf >> stars[i].x[2];
      sbuf >> stars[i].x[3];
      sbuf >> stars[i].v[1];
      sbuf >> stars[i].v[2];
      sbuf >> stars[i].v[3];
      sbuf >> stars[i].t;
      stars[i].work = 0.0;
    }
}




/* 
   =============================================================================
   
   
   
   "ORBIT" I/O: record phase space and an individual time for
   each star. No masses recorded. The individual time could be, 
   e.g., used to record the orbital period for that star.
   
   
   
   =============================================================================
*/


void Ensemble::write_orbits(ostream& fp)
{
  fp << Nstars << '\n';

  for (int i=0; i<Nstars; i++)
    {
      fp << stars[i].t << " "
	 << stars[i].x[1] << " " 
	 << stars[i].x[2] << " " 
	 << stars[i].x[3] << " " 
	 <<	stars[i].v[1] << " " 
	 << stars[i].v[2] << " " 
	 << stars[i].v[3] << '\n';
    }
}




void Ensemble::read_orbits(istream& fp)
{
  char tmp[128];
  
  fp >> tmp;
  fp >> tmp;
  fp >> tmp;
  fp >> tmp;
  fp >> Nstars;
  
  stars.resize(Nstars);
  
  for (int i=0; i<Nstars; i++)
    {
      fp >> stars[i].t;
      fp >> stars[i].x[1];
      fp >> stars[i].x[2];		  
      fp >> stars[i].x[3];
      fp >> stars[i].v[1];
      fp >> stars[i].v[2];
      fp >> stars[i].v[3];
      stars[i].work = 0.0;
    }
}



/*
  =============================================================================
  
  
  BINARY I/O:	Record complete mass, time, and phase information
  in binary format.
  
  
  
  =============================================================================
*/


void Ensemble::binwrite(ostream& fp)
{
  fp.write((char *)&Nstars, sizeof(int));
  for (int i=0; i<Nstars; i++) stars[i].binwrite(fp);
}

void Ensemble::binread(istream& fp)
{
  int n;
  fp.read((char *)&n, sizeof(int));
  
  setsize(n);
  
  for (int i=0; i<Nstars; i++)
    {
      stars[i].binread(fp);
    }
}

