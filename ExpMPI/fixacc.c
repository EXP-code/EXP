/*
  Correct accelerations so that the center of
  mass remains fixed at the origin

  bug fix: frozen particles were being added to mass,
  making subsequent compensating acceleration too low
  when significant mass loss occurred.
  KL 6/29/92

*/

#include "expand.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void fix_acceleration(void)
{
  int i;
  double axcm,aycm,azcm,mtot;
  double axcm1,aycm1,azcm1,mtot1;

  if (myid == 0) return;

  axcm = aycm = azcm = mtot = 0.0;
  axcm1 = aycm1 = azcm1 = mtot1 = 0.0;

  for (i=1; i<=nbodies; i++) {

    if (freeze_particle(i)) continue;

    mtot1=mtot1+mass[i];
    axcm1=axcm1+mass[i]*ax[i];
    aycm1=aycm1+mass[i]*ay[i];
    azcm1=azcm1+mass[i]*az[i];
  }


  MPI_Allreduce(&mtot1, &mtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_SLAVE);
  MPI_Allreduce(&axcm1, &axcm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_SLAVE);
  MPI_Allreduce(&aycm1, &aycm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_SLAVE);
  MPI_Allreduce(&azcm1, &azcm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_SLAVE);

  if (mtot>0.0) {
    axcm=axcm/mtot;
    aycm=aycm/mtot;
    azcm=azcm/mtot;
  }

  for (i=1; i<=nbodies; i++) {
    if (freeze_particle(i)) continue;
    if (!slab) {
      ax[i]=ax[i]-axcm;
      ay[i]=ay[i]-aycm;
    }
    az[i]=az[i]-azcm;
  }

  if (slab) azcm_slab = azcm;

}





/* 	
	insist that the cm remains fixed at origin, by subtracting
cm position from all coordinates.

KL - 6/29/92

*/


void fix_positions_by_component(void)
{
  int i;
  double mcom1, mcom11;
  double mcom2, mcom21;
  double mtot, mtot1;
  double com11[3], com21[3], vcom[3], vcom1[3];
  MPI_Status status;

  if (myid > 0) {

    for (i=0; i<3; i++)
      com1[i] = com11[i] = com2[i] = com21[i] = vcom[i] = vcom1[i] = 0.0;

    mcom1 = mcom11 = mcom2 = mcom21 = 0.0;
    mtot1 = mtot = 0.0;

    for (i=1; i<=nbodies; i++) {

      if (freeze_particle(i)) continue;
    
      if (component[i] == 1) {
	
	mcom11 += mass[i];
	com11[0] += mass[i]*x[i];
	com11[1] += mass[i]*y[i];
	com11[2] += mass[i]*z[i];

      }

      if (component[i] == 2) {

	mcom21 += mass[i];
	com21[0] += mass[i]*x[i];
	com21[1] += mass[i]*y[i];
	com21[2] += mass[i]*z[i];

      }

      if (fixvel) {

	mtot1 += mass[i];
	vcom1[0] += mass[i]*vx[i];
	vcom1[1] += mass[i]*vy[i];
	vcom1[2] += mass[i]*vz[i];
      
      }
    
    }

    MPI_Allreduce(&mcom11, &mcom1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_SLAVE);
    MPI_Allreduce(&mcom21, &mcom2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_SLAVE);
    MPI_Allreduce(com11, com1, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_SLAVE);
    MPI_Allreduce(com21, com2, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_SLAVE);

    if (mcom1 > 0.0) {
      for (i=0; i<3; i++) com1[i] /= mcom1;
    }
    if (mcom2 > 0.0) {
      for (i=0; i<3; i++) com2[i] /= mcom2;
    }

    if (fixvel) {
      MPI_Allreduce(&mtot1, &mtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_SLAVE);
      MPI_Allreduce(vcom1, vcom, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_SLAVE);
      if (mtot > 0.0) {
	for (i=0; i<3; i++) vcom[i] /= mtot;
      }

      for (i=1; i<=nbodies; i++) {

	if (freeze_particle(i)) continue;

	vx[i] -= vcom[0];
	vy[i] -= vcom[1];
	vz[i] -= vcom[2];
	
      }
    }

				/* First slave return com to master */
    if (myid==1) {
      MPI_Send(com1, 3, MPI_DOUBLE, 0, 40, MPI_COMM_WORLD);
      MPI_Send(com2, 3, MPI_DOUBLE, 0, 40, MPI_COMM_WORLD);
    }

  }
  else {
    MPI_Recv(com1, 3, MPI_DOUBLE, MPI_ANY_SOURCE, 40, MPI_COMM_WORLD,
	     &status);
    MPI_Recv(com2, 3, MPI_DOUBLE, MPI_ANY_SOURCE, 40, MPI_COMM_WORLD,
	     &status);

  }

      
}


void fix_positions(void)
{
  int i;
  double xcm, ycm, zcm, mtot;
  double vxcm, vycm, vzcm;

  double xcm1, ycm1, zcm1, mtot1;
  double vxcm1, vycm1, vzcm1;

  if (myid == 0) return;

  xcm = ycm = zcm = mtot = 0.0;
  vxcm = vycm = vzcm = 0.0;

  xcm1 = ycm1 = zcm1 = mtot1 = 0.0;
  vxcm1 = vycm1 = vzcm1 = 0.0;


  for (i=1; i<=nbodies; i++) {
    if (freeze_particle(i)) continue;
    mtot1=mtot1+mass[i];
    xcm1=xcm1+mass[i]*x[i];
    ycm1=ycm1+mass[i]*y[i];
    zcm1=zcm1+mass[i]*z[i];
    vxcm1=vxcm1+mass[i]*vx[i];
    vycm1=vycm1+mass[i]*vy[i];
    vzcm1=vzcm1+mass[i]*vz[i];
  }


  MPI_Allreduce(&mtot1, &mtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_SLAVE);

  MPI_Allreduce(&xcm1, &xcm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_SLAVE);
  MPI_Allreduce(&ycm1, &ycm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_SLAVE);
  MPI_Allreduce(&zcm1, &zcm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_SLAVE);
  
  MPI_Allreduce(&vxcm1, &vxcm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_SLAVE);
  MPI_Allreduce(&vycm1, &vycm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_SLAVE);
  MPI_Allreduce(&vzcm1, &vzcm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_SLAVE);

  if (mtot > 0.0) {
    vxcm /= mtot;
    vycm /= mtot;
    vzcm /= mtot;
    xcm /= mtot;
    ycm /= mtot;
    zcm /= mtot;
  }

  for (i=1; i<=nbodies; i++) {

    if (freeze_particle(i)) continue;

    if (!slab) {
      x[i] -= xcm;
      y[i] -= ycm;
      vx[i] -= vxcm;
      vy[i] -= vycm;
    }
    z[i] -= zcm;
    vz[i] -= vzcm;
    
  }
  
  if (slab) {
    zcm_slab = zcm;
    vzcm_slab = vzcm;
  }

}



void set_global_com(void)
{
  int i;
  double m1, m2, m0;
  double ps0[2], ps1[6], ps2[6];

  m0 = m1 = m2 = 0.0;
  for (i=0; i<6; i++) ps1[i] = ps2[i] = 0.0;
  ps0[0] = ps0[1] = 0.0;

  for (i=1; i<=nbodies; i++) {

    if (component[i] == -1) {
	
      m0 += mass[i];

      ps0[0] += mass[i]*z[i];
      ps0[1] += mass[i]*vz[i];

    }

    if (component[i] == 1) {
	
      m1 += mass[i];

      ps1[0] += mass[i]*x[i];
      ps1[1] += mass[i]*y[i];
      ps1[2] += mass[i]*z[i];
      ps1[3] += mass[i]*vx[i];
      ps1[4] += mass[i]*vy[i];
      ps1[5] += mass[i]*vz[i];

    }

    if (component[i] == 2) {

      m2 += mass[i];

      ps2[0] += mass[i]*x[i];
      ps2[1] += mass[i]*y[i];
      ps2[2] += mass[i]*z[i];
      ps2[3] += mass[i]*vx[i];
      ps2[4] += mass[i]*vy[i];
      ps2[5] += mass[i]*vz[i];

    }

  }

  for (i=0; i<6; i++) {
    if (m1>0.0) ps1[i] /= m1;
    if (m2>0.0) ps2[i] /= m2;
  }
  
  if (m0>0.0) {
    ps0[0] /= m0;
    ps0[1] /= m0;
  }

  for (i=1; i<=nbodies; i++) {

    if (component[i] == -1) {
      
      if (zerocom) z[i] -= ps0[0];
      if (zerovel) vz[i] -= ps0[1];

    }
	
    if (component[i] == 1) {
      
      if (zerocom) {
	x[i] -= ps1[0];
	y[i] -= ps1[1];
	z[i] -= ps1[2];
      }
      if (zerovel) {
	vx[i] -= ps1[3];
	vy[i] -= ps1[4];
	vz[i] -= ps1[5];
      }

    }
	
    if (component[i] == 2) {
      
      if (zerocom) {
	x[i] -= ps2[0];
	y[i] -= ps2[1];
	z[i] -= ps2[2];
      }
      if (zerovel) {
	vx[i] -= ps2[3];
	vy[i] -= ps2[4];
	vz[i] -= ps2[5];
      }
    }
	
  }

}

