#ifndef _gadget_H
#define _gadget_H

struct gadget_header
{
  //! number of particles of each type in this file 
  int npart[6];	 

  //! mass of particles of each type. If 0, then the masses are
  //! explicitly stored in the mass-block of the snapshot file,
  //! otherwise they are omitted
  double mass[6];  

  //! time of snapshot file
  double time;

  //! redshift of snapshot file */
  double redshift;

  //! flags whether the simulation was including star formation */
  int flag_sfr;

  //! flags whether feedback was included (obsolete)
  int flag_feedback;

  //!  total number of particles of each type in this snapshot. This
  //!  can be different from npart if one is dealing with a multi-file
  //!  snapshot.
  unsigned int npartTotal[6];
				     
  //! flags whether cooling was included
  int flag_cooling;

  //! number of files in multi-file snapshot
  int num_files;

  //! box-size of simulation in case periodic boundaries were used
  double BoxSize;

  //! matter density in units of critical density
  double Omega0;

  //! cosmological constant parameter
  double OmegaLambda;

  //! Hubble parameter in units of 100 km/sec/Mpc
  double HubbleParam;          

  //! flags whether the file contains formation times of star particles
  int flag_stellarage;

  //! flags whether the file contains metallicity values for gas and star particles
  int flag_metals;

  //! High word of the total number of particles of each type
  unsigned int npartTotalHighWord[6];   

  char fill[64];
};

#endif
