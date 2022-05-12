#ifndef _tipsy_H
#define _tipsy_H

#include <rpc/types.h>
#include <rpc/xdr.h>

#include <ParticleReader.H>

namespace TipsyXDR
{
  const int MAXDIM = 3;

  using Real = float;

  struct gas_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real rho;
    Real temp;
    Real hsmooth;
    Real metals;
    Real phi;
  } ;

  struct dark_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real eps;
    Real phi;
  } ;

  struct star_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real metals ;
    Real tform  ;
    Real eps ;
    Real phi ;
  } ;
  
  struct dump {
    double time ;
    int nbodies ;
    int ndim  ;
    int nsph  ;
    int ndark ;
    int nstar ;
  } ;
  
  class TipsyFile
  {
  public:

    std::vector<gas_particle>  gas_particles;
    std::vector<dark_particle> dark_particles;
    std::vector<star_particle> star_particles;
    dump header ;

  private:

    FILE* input;
    XDR xdrs;

    int xdr_header()
    {
      int pad;
      
      if (xdr_double(&xdrs, &header.time) != TRUE) return 0;
      if (xdr_int(&xdrs, &header.nbodies) != TRUE) return 0;
      if (xdr_int(&xdrs, &header.ndim)    != TRUE) return 0;
      if (xdr_int(&xdrs, &header.nsph)    != TRUE) return 0;
      if (xdr_int(&xdrs, &header.ndark)   != TRUE) return 0;
      if (xdr_int(&xdrs, &header.nstar)   != TRUE) return 0;
      if (xdr_int(&xdrs, &pad)            != TRUE) return 0;
      
      return 1;
    }
    
    void xdr_gas()
    {
      if (sizeof(Real) == sizeof(float)) {
	xdr_vector(&xdrs, (char *) &gas_particles[0],
		   header.nsph*(sizeof(gas_particle)/sizeof(Real)),
		   sizeof(Real), (xdrproc_t) xdr_float);
      }
    }  
    
    void xdr_dark()
    {
      if (sizeof(Real) == sizeof(float)) {
	xdr_vector(&xdrs, (char *) &dark_particles[0],
		   header.ndark*(sizeof(dark_particle)/sizeof(Real)),
		   sizeof(Real), (xdrproc_t) xdr_float);
      }
    }  
  
    void xdr_star()
    {
      if (sizeof(Real) == sizeof(float)) {
	xdr_vector(&xdrs, (char *) &star_particles[0],
		   header.nstar*(sizeof(star_particle)/sizeof(Real)),
		   sizeof(Real), (xdrproc_t) xdr_float);
      }
    }  
    
    int xdr_init(const std::string& filename)
    {
      input = fopen(filename.c_str(), "r");
      if (fopen==NULL) {
	std::cerr << "TipsyFile: error opening file <" << filename << ">"
		  << std::endl;
	exit(-1);
      }
	
      xdrstdio_create(&xdrs, input, XDR_DECODE);
      
      if (xdr_header() != 1) {
	std::cerr << "TipsyFile: could not read a valid header" << std::endl;
	exit(-1);
      }
    
      int N=0;
    
      if (header.nsph != 0)  gas_particles.resize(header.nsph);
      N++;
    
      if (header.ndark != 0) dark_particles.resize(header.ndark);
      N++;

      if (header.nstar != 0) star_particles.resize(header.nstar);
      N++;
      
      return N;
    }


  public:

    TipsyFile(const std::string filename) {
      int ret = xdr_init(filename);
      if (ret==0) {
	std::cerr << "TipsyFile: no components found??" << std::endl;
	exit(-1);
      }
    }

    ~TipsyFile()
    {
      xdr_destroy(&xdrs);
      fclose(input);
    }

  };

};



#endif
