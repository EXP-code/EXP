#ifndef _tipsy_H
#define _tipsy_H

#include <rpc/types.h>
#include <rpc/xdr.h>

// Forward declare PR::Tipsy so it can friend TipsyFile
//
namespace PR {
  class Tipsy;
};

namespace TipsyReader
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
  public:
    int ID() const {
      union id {Real v; int i;} u;
      u.v = phi; return u.i;
    }
  } ;

  struct dark_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real eps;
    Real phi;
  public:
    int ID() const {
      union id {Real v; int i;} u;
      u.v = phi; return u.i;
    }
  } ;

  struct star_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real metals ;
    Real tform  ;
    Real eps ;
    Real phi ;
  public:
    int ID() const {
      union id {Real v; int i;} u;
      u.v = phi; return u.i;
    }
  } ;
  
  struct dump {
    double time ;
    int nbodies ;
    int ndim  ;
    int nsph  ;
    int ndark ;
    int nstar ;
    int version ;
  } ;
  
  class TipsyFile
  {
  protected:

    virtual void read_gas () = 0;
    virtual void read_dark() = 0;
    virtual void read_star() = 0;


  public:

    std::vector<gas_particle>  gas_particles;
    std::vector<dark_particle> dark_particles;
    std::vector<star_particle> star_particles;
    dump header;

    virtual int readParticles() = 0;

    virtual ~TipsyFile() {}
  };

  class TipsyXDR : public TipsyFile
  {
  private:

    FILE* input;
    XDR xdrs;

    int xdr_header()
    {
      if (xdr_double(&xdrs, &header.time) != TRUE) return 0;
      if (xdr_int(&xdrs, &header.nbodies) != TRUE) return 0;
      if (xdr_int(&xdrs, &header.ndim)    != TRUE) return 0;
      if (xdr_int(&xdrs, &header.nsph)    != TRUE) return 0;
      if (xdr_int(&xdrs, &header.ndark)   != TRUE) return 0;
      if (xdr_int(&xdrs, &header.nstar)   != TRUE) return 0;
      if (xdr_int(&xdrs, &header.version) != TRUE) return 0;
      
      return 1;
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

      return header.nbodies;
    }

    int xdr_read()
    {
      int N=0;
      
      if (header.nsph != 0)  {
	gas_particles.resize(header.nsph);
	read_gas();
	N++;
      }
    
      if (header.ndark != 0) {
	dark_particles.resize(header.ndark);
	read_dark();
	N++;
      }

      if (header.nstar != 0) {
	star_particles.resize(header.nstar);
	read_star();
	N++;
      }
	
      return N;
    }

    void read_gas()
    {
      if (sizeof(Real) == sizeof(float)) {
	xdr_vector(&xdrs, (char *) &gas_particles[0],
		   header.nsph*(sizeof(gas_particle)/sizeof(Real)),
		   sizeof(Real), (xdrproc_t) xdr_float);
      }
    }  
    
    void read_dark()
    {
      if (sizeof(Real) == sizeof(float)) {
	xdr_vector(&xdrs, (char *) &dark_particles[0],
		   header.ndark*(sizeof(dark_particle)/sizeof(Real)),
		   sizeof(Real), (xdrproc_t) xdr_float);
      }
    }  
  
    void read_star()
    {
      if (sizeof(Real) == sizeof(float)) {
	xdr_vector(&xdrs, (char *) &star_particles[0],
		   header.nstar*(sizeof(star_particle)/sizeof(Real)),
		   sizeof(Real), (xdrproc_t) xdr_float);
      }
    }  
    
    friend class PR::Tipsy;

  public:

    TipsyXDR(const std::string filename)
    {
      // Attempt to read tipsy file
      int ret = xdr_init(filename);
      if (ret==0) {
	std::cerr << "TipsyFile: no bodies found??" << std::endl;
	exit(-1);
      }
    }

    int readParticles() { return xdr_read(); }

    ~TipsyXDR()
    {
      xdr_destroy(&xdrs);
      fclose(input);
    }
    
  };

  class TipsyNative : public TipsyFile
  {
  private:

    //! The input file stream
    std::ifstream input;

    //! Convert phi to IDs
    bool convert;

    //! Read the header
    int native_header()
    {
      try {
	input.read((char *)&header, sizeof(header));
      } catch (std::exception& e) {
	std::cerr << "TipsyFile native_header error: " << e.what() << std::endl;
	return 0;
      }
      return 1;
    }
    
    int native_init(const std::string& filename)
    {
      try {
	input.open(filename);
	input.exceptions(ifstream::eofbit | ifstream::failbit | ifstream::badbit);
      } catch (std::exception& e) {
	std::cerr << "TipsyFile native error opening file <" << filename << ">, "
		  << e.what() << std::endl;
	exit(-1);
      }
	
      if (native_header() != 1) {
	std::cerr << "TipsyFile native: could not read a valid header" << std::endl;
	exit(-1);
      }

      return header.nbodies;
    }

    int native_read()
    {
      int N=0;
      
      if (header.nsph != 0)  {
	gas_particles.resize(header.nsph);
	read_gas();
	N++;
      }
    
      if (header.ndark != 0) {
	dark_particles.resize(header.ndark);
	read_dark();
	N++;
      }

      if (header.nstar != 0) {
	star_particles.resize(header.nstar);
	read_star();
	N++;
      }
	
      return N;
    }

    void read_gas()
    {
      try {
	input.read((char *) &gas_particles[0],
		   header.nsph*(sizeof(gas_particle)/sizeof(Real)));
      }
      catch (std::exception& e) {
	std::cerr << "TipsyFile native error reading sph particles: "
		  << e.what() << std::endl;
	exit(-1);
      }
    }
    
    void read_dark()
    {
      try {
	input.read((char *) &dark_particles[0],
		   header.ndark*(sizeof(dark_particle)/sizeof(Real)));
      }
      catch (std::exception& e) {
	std::cerr << "TipsyFile native error reading dark particles: "
		  << e.what() << std::endl;
	exit(-1);
      }
    }  
  
    void read_star()
    {
      try {
	input.read((char *) &star_particles[0],
		   header.nstar*(sizeof(star_particle)/sizeof(Real)));
      }
      catch (std::exception& e) {
	std::cerr << "TipsyFile native error reading star particles: "
		  << e.what() << std::endl;
	exit(-1);
      }
    }  
    
    friend class PR::Tipsy;

  public:

    TipsyNative(const std::string filename, bool ids=false)
    {
      convert = ids;

      // Attempt to read tipsy file
      int ret = native_init(filename);
      if (ret==0) {
	std::cerr << "TipsyFile: no bodies found??" << std::endl;
	exit(-1);
      }
    }

    int readParticles() { return native_read(); }

    ~TipsyNative()
    {
      // Nothing
    }
    
  };

};

#endif
