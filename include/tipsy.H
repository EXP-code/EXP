#ifndef _tipsy_H
#define _tipsy_H

#include <map>

// For MPI awareness
//
#include <localmpi.H>


// Forward declare PR::Tipsy so it can friend TipsyFile
//
namespace PR {
  class Tipsy;
};

// Uncommenting the define below is an old fix for compilers that
// promoted the 28 byte Tipsy header to 32 bytes to put it on a word
// boundary.  I don't this that this is an issue in modern compilers.
//
// #define TIPSY_32BYTE_PAD

namespace TipsyReader
{

  const int MAXDIM = 3;

  using Real = float;

  struct gas_particle
  {
    //@{
    //! Particle data
    Real mass   ;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real rho    ;
    Real temp   ;
    Real hsmooth;
    Real metals ;
    Real phi    ;
    //@}
  } ;

  struct dark_particle
  {
    //@{
    //! Particle data
    Real mass   ;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real eps    ;
    Real phi    ;
    //@}
    
    //! Convert phi to 32-bit index for Bonsai v1
    int ID() const {
      union id {Real v; int i;} u;
      u.v = phi; return u.i;
    }

    //! Convert phi to 64-bit index for Bonsai v2
    unsigned long ID2() const {
      union id {Real v[2]; uint64_t i;} u;
      u.v[0] = eps; u.v[1] = phi; return u.i;
    }

  } ;

  struct star_particle
  {
    //@{
    //! Particle data
    Real mass   ;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real metals ;
    Real tform  ;
    Real eps    ;
    Real phi    ;
    //@}

    //! Convert phi to 32-bit index for Bonsai v1
    int ID() const {
      union id {Real v; int i;} u;
      u.v = phi; return u.i;
    }

    //! Convert exp, phi to 64bit index for Bonsai v2
    unsigned long ID2() const {
      union id {Real v[2]; uint64_t i;} u;
      u.v[0] = eps; u.v[1] = phi; return u.i;
    }
  } ;
  
  struct Header
  {
    double time ;
    int nbodies ;
    int ndim    ;
    int nsph    ;
    int ndark   ;
    int nstar   ;
#ifdef TIPSY_32BYTE_PAD
    int version ;
#endif
  } ;
  
  //! Used to identify particle type for indexing
  enum class Ptype { gas, dark, star };

  //! Base class for Tipsy file reading
  class TipsyFile
  {
  protected:

    virtual void read_gas () = 0;
    virtual void read_dark() = 0;
    virtual void read_star() = 0;

    //! Offset for indexing particles in the file
    std::map<Ptype, unsigned long> index_offset_map
      { {Ptype::gas,  0}, {Ptype::dark, 0}, {Ptype::star, 0} };

  public:

    std::vector<gas_particle>  gas_particles;
    std::vector<dark_particle> dark_particles;
    std::vector<star_particle> star_particles;
    Header header;

    virtual int readParticles() = 0;

    virtual ~TipsyFile() {}

    //! Get the index offset for a given particle type
    unsigned long getIndexOffset(const Ptype ptype)
    {
      return index_offset_map[ptype];
    }
  };

#ifdef HAVE_XDR

#include <rpc/types.h>
#include <rpc/xdr.h>

  //! Tipsy XDR file reader
  class TipsyXDR : public TipsyFile
  {
  private:

    //! The input file stream
    FILE* input;

    //! The XDR stream
    XDR xdrs;

    //! Cache the XDR stream position at beginning of particle group
    unsigned long curpos = 0;

    int xdr_header()
    {
      if (xdr_double(&xdrs, &header.time) != TRUE) return 0;
      if (xdr_int(&xdrs, &header.nbodies) != TRUE) return 0;
      if (xdr_int(&xdrs, &header.ndim)    != TRUE) return 0;
      if (xdr_int(&xdrs, &header.nsph)    != TRUE) return 0;
      if (xdr_int(&xdrs, &header.ndark)   != TRUE) return 0;
      if (xdr_int(&xdrs, &header.nstar)   != TRUE) return 0;
#ifdef TIPSY_32BYTE_PAD
      if (xdr_int(&xdrs, &header.version) != TRUE) return 0;
#endif
      return 1;
    }
    
    int xdr_init(const std::string& filename)
    {
      input = fopen(filename.c_str(), "r");
      if (fopen==NULL) {
	std::ostringstream sout;
	sout << "TipsyFile: error opening file <" << filename << ">";
	throw std::runtime_error(sout.str());
      }
	
      xdrstdio_create(&xdrs, input, XDR_DECODE);
      
      if (xdr_header() != 1) {
	throw std::runtime_error("TipsyFile: could not read a valid header");
      }

      return header.nbodies;
    }

    //! Generate stream position for parallel read
    //! \param nsize The number of particles to read
    //! \param ssize The size of each particle structure in bytes
    //! \note Using unsigned long to avoid wrapping issues
    //! \return The number of particles to read for this process
    unsigned long xdr_psize(unsigned long nsize, unsigned long ssize, Ptype ptype)
    {
      // Number of particles per process
      unsigned long psize = nsize/numprocs;

      // Get current stream position
      curpos = xdr_getpos(&xdrs);
      auto pos = curpos;

      // Set the index offset for indexing particles
      index_offset_map[ptype] = psize * myid;

      // Move the stream position to the start of this process' data
      pos += index_offset_map[ptype] * ssize;
      if (xdr_setpos(&xdrs, pos) == 0) {
	throw std::runtime_error("TipsyFile: xdr_setpos failed");
      }

      // If this is the last process, read the remainder
      if (myid == numprocs-1) {
	index_offset_map[ptype] = (numprocs-1) * psize;
	psize = nsize - index_offset;
      }}

      return psize;
    }

    //! Set stream position to next group of particles
    //! \param nsize The number of particles to read
    //! \param ssize The size of each particle structure in bytes
    //! \note Using unsigned long to avoid wrapping issues
    void xdr_pnext(unsigned long nsize, unsigned long ssize)
    {
      // Move the stream position to the start of the next particle
      // group
      if (xdr_setpos(&xdrs, curpos + nsize * ssize) == 0) {
	throw std::runtime_error("TipsyFile: xdr_setpos failed");
      }
    }

    int xdr_read()
    {
      int N=0;
      
      if (header.nsph != 0)  {
	// Get the number of gas particles for this process and set
	// stream position
	auto psize = xdr_psize(header.nsph, sizeof(gas_particle), Ptype::gas);

	// Do the read
	gas_particles.resize(psize);
	read_gas();
	N++;

	// Set stream position to next group of particles
	xdr_pnext(header.nsph, sizeof(gas_particle));
      }
    
      if (header.ndark != 0) {
	// Get the number of dark particles for this process and set
	// stream position
	auto psize = xdr_psize(header.ndark, sizeof(dark_particle), Ptype::dark);

	// Do the read
	dark_particles.resize(psize);
	read_dark();
	N++;

	// Set stream position to next group of particles
	xdr_pnext(header.ndark, sizeof(dark_particle));
      }

      if (header.nstar != 0) {
	// Get the number of star particles for this process and set
	// stream position
	auto psize = xdr_psize(header.nstar, sizeof(star_particle), Ptype::star);

	// Do the read
	star_particles.resize(psize);
	read_star();
	N++;

	// Set stream position to next group of particles
	xdr_pnext(header.nstar, sizeof(star_particle));
      }
	
      return N;
    }

    void read_gas()
    {
      if (sizeof(Real) == sizeof(float)) {
	xdr_vector(&xdrs, (char *) &gas_particles[0],
		   gas_particles.size()*(sizeof(gas_particle)/sizeof(Real)),
		   sizeof(Real), (xdrproc_t) xdr_float);
      }
    }  
    
    void read_dark()
    {
      if (sizeof(Real) == sizeof(float)) {
	xdr_vector(&xdrs, (char *) &dark_particles[0],
		   dark_particles.size()*(sizeof(dark_particle)/sizeof(Real)),
		   sizeof(Real), (xdrproc_t) xdr_float);
      }
    }  
  
    void read_star()
    {
      if (sizeof(Real) == sizeof(float)) {
	xdr_vector(&xdrs, (char *) &star_particles[0],
		   star_particles.size()*(sizeof(star_particle)/sizeof(Real)),
		   sizeof(Real), (xdrproc_t) xdr_float);
      }
    }  
    
    friend class PR::Tipsy;

  public:

    //! Constructor that opens the file and reads the header
    TipsyXDR(const std::string filename)
    {
      // Attempt to read tipsy file
      int ret = xdr_init(filename);
      if (ret==0) {
	throw std::runtime_error("TipsyFile: no bodies found??");
      }
    }

    //! Read the particles from the file
    int readParticles() { return xdr_read(); }

    //! Destructor that removes the stream and closes the file
    ~TipsyXDR()
    {
      xdr_destroy(&xdrs);
      fclose(input);
    }
    
  };
#endif


  //! Tipsy native file reader
  class TipsyNative : public TipsyFile
  {
  private:

    //! The input file stream
    std::ifstream input;

    //! Cache stream position at beginning of partcile group
    std::istream::pos_type curpos = 0; 

    //! Read the header
    int native_header()
    {
      try {
	input.read((char *)&header, sizeof(header));
      } catch (std::exception& e) {
	std::ostringstream s;
	s << "TipsyFile native error reading header: " << e.what();
	throw std::runtime_error(s.str());
      }
      return 1;
    }
    
    //! Open the file and read the header
    //! \param filename The name of the file to read
    int native_init(const std::string& filename)
    {
      try {
	input.open(filename);
	input.exceptions(ifstream::eofbit | ifstream::failbit | ifstream::badbit);
      } catch (std::exception& e) {
	std::ostringstream s;
	s << "TipsyFile native error opening file <" << filename << ">: "
	  << e.what();
	throw std::runtime_error(s.str());
      }
	
      if (native_header() != 1) {
	throw std::runtime_error("TipsyFile native: could not read a valid header");
      }

      return header.nbodies;
    }

    //! Generate stream position for parallel read
    //! \param nsize The number of particles to read
    //! \param ssize The size of each particle structure in bytes
    //! \note Using unsigned long to avoid wrapping issues
    //! \return The number of particles to read for this process
    unsigned long ios_psize(unsigned long nsize, unsigned long ssize, Ptype ptype)
    {
      // Number of particles per process
      unsigned long psize = nsize/numprocs;

      // Get current stream position
      curpos = input.tellg();
      std::streampos pos(curpos);

      // Set the index offset for indexing particles
      index_offset_map[ptype] = psize * myid;

      // Move the stream position to the start of this process' data
      pos += index_offset_map[ptype] * ssize;

      if (input.seekg(pos).fail()) {
	throw std::runtime_error("TipsyFile: seekg failed");
      }

      std::streampos test = input.tellg();

      // If this is the last process, read the remainder
      if (myid == numprocs-1) {
	index_offset_map[ptype] = (numprocs-1) * psize;
	psize = nsize - index_offset_map[ptype];
      }

      return psize;
    }

    //! Move to next particle group
    //! \param nsize The number of particles to read
    //! \param ssize The size of each particle structure in bytes
    //! \note Using unsigned long to avoid wrapping issues
    void ios_pnext(unsigned long nsize, unsigned long ssize)
    {
      // Move the stream position to the start of the next particle
      // group
      std::streampos pos(curpos); // Make a copy of the current position
      pos += nsize * ssize;	  // Augment the position

      if (input.seekg(pos).fail()) {
	throw std::runtime_error("TipsyFile: seekg failed");
      }
    }

    //! Read the particles from the file
    //! \return The number of particle groups read
    int native_read()
    {
      int N=0;
      
      if (header.nsph != 0)  {
	// Get the number of gas particles for this process and set
	// stream position
	auto psize = ios_psize(header.nsph, sizeof(gas_particle), Ptype::gas);

	// Do the read
	gas_particles.resize(psize);
	read_gas();
	N++;

	// Set stream position to next group of particles
	ios_pnext(header.nsph, sizeof(gas_particle));
      }
    
      if (header.ndark != 0) {
	// Get the number of dark particles for this process and set
	// stream position
	auto psize = ios_psize(header.ndark, sizeof(dark_particle), Ptype::dark);

	// Do the read
	dark_particles.resize(psize);
	read_dark();
	// Stream position debug
	std::streampos test = input.tellg();
	N++;

	// Set stream position to next group of particles
	ios_pnext(header.ndark, sizeof(dark_particle));
      }

      if (header.nstar != 0) {
	// Get the number of star particles for this process and set
	// stream position
	auto psize = ios_psize(header.nstar, sizeof(star_particle), Ptype::star);

	// Do the read
	star_particles.resize(psize);
	read_star();
	N++;

	// Set stream position to next group of particles
	ios_pnext(header.nstar, sizeof(star_particle));
      }
	
      return N;
    }

    //! Read the gas particles from the file
    void read_gas()
    {
      try {
	input.read((char *) &gas_particles[0],
		   gas_particles.size()*sizeof(gas_particle));
      }
      catch (std::exception& e) {
	std::ostringstream s;
	s << "TipsyFile native error reading gas particles: " << e.what();
	throw std::runtime_error(s.str());
      }
    }
    
    //! Read the dark particles from the file
    void read_dark()
    {
      try {
	input.read((char *) &dark_particles[0],
		   dark_particles.size()*sizeof(dark_particle));
      }
      catch (std::exception& e) {
	std::ostringstream s;
	s << "TipsyFile native error reading dark particles: " << e.what();
	throw std::runtime_error(s.str());
      }
    }  
  
    //! Read the star particles from the file
    void read_star()
    {
      try {
	input.read((char *) &star_particles[0],
		   star_particles.size()*sizeof(star_particle));
      }
      catch (std::exception& e) {
	std::ostringstream s;
	s << "TipsyFile native error reading star particles: " << e.what();
	throw std::runtime_error(s.str());
      }
    }  
    
    friend class PR::Tipsy;

  public:

    //! Constructor that opens the file and reads the header
    TipsyNative(const std::string filename)
    {
      // Attempt to read tipsy file
      int ret = native_init(filename);
      if (ret==0) {
	throw std::runtime_error("TipsyFile: no bodies found??");
      }
    }

    //! Read the particles from the file
    int readParticles() { return native_read(); }

    //! Destructor that closes the file by scoping
    ~TipsyNative()
    {
      // Nothing
    }
    
  };

};

#endif
