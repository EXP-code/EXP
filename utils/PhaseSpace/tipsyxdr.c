#include <unistd.h>
#include "tipsydefs.h"
#include <rpc/types.h>
#include <rpc/xdr.h>


static XDR xdrs;
struct gas_particle *gas_particles;
struct dark_particle *dark_particles;
struct star_particle *star_particles;
struct dump header ;

int xdr_header()
{
  int pad;
  
  if(xdr_double(&xdrs, &header.time) != TRUE)
    return 0;
  if(xdr_int(&xdrs, &header.nbodies) != TRUE)
    return 0;
  if(xdr_int(&xdrs, &header.ndim) != TRUE)
    return 0;
  if(xdr_int(&xdrs, &header.nsph) != TRUE)
    return 0;
  if(xdr_int(&xdrs, &header.ndark) != TRUE)
    return 0;
  if(xdr_int(&xdrs, &header.nstar) != TRUE)
    return 0;
  if(xdr_int(&xdrs, &pad) != TRUE)
    return 0;

  return 1;
}

void xdr_gas()
{
  if (sizeof(Real) == sizeof(float))
    {
      xdr_vector(&xdrs, (char *) gas_particles,
		 header.nsph*(sizeof(*gas_particles)/sizeof(Real)),
		 sizeof(Real), xdr_float);
    }
}  

void xdr_dark()
{
  if (sizeof(Real) == sizeof(float))
    {
      xdr_vector(&xdrs, (char *) dark_particles,
		 header.ndark*(sizeof(*dark_particles)/sizeof(Real)),
		 sizeof(Real), xdr_float);
    }
}  

void xdr_star()
{
  if(sizeof(Real) == sizeof(float))
    {
      xdr_vector(&xdrs, (char *) star_particles,
		 header.nstar*(sizeof(*star_particles)/sizeof(Real)),
		 sizeof(Real), xdr_float);
    }
}  

int xdr_init()
{
  xdrstdio_create(&xdrs, stdin, XDR_DECODE);

  if(xdr_header() != 1) {
    fprintf(stderr, "Could not read a valid header\n");
    exit(-1);
  }

  int N=0;

  if(header.nsph != 0) {
    gas_particles = (struct gas_particle *)
      malloc(header.nsph*sizeof(*gas_particles));
    if(gas_particles == NULL) {
      fprintf(stderr, "Out of memory for gas particles\n");
      exit(-1);
    }
    N++;
  }
  else
    gas_particles = NULL;

  if(header.ndark != 0) {
    dark_particles = (struct dark_particle *)
      malloc(header.ndark*sizeof(*dark_particles));
    if(dark_particles == NULL) {
      fprintf(stderr, "Out of memory for dark particles\n");
      exit(-1);
    }
    N++;
  }
  else
    dark_particles = NULL;

  if(header.nstar != 0) {
    star_particles = (struct star_particle *)
      malloc(header.nstar*sizeof(*star_particles));
    if(star_particles == NULL) {
      fprintf(stderr, "Out of memory for star particles\n");
      exit(-1);
    }
    N++;
  }
  else
    star_particles = NULL;

  return N;
}

void xdr_quit()
{
  xdr_destroy(&xdrs);
}
