#include "SParticle.H"

SPtype::SPtype()
{
  // Make MPI datatype
  //
  MPI_Datatype type[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  //                      ^           ^           ^
  //                      |           |           |
  //  mass----------------+           |           |
  //  position------------------------+           |
  //  velocity------------------------------------+
  //  [and that's it]
  
  // Get displacements
  //
  SParticle buf;

  MPI_Aint disp[3];
  MPI_Get_address(&buf.mass,	&disp[0]);
  MPI_Get_address(&buf.pos[0],	&disp[1]);
  MPI_Get_address(&buf.vel[0],	&disp[2]);
  
  for (int i=2; i>=0; i--) disp[i] -= disp[0];
  
  // Block offsets
  //
  int blocklen[3] = {1, 3, 3};
  
  // Make and register the new type
  //
  MPI_Type_create_struct(3, blocklen, disp, type, &Particletype);
  MPI_Type_commit(&Particletype);
}

SPtype::~SPtype()
{
  // Free the registered type
  //
  MPI_Type_free(&Particletype);
}
