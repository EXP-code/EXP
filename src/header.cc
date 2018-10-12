#include <iostream>
#include <iomanip>
#include <header.H>

int ComponentHeader::defaultInfoSize = 1024;

bool ComponentHeader::write(ostream *out)
{
  out->write((char *)&nbod,  sizeof(int));
  out->write((char *)&niatr, sizeof(int));
  out->write((char *)&ndatr, sizeof(int));
  out->write((char *)&ninfochar, sizeof(int));
  out->write((char *)info.get(), ninfochar*sizeof(char));

  if (*out)
    return true;
  else
    return false;
}

bool ComponentHeader::write_mpi(MPI_File& out, MPI_Offset& offset)
{
  MPI_Status status;

  MPI_File_write_at(out, offset, &nbod,      1, MPI_INT, &status);
  offset += sizeof(int);

  MPI_File_write_at(out, offset, &niatr,     1, MPI_INT, &status);
  offset += sizeof(int);

  MPI_File_write_at(out, offset, &ndatr,     1, MPI_INT, &status);
  offset += sizeof(int);

  MPI_File_write_at(out, offset, &ninfochar, 1, MPI_INT, &status);
  offset += sizeof(int);

  MPI_File_write_at(out, offset, info.get(), ninfochar, MPI_CHAR, &status);
  offset += ninfochar;

  if (status.MPI_ERROR == MPI_SUCCESS)
    return true;
  else
    return false;
}

bool ComponentHeader::read(istream *in)
{
  int ninfo;

  in->read((char *)&nbod,  sizeof(int));		if (!*in) return false;
  in->read((char *)&niatr, sizeof(int));		if (!*in) return false;
  in->read((char *)&ndatr, sizeof(int));		if (!*in) return false;
  in->read((char *)&ninfo, sizeof(int));		if (!*in) return false;

  if (ninfo != ninfochar) {
    ninfochar = ninfo;
    info = boost::shared_array<char>(new char [ninfochar]);
  }
  in->read((char *)info.get(), ninfochar*sizeof(char));	if (!*in) return false;

  return true;
}

