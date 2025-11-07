#include <iostream>
#include <iomanip>
#include "header.H"

int ComponentHeader::defaultInfoSize = 1024;

bool ComponentHeader::write(ostream *out)
{
  out->write((const char *)&nbod,  sizeof(int));
  out->write((const char *)&niatr, sizeof(int));
  out->write((const char *)&ndatr, sizeof(int));
  out->write((const char *)&ninfochar, sizeof(int));
  out->write((const char *)info.get(), ninfochar*sizeof(char));

  if (*out)
    return true;
  else
    return false;
}

bool ComponentHeader::write_mpi(MPI_File& out, MPI_Offset& offset)
{
  MPI_Status status;
  char err[MPI_MAX_ERROR_STRING];
  int len, ret;

  ret = MPI_File_write_at(out, offset, &nbod, 1, MPI_INT, &status);

  if (ret != MPI_SUCCESS) {
    MPI_Error_string(ret, err, &len);
    std::cout << "ComponentHeader::write_mpi: " << err
	      << " at line " << __LINE__ << std::endl;
    return false;
  }

  offset += sizeof(int);

  ret = MPI_File_write_at(out, offset, &niatr, 1, MPI_INT, &status);

  if (ret != MPI_SUCCESS) {
    MPI_Error_string(ret, err, &len);
    std::cout << "ComponentHeader::write_mpi: " << err
	      << " at line " << __LINE__ << std::endl;
    return false;
  }

  offset += sizeof(int);

  ret = MPI_File_write_at(out, offset, &ndatr, 1, MPI_INT, &status);

  if (ret != MPI_SUCCESS) {
    MPI_Error_string(ret, err, &len);
    std::cout << "ComponentHeader::write_mpi: " << err
	      << " at line " << __LINE__ << std::endl;
    return false;
  }

  offset += sizeof(int);

  ret = MPI_File_write_at(out, offset, &ninfochar, 1, MPI_INT, &status);

  if (ret != MPI_SUCCESS) {
    MPI_Error_string(ret, err, &len);
    std::cout << "ComponentHeader::write_mpi: " << err
	      << " at line " << __LINE__ << std::endl;
    return false;
  }

  offset += sizeof(int);

  ret = MPI_File_write_at(out, offset, info.get(), ninfochar, MPI_CHAR, &status);

  if (ret != MPI_SUCCESS) {
    MPI_Error_string(ret, err, &len);
    std::cout << "ComponentHeader::write_mpi: " << err
	      << " at line " << __LINE__ << std::endl;
    return false;
  }

  offset += ninfochar;

  return true;
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
    // Use this as of C++17
    // info = std::make_shared<char[]>(ninfochar+1);

    // C++14 workaround:
    info = std::shared_ptr<char>(new char[ninfochar+1],
				 std::default_delete<char[]>());
    // This ensures that info is null terminated
    std::fill(info.get(), info.get()+ninfochar+1, '\0');
  }
  
  in->read((char *)info.get(), ninfochar*sizeof(char));
  if (!*in) return false;

  return true;
}

