#ifndef _ParamParseMPI_H
#define _ParamParseMPI_H

#include <list>
#include <string>
#include <pair.h>
#include <mpi.h>

#include <ParamParse.H>

extern int numprocs, myid;

class ParamParseMPI : public ParamParse
{
public:

  ParamParseMPI(istream* in, string Delim);
};

#endif
