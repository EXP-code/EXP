#include <EXPException.H>
#include <errno.h>

string EXPException::getErrorMessage()
{
  // This combines the error name, error message, and the throw locations.
  // (*errormessage) << ends;  removed for conversion to stringstream
  //  char * buffer = (*errormessage).str();
  string buffer = (*errormessage).str();

  ostringstream wholemessage;
  wholemessage << exceptionname << ": " << buffer << endl;
  wholemessage << "Thrown from " << sourcefilename << ":" << sourcelinenumber;
  return wholemessage.str();  // stringstream returns a string
}

EXPException::EXPException
(string exceptionname, string message, string sourcefilename, int sourcelinenumber)
{
  errormessage = new ostringstream;
  this->sourcefilename   = sourcefilename;
  this->sourcelinenumber = sourcelinenumber; 
  this->exceptionname = exceptionname;
  (*errormessage) << message;
}

EXPException::EXPException(string sourcefilename, int sourcelinenumber)
{
  errormessage = new ostringstream;
  this->sourcefilename   = sourcefilename;
  this->sourcelinenumber = sourcelinenumber; 
}

EXPException::~EXPException()
{
  delete errormessage;
  errormessage = 0;
}

InternalError::InternalError
(string sourcefilename, int sourcelinenumber)
: EXPException(sourcefilename, sourcelinenumber)
{
  exceptionname = "Internal Error Exception";
  (*errormessage) << "Execution reached a state that should not reachable.";
}

InternalError::InternalError
(string msg, string sourcefilename, int sourcelinenumber)
: EXPException(sourcefilename, sourcelinenumber)
{
  exceptionname = "Internal Error Exception";
  (*errormessage) << "Execution reached a state that should not reachable. "
		  << msg << ".";
}

BadIndexException::BadIndexException
(int index, int num,
 string sourcefilename, int sourcelinenumber)
: EXPException(sourcefilename, sourcelinenumber)
{
  exceptionname = "Requested index is not in PartMap: ";
  (*errormessage) << "Invalid index: " << index 
		  << " out of " << num << " particles" << endl;
}

