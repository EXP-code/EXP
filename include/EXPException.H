#ifndef EXPException_H
#define EXPException_H

#include <string>
#include <sstream>

/** Defines an error handler base class EXPException to handle exceptions
 */
class EXPException : public std::exception
{
  public:
  //! Constructor for base
  EXPException(const std::string exceptionname, const std::string message,
	       const std::string sourcefilename, int sourcelinenumber)
  {
    errormessage = new std::ostringstream;
    this->sourcefilename   = sourcefilename;
    this->sourcelinenumber = sourcelinenumber; 
    this->exceptionname = exceptionname;
    (*errormessage) << message;
  }

  //! Destructor for base
  virtual ~EXPException()
  {
    delete errormessage;
    errormessage = 0;
  }

  //! Returns an error message suitable for printing to the user.
  const std::string getErrorMessage() const
  {
    // This combines the error name, error message, and the throw locations.
    std::string buffer = (*errormessage).str();

    std::ostringstream wholemessage;
    wholemessage << exceptionname << ": " << buffer;
    wholemessage << "Thrown from " << sourcefilename << ":" << sourcelinenumber;
    return wholemessage.str();  // stringstream returns a const std::string
  }

  //! std::exception member
  virtual const char* what() const throw ()
  {
    return getErrorMessage().c_str();
  }
  
protected:
  //! Protected so it is only called by properly implemented classes
  EXPException(const std::string sourcefile, int linenumber)
  {
    errormessage = new std::ostringstream;
    this->sourcefilename   = sourcefilename;
    this->sourcelinenumber = sourcelinenumber; 
  }

  //! Friendly name of the exception
  std::string exceptionname;

  //! Error message describing the error in more detail
  std::ostringstream * errormessage;
  
private:
  //! Source file where throw occured
  std::string sourcefilename;

  //! Line number of throw
  int sourcelinenumber;
};

//! Used for explicit detailed messages
class GenericError : public EXPException
{
public:
  //! Use this for reporting an error with a specific message
  //@{
  GenericError(const std::string msg, const std::string sourcefilename,
	       int sourcelinenumber)
    : EXPException(sourcefilename, sourcelinenumber)
  {
    exceptionname = "Execution exception";
    (*errormessage) << "Error details: " << msg << ".";
  }
  //@}
};

//! Used when execution reaches a point it should not reach.
class InternalError : public EXPException
{
public:
  //! Use this when you reach an unexpected state.
  //@{
  InternalError
  (const std::string sourcefilename, int sourcelinenumber)
    : EXPException(sourcefilename, sourcelinenumber)
  {
    exceptionname = "Internal Error Exception";
    (*errormessage) << "Execution reached a state that should not reachable.";
  }

  InternalError
  (const std::string msg, const std::string sourcefilename,
   int sourcelinenumber)
    : EXPException(sourcefilename, sourcelinenumber)
  {
    exceptionname = "Internal Error Exception";
    (*errormessage) << "Execution reached a state that should not reachable. "
		    << msg << ".";
  }
  //@}
};

//! Handle bad range related exception PartMap
class BadIndexException : public EXPException
{ 
public:
  //! Constructor
  BadIndexException(int index, int num,
		    const std::string sourcefilename, int sourcelinenumber)
    : EXPException(sourcefilename, sourcelinenumber)
  {
    exceptionname = "Requested index is not in PartMap: ";
    (*errormessage) << "Invalid index: " << index 
		    << " out of " << num << " particles";
  }

};

//! File creation error
class FileCreateError : public EXPException
{ 
public:
  //! Constructor
  FileCreateError(const std::string filename,
		  std::string sourcefilename, int sourcelinenumber)
  : EXPException(sourcefilename, sourcelinenumber)
  {
    exceptionname = "Cannot create new file";
    (*errormessage) << "File name is <" << filename << ">";
  }

  //! Constructor with method string
  FileCreateError(const std::string filename, const std::string method,
		  const std::string sourcefilename, int sourcelinenumber)
  : EXPException(sourcefilename, sourcelinenumber)
  {
    exceptionname = "Cannot create new file in " + method;
    (*errormessage) << "File name is <" << filename << ">";
  }

};

//! File open error
class FileOpenError : public EXPException
{ 
public:
  //! Constructor
  FileOpenError(const std::string filename,
		const std::string sourcefilename, int sourcelinenumber)
    : EXPException(sourcefilename, sourcelinenumber)
  {
    exceptionname = "Cannot open existing file";
    (*errormessage) << "File name is <" << filename << ">";
  }

  //! Constructor with method string
  FileOpenError(const std::string filename, const std::string method,
		const std::string sourcefilename, int sourcelinenumber)
    : EXPException(sourcefilename, sourcelinenumber)
  {
    exceptionname = "Cannot open existing file in " + method;
    (*errormessage) << "File name is <" << filename << ">";
  }
};

#endif
