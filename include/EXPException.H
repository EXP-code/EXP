#ifndef EXPException_H
#define EXPException_H

#include <set>
#include <vector>
#include <string>
#include <sstream>

#include "libvars.H"

/** Defines an error handler base class EXPException to handle exceptions
 */
class EXPException : public std::exception
{
  public:
  /** Constructor for base
      @param exceptionname     Human-readable identifier for the error
      @param message           Error details for the user
      @param sourcefilename    File where the thrower resides
      @param sourcelinenumber  Line number in the source for the thrower
      @param errorcode         Error code to report at top-level catch
      @param deadlock          Exception may break MPI synchronization
  */
  EXPException(const std::string exceptionname, const std::string message,
	       const std::string sourcefilename, int sourcelinenumber,
	       int errorcode=-1, bool deadlock=false)
  {
    this->sourcefilename   = sourcefilename;
    this->sourcelinenumber = sourcelinenumber; 
    this->exceptionname    = exceptionname;
    this->errormessage     = message;
    this->errorcode        = errorcode;
    this->deadlock         = deadlock;

    msg_ = getErrorMessage();
  }

  //! Destructor for base.  Virtual to permit derived classes.
  virtual ~EXPException() {}

  /** 
      Returns an error message suitable for printing to the user.
      Also called by derived classes to assemble the error message for
      std::exception::what() 
  */
  std::string getErrorMessage()
  {
    // This combines the error name, error message, and the throw locations
    //
    std::ostringstream wholemessage;
    wholemessage << exceptionname << ": " << errormessage;
    if (__EXP__::sourceline)
    wholemessage << "Thrown from " << sourcefilename << ":" << sourcelinenumber;

    return wholemessage.str();
  }

  /**  std::exception member returns a pointer to the error message

       @return A pointer to a const char*. 

       The std::string "buffer" used for assembling and passing the
       message to std::exception will be deleted when the exception is
       deleted.
  */
  virtual const char* what() const throw ()
  {
    return msg_.c_str();
  }
  
  //! Get errorcode
  int getErrorcode() { return errorcode; }

  //! Get deadlock status
  bool getDeadlock() { return deadlock; }

protected:

  //! Protected so it is only called by properly implemented classes
  EXPException(const std::string sourcefile, int linenumber,
	       int errorcode=-1, bool deadlock=false)
  {
    this->sourcefilename   = sourcefile;
    this->sourcelinenumber = linenumber; 
    this->errorcode        = errorcode;
    this->deadlock         = deadlock;

    msg_ = getErrorMessage();
  }

  //! Friendly name of the exception
  std::string exceptionname;

  //! Error message describing the error in more detail
  std::string errormessage;

  //! Message buffer
  std::string msg_;

  //! Possible deadlock?
  bool deadlock;

  //! Return code
  int errorcode;

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
	       int sourcelinenumber, int errorcode=-1, bool deadlock=false)
    : EXPException(sourcefilename, sourcelinenumber, errorcode, deadlock)
  {
    exceptionname = "Execution Exception";
    std::ostringstream ostr; ostr << msg << ".";
    errormessage = ostr.str();

    // Assemble and store the error message
    msg_ = getErrorMessage();
  }
  //@}
};

//! Used when execution reaches a point it should not reach
class InternalError : public EXPException
{
public:
  //! Use this when you reach an unexpected state.
  //@{
  InternalError
  (const std::string sourcefilename, int sourcelinenumber,
   int errorcode=-1, bool deadlock=true)
    : EXPException(sourcefilename, sourcelinenumber, errorcode, deadlock)
  {
    exceptionname = "Internal Error Exception";
    errormessage  = "Execution reached a state that should not reachable.";

    msg_ = getErrorMessage();
  }

  InternalError
  (const std::string msg,
   const std::string sourcefilename, int sourcelinenumber,
    int errorcode=-1, bool deadlock=true)
    : EXPException(sourcefilename, sourcelinenumber, errorcode, deadlock)
  {
    exceptionname = "Internal Error Exception";
    errormessage  = "Execution reached a state that should not reachable. "
      +  msg + ".";

    msg_ = getErrorMessage();
  }
  //@}
};

//! Handle bad range related exception PartMap
class BadIndexException : public EXPException
{ 
public:
  //! Constructor
  BadIndexException(int index, int num,
		    const std::string sourcefilename, int sourcelinenumber,
		    int errorcode=-1, bool deadlock=true)
    : EXPException(sourcefilename, sourcelinenumber, errorcode, deadlock)
  {
    exceptionname = "Requested index is not in PartMap: ";
    std::ostringstream sout;
    sout << "Invalid index: " << index << " out of " << num << " particles";
    errormessage = sout.str();

    msg_ = getErrorMessage();
  }

};

//! File creation error
class FileCreateError : public EXPException
{ 
public:
  //! Constructor
  FileCreateError(const std::string filename,
		  std::string sourcefilename, int sourcelinenumber,
		  int errorcode=-1, bool deadlock=false)
    : EXPException(sourcefilename, sourcelinenumber, errorcode, deadlock)
  {
    exceptionname = "Cannot create new file";
    errormessage  = "File name is <" + filename + ">";

    msg_ = getErrorMessage();
  }

  //! Constructor with method string
  FileCreateError(const std::string filename, const std::string method,
		  const std::string sourcefilename, int sourcelinenumber,
		  int errorcode=-1, bool deadlock=false)
    : EXPException(sourcefilename, sourcelinenumber, errorcode, deadlock)
  {
    exceptionname = "Cannot create new file in " + method;
    errormessage  = "File name is <" + filename + ">";

    msg_ = getErrorMessage();
  }

};

//! File open error
class FileOpenError : public EXPException
{ 
public:
  //! Constructor
  FileOpenError(const std::string filename,
		const std::string sourcefilename, int sourcelinenumber,
		int errorcode=-1, bool deadlock=false)
    : EXPException(sourcefilename, sourcelinenumber, errorcode, deadlock)
  {
    exceptionname = "Cannot open existing file";
    errormessage  = "File name is <" + filename + ">";

    msg_ = getErrorMessage();
  }

  //! Constructor with method string
  FileOpenError(const std::string filename, const std::string method,
		const std::string sourcefilename, int sourcelinenumber,
		int errorcode=-1, bool deadlock=false)
    : EXPException(sourcefilename, sourcelinenumber, errorcode, deadlock)
  {
    exceptionname = "Cannot open existing file in " + method;
    errormessage  = "File name is <" + filename + ">";

    msg_ = getErrorMessage();
  }
};

//! Configuration error
class YamlConfigError : public EXPException
{ 
public:
  //! Constructor
  YamlConfigError(const std::string classname, const std::string type,
		  const std::set<std::string> keys,
		  const std::string sourcefilename, int sourcelinenumber,
		  int errorcode=-1, bool deadlock=false)
    : EXPException(sourcefilename, sourcelinenumber, errorcode, deadlock)
  {
    exceptionname = classname + " unmatched " + type + " configuration keys";
    errormessage = "<";
    for (auto s : keys) {
      if (errormessage.size()>1) errormessage += " ";
      errormessage += s;
    }
    errormessage += ">";

    msg_ = getErrorMessage();
  }
};

#endif
