// -*- C++ -*-

#ifndef EXPException_h
#define EXPException_h

#include <vector>
#include <string>
#include <sstream>


/** Defines an Error Handler base class EXPException to handle
    exceptions
*/
class Exception : public std::exception
{
public:
  //! Generic exception
  explicit Exception(std::string exceptionname, std::string message,
			std::string sourcefilename, int sourcelinenumber);
  
  //! Copy constructor
  Exception(const Exception& e);
  
  //! Destructor.
  virtual ~Exception() {}
  
  //! Returns an error message suitable for printing to the user.
  std::string getErrorMessage();
  
  /** Returns a pointer to the (constant) error description.
   *  @return A pointer to a const char*. The underlying memory
   *          is in posession of the Exception object. Callers must
   *          not attempt to free the memory.
   */
  virtual const char* what() const throw (){
    return getErrorMessage().c_str();
  }


protected:
  //! Protected so it is only called by properly implemented classes.
  Exception(string sourcefile, int linenumber);
  
  //! Friendly name of the exception.
  std::string exceptionname;
  
  //! Error message describing the error in more detail.
  std::ostringstream errormessage;
  
protected:
  //! Source file where throw occured.
  std::string sourcefilename;
  
  //! Line number of throw.
  int sourcelinenumber;
};


//! Used when execution reaches a point it should not reach.
class InternalError : public Exception 
{
public:
  //! Use this when you reach an unexpected state.
  //@{
  InternalError(std::string sourcefilename, int sourcelinenumber);
  
  InternalError(std::string msg, 
		std::string sourcefilename, int sourcelinenumber);
  
  InternalError(int err, std::string msg, 
		std::string sourcefilename, int sourcelinenumber);
  //@}
};

    
//! Handles open file related exceptions.
class FileOpenException : public Exception 
{
public: 
  //! Reports filename, errno error, and location of throw
  FileOpenException(std::string filename, int errno_, 
		    std::string sourcefilename, int sourcelinenumber);
};

//! Handle file creation related exceptions.
class FileCreateException : public Exception 
{
public: 
  //! Reports filename, errno error, and location of throw
  FileCreateException(std::string filename, int errno_, 
		      std::string sourcefilename, int sourcelinenumber);
};

//! Handle file format related exceptions for DataStreams
class FileFormatException : public Exception 
{
public:
  //! Constructor: supply diagnostic message
  FileFormatException(std::string message, 
		      std::string sourcefilename, int sourcelinenumber);
};
  
  
#endif
