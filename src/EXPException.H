#ifndef EXPException_H
#define EXPException_H

#include <string>
#include <sstream>

using namespace std;

/** Defines an Error Handler base class EXPException to handle exceptions
 */
class EXPException {
  public:
  //! Generic exception
  EXPException(string exceptionname, string message,
	       string sourcefilename, int sourcelinenumber);
  //! Destructor.
  virtual ~EXPException();
  //! Returns an error message suitable for printing to the user.
  string getErrorMessage();
    
protected:
  //! Protected so it is only called by properly implemented classes.
  EXPException(string sourcefile, int linenumber);
  //! Friendly name of the exception.
  string exceptionname;
  //! Error message describing the error in more detail.
  ostringstream * errormessage;
  
private:
  //! Source file where throw occured.
  string sourcefilename;
  //! Line number of throw.
  int sourcelinenumber;
};

//! Used when execution reaches a point it should not reach.
class InternalError : public EXPException {
public:
  //! Use this when you reach an unexpected state.
  //@{
  InternalError(string sourcefilename, int sourcelinenumber);
  InternalError(string msg, string sourcefilename, int sourcelinenumber);
  //@}
};

//! handle Bad range related exception PartMap
class BadIndexException : public EXPException { 
public:
  //! Constructor
  BadIndexException(int index, int num,
		    string sourcefilename, int sourcelinenumber);
};

#endif
