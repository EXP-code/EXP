#ifndef _OutPCA_H
#define _OutPCA_H

#include <OutPCA.H>

//! Print PCA data; launched by Basis when pca is selected
class OutPCA : public Output
{

private:

  string pcaname;
  void zeroth_order(void);

public:

  //! Constructor: <code>param line</code> is the filename
  OutPCA(string& line);

  //! Make the output
  /*!
    \param nstep is the current time step used to decide whether or not
    to dump
  */
  void Run(int nstep);

};

#endif
