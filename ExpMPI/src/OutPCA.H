#ifndef _OutPCA_H
#define _OutPCA_H

#include <OutPCA.H>

//! Print PCA data; launched by Basis when pca is selected
class OutPCA : public Output
{

private:

  string pcaname;
  void zeroth_order(void)

public:

  OutPCA(string& line);

  void Run(int nstep);

};

#endif
