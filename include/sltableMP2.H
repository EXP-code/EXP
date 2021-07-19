#ifndef _sltableMP2_H
#define _sltableMP2_H

#include <Eigen/Eigen>

class TableCyl 
{
 public:
  int m;
  int k;

  Eigen::VectorXd ev;
  Eigen::MatrixXd ef;
};

class TableSph 
{
 public:
  int l;

  Eigen::VectorXd ev;
  Eigen::MatrixXd ef;
};

class TableSlab 
{
public:
  int kx;
  int ky;

  Eigen::VectorXd ev;
  Eigen::MatrixXd ef;
};

#endif
