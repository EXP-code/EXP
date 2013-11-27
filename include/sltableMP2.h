#include <Vector.h>

class TableCyl 
{
 public:
  int m;
  int k;

  Vector ev;
  Matrix ef;
};

class TableSph 
{
 public:
  int l;

  Vector ev;
  Matrix ef;
};

class TableSlab 
{
public:
  int kx;
  int ky;

  Vector ev;
  Matrix ef;
};

