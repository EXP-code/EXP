#include <Vector.h>

struct TableCyl {

  int m;
  int k;

  Vector ev;
  Matrix ef;
};

struct TableSph {

  int l;

  Vector ev;
  Matrix ef;
};

struct TableSlab {

  int kx;
  int ky;

  Vector ev;
  Matrix ef;
};

