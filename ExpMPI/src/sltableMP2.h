
//! Holds the tabulated Sturm-Liouville basis for a cylindrical basis
struct TableCyl {

  int m;
  int k;

  Vector ev;
  Matrix ef;
};


//! Holds the tabulated Sturm-Liouville basis for a spherical basis
struct TableSph {

  int l;

  Vector ev;
  Matrix ef;
};

//! Holds the tabulated Sturm-Liouville basis for a slab basis
struct TableSlab {

  int kx;
  int ky;

  Vector ev;
  Matrix ef;
};

