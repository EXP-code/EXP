
//! Holds the tabulated Sturm-Liouville basis for a cylindrical basis
struct TableCyl {

  //! The azimuthal order
  int m;
  //! The vertical order
  int k;

  //! Eigenvalues
  Vector ev;

  //! Eigenvectors
  Matrix ef;
};


//! Holds the tabulated Sturm-Liouville basis for a spherical basis
struct TableSph {

  //! Harmonic order
  int l;

  //! Eigenvalues
  Vector ev;

  //! Eigenvectors
  Matrix ef;
};

//! Holds the tabulated Sturm-Liouville basis for a slab basis
struct TableSlab {

  //! Number of x wavenumber
  int kx;

  //! Number of y wavenumber
  int ky;

  //! Eigenvalues
  Vector ev;

  //! Eigenvectors
  Matrix ef;
};

