#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <FieldGenerator.H>

namespace py = pybind11;

#include "TensorToArray.H"

void FieldGeneratorClasses(py::module &m) {

  m.doc() =
    R"(
    FieldGenerator class bindings

    FieldGenerator
    ==============
    This class computes surfaces and volumes for visualizing the physical
    quantities implied by your basis and associated coefficients.

    The generator is constructed by passing a vector of desired times
    that must be in the coefficient object and list of lower bounds,
    a list of upper bounds, and a list of knots per dimension.  These
    lists all have rank 3 for (x, y, z).  For a two-dimensional surface,
    one of the knot array values must be zero.  The member functions
    lines, slices and volumes, called with the basis and coefficient
    objects, return a numpy.ndarray containing the field evaluations.
    Each of these functions returns a dictionary of times to a dictionary
    of field names to numpy.ndarrays at each time.  There are also members
    which will write these generated fields to files. The linear probe
    members, 'lines' and 'file_lines', evaluate 'num' field points along
    a user-specified segment between the 3d points 'beg' and 'end'.  See
    help(pyEXP.basis) and help(pyEXP.coefs) for info on the basis and
    coefficient objects.

    Data packing
    ------------
    All slices and volumes are returned as numpy.ndarrays in row-major
    order.  That is, the first index is x, the second is y, and the
    third is z.  The ranks are specified by the 'gridsize' array with
    (nx, ny, nz) as input to the FieldGenerator constructor.

    Coordinate systems
    ------------------
    The FieldGenerator class supports spherical, cylindrical, and 
    Cartesian for force field components.  These are selected by your
    basis instance, consistent with the natural coordinate system for
    that basis.  You may change the default coorindate system for a basis
    using the 'setFieldType()' member function.

    Field names
    -----------
    The data fields are as follows:
    * 'dens'        the total density
    * 'dens m>0'    the non-axisymmetric component of the density
    * 'dens m=0'    the axisymmetric component of the density
    * 'rad force'   the radial force
    * 'mer force'   the meridional force
    * 'azi force'   the azimuthal force
    * 'potl'        the total potential
    * 'potl m>0'    the non-axisymmetric component of the potential
    * 'potl m=0'    the axisymmetric component of the potential

    For spherical coordinates (coord="Spherical"):
    * 'rad force'   the radial force
    * 'mer force'   the meridional force
    * 'azi force'   the azimuthal force

    For cylindrical coordinates (coord="Cylindrical"):
    * 'rad force'   the radial force
    * 'ver force'   the meridional force
    * 'azi force'   the azimuthal force

    For Cartesian coordinates (coord="Cartesian"):
    * 'x force'     the radial force
    * 'y force'     the meridional force
    * 'z force'     the azimuthal force

    Notes
    -----
    Note that the 'dens' field is the sum of the 'dens m=0' and 'dens m>0'
    fields and, similarly, the 'potl' field is the sum of 'potl m=0' and 
    'potl m>0' fields. These redundant entries are provided for convenience 
    and conceptual clarity.  For spherical bases, the 'm' index should be
    interpreted as the 'l' index.
    )";
    

  using namespace Field;

  py::class_<Field::FieldGenerator, std::shared_ptr<Field::FieldGenerator>>
    f(m, "FieldGenerator");

  f.def(py::init<const std::vector<double>, const std::vector<double>,
	const std::vector<double>, const std::vector<int>>(),
	R"(
        Create fields for given times and lower and upper bounds and grid sizes

        Parameters
        ----------
        times : list(float,...)
            list of evaluation times
        lower : list(float, float, float)
            lower grid boundaries for x, y, z
        upper : list(float, float, float)
            upper grid boundaries for x, y, z
        gridsize : list(int, int, int)
            grid sizes for each dimension; use zero in one dimension for a plane

        Returns
        -------
        FieldGenerator
            new object
        )",
	py::arg("times"), py::arg("lower"), py::arg("upper"),
	py::arg("gridsize"));

  f.def("setMidplane", &Field::FieldGenerator::setMidplane,
	R"(
        Set the field generator to generate midplane fields

        Parameters
        ----------
        on : bool
           True to generate midplane fields
        )", py::arg("on"));

  f.def("setColumnHeight", &Field::FieldGenerator::setColumnHeight,
	R"(
        Set the column extent for midplane position search

        Parameters
        ----------
        colheight : double
           Number of scale heights above and below plane for search
        )", py::arg("colheight"));

  f.def("slices", &Field::FieldGenerator::slices,
	R"(
        Return a dictionary of grids (2d numpy arrays) indexed by time and field type

        Parameters
        ----------
        basis : Basis
            basis instance of any geometry; geometry will be deduced by the generator
        coefs : Coefs
            coefficient container instance

        Returns
        -------
        dict({time: {field-name: numpy.ndarray})
            dictionary of times, field names, and data arrays

        Notes
        -----
        Each data array in the dictionary is a surface with the range given
        by the constructor.  This can be visualized using the standard matplotlib
        routines.  See pyEXP-examples for a tutorial.

        See also
        --------
        lines : generate fields along a line given by its end points
        volumes : generate fields in volume given by the initializtion grid
       )", py::arg("basis"), py::arg("coefs"));
  
  f.def("lines", &Field::FieldGenerator::lines,
	R"(
        Return a dictionary of arrays indexed by time and field type

        Parameters
        ----------
        basis : Basis
            basis instance of any geometry; geometry will be deduced by the generator
        coefs : Coefs
            coefficient container instance

        beg : list(float, float, float)
            initial evaluation point

        end : list(float, float, float)
            final evaluation point

        num : int
            number of evaluations

        Returns
        -------
        dict({time: {field-name: numpy.ndarray}})
            dictionary of times, field names, and data arrays

        Notes
        -----
        Each data array in the dictionary is 1d Numpy array with the range given
        by the beg and end position lists.

        See also
        --------
        slices : generate fields in a surface slice given by the initializtion grid
        volumes : generate fields in volume given by the initializtion grid
        )",
	py::arg("basis"), py::arg("coefs"),
	py::arg("beg"), py::arg("end"), py::arg("num"));
  
  f.def("histo2d", &Field::FieldGenerator::histogram2d,
	R"(
        Compute a surface histogram from particlesReturn a density histogram

        Parameters
        ----------
        reader : ParticleReader
            particle reader instance
        center : list(float, float, float), default=[0, 0, 0]
            origin for computing the histogram

        Returns
        -------
        dict({time: numpy.ndarray})

        Notes
        -----
        Range for histogram is taken from the grid ranges in the constructor
        )",
	py::arg("reader"),
	py::arg("center") = std::vector<double>(3, 0.0));

  f.def("histo1d", &Field::FieldGenerator::histogram1d,
	R"(
        Make a 1d density histogram (array) for a chosen projection

        Parameters
        ----------
        reader : ParticleReader
            particle reader instance
        rmax : float
            linear extent of the histogram
        nbins : int
            number of bins
        projection : str
	    projection indicated by one of the string: \"xy\", \"xz\", \"yz\", \"r\"",
        center : list(float, float, float), default=[0, 0, 0]
            origin for computing the histogram

        Returns
        -------
        numpy.ndarray
            the computed 1d histogram
        )",
	py::arg("reader"), py::arg("rmax"), py::arg("nbins"),
	py::arg("projection"),
	py::arg("center") = std::vector<double>(3, 0.0));

  f.def("file_lines", &Field::FieldGenerator::file_lines,
	R"(
        Write field arrays to files using the supplied string prefix.

        Parameters
        ----------
        basis : Basis
            basis instance of any geometry; geometry will be deduced by the generator
        coefs : Coefs
            coefficient container instance
        beg : list(float, float, float)
            initial evaluation point
        end : list(float, float, float)
            final evaluation point
        num : int
            number of evaluations
        filename : str
            file name for output
        dir : str, default='.'
            directory to write files

        Returns
        -------
        None

        Notes
        -----
	The files are ascii tables with column headers.

        See also
        --------
        slices : generate fields in a surface slice given by the initializtion grid
        volumes : generate fields in volume given by the initializtion grid
        lines : generate fields along a line given by its end points
        )", 
	py::arg("basis"),
	py::arg("coefs"), py::arg("beg"), py::arg("end"),
	py::arg("num")=1000, py::arg("filename"), py::arg("dir")=".");

  f.def("file_slices", &Field::FieldGenerator::file_slices,
	R"(
        Write 2d field grids to files using the supplied string prefix. "

        Parameters
        ----------
        basis : Basis
            basis instance of any geometry; geometry will be deduced by the generator
        coefs : Coefs
            coefficient container instance
        filename : str
            file name for output
        dir : str, default='.'
            directory to write files

        Returns
        -------
        None

        Notes
        -----
	The files will be ascii or VTK rectangular grid files (if you compiled with VTK)

        See also
        --------
        slices : generate fields in a surface slice given by the initializtion grid
        volumes : generate fields in volume given by the initializtion grid
        lines : generate fields along a line given by its end points
        )",
	py::arg("basis"), py::arg("coefs"), py::arg("filename"),
	py::arg("dir")=".");


  f.def("volumes", [](FieldGenerator& A,
		      BasisClasses::BasisPtr basis, CoefClasses::CoefsPtr coefs)
  {
    std::map<double, std::map<std::string, py::array_t<float>>> ret;
    auto vols = A.volumes(basis, coefs);
    for (auto & v : vols) {
      for (auto & u : v.second) {
	ret[v.first][u.first] = make_ndarray<float>(u.second);
      }
    }

    return ret;
  },
	R"(
        Volume grids (3d numpy arrays) indexed by time and field type

        Parameters
        ----------
        basis : Basis
            basis instance of any geometry; geometry will be deduced by the generator
        coefs : Coefs
            coefficient container instance

        Returns
        -------
        dict({time: {field-name: numpy.ndarray})
            dictionary of times, field names, and data arrays

        Notes
        -----
        Each data array in the dictionary is a volume with the range given
        by the constructor.

        See also
        --------
        lines : generate fields along a line given by its end points
        slices : generate fields in a slice given by the initializtion grid
        )");

  f.def("file_volumes", &Field::FieldGenerator::file_volumes,
	R"(
	Write 3d field grids to files using the supplied string prefix

        Parameters
        ----------
        basis : Basis
            basis instance of any geometry; geometry will be deduced by the generator
        coefs : Coefs
            coefficient container instance
        filename : str
            file name for output
        dir : str, default='.'
            directory to write files

        Returns
        -------
        None

        Notes
        -----
	The files will be ascii or VTK rectangular grid files (if you ompiled with VTK

        See also
        --------
        file_lines : generate files with fields along a line given by its end points
        file_slices : generate files with fields along surfaces
	)",
	py::arg("basis"), py::arg("coefs"), py::arg("filename"),
	py::arg("dir")=".");
}
