#include <functional>

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <BiorthBasis.H>
#include <FieldBasis.H>

namespace py = pybind11;
#include <TensorToArray.H>

void BasisFactoryClasses(py::module &m)
{
  m.doc() =
    R"(
    BasisFactory class bindings
    
    This module provides a factory class that will create biorthogonal
    bases from input YAML configuration files.  Each basis can then be
    used to compute coefficients, provide field quantities such as
    forces and, together with the FieldGenerator, surfaces and fields for
    visualization.

    Five bases are currently implemented:
     1. SphericalSL, the Sturm-Liouiville spherical basis;
     2. Cylindrical, created created by computing empirical orthogonal functions
        over a densely sampled SphericalSL basis;
     3. FlatDisk, an EOF rotation of the finite Bessel basis; and
     4. Cube, a periodic cube basis whose functions are the Cartesian eigenfunctions
        of the Cartesian Laplacian: sines and cosines.
     5. VelocityBasis, a basis for computing the mean field velocity fields from
        a phase-space snapshot.

    Each of these bases take a YAML configuration file as input. These parameter
    lists are as subset of and have the same structure as thosed used by EXP. The
    factory and the individual constructors will check the parameters keys
    and warn of mismatches for safety.  See the EXP documentation and
    the pyEXP examples for more detail.  Other bases in EXP but not in
    pyEXP include those for cubic and slab geometries and other special-
    purpose bases such as the Hernquist, Clutton-Brock sphere and two-
    dimensional disk basis.  These will be made available in a future
    release if there is demand.  Note however that the Hernquist and
    Clutton-Brock spheres can be constructed using SphericalSL with a
    Hernquist of modified Plummer model as input.

    The primary functions of these basis classes are:
      1. To compute BFE coefficients from phase-space snapshots
         using the ParticleReader class. See help(pyEXP.read).
      2. To evaluate the fields from the basis and a coefficient
         object. See help(pyEXP.coefs) and help(pyEXP.field).

    Introspection
    -------------
    The two bases have a 'cacheInfo(str)' member that reports the
    parameters used to create the cached basis.  This may be used
    grab the parameters for creating a basis.  At this point, you
    must create the YAML configuration for the basis even if the
    basis is cached.  This is a safety and consistency feature that
    may be relaxed in a future version.

    Coefficient creation
    --------------------
    The Basis class creates coefficients from phase space with two
    methods: 'createFromReader()' and 'createFromArray()'.  The first
    uses a ParticleReader, see help(pyEXP.read), and the second uses
    arrays of mass and 3d position vectors.  Both methods take an
    optional center vector (default: 0, 0, 0).  You may also register
    and an optional boolean functor used to select which particles to
    using the 'setSelector(functor)' member.  An example functor
    would be defined in Python as follows:
       def myFunctor(m, pos, vel, index):
          ret = False  # Default return value
          # some caculation with scalar mass, pos array, vel array and
          # integer index that sets ret to True if desired . . . 
          return ret
    If you are using 'createFromArray()', you will only have access to
    the mass and position vector.   You may clear and turn off the
    selector using the 'clrSelector()' member.

    Scalablility
    ------------
    createFromArray() is a convenience method allows you to transform
    coordinates and preprocess phase space using your own methods and
    readers.  Inside this method are three member functions calls that
    separately initialize, accumulate the coefficient contributions from
    the provided vectors, and finally construct and return the new coeffi-
    cient instance (Coefs).  For scalability, we provide access to each 
    of these three methods so that the phase space may be partitioned into
    any number of smaller pieces.  These three members are: initFromArray(),
    addFromArray(), makeFromArray().  The initFromArray() is called once to
    begin the creation and the makeFromArray() method is called once to
    build the final set of coefficients.  The addFromArray() may be called
    any number of times in between.  For example, the addFromArray() call
    can be inside of a loop that iterates over any partition of phase space
    from your own pipeline.  The underlying computation is identical to
    createFromArray().  However, access to the three underlying steps allows
    you to scale your phase-space processing to snapshots of any size.
    For reference, the createFromReader() method uses a producer-consumer
    pattern internally to provide scalability.  These three methods allow
    you to provide the same pattern in your own pipeline.

    Coordinate systems
    -------------------
    Each basis is assigned a natural coordinate system for field evaluation
    as follows:
     1. SphericalSL uses spherical coordinates
     2. Cylindrical uses cylindrical coordinates
     3. FlatDisk uses cylindrical coordinates
     4. Cube uses Cartesian coordinates
    These default choices may be overridden by passing a string argument
    to the 'setFieldType()' member. The argument is case insensitive and only 
    distinguishing characters are necessary.  E.g. for 'Cylindrical', the 
    argument 'cyl' or even 'cy' is sufficient.  The argument 'c' is clearly 
    not enough.

    Orbit integration
    -----------------
    The IntegrateOrbits routine uses a fixed time step leap frog integrator
    to advance orbits from tinit to tfinal with time step h.  The initial
    positions and velocities are supplied in an nx6 NumPy array.  Tuples
    of the basis (a Basis instance) and coefficient database (a Coefs
    instance) for each component is supplied to IntegrateOrbtis as a list.
    Finally, the type of acceleration is an instance of the AccelFunc class.
    The acceleration at each time step is computed by setting a coefficient
    set in Basis and evaluating and accumulating the acceleration for each
    phase-space point.  The coefficient are handled by implementing the
    evalcoefs() method of AccelFunc. We supply two implemented derived
    classes, AllTimeFunc and SingleTimeFunc.  The first interpolates on the
    Coefs data base and installs the interpolated coefficients for the
    current time in the basis instance.  The SingleTimeFunc interpolates on
    the Coefs data base for a single fixed time and sets the interpolated
    coefficients once at the beginning of the integration.  This implementes
    a fixed potential model.  AccelFunc can be inherited by a native Python
    class and the evalcoefs() may be implemented in Python and passed to
    IntegrateOrbits in the same way as a native C++ class.
    )";

  using namespace BasisClasses;

  //! Need an alias to prevent the pybind11 macro from expanding the
  //! STL signature
  using DictMapStr = std::map<std::string, std::string>;

  class PyBasis : public Basis
  {
  protected:
    std::vector<double> sph_eval(double r, double costh, double phi) override
    {
      PYBIND11_OVERRIDE_PURE(std::vector<double>, Basis,
			     sph_eval, r, costh, phi);
    }
    
    std::vector<double> cyl_eval(double R, double z, double phi) override
    {
      PYBIND11_OVERRIDE_PURE(std::vector<double>, Basis,
			     cyl_eval, R, z, phi);
    }
    
    std::vector<double> crt_eval(double x, double y, double z) override
    {
      PYBIND11_OVERRIDE_PURE(std::vector<double>, Basis,
			     crt_eval, x, y, z);
    }
    
    void load_coefs(CoefClasses::CoefStrPtr coefs, double time) override {
      PYBIND11_OVERRIDE_PURE(void, Basis, load_coefs, coefs, time);
    }
    
    void set_coefs(CoefClasses::CoefStrPtr coefs) override {
      PYBIND11_OVERRIDE_PURE(void, Basis, set_coefs, coefs);
    }

    const std::string classname() override {
      PYBIND11_OVERRIDE_PURE(std::string, Basis, classname);
    }

    const std::string harmonic() override {
      PYBIND11_OVERRIDE_PURE(std::string, Basis, harmonic);
    }

    std::vector<std::string> getFieldLabels(const Coord ctype) override {
      PYBIND11_OVERRIDE_PURE(std::vector<std::string>, Basis, getFieldLabels, ctype);
    }


  public:
    // Inherit the constructors
    using BasisClasses::Basis::Basis;

    std::vector<double> getFields(double x, double y, double z) override
    {
      PYBIND11_OVERRIDE_PURE(std::vector<double>, Basis, getFields,
			     x, y, z);
    }

    void accumulate(double x, double y, double z, double mass) override {
      PYBIND11_OVERRIDE_PURE(void, Basis, accumulate, mass, x, y, z, mass);
    }

    void reset_coefs(void) override {
      PYBIND11_OVERRIDE_PURE(void, Basis, reset_coefs,);
    }

    void make_coefs(void) override {
      PYBIND11_OVERRIDE_PURE(void, Basis, make_coefs,);
    }

    virtual CoefClasses::CoefStrPtr
    createFromReader(PR::PRptr reader, std::vector<double> ctr) override {
      PYBIND11_OVERRIDE_PURE(CoefClasses::CoefStrPtr, Basis, createFromReader, reader, ctr);
    }


    virtual void initFromArray(std::vector<double> ctr) override {
      PYBIND11_OVERRIDE_PURE(void, Basis, initFromArray, ctr);
    }

    virtual void
    addFromArray(Eigen::VectorXd& m, RowMatrixXd& p, bool roundrobin) override {
      PYBIND11_OVERRIDE_PURE(void, Basis, addFromArray, m, p, roundrobin);
    }
  };

  class PyFieldBasis : public FieldBasis
  {
  protected:
    std::vector<double> sph_eval(double r, double costh, double phi) override
    {
      PYBIND11_OVERRIDE(std::vector<double>, FieldBasis,
			sph_eval, r, costh, phi);
    }
    
    std::vector<double> cyl_eval(double R, double z, double phi) override
    {
      PYBIND11_OVERRIDE(std::vector<double>, FieldBasis,
			cyl_eval, R, z, phi);
    }
    
    std::vector<double> crt_eval(double x, double y, double z) override
    {
      PYBIND11_OVERRIDE(std::vector<double> , FieldBasis,
			     crt_eval, x, y, z);
    }
    
    void load_coefs(CoefClasses::CoefStrPtr coefs, double time) override
    {
      PYBIND11_OVERRIDE(void, FieldBasis, load_coefs, coefs, time);
    }
    
    void set_coefs(CoefClasses::CoefStrPtr coefs) override
    {
      PYBIND11_OVERRIDE(void, FieldBasis, set_coefs, coefs);
    }

    const std::string classname() override
    {
      PYBIND11_OVERRIDE(std::string, FieldBasis, classname);
    }

    const std::string harmonic() override
    {
      PYBIND11_OVERRIDE(std::string, FieldBasis, harmonic);
    }

  public:
    // Inherit the constructors
    using FieldBasis::FieldBasis;

    std::vector<double> getFields(double x, double y, double z) override
    {
      PYBIND11_OVERRIDE(std::vector<double>, FieldBasis, getFields,
			x, y, z);
    }

    void accumulate(double m, double x, double y, double z,
		    double u, double v, double w) override
    {
      PYBIND11_OVERRIDE(void, FieldBasis, accumulate, m, x, y, z, u, v, w);
    }
    
    void reset_coefs(void) override {
      PYBIND11_OVERRIDE(void, FieldBasis, reset_coefs,);
    }

    void make_coefs(void) override {
      PYBIND11_OVERRIDE(void, FieldBasis, make_coefs,);
    }

  };

  class PyBiorthBasis : public BiorthBasis
  {
  protected:
    std::vector<double> sph_eval(double r, double costh, double phi) override
    {
      PYBIND11_OVERRIDE_PURE(std::vector<double>, BiorthBasis,
			sph_eval, r, costh, phi);
    }
    
    std::vector<double> cyl_eval(double R, double z, double phi) override
    {
      PYBIND11_OVERRIDE_PURE(std::vector<double>, BiorthBasis,
			cyl_eval, R, z, phi);
    }
    
    std::vector<double> crt_eval(double x, double y, double z) override
    {
      PYBIND11_OVERRIDE_PURE(std::vector<double> , BiorthBasis,
			crt_eval, x, y, z);
    }
    
    void load_coefs(CoefClasses::CoefStrPtr coefs, double time) override {
      PYBIND11_OVERRIDE_PURE(void, BiorthBasis, load_coefs, coefs, time);
    }
    
    const std::string classname() override {
      PYBIND11_OVERRIDE_PURE(std::string, BiorthBasis, classname);
    }

    const std::string harmonic() override {
      PYBIND11_OVERRIDE_PURE(std::string, BiorthBasis, harmonic);
    }

  public:
    // Inherit the constructors
    using BasisClasses::BiorthBasis::BiorthBasis;

    /*
    std::vector<double> getFields(double x, double y, double z) override
    {
      PYBIND11_OVERRIDE(std::vector<double>, BiorthBasis, getFields,
			     x, y, z);
    }
    */

    void accumulate(double x, double y, double z, double mass) override {
      PYBIND11_OVERRIDE_PURE(void, BiorthBasis, accumulate, x, y, z, mass);
    }

    void reset_coefs(void) override {
      PYBIND11_OVERRIDE_PURE(void, BiorthBasis, reset_coefs,);
    }

    void make_coefs(void) override {
      PYBIND11_OVERRIDE_PURE(void, BiorthBasis, make_coefs,);
    }

    void set_coefs(CoefClasses::CoefStrPtr coefs) override {
      PYBIND11_OVERRIDE_PURE(void, BiorthBasis, set_coefs, coefs);
    }

  };

  class PySphericalSL : public SphericalSL
  {
  protected:

    void load_coefs(CoefClasses::CoefStrPtr coefs, double time) override {
      PYBIND11_OVERRIDE(void, SphericalSL, load_coefs, coefs, time);
    }
    
    void set_coefs(CoefClasses::CoefStrPtr coefs) override {
      PYBIND11_OVERRIDE(void, SphericalSL, set_coefs, coefs);
    }

    const std::string classname() override {
      PYBIND11_OVERRIDE(std::string, SphericalSL, classname);
    }

    const std::string harmonic() override {
      PYBIND11_OVERRIDE(std::string, SphericalSL, harmonic);
    }

  public:

    // Inherit the constructors
    using SphericalSL::SphericalSL;

    std::vector<double> getFields(double x, double y, double z) override {
      PYBIND11_OVERRIDE(std::vector<double>, SphericalSL, getFields,
			x, y, z);
    }

    void accumulate(double x, double y, double z, double mass) override {
      PYBIND11_OVERRIDE(void, SphericalSL, accumulate, x, y, z, mass);
    }

    void reset_coefs(void) override {
      PYBIND11_OVERRIDE(void, SphericalSL, reset_coefs,);
    }

    void make_coefs(void) override {
      PYBIND11_OVERRIDE(void, SphericalSL, make_coefs,);
    }

  };

  class PyCylindrical : public Cylindrical
  {
  protected:

    void load_coefs(CoefClasses::CoefStrPtr coefs, double time) override {
      PYBIND11_OVERRIDE(void, Cylindrical, load_coefs, coefs, time);
    }
    
    void set_coefs(CoefClasses::CoefStrPtr coefs) override {
      PYBIND11_OVERRIDE(void, Cylindrical, set_coefs, coefs);
    }

    const std::string classname() override {
      PYBIND11_OVERRIDE(std::string, Cylindrical, classname);
    }

    const std::string harmonic() override {
      PYBIND11_OVERRIDE(std::string, Cylindrical, harmonic);
    }

  public:

    // Inherit the constructors
    using Cylindrical::Cylindrical;

    std::vector<double> getFields(double x, double y, double z) override {
      PYBIND11_OVERRIDE(std::vector<double>, Cylindrical, getFields,
			x, y, z);
    }

    void accumulate(double x, double y, double z, double mass) override {
      PYBIND11_OVERRIDE(void, Cylindrical, accumulate, x, y, z, mass);
    }

    void reset_coefs(void) override {
      PYBIND11_OVERRIDE(void, Cylindrical, reset_coefs,);
    }

    void make_coefs(void) override {
      PYBIND11_OVERRIDE(void, Cylindrical, make_coefs,);
    }

  };


  class PyFlatDisk : public FlatDisk
  {
  protected:

    std::vector<double> sph_eval(double r, double costh, double phi) override
    {
      PYBIND11_OVERRIDE(std::vector<double>, FlatDisk, sph_eval, r, costh, phi);
    }
    
    std::vector<double> cyl_eval(double R, double z, double phi) override
    {
      PYBIND11_OVERRIDE(std::vector<double>, FlatDisk, cyl_eval, R, z, phi);
    }
    
    std::vector<double> crt_eval(double x, double y, double z) override
    {
      PYBIND11_OVERRIDE(std::vector<double>, FlatDisk, crt_eval, x, y, z);
    }
    
    void load_coefs(CoefClasses::CoefStrPtr coefs, double time) override
    {
      PYBIND11_OVERRIDE(void, FlatDisk, load_coefs, coefs, time);
    }
    
    void set_coefs(CoefClasses::CoefStrPtr coefs) override
    {
      PYBIND11_OVERRIDE(void, FlatDisk, set_coefs, coefs);
    }

    const std::string classname() override {
      PYBIND11_OVERRIDE(std::string, FlatDisk, classname);
    }

    const std::string harmonic() override {
      PYBIND11_OVERRIDE(std::string, FlatDisk, harmonic);
    }

  public:

    // Inherit the constructors
    using FlatDisk::FlatDisk;

    std::vector<double> getFields(double x, double y, double z) override
    {
      PYBIND11_OVERRIDE(std::vector<double>, FlatDisk, getFields,
			x, y, z);
    }

    void accumulate(double x, double y, double z, double mass) override
    {
      PYBIND11_OVERRIDE(void, FlatDisk, accumulate, x, y, z, mass);
    }

    void reset_coefs(void) override
    {
      PYBIND11_OVERRIDE(void, FlatDisk, reset_coefs,);
    }

    void make_coefs(void) override
    {
      PYBIND11_OVERRIDE(void, FlatDisk, make_coefs,);
    }

  };


  class PyCube : public Cube
  {
  protected:

    std::vector<double> sph_eval(double r, double costh, double phi) override
    {
      PYBIND11_OVERRIDE(std::vector<double>, Cube, sph_eval, r, costh, phi);
    }
    

    std::vector<double> cyl_eval(double R, double z, double phi) override
    {
      PYBIND11_OVERRIDE(std::vector<double>, Cube, cyl_eval, R, z, phi);
    }
    
    std::vector<double> crt_eval(double x, double y, double z) override
    {
      PYBIND11_OVERRIDE(std::vector<double>, Cube, crt_eval, x, y, z);
    }
    
    void load_coefs(CoefClasses::CoefStrPtr coefs, double time) override
    {
      PYBIND11_OVERRIDE(void, Cube, load_coefs, coefs, time);
    }
    
    void set_coefs(CoefClasses::CoefStrPtr coefs) override
    {
      PYBIND11_OVERRIDE(void, Cube, set_coefs, coefs);
    }

    const std::string classname() override
    {
      PYBIND11_OVERRIDE(std::string, Cube, classname);
    }

    const std::string harmonic() override
    {
      PYBIND11_OVERRIDE(std::string, Cube, harmonic);
    }

  public:

    // Inherit the constructors
    using Cube::Cube;

    std::vector<double> getFields(double x, double y, double z) override
    {
      PYBIND11_OVERRIDE(std::vector<double>, Cube, getFields, x, y, z);
    }

    void accumulate(double x, double y, double z, double mass) override
    {
      PYBIND11_OVERRIDE(void, Cube, accumulate, x, y, z, mass);
    }

    void reset_coefs(void) override
    {
      PYBIND11_OVERRIDE(void, Cube, reset_coefs,);
    }

    void make_coefs(void) override
    {
      PYBIND11_OVERRIDE(void, Cube, make_coefs,);
    }

  };


  class PyAccelFunc : public AccelFunc
  {
  public:
    // Inherit the constructors
    using BasisClasses::AccelFunc::AccelFunc;

    // The coefficient evaluation member
    void evalcoefs(double t, BasisCoef mod) override {
      PYBIND11_OVERRIDE_PURE(void, AccelFunc, evalcoefs, t, mod);
    }
  };


  py::class_<BasisClasses::Basis, std::shared_ptr<BasisClasses::Basis>, PyBasis>
    (m, "Basis")
    .def("factory",            &BasisClasses::BiorthBasis::factory_string,
	 R"(
         Generate a basis from a YAML configuration supplied as a string

         Parameters
         ----------
         config : str
             the YAML config string

         Returns
         -------
         None
         )")
    .def("__call__", [](BasisClasses::Basis& A,	double x, double y, double z)
    {
      return A.evaluate(x, y, z);
    },
      R"(
         Evaluate the field at the given point, returning the tuple of
         field values and field labels as a pair of lists.

         Parameters
         ----------
         x, y, z : float values
             desired position for field evaluation

         Returns
         -------
         tuple of lists
             the field array and label array

         Note
         ----
         This is an experimental feature
         )"
      )
    .def("createFromArray",
	 [](BasisClasses::Basis& A, Eigen::VectorXd& mass, RowMatrixXd& ps,
	    double time, std::vector<double> center, bool roundrobin)
	 {
	   return A.createFromArray(mass, ps, time, center, roundrobin);
	 },
	 R"(
         Generate the coefficients from a mass and position array or,
	 phase-space array, time, and an optional expansion center location. 

         Parameters
         ----------
         mass : list
             vector containing the masses for the n particles
         ps   : numpy.ndarray
             an array with n rows and 6 columns (x, y, z, u, v, w)
             or 3 columns (x, y, z) for a biorthogonal basis
         roundrobin : bool
             the particles will be accumulated for each process 
             round-robin style with MPI by default.  This may be 
             disabled with 'roundrobin=false'.

         Returns
         -------
         CoefStruct
             the coefficient structure derived from the particles

         See also
         --------
         initFromArray : initialize for coefficient contributions
         addFromArray : add contribution for particles
         makeFromArray: create coefficients contributions
         )",
	 py::arg("mass"), py::arg("pos"), py::arg("time"),
	 py::arg("center") = std::vector<double>(3, 0.0),
	 py::arg("roundrobin") = true)
    .def("makeFromArray",
	 [](BasisClasses::Basis& A, double time)
	 {
	   return A.makeFromArray(time);
	 },
	 R"(
         Make the coefficients

         This is the final call in the initFromArray(), addFromArray()...
	 addFromArray()...makeFromArray() call sequence.

         Parameters
         ----------
         time : float
             snapshot time

         Returns
         -------
         CoefStructure
             the coefficient structure created from the particles

         See also
         --------
         initFromArray : initialize for coefficient contributions
         addFromArray : add contribution for particles
         )",
	 py::arg("time")
	 )
    .def("setSelector", &BasisClasses::Basis::setSelector,
	 R"(
         Register a Python particle selection functor. 

         Returns
         -------
         None

         See also
         --------
         clrSelector : clears the selection set here
         )")
    .def("clrSelector", &BasisClasses::Basis::clrSelector,
	 R"(
         Clear the previously registered particle selection functor

         Returns
         -------
         None
         )");


  py::class_<BasisClasses::BiorthBasis, std::shared_ptr<BasisClasses::BiorthBasis>, PyBiorthBasis, BasisClasses::Basis>
    (m, "BiorthBasis")
    .def(py::init<const std::string&>(),
	 R"(
         Initialize a biorthogonal basis

         Parameters
         ----------
         YAMLstring : str
             the serialized YAML configuration

         Returns
         -------
         Basis
             the Basis object
        )", py::arg("YAMLstring"))
    .def("createFromReader", &BasisClasses::BiorthBasis::createFromReader,
	 R"(
         Generate the coefficients from the supplied ParticleReader

         Parameters
         ----------
         reader : Particle reader
             the ParticleReader instance
         center : list, default=[0, 0, 0]
	     an optional expansion center location

         Returns
         -------
         CoefStruct
             the basis coefficients computed from the particles
         )",
	 py::arg("reader"), 
	 py::arg("center") = std::vector<double>(3, 0.0))
    .def("createFromArray",
	 [](BasisClasses::BiorthBasis& A, Eigen::VectorXd& mass, RowMatrixXd& pos,
	    double time, std::vector<double> center, bool roundrobin)
	 {
	   return A.createFromArray(mass, pos, time, center, roundrobin);
	 },
	 R"(
         Generate the coefficients from a mass and position array,
	 time, and an optional expansion center location. 

         Parameters
         ----------
         mass : list
             vector containing the masses for the n particles
         pos  : numpy.ndarray
             an array with n rows and 3 columns (x, y, z)
         roundrobin : bool
             the particles will be accumulated for each process 
             round-robin style with MPI by default.  This may be 
             disabled with 'roundrobin=false'.

         Returns
         -------
         CoefStruct
             the coefficient structure derived from the particles

         See also
         --------
         initFromArray : initialize for coefficient contributions
         addFromArray : add contribution for particles
         makeFromArray: create coefficients contributions
         )",
	 py::arg("mass"), py::arg("pos"), py::arg("time"),
	 py::arg("center") = std::vector<double>(3, 0.0),
	 py::arg("roundrobin") = true)
    .def("initFromArray",
	 [](BasisClasses::BiorthBasis& A, std::vector<double> center)
	 {
	   return A.initFromArray(center);
	 },
	 R"(
         Initialize coefficient accumulation

         Parameters
         ----------
         center : list, default=[0, 0, 0]
             vector of center positions

         Returns
         -------
         None

         Notes
         -----
	 After initialization, phase-space data is then added with 
         addFromArray() call.  addFromArray() may be called multiple times 
         with any unique partition of phase space. The final generation is 
         finished with a call to makeFromArray() with the snapshot time.  
         This final call returns the coefficient set. This sequence of 
         calls is identical to createFromArray() for a single set of 
         phase space arrays but allows for generation from very large 
         phase-space sets that can not be stored in physical memory.
         )",
	 py::arg("center") = std::vector<double>(3, 0.0))
    .def("addFromArray",
	 [](BasisClasses::BiorthBasis& A, Eigen::VectorXd& mass, RowMatrixXd& pos)
	 {
	   return A.addFromArray(mass, pos);
	 },
	 R"(
         Add particle contributions to coefficients

         Parameters
         ----------
         mass : list
             vector containing the masses for the n particles
         pos : numpy.ndarray
             an array with n rows and 3 columns (x, y, z)

         Returns
         -------
         None

         See also
         --------
         initFromArray : initialize for coefficient contributions
         makeFromArray: create coefficients contributions
         )",
	 py::arg("mass"), py::arg("pos"))
    .def("getFields", &BasisClasses::BiorthBasis::getFields,
	 R"(
         Return the density, potential, and forces for a cartesian position.

	 Field order is: dens0, potl0, dens, potl, fx, fy, fz. Dens0 and
	 potl0 are the fields evaluated for l=0 or m=0 and dens and potl
	 are evaluated for l>0 or m>0

         Parameters
         ----------
         x : float
             x-axis position
         y : float
             y-axis position
         z : float
             z-axis position

         Returns
         -------
         None
         )",
	 py::arg("x"), py::arg("y"), py::arg("z"))
    .def("setFieldType",       &BasisClasses::BiorthBasis::setFieldType,
         R"(
         Set the coordinate system for force evaluations.  The natural 
         coordinates for the basis class are the default; spherical
         coordinates for SphericalSL, cylindrical coordinates for
         Cylindrical and FlatDisk, and Cartesian coordinates for Cube.
         This member function can be used to override the default.  The
         available coorindates are: 'spherical', 'cylindrical', 'cartesian'.

         Parameters
         ----------
         coord : str
             the coordinate system

         Returns
         -------
         None
         )",
	 py::arg("coord"))
    .def("getFieldType",       &BasisClasses::BiorthBasis::getFieldType,
         R"(
         Get the coordinate system for force evaluations for inspection.

         Parameters
         ----------
         None

         Returns
         -------
         None
         )")
    .def("accumulate", [](BasisClasses::BiorthBasis& A, double x, double y, double z, double mass)
    {
      return A.accumulate(x, y, z, mass);
    },
      R"(
         Add the contribution of a single particle to the coefficients

         Parameters
         ----------
         x : float
             x-axis position
         y : float
             y-axis position
         z : float
             z-axis position
         mass : float
             particle mass

         Returns
         -------
         None
        )", 
	 py::arg("x"), py::arg("y"), py::arg("z"), py::arg("mass"))
    .def("getMass",            &BasisClasses::BiorthBasis::getMass,
	 R"(
         Return the total mass of particles contributing the current coefficient set

         Returns
         -------
         out : float
            total mass value
         )")
    .def("reset_coefs",        &BasisClasses::BiorthBasis::reset_coefs,
	 R"(
         Reset the coefficients to begin a generating a new set

         Returns
         -------
         None
         )")
    .def("make_coefs",         &BasisClasses::BiorthBasis::make_coefs,
	 R"(
         Create the coefficients after particle accumuluation is complete

         Returns
         -------
         None
         )")
    .def("set_coefs",          &BasisClasses::BiorthBasis::set_coefs,
	 R"(
         Install a new set of coefficients from a CoefStruct

         Parameters
         ----------
         coefs : CoefStruct

         Returns
         -------
         None
         )", py::arg("coefs"));

    py::class_<BasisClasses::SphericalSL, std::shared_ptr<BasisClasses::SphericalSL>, PySphericalSL, BasisClasses::BiorthBasis>(m, "SphericalSL")
      .def(py::init<const std::string&>(),
	 R"(
         Create a spherical Sturm-Liouville basis

         Parameters
         ----------
         YAMLstring : str
             The YAML configuration for the spherical basis

         Returns
         -------
         SphericalSL
              the new instance
         )", py::arg("YAMLstring"))

      .def("getBasis", &BasisClasses::SphericalSL::getBasis,
	   R"(
           Get basis functions

	   Evaluate the potential-density basis functions on a logarithmically
	   spaced grid for inspection. The structure is a two-grid of dimension
	   lmax by nmax each pointing to a dictionary of 1-d arrays ('potential',
	   'density', 'rforce') of dimension numr.

           Parameters
           ----------
           logxmin : float, default=-3.0
                minimum mapped radius in log10 units
           logxmax : float, default=0.5
                maximum mapped radius in log10 units
           numr : int, default=400
                number of equally spaced output points

           Returns
           -------
           list(list(dict))
               dictionaries of basis functions as lists indexed by l, n
           )",
	   py::arg("logxmin")=-3.0,
	   py::arg("logxmax")=0.5,
	   py::arg("numr")=400)
      // The following member needs to be a lambda capture because
      // orthoCheck is not in the base class and needs to have
      // different parameters depending on the basis type.  Here the
      // user can and will often need to specify a quadrature value.
      .def("orthoCheck", [](BasisClasses::SphericalSL& A, int knots)
      {
	return A.orthoCheck(knots);
      },
	R"(
        Check orthgonality of basis functions by quadrature

        Inner-product matrix of Sturm-Liouville solutions indexed by
        harmonic order used to assess fidelity.

        Parameters
        ----------
        knots : int, default=40
            Number of quadrature knots

        Returns
        -------
        list(numpy.ndarray)
	    list of numpy.ndarrays from [0, ... , Lmax]
        )",
	py::arg("knots")=40)
      .def_static("cacheInfo", [](std::string cachefile)
      {
	return BasisClasses::SphericalSL::cacheInfo(cachefile);
      },
	R"(
        Report the parameters in a basis cache file and return a dictionary

        Parameters
        ----------
        cachefile : str
            name of cache file

        Returns
        -------
        dict({tag: value},...)
            cache parameters
        )",
	py::arg("cachefile"));

  py::class_<BasisClasses::Cylindrical, std::shared_ptr<BasisClasses::Cylindrical>, PyCylindrical, BasisClasses::BiorthBasis>(m, "Cylindrical")
    .def(py::init<const std::string&>(),
	 R"(
	 Create a cylindrical EOF basis

         Parameters
         ----------
         YAMLstring : str
             The YAML configuration for the cylindrical basis

         Returns
         -------
         Cylindrical
              the new instance
         )", py::arg("YAMLstring"))
    .def("getBasis", &BasisClasses::Cylindrical::getBasis,
	 R"(

         Evaluate basis on grid for visualization

         Evaluate the potential-density basis functions on a linearly spaced
	 2d-grid for inspection.  The structure is a two-grid of dimension
	 lmax by nmax each pointing to a dictionary of 2-d arrays ('potential',
	 'density', 'rforce', 'zforce') of dimension numr X numz.

         Parameters
         ----------
         xmin : float, default=0.0
              minimum value in mapped radius
         xmax : float, default=1.0
              maximum value in mapped radius
         numr : int, default=40
              number of linearly-space evaluation points in radius
         zmin : float, default=-0.1
              minimum value in vertical height
         zmax : float, default=0.1
              maximum value in vertical height
         numz : int, default=40
              number of linearly-space evaluation points in height

         Returns
         -------
         list(list(dict))
             dictionaries of basis functions as lists indexed by m, n
         )",
	 py::arg("xmin")=0.0,
	 py::arg("xmax")=1.0,
	 py::arg("numr")=40,
	 py::arg("zmin")=-0.1,
	 py::arg("zmax")=0.1,
	 py::arg("numz")=40 )
    // The following member needs to be a lambda capture because
    // orthoCheck is not in the base class and needs to have different
    // parameters depending on the basis type.  Here, the quadrature
    // is determined by the scale of the meridional grid.
    .def("orthoCheck", [](BasisClasses::Cylindrical& A)
	 {
	   return A.orthoCheck();
	 },
	R"(
        Check orthgonality of basis functions by quadrature

        Inner-product matrix of Sturm-Liouville solutions indexed by
        harmonic order used to assess fidelity.

        Parameters
        ----------
        knots : int
            Number of quadrature knots

        Returns
        -------
        list(numpy.ndarray)
	    list of numpy.ndarrays from [0, ... , Mmax]
        )")
    .def_static("cacheInfo", [](std::string cachefile)
    {
      return BasisClasses::Cylindrical::cacheInfo(cachefile);
    },
      R"(
      Report the parameters in a basis cache file and return a dictionary

      Parameters
      ----------
      cachefile : str
          name of cache file

      Returns
      -------
      dict({tag: value},...)
          cache parameters
      )",
      py::arg("cachefile"));

  py::class_<BasisClasses::FlatDisk, std::shared_ptr<BasisClasses::FlatDisk>, PyFlatDisk, BasisClasses::BiorthBasis>(m, "FlatDisk")
    .def(py::init<const std::string&>(),
	 R"(
         Create a 2d disk basis

         Parameters
         ----------
         YAMLstring : str
             The YAML configuration for the razor-thin EOF basis.  The default 
             parameters will give an exponential disk with scale length of
             0.01 units. Set the disk scale length using the 'scale'  parameter.

         Returns
         -------
         FlatDisk
             the new instance
         )", py::arg("YAMLstring"))
    .def("getBasis", &BasisClasses::FlatDisk::getBasis,
	 R"(
         Evaluate the potential-density basis functions

         Returned functions will linearly spaced 2d-grid for inspection. The 
         min/max radii are given in log_10 units.  The structure is a two-grid of 
         dimension mmax by nmax each pointing to a dictionary of 1-d arrays 
         ('potential', 'density', 'rforce') of dimension numr.

         Parameters
         ----------
         logxmin : float, default=-4.0
             the minimum radius in log10 scaled units
         logxmax : float, default=-1.0
             the maximum radius in log10 scaled units
         numr : int
             the number of output evaluations

         Returns
         -------
         list(dict{str: numpy.ndarray})
             list of lists of dictionaries in harmonic and radial pointing to 
             density and potential basis functions
         )",
	 py::arg("logxmin")=-4.0,
	 py::arg("logxmax")=-1.0,
	 py::arg("numr")=400)
    // The following member needs to be a lambda capture because
    // orthoCheck is not in the base class and needs to have different
    // parameters depending on the basis type.  Here, the quadrature
    // is determined by the scale of the meridional grid.
    .def("orthoCheck", [](BasisClasses::FlatDisk& A)
    {
      return A.orthoCheck();
    },
      R"(
      Check orthgonality of basis functions by quadrature

      Inner-product matrix of Sturm-Liouville solutions indexed by
      harmonic order used to assess fidelity.

      Parameters
      ----------
      knots : int, default=40
          Number of quadrature knots

      Returns
      -------
      list(numpy.ndarray)
          list of numpy.ndarrays from [0, ... , Mmax]
       )"
      )
    .def_static("cacheInfo", [](std::string cachefile)
    {
      return BasisClasses::FlatDisk::cacheInfo(cachefile);
    },
      R"(
      Report the parameters in a basis cache file and return a dictionary

      Parameters
      ----------
      cachefile : str
          name of cache file

      Returns
      -------
      out : dict({tag: value})
          cache parameters
      )",
      py::arg("cachefile"));

  py::class_<BasisClasses::Cube, std::shared_ptr<BasisClasses::Cube>, PyCube, BasisClasses::BiorthBasis>(m, "Cube")
    .def(py::init<const std::string&>(),
	 R"(
         Create a 3d periodic cube basis

         Parameters
         ----------
         YAMLstring : str
             The YAML configuration for the periodic cube basis.  The coordinates
             are the unit cube with origin at (0, 0, 0) and maximum extent (1, 1, 1).
             The default  parameters will wave numbers between [-6,...,6] in each
             dimension.

         Returns
         -------
         Cube
             the new instance
         )", py::arg("YAMLstring"))
    .def("orthoCheck", [](BasisClasses::Cube& A)
    {
      return A.orthoCheck();
    },
      R"(
      Check orthgonality of basis functions by quadrature

      Inner-product matrix of indexed by flattened wave number (nx, ny, nz) where
      each of nx is in [-nmaxx, nmaxx], and so on for ny and nz.  Each dimension 
      has dx=2*nmaxx+1 wave numbers and similarly for dy and dz.  The index into the
      array is index=(nx+nmaxx)*dx*dy + (ny+nmaxy)*dy + (nz+nmaxz).   This is an 
      analyic basis so the orthogonality matrix is not a check of andy numerical
      computation other than the quadature itself.  It is included for completeness.

      Parameters
      ----------
      None

      Returns
      -------
      numpy.ndarray)
          list of numpy.ndarrays from [0, ... , dx*dy*dz]
       )"
      );


  py::class_<BasisClasses::FieldBasis, std::shared_ptr<BasisClasses::FieldBasis>, PyFieldBasis, BasisClasses::Basis>(m, "FieldBasis")
      .def(py::init<const std::string&>(),
	 R"(
         Create a orthogonal basis for representing general fields

         Parameters
         ----------
         YAMLstring : str
             The YAML configuration for the field basis

         Returns
         -------
         FieldBasis
              the new instance
         )", py::arg("YAMLstring"))
    .def("addPSFunction", &BasisClasses::FieldBasis::addPSFunction,
	 R"(
         Register a functor that returns a list of derived phase-space
         fields and a list of labels desribing them

         Parameters
         ----------
         function : [float,...] = f(float, [float, float, float], [float, float, float]
             Returns a list of n float values derived from mass, position,
             and velocity of the following form:
             [a, b, c] = function(mass, [x, y, z], [u, v, w])
         labels : list of str
	     labels for the fields, e.g. ['a', 'b', 'c']

         Returns
         -------
         None
         )",
	 py::arg("function"), py::arg("labels"))
    .def("initFromArray",
	 [](BasisClasses::FieldBasis& A, std::vector<double> center)
	 {
	   return A.initFromArray(center);
	 },
	 R"(
         Initialize coefficient accumulation

         Parameters
         ----------
         center : list, default=[0, 0, 0]
             vector of center positions

         Returns
         -------
         None

         Notes
         -----
	 After initialization, phase-space data is then added with 
         addFromArray() call.  addFromArray() may be called multiple times 
         with any unique partition of phase space. The final generation is 
         finished with a call to makeFromArray() with the snapshot time.  
         This final call returns the coefficient set. This sequence of 
         calls is identical to createFromArray() for a single set of 
         phase space arrays but allows for generation from very large 
         phase-space sets that can not be stored in physical memory.
         )",
	 py::arg("center") = std::vector<double>(3, 0.0))
    .def("addFromArray",
	 [](BasisClasses::FieldBasis& A, Eigen::VectorXd& mass, RowMatrixXd& ps)
	 {
	   return A.addFromArray(mass, ps);
	 },
	 R"(
         Add particle contributions to coefficients

         Parameters
         ----------
         mass : list
             vector containing the masses for the n particles
         pos : numpy.ndarray
             an array with n rows and 6 columns (x, y, z, u, v, w)

         Returns
         -------
         None

         See also
         --------
         initFromArray : initialize for coefficient contributions
         makeFromArray: create coefficients contributions
         )",
	 py::arg("mass"), py::arg("pos"))
    .def("makeFromArray",
	 [](BasisClasses::FieldBasis& A, double time)
	 {
	   return A.makeFromArray(time);
	 },
	 R"(
         Make the coefficients

         This is the final call in the initFromArray(), addFromArray()...
	 addFromArray()...makeFromArray() call sequence.

         Parameters
         ----------
         time : float
             snapshot time

         Returns
         -------
         CoefStructure
             the coefficient structure created from the particles

         See also
         --------
         initFromArray : initialize for coefficient contributions
         addFromArray : add contribution for particles
         )",
	 py::arg("time")
	 )
      .def("getBasis", &BasisClasses::FieldBasis::getBasis,
	   R"(
           Get basis functions

	   Evaluate the orthogonal basis functions on a logarithmically
	   spaced grid for inspection. The structure is a two-grid of dimension
	   lmax by nmax each pointing to a dictionary of 1-d arrays for the 
           velocity fields depending on the geometry (dof=2->cylindrical, 
           dof=3->spherical)

           Parameters
           ----------
           logxmin : float, default=-3.0
                minimum mapped radius in log10 units
           logxmax : float, default=0.5
                maximum mapped radius in log10 units
           numr : int, default=400
                number of equally spaced output points

           Returns
           -------
           list(list(dict))
               dictionaries of basis functions as lists indexed by l, n
           )",
	   py::arg("logxmin")=-3.0,
	   py::arg("logxmax")=0.5,
	   py::arg("numr")=400)
      // The following member needs to be a lambda capture because
      // orthoCheck is not in the base class and needs to have
      // different parameters depending on the basis type.  Here the
      // user can and will often need to specify a quadrature value.
      .def("orthoCheck", [](BasisClasses::FieldBasis& A)
      {
	return A.orthoCheck();
      },
	R"(
        Check orthgonality of basis functions by quadrature

        Inner-product matrix of orthogonal functions

        Parameters
        ----------
        None

        Returns
        -------
        numpy.ndarray
	    orthogonality matrix
        )"
	);

  py::class_<BasisClasses::VelocityBasis, std::shared_ptr<BasisClasses::VelocityBasis>, BasisClasses::FieldBasis>(m, "VelocityBasis")
      .def(py::init<const std::string&>(),
	 R"(
         Create a orthogonal velocity-field basis

         Parameters
         ----------
         YAMLstring : str
             The YAML configuration for the velocity basis

         Returns
         -------
         VelocityBasis
              the new instance

         Notes
         -----
         This is a FieldBasis specialized to return velocity fields in
         spherical, cylindrical, or Cartesian systems.
         )", py::arg("YAMLstring"));

  py::class_<BasisClasses::AccelFunc, std::shared_ptr<BasisClasses::AccelFunc>, PyAccelFunc>(m, "AccelFunc")
    .def(py::init<>(),
	 R"(
         Create a acceleration functor (AccelFunc instance) for a basis component

         Returns
         -------
         AccelFunc
         )")
    .def("F", &BasisClasses::AccelFunc::F,
	 R"(
         Computes and returns the acceleration array

         Parameters
         ----------
         time : float
             evaluation time
         ps : numpy.ndarray
             n x 6 phase space array
         accel : numpy.ndarray
             n x 6 array of accelerations
         mod : BasisCoef = tuple(Basis, Coefs)
             model description

         Returns
         -------
         accel : numpy.ndarray
             n x 6 array of accelerations
         )",
	 py::arg("time"), py::arg("ps"), py::arg("accel"), py::arg("mod"));

  py::class_<BasisClasses::AllTimeAccel, std::shared_ptr<BasisClasses::AllTimeAccel>, BasisClasses::AccelFunc>(m, "AllTimeAccel")
    .def(py::init<>(),
	 R"(
         AccelFunc instance that interpolates coefficients from the Coefs database for every time

         Returns
         -------
         AllTimeAccel : AccelFunc

         See also
         --------
         AccelFunc
         )");

  py::class_<BasisClasses::SingleTimeAccel, std::shared_ptr<BasisClasses::SingleTimeAccel>, BasisClasses::AccelFunc>(m, "SingleTimeAccel")
    .def(py::init<double, std::vector<BasisClasses::BasisCoef>>(),
	 R"(
         AccelFunc instance that uses a single time Coefs database

         Parameters
         ----------
         time : float
             evaluation time
         mod : BasisCoef = tuple(Basis, Coefs)
             model description

         Returns
         -------
         SingleTimeAccel : AccelFunc

         See also
         --------
         AccelFunc
         AllTimeAccel
         )", py::arg("time"), py::arg("mod"));
  
  m.def("IntegrateOrbits", 
	[](double tinit, double tfinal, double h, Eigen::MatrixXd ps,
	   std::vector<BasisClasses::BasisCoef> bfe,
	   BasisClasses::AccelFunc& func, int stride)
	{
	  Eigen::VectorXd T;
	  Eigen::Tensor<float, 3> O;

	  AccelFunctor F = [&func](double t, Eigen::MatrixXd& ps, Eigen::MatrixXd& accel, BasisCoef mod)->Eigen::MatrixXd& { return func.F(t, ps, accel, mod);};

	  std::tie(T, O) =
	    BasisClasses::IntegrateOrbits(tinit, tfinal, h, ps, bfe, F, stride);

	  py::array_t<float> ret = make_ndarray<float>(O);
	  return std::tuple<Eigen::VectorXd, py::array_t<float>>(T, ret);
	},
	R"(
        Compute particle orbits in gravitational field from the bases

        Integrate a list of initial conditions from tinit to tfinal with a
	step size of h using the list of basis and coefficient pairs. Every
	step will be included in return unless you provide an explicit
	value for 'nout', the number of desired output steps.  This will
	choose the 'nout' points closed to the desired time.

        Parameters
        ----------
        tinit : float
            the intial time
        tfinal : float
            the final time
        h : float
            the integration step size
        ps : numpy.ndarray
            an n x 6 table of phase-space initial conditions
        bfe : list(BasisCoef)
            a list of BFE coefficients used to generate the gravitational 
            field
        func : AccelFunctor
            the force function
        nout : int 
            the number of output intervals

        Returns
        -------
        tuple(numpy.array, numpy.ndarray)
            time and phase-space arrays
        )",
	py::arg("tinit"), py::arg("tfinal"), py::arg("h"),
	py::arg("ps"), py::arg("basiscoef"), py::arg("func"),
	py::arg("nout")=std::numeric_limits<int>::max());

}
