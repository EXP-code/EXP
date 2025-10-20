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
    forces and, together with the FieldGenerator, surfaces and fields
    for visualization.

    Eight bases are currently implemented:

     1. SphericalSL, the Sturm-Liouiville spherical basis;

     2. Bessel, the classic spherical biorthogonal constructed from
        the eigenfunctions of the spherical Laplacian;

     3. Cylindrical, created created by computing empirical orthogonal functions
        over a densely sampled SphericalSL basis;

     4. FlatDisk, an EOF rotation of the finite Bessel basis;

     5. CBDisk, the Clutton-Brock disk basis for testing;

     6. Slab, a biorthogonal basis for a slab geometry with a finite
        finite vertical extent.  The basis is constructed from direct
        solution of the Sturm-Liouville equation.

     7. Cube, a periodic cube basis whose functions are the Cartesian
        eigenfunctions of the Cartesian Laplacian: sines and cosines.

     8. FieldBasis, for computing user-provided quantities from a
        phase-space snapshot.

     9. VelocityBasis, for computing the mean field velocity fields from
        a phase-space snapshot.  This is a specialized version of FieldBasis.

    Each of these bases take a YAML configuration file as input. These parameter
    lists are as subset of and have the same structure as thosed used by EXP.
    The factory and the individual constructors will check the parameters keys
    and warn of mismatches for safety.  See the EXP documentation and the pyEXP
    examples for more detail.  The first four bases are the most often used bi-
    orthogonal basis types used for computing the potential and forces from
    density distributions.  Other biorthgonal bases in EXP but not in pyEXP
    include those for cubic and slab geometries and other special-purpose bases
    such as the Hernquist, Clutton-Brock sphere and two-dimensional disk basis.
    These will be made available in a future release if there is demand.  Note
    that the Hernquist and Clutton-Brock spheres can be constructed using
    SphericalSL with a Hernquist of modified Plummer model as input.  The
    FieldBasis and VelocityBasis are designed for producing summary data for
    post-production analysis (using mSSA or eDMD, for example) and for
    simulation cross-comparison.

    The primary functions of these basis classes are:

      1. To compute BFE coefficients from phase-space snapshots
         using the ParticleReader class. See help(pyEXP.read).

      2. To evaluate the fields from the basis and a coefficient
         object. See help(pyEXP.coefs) and help(pyEXP.field).

      3. To provide compact summary field data for post-production
         analysis.  See help(pyEXP.basis.FieldBasis) and
         help(pyEXP.basis.VelocityBasis).

    Introspection
    -------------
    The first two bases have a 'cacheInfo(str)' member that reports the
    parameters used to create the cached basis.  This may be used to
    grab the parameters for creating a basis.  Cache use ensures that
    your analyses are computed with the same bases used in a simulation
    or with the same basis used on previous pyEXP invocations.  At this
    point, you must create the YAML configuration for the basis even if
    the basis is cached.  This is a safety and consistency feature that
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

    The FieldBasis class requires a user-specified phase-space field
    functor that produces an list of quantities derived from the
    phase space for each particle.  For example, to get a total
    velocity field, we could use:

       def totalVelocity(m, pos, vel):
          # Some caculation with scalar mass, pos array, vel array.
          # Total velocity for this example...
          return [(vel[0]**2 + vel[1]**2 + vel[2]**2)**0.5]

    This function is registered with the FieldBasis using:

       basis->addPSFunction(totalVelocity, ['total velocity'])

    The VelocityBasis is a FieldBasis that automatically sets the
    phase-space field functor to cylindrical or spherical velocities
    based on the 'dof' parameter.  More on 'dof' below.

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

     4. CBDisk uses cylindrical coordinates

     5. Slab uses Cartesian coordinates

     6. Cube uses Cartesian coordinates

     7. FieldBasis and VelocityBasis provides two natural geometries for
        field evaluation: a two-dimensional (dof=2) polar disk and a
        three-dimensional (dof=3) spherical geometry that are chosen using
        the 'dof' parameter.  These use cylindrical and spherical
        coordinates, respectively, by default.

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

    Non-inertial frames of reference
    --------------------------------
    Each component of a multiple component simulation may have its own expansion
    center. Orbit integration in the frame of reference of the expansion is
    accomplished by defining a moving frame of reference using the setNonInertial()
    call with either an array of n times and center positions (as an nx3 array)
    or by initializing with an EXP orient file.

    We provide a member function, setNonInertialAccel(t), to estimate the frame
    acceleration at a given time.  This may be useful for user-defined acceleration
    routines.  This is automatically called default C++ evalcoefs() routine.
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
      PYBIND11_OVERRIDE(std::vector<double>, Basis, getFields, x, y, z);
    }

    using FCReturn = std::tuple<std::map<std::string, Eigen::VectorXd>,
				Eigen::VectorXd>;

    FCReturn getFieldsCoefs
    (double x, double y, double z, CoefClasses::CoefsPtr coefs) override
    {
      PYBIND11_OVERRIDE(FCReturn, Basis, getFieldsCoefs, x, y, z, coefs);
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
    createFromReader(PR::PRptr reader, Eigen::Vector3d ctr,
		     RowMatrix3d rot) override {
      PYBIND11_OVERRIDE_PURE(CoefClasses::CoefStrPtr, Basis, createFromReader, reader, ctr, rot);
    }


    void initFromArray
    (Eigen::Vector3d ctr, RowMatrix3d rot) override {
      PYBIND11_OVERRIDE_PURE(void, Basis, initFromArray, ctr, rot);
    }

    virtual void
    addFromArray(Eigen::VectorXd& m, RowMatrixXd& p, bool roundrobin, bool posvelrows) override {
      PYBIND11_OVERRIDE_PURE(void, Basis, addFromArray, m, p, roundrobin, posvelrows);
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
      PYBIND11_OVERRIDE(std::vector<double>, FieldBasis, getFields, x, y, z);
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
    
    void load_coefs(CoefClasses::CoefStrPtr coefs, double time) override
    {
      PYBIND11_OVERRIDE_PURE(void, BiorthBasis, load_coefs, coefs, time);
    }
    
    const std::string classname() override {
      PYBIND11_OVERRIDE_PURE(std::string, BiorthBasis, classname);
    }

    const std::string harmonic() override {
      PYBIND11_OVERRIDE_PURE(std::string, BiorthBasis, harmonic);
    }

    void computeAccel(double x, double y, double z,
		      Eigen::Ref<Eigen::Vector3d> v) override
    {
      PYBIND11_OVERRIDE_PURE(void, BiorthBasis, computeAccel, x, y, z, v);
    }

  public:
    // Inherit the constructors
    using BasisClasses::BiorthBasis::BiorthBasis;

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

  class PySpherical : public Spherical
  {
  protected:

    void load_coefs(CoefClasses::CoefStrPtr coefs, double time) override {
      PYBIND11_OVERRIDE(void, Spherical, load_coefs, coefs, time);
    }
    
    void set_coefs(CoefClasses::CoefStrPtr coefs) override {
      PYBIND11_OVERRIDE(void, Spherical, set_coefs, coefs);
    }

    const std::string classname() override {
      PYBIND11_OVERRIDE(std::string, Spherical, classname);
    }

    const std::string harmonic() override {
      PYBIND11_OVERRIDE(std::string, Spherical, harmonic);
    }

    void get_pot(Eigen::MatrixXd& tab, double x) override {
      PYBIND11_OVERRIDE_PURE(void, Spherical, get_pot, tab, x);
    }
    
    void get_dens(Eigen::MatrixXd& tab, double x) override {
      PYBIND11_OVERRIDE_PURE(void, Spherical, get_dens, tab, x);
    }

    void get_force(Eigen::MatrixXd& tab, double x) override {
      PYBIND11_OVERRIDE_PURE(void, Spherical, get_force, tab, x);
    }

    void computeAccel(double x, double y, double z,
		      Eigen::Ref<Eigen::Vector3d> v) override
    {
      PYBIND11_OVERRIDE(void, Spherical, computeAccel, x, y, z, v);
    }

  public:

    // Inherit the constructors
    using Spherical::Spherical;

    std::vector<double> getFields(double x, double y, double z) override {
      PYBIND11_OVERRIDE(std::vector<double>, Spherical, getFields, x, y, z);
    }

    void accumulate(double x, double y, double z, double mass) override {
      PYBIND11_OVERRIDE(void, Spherical, accumulate, x, y, z, mass);
    }

    void reset_coefs(void) override {
      PYBIND11_OVERRIDE(void, Spherical, reset_coefs,);
    }

    void make_coefs(void) override {
      PYBIND11_OVERRIDE(void, Spherical, make_coefs,);
    }

    std::vector<Eigen::MatrixXd> orthoCheck(int knots) override {
      PYBIND11_OVERRIDE_PURE(std::vector<Eigen::MatrixXd>, Spherical, orthoCheck, knots);
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

    void computeAccel(double x, double y, double z,
		      Eigen::Ref<Eigen::Vector3d> v) override
    {
      PYBIND11_OVERRIDE(void, Cylindrical, computeAccel, x, y, z, v);
    }


  public:

    // Inherit the constructors
    using Cylindrical::Cylindrical;

    std::vector<double> getFields(double x, double y, double z) override {
      PYBIND11_OVERRIDE(std::vector<double>, Cylindrical, getFields, x, y, z);
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

    void computeAccel(double x, double y, double z,
		      Eigen::Ref<Eigen::Vector3d> v) override
    {
      PYBIND11_OVERRIDE(void, FlatDisk, computeAccel, x, y, z, v);
    }


  public:

    // Inherit the constructors
    using FlatDisk::FlatDisk;

    std::vector<double> getFields(double x, double y, double z) override
    {
      PYBIND11_OVERRIDE(std::vector<double>, FlatDisk, getFields, x, y, z);
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


  class PyCBDisk : public CBDisk
  {
  protected:

    std::vector<double> sph_eval(double r, double costh, double phi) override
    {
      PYBIND11_OVERRIDE(std::vector<double>, CBDisk, sph_eval, r, costh, phi);
    }
    
    std::vector<double> cyl_eval(double R, double z, double phi) override
    {
      PYBIND11_OVERRIDE(std::vector<double>, CBDisk, cyl_eval, R, z, phi);
    }
    
    std::vector<double> crt_eval(double x, double y, double z) override
    {
      PYBIND11_OVERRIDE(std::vector<double>, CBDisk, crt_eval, x, y, z);
    }
    
    void load_coefs(CoefClasses::CoefStrPtr coefs, double time) override
    {
      PYBIND11_OVERRIDE(void, CBDisk, load_coefs, coefs, time);
    }
    
    void set_coefs(CoefClasses::CoefStrPtr coefs) override
    {
      PYBIND11_OVERRIDE(void, CBDisk, set_coefs, coefs);
    }

    const std::string classname() override {
      PYBIND11_OVERRIDE(std::string, CBDisk, classname);
    }

    const std::string harmonic() override {
      PYBIND11_OVERRIDE(std::string, CBDisk, harmonic);
    }

    void computeAccel(double x, double y, double z,
		      Eigen::Ref<Eigen::Vector3d> v) override
    {
      PYBIND11_OVERRIDE(void, CBDisk, computeAccel, x, y, z, v);
    }

  public:

    // Inherit the constructors
    using CBDisk::CBDisk;

    std::vector<double> getFields(double x, double y, double z) override
    {
      PYBIND11_OVERRIDE(std::vector<double>, CBDisk, getFields, x, y, z);
    }

    void accumulate(double x, double y, double z, double mass) override
    {
      PYBIND11_OVERRIDE(void, CBDisk, accumulate, x, y, z, mass);
    }

    void reset_coefs(void) override
    {
      PYBIND11_OVERRIDE(void, CBDisk, reset_coefs,);
    }

    void make_coefs(void) override
    {
      PYBIND11_OVERRIDE(void, CBDisk, make_coefs,);
    }

  };


  class PySlab : public Slab
  {
  protected:

    std::vector<double> sph_eval(double r, double costh, double phi) override
    {
      PYBIND11_OVERRIDE(std::vector<double>, Slab, sph_eval, r, costh, phi);
    }
    

    std::vector<double> cyl_eval(double R, double z, double phi) override
    {
      PYBIND11_OVERRIDE(std::vector<double>, Slab, cyl_eval, R, z, phi);
    }
    
    std::vector<double> crt_eval(double x, double y, double z) override
    {
      PYBIND11_OVERRIDE(std::vector<double>, Slab, crt_eval, x, y, z);
    }
    
    void load_coefs(CoefClasses::CoefStrPtr coefs, double time) override
    {
      PYBIND11_OVERRIDE(void, Slab, load_coefs, coefs, time);
    }
    
    void set_coefs(CoefClasses::CoefStrPtr coefs) override
    {
      PYBIND11_OVERRIDE(void, Slab, set_coefs, coefs);
    }

    const std::string classname() override
    {
      PYBIND11_OVERRIDE(std::string, Slab, classname);
    }

    const std::string harmonic() override
    {
      PYBIND11_OVERRIDE(std::string, Slab, harmonic);
    }

    void computeAccel(double x, double y, double z,
		      Eigen::Ref<Eigen::Vector3d> v) override
    {
      PYBIND11_OVERRIDE(void, Slab, computeAccel, x, y, z, v);
    }

  public:

    // Inherit the constructors
    using Slab::Slab;

    std::vector<double> getFields(double x, double y, double z) override
    {
      PYBIND11_OVERRIDE(std::vector<double>, Slab, getFields, x, y, z);
    }

    void accumulate(double x, double y, double z, double mass) override
    {
      PYBIND11_OVERRIDE(void, Slab, accumulate, x, y, z, mass);
    }

    void reset_coefs(void) override
    {
      PYBIND11_OVERRIDE(void, Slab, reset_coefs,);
    }

    void make_coefs(void) override
    {
      PYBIND11_OVERRIDE(void, Slab, make_coefs,);
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

    void computeAccel(double x, double y, double z,
		      Eigen::Ref<Eigen::Vector3d> v) override
    {
      PYBIND11_OVERRIDE(void, Cube, computeAccel, x, y, z, v);
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
    .def("factory", &BasisClasses::BiorthBasis::factory_string,
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
         tuple of numpy.ndarray and list of labels
             the field array and label array

         Note
         ----
         This is an experimental feature
         )"
      )
    .def("createFromArray", &BasisClasses::Basis::createFromArray,
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
         center : numpy.ndarray
             the center of expansion for the basis functions
         rotation : numpy.ndarray
             the rotation matrix for the basis functions
             (default is identity matrix)
         roundrobin : bool
             the particles will be accumulated for each process 
             round-robin style with MPI by default.  This may be 
             disabled with 'roundrobin=False'.
         posvelrows : bool
             positions (and optionally velocities) will be packed
             in rows instead of columns.  This accommodates the numpy
             construction [xpos, ypos, zpos] where xpos, ypos, zpos are
             arrays.  Defaults to True.

         Returns
         -------
         CoefStruct
             the coefficient structure derived from the particles

         See also
         --------
         initFromArray : initialize for coefficient contributions
         addFromArray  : add contribution for particles
         makeFromArray : create coefficients contributions
         )",
	 py::arg("mass"), py::arg("pos"), py::arg("time"),
	 py::arg("center") = Eigen::Vector3d::Zero(),
	 py::arg("rotation") = RowMatrix3d::Identity(),
	 py::arg("roundrobin") = true, py::arg("posvelrows") = true)
    .def("createFromArray",
	 [](BasisClasses::Basis& A, Eigen::VectorXd& mass, RowMatrixXd& ps,
	    double time, Eigen::Vector3d center, bool roundrobin, bool posvelrows)
	 {
	   return A.createFromArray(mass, ps, time, center, RowMatrix3d::Identity(),
				    roundrobin, posvelrows);
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
         center : numpy.ndarray
             the center of expansion for the basis functions
         roundrobin : bool
             the particles will be accumulated for each process 
             round-robin style with MPI by default.  This may be 
             disabled with 'roundrobin=False'.
         posvelrows : bool
             positions (and optionally velocities) will be packed
             in rows instead of columns.  This accommodates the numpy
             construction [xpos, ypos, zpos] where xpos, ypos, zpos are
             arrays.  Defaults to True.

         Returns
         -------
         CoefStruct
             the coefficient structure derived from the particles

         See also
         --------
         initFromArray : initialize for coefficient contributions
         addFromArray  : add contribution for particles
         makeFromArray : create coefficients contributions
         )",
	 py::arg("mass"), py::arg("pos"), py::arg("time"),
	 py::arg("center") = Eigen::Vector3d::Zero(),
	 py::arg("roundrobin") = true, py::arg("posvelrows") = true)
    .def("makeFromArray", &BasisClasses::Basis::makeFromArray,
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
         CoefStruct
             the coefficient structure created from the particles

         See also
         --------
         initFromArray : initialize for coefficient contributions
         addFromArray  : add contribution for particles
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
         )")
    .def("getFieldLabels",
	 [](BasisClasses::Basis& A)
	 {
	   return A.getFieldLabels();
	 },
	 R"(
         Provide the field labels for the basis functions

         Parameters
         ----------
         None

         Returns
         -------
         list: str
           list of basis function labels
         )"
	 )
    .def("setInertial", &BasisClasses::Basis::setInertial,
	 R"(
         Reset to inertial coordinates

         Parameters
         ----------
         None

         Returns
         -------
         None

         See also
         --------
         setNonInertial : set non-inertial data
         setNonInertialAccel : set the non-inertial acceleration
         )"
	 )
    .def("setNonInertial",
	 [](BasisClasses::Basis& A,
	    int N, const Eigen::VectorXd& times, const Eigen::MatrixXd& pos) {
	   A.setNonInertial(N, times, pos);
	 },
	 R"(
         Initialize for pseudo-force computation with a time series of positions
         using (1) a time vector and (2) a center position matrix with rows of three
         vectors

         Parameters
         ----------
         N : int
           number of previous positions to use for quadratic fit
         times : list or numpy.ndarray
           list of time points
         pos : numpy.ndarray
           an array with N rows and 3 columns of center positions

         Returns
         -------
         None

         See also
         --------
         setNonInertial : set non-inertial data from an Orient file
         setNonInertialAccel : set the non-inertial acceration
         )",
	 py::arg("N"), py::arg("times"), py::arg("pos")
         )
    .def("setNonInertial",
	 [](BasisClasses::Basis& A, int N, const std::string orient)
	 {
	   A.setNonInertial(N, orient);
	 },
	 R"(
         Initialize for pseudo-force computation with a time series of positions
         using a EXP orient file

         Parameters
         ----------
         N : int
           number of previous positions to use for quadratic fit
         orient : str
           name of the orient file

         Returns
         -------
         None

         See also
         --------
         setNonInertial : set non-inertial data from a time series of values
         setNonInertialAccel : set the non-inertial acceration
         )",
	 py::arg("N"), py::arg("orient")
         )
    .def("setNonInertialAccel", &BasisClasses::Basis::setNonInertialAccel,
	 R"(
         Set the pseudo acceleration for the non-inertial data at a given time

         Parameters
         ----------
         time : float
           evaluation time

         Returns
         -------
         None

         See also
         --------
         setNonInertial : set non-inertial data from a time series of values
         setNonInertial : set non-inertial data from an EXP orient file
         )",
	 py::arg("time")
         );
    
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
         BiorthBasis
             the BiorthBasis object

         Notes
         -----
         Needed for copying objects in the Python interpreter.  This can not be
         instantiated directly.
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
         rotation : numpy.ndarray, default=Identity
             rotation matrix to apply to the phase-space coordinates

         Returns
         -------
         CoefStruct
             the basis coefficients computed from the particles
         )",
	 py::arg("reader"), 
	 py::arg("center") = Eigen::Vector3d::Zero(),
	 py::arg("rotation") = RowMatrix3d::Identity())
    .def("createFromArray", &BasisClasses::BiorthBasis::createFromArray,
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
         posvelrows : bool
             positions (and optionally velocities) will be packed
             in rows instead of columns.  This accommodates the numpy
             construction [xpos, ypos, zpos] where xpos, ypos, zpos are
             arrays.  Defaults to True.

         Returns
         -------
         CoefStruct
             the coefficient structure derived from the particles

         See also
         --------
         initFromArray : initialize for coefficient contributions
         addFromArray  : add contribution for particles
         makeFromArray : create coefficients contributions
         )",
	 py::arg("mass"), py::arg("pos"), py::arg("time"),
	 py::arg("center") = Eigen::Vector3d::Zero(),
	 py::arg("rotation") = RowMatrix3d::Identity(),
	 py::arg("roundrobin") = true, py::arg("posvelrows") = true)
    .def("createFromArray",
	 [](BasisClasses::BiorthBasis& A, Eigen::VectorXd& mass, RowMatrixXd& pos,
	    double time, Eigen::Vector3d center, bool roundrobin, bool posvelrows)
	 {
	   return A.createFromArray(mass, pos, time, center, RowMatrix3d::Identity(),
				    roundrobin, posvelrows);
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
      	 center : numpy.ndarray
             the center of expansion for the basis functions
         roundrobin : bool
             the particles will be accumulated for each process 
             round-robin style with MPI by default.  This may be 
             disabled with 'roundrobin=false'.
         posvelrows : bool
             positions (and optionally velocities) will be packed
             in rows instead of columns.  This accommodates the numpy
             construction [xpos, ypos, zpos] where xpos, ypos, zpos are
             arrays.  Defaults to True.

         Returns
         -------
         CoefStruct
             the coefficient structure derived from the particles

         See also
         --------
         initFromArray : initialize for coefficient contributions
         addFromArray  : add contribution for particles
         makeFromArray : create coefficients contributions
         )",
	 py::arg("mass"), py::arg("pos"), py::arg("time"),
	 py::arg("center") = Eigen::Vector3d::Zero(),
	 py::arg("roundrobin") = true, py::arg("posvelrows") = true)
    .def("initFromArray", &BasisClasses::BiorthBasis::initFromArray,
	 R"(
         Initialize coefficient accumulation

         Parameters
         ----------
         center : list, default=[0, 0, 0]
             vector of center positions
         rotation : numpy.ndarray
             the rotation matrix for the basis functions
             (default is identity matrix)

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
	 py::arg("center") = Eigen::Vector3d::Zero(),
	 py::arg("rotation") = RowMatrix3d::Identity())
    .def("addFromArray",
	 [](BasisClasses::BiorthBasis& A,
	    Eigen::VectorXd& mass, RowMatrixXd& pos)
	 {
	   A.addFromArray(mass, pos);
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
         makeFromArray : create coefficients contributions
         )",
	 py::arg("mass"), py::arg("pos"))
    .def("getFields", &BasisClasses::BiorthBasis::getFields,
	 R"(
         Return the field evaluations for a given cartesian position. The
         fields include density, potential, and force.  The density and
         potential evaluations are separated into full, axisymmetric and
         non-axisymmetric contributions.

         You can get the fields labels by using the __call__ method of the
         basis object.  This is equilevalent to a tuple of the getFields()
         output with a list of field labels.

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
         fields: numpy.ndarray

         See also
         --------
         getFieldsCoefs : get fields for each coefficient set
         __call__       : same as getFields() but provides field labels in a tuple
         )",
	 py::arg("x"), py::arg("y"), py::arg("z"))
    .def("getAccel", static_cast<Eigen::Vector3d& (BasisClasses::BiorthBasis*)(double, double, double)>(&BasisClasses::BiorthBasis::getAccel)
	 R"(
         Return the acceleration for a given cartesian position

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
         fields: numpy.ndarray

         See also
         --------
         getFields      : returns density, potential and acceleration
         getFieldsCoefs : get fields for each coefficient set
         __call__       : same as getFields() but provides field labels in a tuple
         )",
	 py::arg("x"), py::arg("y"), py::arg("z"))
    .def("getAccelArray", static_cast<RowMatrixXd& (BasisClasses::BiorthBasis*)(Eigen::VectorXd&, Eigen::VectorXd&, Eigen::VectorXd&)>(&BasisClasses::BiorthBasis::getAccel)
	 R"(
         Return the acceleration for a given cartesian position

         Parameters
         ----------
         x : ndarray
             x-axis positions
         y : ndarray
             y-axis positions
         z : ndarray
             z-axis positions

         Returns
         -------
         accel: numpy.ndarray

         See also
         --------
         getFields      : returns density, potential and acceleration
         getFieldsCoefs : get fields for each coefficient set
         __call__       : same as getFields() but provides field labels in a tuple
         )",
	 py::arg("x"), py::arg("y"), py::arg("z"))
    .def("getAccelArray", &BasisClasses::BiorthBasis::getAccelArray,
	 R"(
         Return the acceleration for a given cartesian position

         Parameters
         ----------
         x : ndarray
             x-axis positions
         y : ndarray
             y-axis positions
         z : ndarray
             z-axis positions

         Returns
         -------
         accel: numpy.ndarray

         See also
         --------
         getFields      : returns density, potential and acceleration
         getFieldsCoefs : get fields for each coefficient set
         __call__       : same as getFields() but provides field labels in a tuple
         )",
	 py::arg("x"), py::arg("y"), py::arg("z"))
    .def("getAccelArray", &BasisClasses::BiorthBasis::getAccelArray,
	 R"(
         Return the acceleration for a given cartesian position

         Parameters
         ----------
         x : ndarray
             x-axis positions
         y : ndarray
             y-axis positions
         z : ndarray
             z-axis positions

         Returns
         -------
         accel: numpy.ndarray

         See also
         --------
         getFields      : returns density, potential and acceleration
         getFieldsCoefs : get fields for each coefficient set
         __call__       : same as getFields() but provides field labels in a tuple
         )",
	 py::arg("x"), py::arg("y"), py::arg("z"))

py::class_<Pet>(m, "Pet")
   .def(py::init<const std::string &, int>())
   .def("set", static_cast<void (Pet::*)(int)>(&Pet::set), "Set the pet's age")
   .def("set", static_cast<void (Pet::*)(const std::string &)>(&Pet::set), "Set the pet's name");

    .def("getAccel", &BasisClasses::BiorthBasis::getAccelArray,
	 R"(
         This is a alias/overload of getAccelArray to provide a vector
         version of getAccel.

         Parameters
         ----------
         x : ndarray
             x-axis positions
         y : ndarray
             y-axis positions
         z : ndarray
             z-axis positions

         Returns
         -------
         accel: numpy.ndarray

         See also
         --------
         getAccelArray  : return the acceleration vectors given position vectors
         getFields      : returns density, potential and acceleration
         getFieldsCoefs : get fields for each coefficient set
         __call__       : same as getFields() but provides field labels in a tuple
         )",
	 py::arg("x"), py::arg("y"), py::arg("z"))
    .def("getFieldsCoefs", &BasisClasses::BiorthBasis::getFieldsCoefs,
	 R"(
         Return the field evaluations for a given cartesian position
         for every frame in a coefficient set.  The field evaluations are
         produced by a call to getFields().

         You get a dictionary of fields keyed by field name and an array
         of evaluation times for convenience.  These times will be the same
         as Times() for the coefficient object.

         Parameters
         ----------
         x : float
             x-axis position
         y : float
             y-axis position
         z : float
             z-axis position
         coefs: CoefClasses::Coefs
             the coefficient set

         Returns
         -------
         tuple of a dictionary of fields of array values, and an 
             array of evaluation times

         See also
         --------
         getFields  : get fields for the currently assigned coefficients
         __call__   : same getFields() but provides field labels in a tuple
         )",
	 py::arg("x"), py::arg("y"), py::arg("z"), py::arg("coefs"))
    .def("setFieldType",       &BasisClasses::BiorthBasis::setFieldType,
         R"(
         Set the coordinate system for force evaluations.  The natural 
         coordinates for the basis class are the default; spherical
         coordinates for SphericalSL and Bessel, cylindrical coordinates for
         Cylindrical, FlatDisk, and CBDisk, and Cartesian coordinates for the 
         Slab and Cube. This member function can be used to override the
         default. The available coorindates are: 'spherical', 'cylindrical', 
         'cartesian'.

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
         linear : bool, default=True
              use linear spacing

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
	 py::arg("numz")=40,
	 py::arg("linear")=true)
    // The following member needs to be a lambda capture because
    // orthoCheck is not in the base class and needs to have different
    // parameters depending on the basis type.  Here, the quadrature
    // is determined by the scale of the meridional grid.
    .def("orthoCheck", [](BasisClasses::Cylindrical& A, int knots)
	 {
	   return A.orthoCheck(knots);
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
        )", py::arg("knots")=400)
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

    py::class_<BasisClasses::Spherical, std::shared_ptr<BasisClasses::Spherical>, PySpherical, BasisClasses::BiorthBasis>(m, "Spherical")
      .def(py::init<const std::string&, const std::string&>(),
	 R"(
         Create a spherical basis

         Parameters
         ----------
         YAMLstring : str
             The YAML configuration for the spherical basis
         ForceID : str
             The string identifier for this force type

         Returns
         -------
         Spherical
              the new instance
         )", py::arg("YAMLstring"), py::arg("ForceID"))
      .def("getBasis", &BasisClasses::Spherical::getBasis,
	   R"(
           Get basis functions

	   Evaluate the potential-density basis functions on a linearly
	   spaced grid for inspection. The structure is a two-grid of dimension
	   lmax by nmax each pointing to a dictionary of 1-d arrays ('potential',
	   'density', 'rforce') of dimension numr.

           Parameters
           ----------
           rmin : float, default=0.0
                minimum radius
           rmax : float, default=1.0
                maximum radius
           numr : int, default=400
                number of equally spaced output points

           Returns
           -------
           list(list(dict))
               dictionaries of basis functions as lists indexed by l, n
           )",
	   py::arg("rmin")=0.0,
	   py::arg("rmax")=1.0,
	   py::arg("numr")=400)
      // The following member needs to be a lambda capture because
      // orthoCheck is not in the base class and needs to have
      // different parameters depending on the basis type.  Here the
      // user can and will often need to specify a quadrature value.
      .def("orthoCheck", [](BasisClasses::Spherical& A, int knots)
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
	return BasisClasses::Spherical::cacheInfo(cachefile);
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
	py::arg("cachefile"))
      .def_static("I", [](int l, int m)
      {
	if (l<0) throw std::runtime_error("l must be greater than 0");
	if (m<0) throw std::runtime_error("m must be greater than 0");
	if (abs(m)>l) throw std::runtime_error("m must be less than or equal to l");
	return (l * (l + 1) / 2) + m;
      },
	R"(
        Calculate the index of a spherical harmonic element given the angular numbers l and m .

        Parameters
        ----------
        l : int
            spherical harmonic order l
        m : int
            azimuthal order m

        Returns
        -------
        I : int
            index array packing index
      )",
	py::arg("l"), py::arg("m"))
      .def_static("invI", [](int I)
      {
	if (I<0) std::runtime_error("I must be an interger greater than or equal to 0");
	int l = std::floor(0.5*(-1.0 + std::sqrt(1.0 + 8.0 * I)));
	int m = I - int(l * (l + 1) / 2);
	return std::tuple<int, int>(l, m);
      },
	R"(
        Calculate the spherical harmonic indices l and m from the coefficient array packing index I

        Parameters
        ----------
        I : int
            the spherical coefficient array index

        Returns
        -------
        (l, m) : tuple
            the harmonic indices (l, m).
      )", py::arg("I"));

  
    py::class_<BasisClasses::SphericalSL, std::shared_ptr<BasisClasses::SphericalSL>, BasisClasses::Spherical>(m, "SphericalSL")
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
	py::arg("knots")=40);

  
    py::class_<BasisClasses::Bessel, std::shared_ptr<BasisClasses::Bessel>, BasisClasses::Spherical>(m, "Bessel")
      .def(py::init<const std::string&>(),
	 R"(
         Create a spherical Bessel-function basis

         Parameters
         ----------
         YAMLstring : str
             The YAML configuration for the spherical Bessel-function basis

         Returns
         -------
         Bessel
              the new instance
         )", py::arg("YAMLstring"))
      // The following member needs to be a lambda capture because
      // orthoCheck is not in the base class and needs to have
      // different parameters depending on the basis type.  Here the
      // user can and will often need to specify a quadrature value.
      .def("orthoCheck", [](BasisClasses::Bessel& A, int knots)
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
	py::arg("knots")=40);

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

  py::class_<BasisClasses::CBDisk, std::shared_ptr<BasisClasses::CBDisk>, PyCBDisk, BasisClasses::BiorthBasis>(m, "CBDisk")
    .def(py::init<const std::string&>(),
	 R"(
         Create a Clutton-Brock 2d disk basis

         Parameters
         ----------
         YAMLstring : str
             The YAML configuration for the razor-thin EOF basis.  The default 
             parameters will give an exponential disk with scale length of
             0.01 units. Set the disk scale length using the 'scale'  parameter.

         Returns
         -------
         CBDisk
             the new instance
         )", py::arg("YAMLstring"))
    .def("getBasis", &BasisClasses::CBDisk::getBasis,
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
    .def("orthoCheck", [](BasisClasses::CBDisk& A)
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
      return BasisClasses::CBDisk::cacheInfo(cachefile);
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

  py::class_<BasisClasses::Slab, std::shared_ptr<BasisClasses::Slab>, PySlab, BasisClasses::BiorthBasis>(m, "Slab")
    .def(py::init<const std::string&>(),
	 R"(
         Create a slab basis, periodic on the unit square and finite in vertical extent

         Parameters
         ----------
         YAMLstring : str
             The YAML configuration for the Slab basis. The coordinates are the
             unit square in x, y with origin at (0, 0) and maximum extent (1, 1)
             and maximum vertical extent of -zmax to zmax. The default 
             parameters are wave numbers between [-6,...,6] in x, y, and order 6
             for the vertical basis.

         Returns
         -------
         Cube
             the new instance
         )", py::arg("YAMLstring"))
    .def("getBasis", &BasisClasses::Slab::getBasis,
	 R"(
         Get basis functions

	 Evaluate the potential-density basis functions on a linearly spaced grid for 
         inspection. The structure is a three-dimensional grid of dimension (nmaxx+1) by 
         (nmaxy+1) by (nmaxz) each pointing to a dictionary of 1-d arrays ('potential',
	 'density', 'rforce') of dimension numz.

         Parameters
         ----------
         zmin : float, default=-1.0
             minimum height
         zmax : float, default=1.0
             maximum height
         numz : int, default=400
             number of equally spaced output points

         Returns
         -------
         list(list(dict))
               dictionaries of basis functions as lists indexed by nx, ny, nz

         Example
         -------
         To plot the nx=ny=0 basis functions, you might use the following code:

         >>> mat = slab_basis.getBasis(-0.5, 0.5, 200)
         >>> z = np.linspace(-0.5, 0.5, 200)
         >>> for n in range(len(mat[0][0])):
         >>> plt.plot(z, mat[0][0][n]['potential'], label=str(n))
         >>> plt.legend()
         >>> plt.xlabel('z')
         >>> plt.ylabel('potential')
         >>> plt.show()
         )",
	 py::arg("zmin")=-1.0,
	 py::arg("zmax")=1.0,
	 py::arg("numz")=400)
    .def("orthoCheck", [](BasisClasses::Cube& A)
    {
      return A.orthoCheck();
    },
      R"(
      Check orthgonality of basis functions by quadrature

      Inner-product matrix of indexed by flattened wave number (nx, ny, nz) where
      each of nx is in [-nmaxx, nmaxx], ny is in [-nmaxy, nmaxy] and nz is in 
      [0, nmaxz-1].  This is an analytic basis so the orthogonality matrix is not a 
      check of any numerical computation other than the quadrature itself.  It is 
      included for completeness.

      Parameters
      ----------
      None

      Returns
      -------
      numpy.ndarray
          list of numpy.ndarrays
      )"
      );


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
      analyic basis so the orthogonality matrix is not a check of any numerical
      computation other than the quadrature itself.  It is included for completeness.

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
         Create a orthogonal basis for representing general phase-space fields

         Parameters
         ----------
         YAMLstring : str
             The YAML configuration for the field basis

         Returns
         -------
         FieldBasis
              the new instance

         Notes
         -----
         A FieldBasis instance may be created directly or by the Basis.factory().
         The fields are defined by a user provided function, addPSFunction().
         VelocityBasis is derived from this class and provides a field function
         for velocity fields by default.

         The evaluation geometry is a polar using the phase-space x-y coordinates
         if the 'dof' parameter is 2 and spherical using x-y-z coordinates if
         'dof' is 3.  The output values and their coordinates are determined by
         the user-suppied field function.

         See also
         --------
         addPSFunction: register the user-supplied field generating function and
                        field names.
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
    .def("createFromReader", &BasisClasses::FieldBasis::createFromReader,
	 R"(
         Generate the coefficients from the supplied ParticleReader

         Parameters
         ----------
         reader : Particle reader
             the ParticleReader instance
         center : list, default=[0, 0, 0]
	     an optional expansion center location
         rotation : numpy.ndarray
             the rotation matrix for the basis functions
             (default is identity matrix)

         Returns
         -------
         CoefStruct
             the basis coefficients computed from the particles
         )",
	 py::arg("reader"), 
	 py::arg("center") = Eigen::Vector3d::Zero(),
	 py::arg("rotation") = RowMatrix3d::Identity())
    .def("initFromArray", &BasisClasses::FieldBasis::initFromArray,
	 R"(
         Initialize coefficient accumulation

         Parameters
         ----------
         center : list, default=[0, 0, 0]
             vector of center positions
         rotation : numpy.ndarray, default=Identity
             rotation matrix to apply to the phase-space coordinates

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
	 py::arg("center") = Eigen::Vector3d::Zero(),
	 py::arg("rotation") = RowMatrix3d::Identity())
    .def("initFromArray",
	 [](BasisClasses::FieldBasis& A, Eigen::Vector3d center)
	 {
	   return A.initFromArray(center, RowMatrix3d::Identity());
	 },
	 R"(
         Initialize coefficient accumulation

         Parameters
         ----------
         center : list, default=[0, 0, 0]
             vector of center positions
         rotation : numpy.ndarray, default=Identity
             rotation matrix to apply to the phase-space coordinates

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
	 py::arg("center") = Eigen::Vector3d::Zero())
    .def("addFromArray",
	 [](BasisClasses::FieldBasis& A,
	    Eigen::VectorXd& mass, RowMatrixXd& pos)
	 {
	   A.addFromArray(mass, pos);
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
         makeFromArray : create coefficients contributions
         )",
	 py::arg("mass"), py::arg("pos"))
    .def("makeFromArray", &BasisClasses::FieldBasis::makeFromArray,
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
         CoefStruct
             the coefficient structure created from the particles

         See also
         --------
         initFromArray : initialize for coefficient contributions
         addFromArray  : add contribution for particles
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
         spherical, cylindrical, or Cartesian systems.  By default, the
         coordinates will be cylindrical if the 'dof' parameter is 2 and
         spherical if the 'dof' parameter is 3. This default can be changed
         using the setFieldType() member.  For example, the following call:
              basis->setFieldType('cart')
         will selection Cartesian velocities u, v, w for output.
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
         AccelFunc instance that interpolates coefficients from the Coefs 
         database for every time

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

	  py::array_t<float> ret = make_ndarray3<float>(O);
	  return std::tuple<Eigen::VectorXd, py::array_t<float>>(T, ret);
	},
	R"(
        Compute particle orbits in gravitational field from the bases

        Integrate a list of initial conditions from 'tinit' to 'tfinal' with
	a step size of 'h' using the list of basis and coefficient pairs. The
        step size will be adjusted to provide uniform sampling.  Every
	step will be returned unless you provide an explicit value for 'nout',
        the number of desired output steps.  In this case, the code will
        choose new set step size equal or smaller to the supplied step size
        with a stride to provide exactly 'nout' output times.

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
            the number of output points, if specified

        Returns
        -------
        tuple(numpy.array, numpy.ndarray)
            time and phase-space arrays
        )",
	py::arg("tinit"), py::arg("tfinal"), py::arg("h"),
	py::arg("ps"), py::arg("basiscoef"), py::arg("func"),
	py::arg("nout")=0);
}
