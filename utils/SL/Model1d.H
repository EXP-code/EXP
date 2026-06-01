#ifndef _Model1d_H_
#define _Model1d_H_

#include <string>

#include "interp.H"

class OneDModel : public std::enable_shared_from_this<OneDModel>
{

protected:

  int Nfgrid = 800;		// Number of points for tabulating
				// distribution function for
				// interpolation in getE(M)

  Linear1d dfgrid;		// Interpolation object for
				// distribution function grid (for
				// getE(M))

  std::vector<double> Egrid, Mgrid;

  void compute_dfgrid(void);

public:
  bool dist_defined;

  virtual double get_mass(const double) = 0;
  virtual double get_density(const double) = 0;
  virtual double get_pot(const double) = 0;
  virtual double get_dpot(const double) = 0;
  virtual double get_dpot2(const double) = 0;
  virtual void   get_pot_dpot(const double, double&, double&) = 0;

  // Required members of mass model

  double get_mass(const double x, const double y, const double z) {
    return get_mass(z);}

  double get_density(const double x, const double y, const double z) {
    return get_density(z);}

  double get_pot(const double x, const double y, const double z) {
    return get_pot(z);}


  // Addiional member functions

  virtual double get_min_radius(void) = 0;
  virtual double get_max_radius(void) = 0;
  virtual double get_scale_height(void) = 0;
  virtual double distf(const double, const double V=0.0) = 0;
  virtual double dfde (const double, const double V=0.0) = 0;
  virtual double dfdv (const double, const double V=0.0) = 0;

  //! Model types
  enum class ModelType {table, lowiso, sech2, sech2u, sech2halo, sech2mu,
			uniform, cosine};


  //! Model names for sysout display (avoid temporaries and dynamic
  //! allocations on every class instantiation)
  inline static const std::map<ModelType, std::string_view> model_names = {
    {ModelType::table,     "file"},
    {ModelType::lowiso,    "LowIso"},
    {ModelType::sech2,     "Sech2"},
    {ModelType::sech2u,    "Sech2u"},
    {ModelType::sech2halo, "Sech2Halo"},
    {ModelType::sech2mu,   "Sech2(mu)"},
    {ModelType::uniform,   "Uniform"},
    {ModelType::cosine,    "Cosine"}
  };

  //! Model types for command line parsing
  //! 'inline static const' allows initialization inside the class definition.
  //! 'std::string_view' keys completely prevent temporary string creations.
  inline static const std::map<std::string_view, ModelType> model_map = {
    {"table",     ModelType::table},
    {"lowiso",    ModelType::lowiso},
    {"sech2",     ModelType::sech2},
    {"sech2u",    ModelType::sech2u},
    {"sech2halo", ModelType::sech2halo},
    {"sech2mu",   ModelType::sech2mu},
    {"uniform",   ModelType::uniform},
    {"cosine",    ModelType::cosine}
  };
  
  //! Model type for sysout display
  static std::string get_model_name(ModelType type) {
    return std::string(model_names.at(type));
  }

  //! Command line parsing for model types
  static ModelType parse_model_type(const std::string& str) {
    std::string MODEL = str;
    std::transform(MODEL.begin(), MODEL.end(), MODEL.begin(), 
		   [](unsigned char c){ return std::tolower(c); });
    
    // A std::string implicitly and cheaply converts to a std::string_view, 
    // so finding it in our map works cleanly and creates zero temporaries.
    auto it = model_map.find(MODEL);
    if (it != model_map.end()) {
      return it->second;
    } else {
      std::ostringstream oss;
      oss << "Invalid model type: " << str << ". Valid types are: ";
      for (const auto& pair : model_map) {
	oss << pair.first << " ";
      }
      throw std::invalid_argument(oss.str());
    }
  }

};

class OneDModelTable : public OneDModel
{
protected:

  Spline1d mass, dens, pot;

  int even;
  int num;
  int numdf;
  double half_height;
  std::vector<double> params;

public:

  OneDModelTable() {};

  OneDModelTable(string filename, int PARM=0);

  OneDModelTable(int num, double *r, double *d, 
		 double *m, double *p, string ID = "" );

  // Required member functions

  double get_mass   (const double);
  double get_density(const double);
  double get_pot    (const double);
  double get_dpot   (const double);
  double get_dpot2  (const double);
  void   get_pot_dpot (const double, double&, double&);
  
  // Additional member functions

  const int get_num_param(void) { return params.size(); }
  const double get_param(int i) { return params[i-1]; }
  double get_scale_height(void) { return half_height; }

  double get_min_radius(void) { return mass.xlo(); }
  double get_max_radius(void) { return mass.xhi(); }
  int    grid_size(void) { return num; }

//  double distf(double E, double V);
//  double dfde(double E, double V);
//  double dfdv(double E, double V);
};

class LowIso : public OneDModelTable
{
private:
  double w0, Bfac, betak, gammak;
  double dispx, normx;

  void setup_model(void);

public:
  LowIso(string filename, double DISPX=0.159154943091895335768) :
    OneDModelTable(filename, 1) { 
      dispx = DISPX;
      setup_model();
    }

  double get_pot(const double);
  double get_dpot(const double);
  double get_dpot2(const double);
  void   get_pot_dpot(const double, double&, double&);

  double distf(const double E, const double V=0.0);
  double dfde(const double E, const double V=0.0);
  double dfdv(const double E, const double V=0.0);
};


class Sech2 : public OneDModel
{
private:
  double h;
  double dispz, dispx;
  double norm;

  static double HMAX;

public:

  Sech2(void)
  {
    dispz = 1.0;
    dispx = 1.0;
    //
    // Units: G = rho_o = 1
    //
    h = sqrt(dispz/(2.0*M_PI));
    norm = 1.0/( sqrt(2.0*M_PI*dispz) );
    dist_defined = true;
  }

  Sech2(const double DISPZ, const double DISPX=1.0)
  {
    dispz = DISPZ;
    dispx = DISPX;
    //
    // Units: G = rho_o = 1
    //
    h = sqrt(dispz/(2.0*M_PI));
    norm = 1.0/( sqrt(2.0*M_PI*dispz) );
    dist_defined = true;
  }
  
  double get_mass(const double z)
  {
    return 2.0*h/(1.0 + exp(-2.0*z/h));
  }

  double get_density(const double z)
  {
    double zz = fabs(z);
    double ret = 2.0*exp(-zz/h)/(1.0 + exp(-2.0*zz/h));
    return ret*ret;
  }

  double get_pot(const double z)
  {
    double zz = fabs(z);
    return 4.0*M_PI*h*(zz + h*log(1.0 + exp(-2.0*zz/h)) - h*M_LN2);
  }

  double get_dpot(const double z)
  {
    double zz = fabs(z);
    double ret = (1.0 - exp(-2.0*zz/h))/(1.0 + exp(-2.0*zz/h));
    return 4.0*M_PI*h* ret * z/(zz+1.0e-18);
  }

  double get_dpot2(const double z)
  {
    double zz = fabs(z);
    double ret = 2.0*exp(-zz/h)/(1.0 + exp(-2*zz/h));
    return 4.0*M_PI*ret*ret;
  }

  void get_pot_dpot(const double z, double& p, double& dp)
  {
    double zz = fabs(z);
    p = 4.0*M_PI*h*(zz + h*log(1.0 + exp(-2.0*zz/h)) - h*M_LN2);
    double ret = (1.0 - exp(-2.0*zz/h))/(1.0 + exp(-2.0*zz/h));
    dp = 4.0*M_PI*h* ret * z/(zz+1.0e-18);
  }

  double get_mass(const double x, const double y, const double z)
  {
    return get_mass(z);
  }

  double get_density(const double x, const double y, const double z)
  {
    return get_density(z);
  }

  double get_pot(const double x, const double y, const double z)
  {
    return get_pot(z);
  }

  double get_min_radius(void) { return 0.0; }
  double get_max_radius(void) { return HMAX*h; }
  double get_scale_height(void) { return h; }

  static void set_hmax(double hmax) { HMAX = hmax; }

  double distf(const double E, const double p)
  {
    return exp(-E/dispz - 0.5*p*p/dispx) * norm;
  }

  double dfde(const double E, const double p)
  {
    return -exp(-E/dispz - 0.5*p*p/dispx)/dispz * norm;
  }

  double dfdv(const double E, const double p)
  {
    return -exp(-E/dispz - 0.5*p*p/dispx)*p/dispx * norm;
  }

};



class Sech2mu : public OneDModel
{
private:
  double mu, h;
  double dispz, dispx;
  double dnorm, knorm, fnorm;

  static double HMAX;

public:

  Sech2mu(void)
  {
    mu    = 1.0;
    dispz = 1.0;
    dispx = 1.0;
    h     = 1.0;
    //
    // Units: G = mu = 1
    //
    dnorm = 0.25*mu/h;
    knorm = 2.0*M_PI*mu*h/dispz;
    fnorm = mu/(4.0*h*knorm*sqrt(2.0*M_PI*dispz));
    dist_defined = true;
  }

  Sech2mu(const double DISPZ, const double H, const double DISPX=1.0)
  {
    mu    = 1.0;
    dispz = DISPZ;
    dispx = DISPX;
    h     = H;
    //
    // Units: G = mu = 1
    //
    dnorm = 0.25*mu/h;
    knorm = 2.0*M_PI*mu*h/dispz;
    fnorm = mu/(4.0*h*knorm*sqrt(2.0*M_PI*dispz));
    dist_defined = true;
  }
  
  void setMu(double Sigma0)
  {
    mu = Sigma0;
    //
    // Units: G = 1, mu = Sigma0
    //
    dnorm = 0.25*mu/h;
    knorm = 2.0*M_PI*mu*h/dispz;
    fnorm = mu/(4.0*h*knorm*sqrt(2.0*M_PI*dispz));
    dist_defined = true;
  }

  double get_mass(const double z) {
    return mu*exp(z/h)/(1.0 + exp(z/h));
  }

  double get_density(const double z) {
    double zz  = fabs(z);
    double fac = exp(zz/h);
    return 4.0*dnorm/(fac + 1.0/fac + 2.0);
  }

  double get_pot(const double z) {
    double zz = fabs(z);
    return 2.0*dispz*(-M_LN2 + 0.5*zz/h + log(1.0 + exp(-zz/h)) ) - dispz*log(knorm);
  }

  double get_dpot(const double z) {
    double zz = fabs(z);
    double ret = (1.0 - exp(-zz/h))/(1.0 + exp(-zz/h));
    return 0.5*dispz*ret/h * z/(zz+1.0e-18);
  }

  double get_dpot2(const double z) {
    double zz = fabs(z);
    double ret = 2.0*exp(-0.5*zz/h)/(1.0 + exp(-zz/h));
    return 0.5*dispz*ret*ret/h/h;
  }

  void get_pot_dpot(const double z, double& p, double& dp) {
    double zz = fabs(z);
    p = 2.0*dispz*(-M_LN2 + 0.5*zz/h + log(1.0 + exp(-zz/h)) ) - dispz*log(knorm);
    double ret = (1.0 - exp(-zz/h))/(1.0 + exp(-zz/h));
    dp = 0.5*dispz*ret/h * z/(zz+1.0e-18);
  }

  double get_mass(const double x, const double y, const double z) {
    return get_mass(z);
  }

  double get_density(const double x, const double y, const double z) {
    return get_density(z);
  }

  double get_pot(const double x, const double y, const double z) {
    return get_pot(z);
  }

  double get_min_radius(void) { return 0.0; }
  double get_max_radius(void) { return HMAX*h; }
  double get_scale_height(void) { return h; }

  static void set_hmax(double hmax) { HMAX = hmax; }

  double distf(const double E, const double p) {
    return exp(-E/dispz - 0.5*p*p/dispx) * fnorm;
  }

  double dfde(const double E, const double p) {
    return -exp(-E/dispz - 0.5*p*p/dispx)/dispz * fnorm;
  }

  double dfdv(const double E, const double p) {
    return -exp(-E/dispz - 0.5*p*p/dispx)*p/dispx * fnorm;
  }

};


class Sech2Halo : public OneDModelTable
{
private:
  double h, rho0;
  double dispz, dispx;
  double dratio, hratio;
  double hh, rho0h;
  double norm, hmax;

  bool model_computed;

  static double HMAX;

  void reset();

public:

  static int NTABLE;
  static double OFFSET;
  static bool MU;

  Sech2Halo(void) {
    dispz = 1.0;
    dispx = 1.0;
    dratio = 0.0;
    hratio = 1.0;

    dist_defined = false;
    model_computed = false;
  }

  Sech2Halo(const double DISPZ, const double DRATIO, const double HRATIO,
	    const double DISPX=1.0);

  double get_pot(const double z) {
    if (!model_computed) reset();
    return OneDModelTable::get_pot(z) + 4.0*M_PI*hh*hh*rho0h*log(cosh(z/hh));
  }

  double get_dpot(const double z) {
    if (!model_computed) reset();
    return OneDModelTable::get_dpot(z) + 4.0*M_PI*hh*rho0h*tanh(z/hh);
  }

  double get_dpot2(const double z) {
    if (!model_computed) reset();
    double sech = 1.0/cosh(z/hh);
    return OneDModelTable::get_dpot2(z) + 4.0*M_PI*rho0h*sech*sech;
  }

  void get_pot_dpot(const double z, double& p, double& dp) {
    if (!model_computed) reset();
    OneDModelTable::get_pot_dpot(z, p, dp);
    p += 4.0*M_PI*hh*hh*rho0h*log(cosh(z/hh));
    dp += 4.0*M_PI*hh*rho0h*tanh(z/hh);
  }


  double get_min_radius(void) { return 0.0; }
  double get_max_radius(void) { return hmax*h; }
  double get_scale_height(void) { return h; }
  double get_scale_height_halo(void) { return hh; }
  double get_rho0(void) { return rho0; }
  double get_rho0_halo(void) { return rho0h; }

  static void set_hmax(double hmax) { HMAX = hmax; }
  
  double distf(const double E, const double p=0.0) {
    if (!model_computed) reset();
    return exp(-E/dispz - 0.5*p*p/dispx) * norm;
  }

  double dfde(const double E, const double p=0.0) {
    if (!model_computed) reset();
    return -exp(-E/dispz - 0.5*p*p/dispx)/dispz * norm;
  }

  double dfdv(const double E, const double p=0.0) {
    if (!model_computed) reset();
    return -exp(-E/dispz - 0.5*p*p/dispx)*p/dispx * norm;
  }
};

class Uniform : public OneDModel
{
private:
  // Half-thickness of the uniform slab
  double H;
  // Velocity widths in x and y directions (for distribution function)
  double Vx, Vy;
  // Vertical frequency for the uniform slab potential
  double Omega;
  // Normalization constant for the distribution function
  double norm;
  
public:
  
  Uniform(void)
  {
    Vx = 1.0;
    Vy = 1.0;
    //
    // Units: G = Sigma_o = 1
    //
    H = 1.0;
    Omega = std::sqrt(2.0*M_PI/H);
    norm  = 1.0/(std::sqrt(2.0)*M_PI*H);
    // ModelID = "Uniform";
  }

  Uniform(double h, double Vx)
  {
    this->H  = h;
    this->Vx = Vx;
    this->Vy = Vx;
    //
    // Units: G = Sigma_o = 1
    //
    Omega = std::sqrt(2.0*M_PI/H);
    norm  = 1.0/(std::sqrt(2.0)*M_PI*H);
    // ModelID = "Uniform";
  }
  
  Uniform(double h, double Vx, double Vy)
  {
    this->H  = h;
    this->Vx = Vx;
    this->Vy = Vy;
    //
    // Units: G = Sigma_o = 1
    //
    Omega = std::sqrt(2.0*M_PI/H);
    norm  = 1.0/(std::sqrt(2.0)*M_PI*H);
    // ModelID = "Uniform";
  }
  
  double get_mass(const double z)
  {
    if (z<-H)
      return 0.0;
    else if (z > H)
      return 1.0;
    else {
      return (z + H)/(2.0*H);
    }
  }

  double get_density(const double z)
  {
    if (fabs(z) > H)
      return 0.0;
    else
      return 1.0/(2.0*H);
  }

  double get_pot(const double z)
  {
    // Offset defined so that potential is zero at |z|=H, and negative
    // for |z|<H inside the slab.
    double zz = fabs(z);
    if (zz > H)
      return 2.0*M_PI*zz - 2.0*M_PI*H;
    else
      return M_PI*zz*zz/H - M_PI*H;
  }

  double get_dpot(const double z)
  {
    double zz = fabs(z);
    if (zz > H)
      return 2.0*M_PI*z/zz;
    else
      return 2.0*M_PI*z/H;
  }

  double get_dpot2(const double z)
  {
    double zz = fabs(z);
    if (zz > H)
      return 0.0;
    else
      return 2.0*M_PI/H;
  }

  void get_pot_dpot(const double z, double& p, double& dp)
  {
    double zz = fabs(z);
    if (zz > H) {
      p  = 2.0*M_PI*zz - 2.0*M_PI*H;
      dp = 2.0*M_PI*z/zz;
    } else {
      p  = M_PI*zz*zz/H - M_PI*H;
      dp = 2.0*M_PI*z/H;
    }
  }

  double get_mass(const double x, const double y, const double z)
  {
    return get_mass(z);
  }

  double get_density(const double x, const double y, const double z)
  {
    return get_density(z);
  }

  double get_pot(const double x, const double y, const double z)
  {
    return get_pot(z);
  }

  double get_min_radius(void)   { return 0.0; }
  double get_max_radius(void)   { return H;   }
  double get_scale_height(void) { return H;   }

  double distf(const double E, const double p=0.0)
  {
    // Convert energy to potential, so that the distribution function
    // is nonzero only.  That is, assume E = Omega*Ja, not zero at +/-
    // H.
    double OJ = fabs(E + M_PI*H);
    return norm / std::sqrt(OJ);
  }

  double dfde(const double E, const double p=0.0)
  {
    // Convert energy to potential, so that the distribution function
    // is nonzero only.  That is, assume E = Omega*Ja, not zero at +/-
    // H.
    double OJ = fabs(E + M_PI*H);
    return -0.5*Omega*norm / std::pow(OJ, 1.5);
  }

  double dfdv(const double E, const double p=0.0)
  {
    // The distribution function is independent of velocity, so this
    // derivative is zero. There are delta functions at v=-Vx and
    // v=+Vx, but those are not captured by this simple functional
    // form.
    return 0.0;
  }

};

//! A simple consine-bell slab model with finite thickness
class Cosine : public OneDModel
{
private:
  //! Half-thickness of the uniform slab
  double H;

  //! Central density of the cosine slab
  double rho0;

  //! Potential zero point is defined at |z|=H, so that potential is
  //! negative inside the slab and zero at the edges
  double Phi0;

  //! Velocity widths in x and y directions (for distribution function)
  double Vx, Vy;

  //! Vertical frequency for the uniform slab potential
  double Omega;

  //! Normalization constant for the distribution function
  double norm;

  //! Regularization parameter for the distribution function (if needed)
  double epsilon;

  //! Number of points in the tabulated distribution function
  int numdf;

  //! Number of quadrature points for computing DF
  int nint;

  //! Tabulated action and distribution function for interpolation
  std::vector<double> Etab, Ftab;

  //! DF table interpolation object
  Linear1d df_interp;

  //! Precompute the distribution function table for interpolation
  void compute_model();

public:

  //! Null constructor with default parameters
  Cosine(void)
  {
    Vx = 1.0;
    Vy = 1.0;
    //
    // Units: G = Sigma_o = 1
    //
    H = 1.0;
    rho0 = 1.0/H;
    Phi0 = M_PI*H*(4.0/(M_PI*M_PI) + 1.0);
    // ModelID = "Cosine";
    numdf = 800;
    nint  = 800;
    dist_defined = true;
    compute_model();
  }

  //! Constructor with specified half-thickness and velocity width
  //! (same for x and y)
  Cosine(double h, double Vx, double epsilon)
  {
    this->H  = h;
    this->Vx = Vx;
    this->Vy = Vx;
    this->epsilon = epsilon;
    //
    // Units: G = Sigma_o = 1
    //
    rho0 = 1.0/H;
    Phi0 = M_PI*H*(4.0/(M_PI*M_PI) + 1.0);
    // ModelID = "Cosine";
    numdf = 800;
    nint  = 800;
    dist_defined = true;
    compute_model();
  }
  
  //! Constructor with specified half-thickness and velocity widths in
  //! x and y
  Cosine(double h, double Vx, double Vy, double epsilon)
  {
    this->H  = h;
    this->Vx = Vx;
    this->Vy = Vy;
    this->epsilon = epsilon;
    //
    // Units: G = Sigma_o = 1
    //
    rho0 = 1.0/H;
    Phi0 = M_PI*H*(4.0/(M_PI*M_PI) + 1.0);
    // ModelID = "Cosine";
    numdf = 800;
    nint  = 800;
    dist_defined = true;
    compute_model();
  }
  
  void setDF(int numdf, int nint=800)
  {
    this->numdf = numdf;
    this->nint  = nint;
    compute_model();
  }

  double get_mass(const double z)
  {
    if (z<-H)
      return 0.0;
    else if (z > H)
      return 1.0;
    else {
      return 0.5*rho0*(z + H + H/M_PI*sin(M_PI*z/H));
    }
  }

  double get_density(const double z)
  {
    if (fabs(z) > H)
      return 0.0;
    else {
      double cosfac = cos(M_PI*z/(2.0*H));
      return rho0*cosfac*cosfac;
    }
  }

  double get_pot(const double z)
  {
    // Offset defined so that potential is zero at |z|=H, and negative
    // for |z|<H inside the slab.
    double zz = fabs(z);
    if (zz > H)
      return 2.0*M_PI*rho0*H*(zz - H);
    else
      return 2.0*M_PI*rho0*(0.5*z*z + H*H/(M_PI*M_PI)*(1 - cos(M_PI*z/H))) - Phi0;
  }

  double get_dpot(const double z)
  {
    double zz = fabs(z);
    if (zz > H)
      return 2.0*M_PI*rho0*H*z/zz;
    else
      return 2.0*M_PI*rho0*(z + H/M_PI*sin(M_PI*z/H));
  }
  
  double get_dpot2(const double z)
  {
    double zz = fabs(z);
    if (zz > H)
      return 0.0;
    else
      return 2.0*M_PI*rho0*(1 + cos(M_PI*z/H));
  }

  void get_pot_dpot(const double z, double& p, double& dp)
  {
    double zz = fabs(z);
    if (zz > H) {
      p  = 2.0*M_PI*rho0*H*(zz - H);
      dp = 2.0*M_PI*rho0*H*z/zz;
    } else {
      p  = 2.0*M_PI*rho0*(0.5*z*z + H*H/(M_PI*M_PI)*(1 - cos(M_PI*z/H))) - Phi0;
      dp = 2.0*M_PI*rho0*(z + H/M_PI*sin(M_PI*z/H));
    }
  }

  double get_mass(const double x, const double y, const double z)
  {
    return get_mass(z);
  }

  double get_density(const double x, const double y, const double z)
  {
    return get_density(z);
  }

  double get_pot(const double x, const double y, const double z)
  {
    return get_pot(z);
  }

  double get_min_radius(void)   { return 0.0; }
  double get_max_radius(void)   { return H;   }
  double get_scale_height(void) { return H;   }

  double distf(const double E, const double p=0.0)
  {
    if (E>=0.0)
      return 0.0;
    else
      return df_interp.eval(E);
  }

  double dfde(const double E, const double p=0.0)
  {
    if (E>=0.0)
      return 0.0;
    else
      return df_interp.deriv(E);
  }

  double dfdv(const double E, const double p=0.0)
  {
    // The distribution function is independent of velocity, so this
    // derivative is zero. There are delta functions at v=-Vx and
    // v=+Vx, but those are not captured by this simple functional
    // form.
    return 0.0;
  }

};

#endif

