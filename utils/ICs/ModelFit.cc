#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

#include <SimAnn.H>
#include <libvars.H>
#include <cxxopts.H>

//! Base class for density fitting functions
class ModelFit : public std::function<double(std::vector<double>&)>
{
protected:
  //! ID string
  std::string ID;

  //! Data
  std::vector<double> r, d;

  //! Model dimension
  int ndim;

  //! Model function
  virtual double fct(double r) = 0;

  //! Maps infinite interval into finite interval using inverse tangent function
  static double scale(double x, double s)
  {
    return s*(atan(x) + M_PI_2)/M_PI;
  }

  //! Maps finite interval into infinite interval using tangent function
  static double invscale(double y, double s)
  {
    return tan(M_PI*y/s - M_PI_2);
  }

public:

  //! Return string ID
  std::string getID() { return ID; }

  //! Return model dimension
  int ModelDim() { return ndim; }

  //! Evaluate the cost function
  virtual double operator()(std::vector<double>& p) = 0;

  //! Evaluate fitting function
  virtual double fit(std::vector<double>& p, double r) = 0;

  //! Convert true parameters to SimAnn scaled parameters in [-inf, inf]
  virtual std::vector<double> ScaleParam(const std::vector<double> p) = 0;

  //! Convert SimAnn scaled parameters to true parameters
  virtual std::vector<double> TrueParam (const std::vector<double> p) = 0;

  //! Label for each dimension for pretty printing
  virtual std::vector<std::string> ParamLabels() = 0;
};



class TwoPowerTrunc : public ModelFit
{
private:
  double m, a, b, c, t, w, expn;

  virtual double fct(double r)
  {
    return m*pow(r, -a)*pow(1.0 + (r/c), -b)*(1.0 + erf(-(r - t)/w))*0.5;
  }
    
public:
  
  /** Constructor
      @param r is the location of the radial bins
      @param d is the density computed from the bins
      @param expn is the density weighting (expn=1 is unweighted)
  */
  TwoPowerTrunc(std::vector<double>& r, std::vector<double>& d,
		double expn=1.0)
  {
    ID = "TwoPowerTrunc";

    this->r    = r;
    this->d    = d;
    this->ndim = 6;
    this->expn = expn;
  }

  virtual double operator()(std::vector<double>& p);

  virtual double fit(std::vector<double>& p, double r)
  {
    m = p[0];
    a = scale(p[1], 4.0);
    b = scale(p[2], 4.0);
    c = scale(p[3], 2.0);
    t = scale(p[4], 2.0);
    w = scale(p[5], 2.0);
    return fct(r);
  }

  std::vector<double> TrueParam(const std::vector<double> p)
  {
    return {p[0], scale(p[1], 4.0), scale(p[2], 4.0), 
	    scale(p[3], 2.0), scale(p[4], 2.0), scale(p[5], 2.0)};
  }

  std::vector<double> ScaleParam(const std::vector<double> p)
  {
    return {p[0], invscale(p[1], 4.0), invscale(p[2], 4.0), 
	    invscale(p[3], 2.0), invscale(p[4], 2.0), invscale(p[5], 2.0)};
  }

  std::vector<std::string> ParamLabels()
  {
    return {"Norm", "Alpha", "Beta", "r_s", "r_trunc", "w_trunc"};
  }

};


double TwoPowerTrunc::operator()(std::vector<double>& params)
{
  m = params[0];
  a = scale(params[1], 4.0);
  b = scale(params[2], 4.0);
  c = scale(params[3], 2.0);
  t = scale(params[4], 2.0);
  w = scale(params[5], 2.0);

  if (t<=0.05) return std::numeric_limits<double>::max();
  if (t >1.0)  return std::numeric_limits<double>::max();
  if (w<0.01)  return std::numeric_limits<double>::max();
  if (w>0.20)  return std::numeric_limits<double>::max();


  int num    = r.size();
  double ans = 0.0;

  for (int i=0; i<num; i++) {
    double val = fct(r[i]);
    double del = (d[i] - val)/pow(val+1.0e-24, expn);
    ans += del*del;
  }
  
  return ans;
}


class OnePowerTrunc : public ModelFit
{
private:
  double m, a, t, w, expn;

  virtual double fct(double r)
  {
    return m*pow(r, -a)*(1.0 + erf(-(r - t)/w))*0.5;
  }
    
public:
  
  /** Constructor
      @param r is the location of the radial bins
      @param d is the density computed from the bins
      @param expn is the density weighting (expn=1 is unweighted)
  */
  OnePowerTrunc(std::vector<double>& r, std::vector<double>& d,
		double expn=1.0)
  {
    ID = "OnePowerTrunc";

    this->r    = r;
    this->d    = d;
    this->ndim = 4;
    this->expn = expn;
  }

  virtual double operator()(std::vector<double>& p);

  virtual double fit(std::vector<double>& p, double r) {
    m = scale(p[0], 100.0);
    a = scale(p[1], 4.0);
    t = scale(p[2], 2.0);
    w = scale(p[3], 2.0);
    return fct(r);
  }

  std::vector<double> TrueParam(const std::vector<double> p)
  {
    return {scale(p[0], 100.0),
	    scale(p[1], 4.0 ),
	    scale(p[2], 2.0 ), 
	    scale(p[3], 2.0 )};
  }

  std::vector<double> ScaleParam(const std::vector<double> p)
  {
    return {invscale(p[0], 40.0),
	    invscale(p[1], 4.0 ),
	    invscale(p[2], 4.0 ), 
	    invscale(p[3], 2.0 )};
  }

  std::vector<std::string> ParamLabels()
  {
    return {"Norm", "Alpha", "r_trunc", "w_trunc"};
  }

};


double OnePowerTrunc::operator()(std::vector<double>& params)
{
  m = scale(params[0], 100.0);
  a = scale(params[1], 4.0  );
  t = scale(params[2], 2.0  );
  w = scale(params[3], 2.0  );

  // Parameter sanity

  if (t < 0.001)  return std::numeric_limits<double>::max();
  if (t > 1.0  )  return std::numeric_limits<double>::max();
  if (w < 0.001)  return std::numeric_limits<double>::max();
  if (w > 0.4  )  return std::numeric_limits<double>::max();

  double val, ans = 0.0;

  int num = r.size();

  for (int i=0; i<num; i++) {
    double val = fct(r[i]);
    double del = (d[i] - val)/pow(val+1.0e-24, expn);
    ans += del*del;
  }
  
  return ans;
}


int main(int argc, char** argv)
{
  std::string binfile, modfile, model, type;
  int maxiter;
  double rate, expon;

  // Option parsing
  //
  cxxopts::Options options(argv[0], "Read a PSP file (either OUT or SPL) and print the metadata\n");

  options.add_options()
    ("h,help", "print this help message")
    ("d,data", "Binned positions and density",
     cxxopts::value<std::string>(binfile)->default_value("data.dat"))
    ("o,model", "Output model file",
     cxxopts::value<std::string>(modfile)->default_value("model.fit"))
    ("t,type", "Model type (currently: OnePowerTrunc, TwoPowerTrunc",
     cxxopts::value<std::string>(type)->default_value("OnePowerTrunc"))
    ("n,iterations", "Maximum number of SA iterations",
     cxxopts::value<int>(maxiter)->default_value("300000"))
    ("r,learnrate", "Simulation learning rate",
     cxxopts::value<double>(rate)->default_value("0.01"))
    ("e,expon", "Density weighting exponent (1.0 is no weighting)",
     cxxopts::value<double>(expon)->default_value("1.0"))
    ;

  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("help")) {
    std::cout << options.help() << std::endl;
    return 1;
  }

  // Read in the model
  std::ifstream in(binfile);
  std::vector<double> r, d;
  while (in) {
    std::string s;
    std::getline(in, s);
    std::istringstream sin(s);
    double radius, density;
    sin >> radius >> density;
    if (sin.good()) {
      r.push_back(radius);
      d.push_back(density);
    }
  }

  // Assign model
  //
  std::string mtype(type);	// Lower-case conversion
  std::transform(mtype.begin(), mtype.end(), mtype.begin(),
		 [](unsigned char c){ return std::tolower(c); });

  // The model
  std::shared_ptr<ModelFit> chi2;

  // SA parameters in [-inf, inf]
  std::vector<double> sa_param;

  if (mtype.find("twopowertrunc") == 0) {
    chi2 = std::make_shared<TwoPowerTrunc>(r, d, expon);
    sa_param = chi2->ScaleParam({1.0, 1.0, 2.0, 0.05, 0.6, 0.1});
  }
  else if (mtype.find("onepowertrunc") == 0) {
    chi2 = std::make_shared<OnePowerTrunc>(r, d, expon);
    sa_param = chi2->ScaleParam({1.0, 1.0, 0.6, 0.1});
  }
  else {
    std::cerr << "ModelFit: no such model <" << type << ">" << std::endl;
    exit(-1);
  }

  // Create the optimization instance
  //
				// There must be a better way than
				// lambda binding, but maybe not . . .
  auto fct = [chi2](std::vector<double>& x) { return (*chi2)(x); };

  SimAnn sa(fct, chi2->ModelDim());

  // Initial values
  //
  sa.learning_rate(rate);
  sa.initial(sa_param);
  sa.melt();

  std::cout << "#" << std::endl
	    << "# Model     [" << chi2->getID() << "]"       << std::endl
	    << "# SimAnneal [" << sa.anneal(maxiter) << "] " << std::endl
	    << "#" << std::endl;

  // Compute and retrieve the result
  //
  sa_param = sa.optimum();
  std::vector<double> param = chi2->TrueParam(sa_param);

  for (int i=0; i<chi2->ModelDim(); i++) {
    std::ostringstream sout; sout << "# " << chi2->ParamLabels()[i];
    std::cout << std::setw(10) << std::left << sout.str() << std::right
	      << std::setw(14) << param[i] << std::endl;
  }

  double z = (*chi2)(sa_param);
  std::cout << "#" << std::endl
	    << "# (energy = " << z << ")" << std::endl
	    << "#" << std::endl << std::endl;

  for (int j=0; j<r.size(); j++) {
    std::cout << std::setw(15) << r[j]
	      << std::setw(15) << d[j]
	      << std::setw(15) << chi2->fit(sa_param, r[j])
	      << std::endl;
  }

  return 0;
}
