//
// Embedded disk model
//

#include <cmath>
#include <string>
#include <memory>

#include "QPDistF.H"
#include "massmodel.H"

EmbeddedDiskModel::EmbeddedDiskModel (std::vector<AxiSymModPtr>& T, 
				      std::vector<double>& M_scale,
				      std::vector<double>& R_scale, 
				      int NUMBER)
{
  number = NUMBER;
  t = T;
  m_scale = M_scale;
  r_scale = R_scale;

  if (t[0]->dof() != 2) bomb ("First model must be a flat disk");
  dim = 2;

  ModelID = "EmbeddedDiskModel(";
  rmin = 0.0;
  rmax = 1.0e30;
  for (int i=0; i<number; i++) {
    if (i==0)
      ModelID += t[i]->ModelID;
    else
      ModelID += "," + t[i]->ModelID;

    rmin = rmin < t[i]->get_min_radius()*r_scale[i] ? t[i]->get_min_radius()*r_scale[i] : rmin;
    rmax = rmax > t[i]->get_max_radius()*r_scale[i] ? t[i]->get_max_radius()*r_scale[i] : rmax;
  }

  dist_defined = false;

}

EmbeddedDiskModel::~EmbeddedDiskModel()
{
  // NADA
}

double EmbeddedDiskModel::get_mass(const double r)
{
  return t[0]->get_mass(r/r_scale[0]) * m_scale[0];
}

double EmbeddedDiskModel::get_density(const double r)
{
  return 
    t[0]->get_density(r/r_scale[0]) * 
      m_scale[0]/(r_scale[0]*r_scale[0]*r_scale[0]);
}

double EmbeddedDiskModel::get_pot(const double r)
{
  double pot=0.0;
  for (int i=0; i<number; i++) pot += t[i]->get_pot(r/r_scale[i])*m_scale[i]/r_scale[i];
  return pot;
}

double EmbeddedDiskModel::get_dpot(const double r)
{
  double dpot=0.0;
  for (int i=0; i<number; i++) dpot += t[i]->get_dpot(r/r_scale[i])*m_scale[i]/(r_scale[i]*r_scale[i]);
  return dpot;
}

double EmbeddedDiskModel::get_dpot2(const double r)
{
  double dpot2=0.0;
  for (int i=0; i<number; i++) dpot2 += t[i]->get_dpot2(r/r_scale[i])*m_scale[i]/(r_scale[i]*r_scale[i]*r_scale[i]);
  return dpot2;
}

void EmbeddedDiskModel::get_pot_dpot(const double r, double& p, double& p2)
{
  double Pot, dPot;
  p = p2 = 0.0;
  for (int i=0; i<number; i++) {
    t[i]->get_pot_dpot(r/r_scale[i], Pot, dPot);
    p += Pot*m_scale[i]/r_scale[i];
    p2 += dPot*m_scale[i]/(r_scale[i]*r_scale[i]);
  }
}

void EmbeddedDiskModel::setup_df(int egrid, int kgrid, int mgrid,
				 double lambda, double alpha, double beta,
				 double gama, double sigma, double rmmax, 
				 double roff, double eoff, double koff, 
				 double kmin, double kmax,
				 int nint, int numt)
{
  if (rmmax<0.0) rmmax = rmax;
  double remax = rmax;

  df = std::make_shared<QPDistF>(shared_from_this(),
				 rmmax, remax, egrid, kgrid, mgrid, 
				   sigma, lambda, alpha, beta, gama,
				 roff, eoff, koff, kmin, kmax, nint, numt);

  dist_defined = true;
}

void EmbeddedDiskModel::setup_df(string &file)
{
  df = std::make_shared<QPDistF>(shared_from_this(), file);
  dist_defined = true;
}

