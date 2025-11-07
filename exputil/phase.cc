#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "numerical.H"
#include "phase.H"

using namespace std;

Eigen::VectorXd binread_vectorX(std::istream& in);
void binwrite_vectorX(Eigen::VectorXd& x, std::ostream& out);

Eigen::Vector3d binread_vector3(std::istream& in);
void binwrite_vector3(Eigen::Vector3d& x, std::ostream& out);

Frame_Rotation Frame;

double Phase::tolerance=1.e-5;
double Phase::xscale=1.0;
double Phase::vscale=1.0;
int Phase::potential_defined = 0;
char Phase::integrator_name[8];
ode_integrator Phase::integrator = rkqc;
symp_integrator Phase::mapping = sia4;
force_func_ptr Phase::Force = Phase::Default_Force;
pot_func_ptr Phase::Potential = Phase::Default_Potential;

//	Eigen::Vector3d Force(double, Eigen::Vector3d &, Eigen::Vector3d &);

Phase &Phase::operator=(const Phase &p)
{
  v = p.v;
  x = p.x;
  t = p.t;
  work = p.work;
  m = p.m;
  
  return *this;
}


Phase::Phase(void)
{
  work = 0.0;
  t=0.0;
  m=1.0;
}


Phase::Phase(const Phase &p)
{
  v = p.v;
  x = p.x;
  t = p.t;
  m = p.m;
  work = p.work;
}

Eigen::Vector3d Phase::Default_Force(double, Eigen::Vector3d &, Eigen::Vector3d &)
{
  std::cerr << "no potential registered" << std::endl;
  std::cerr << "you must call Phase::Register_Potential()" << std::endl;
  exit(0);
}

double Phase::Default_Potential(double, Eigen::Vector3d &)
{
  cerr << "no potential registered" << endl;
  cerr << "you must call Phase::Register_Potential()" << endl;
  exit(0);
}

double Phase::Jacobi(void)
{
  double centrifugal, w;
  
  w = Frame.omega;
  centrifugal = 0.5*m*w*w*(x[0]*x[0] + x[1]*x[1]);
  
  return Energy() - centrifugal;
}


double Phase::Accuracy(Phase &initial)
{
  return (m*work - 0.5*m*(v.dot(v) - initial.v.dot(initial.v)))/
    (0.5*m*initial.v.dot(initial.v) + TINY);
}



void Phase::print(ostream& out)
{
  out << t << " "
      << m << " "
      << x[0] << " "
      << x[1] << " "
      << x[2] << " "
      << v[0] << " "
      << v[1] << " "
      << v[2] << " "
      << '\n';
}

void Newton(double t, Eigen::VectorXd& u, Eigen::VectorXd& dudt)
{
  static Eigen::Vector3d f, x, v;
  
  for (int i=0; i<3; i++)
    {
      dudt[i] = u[i+3];
      v[i] = u[i+3];
      x[i] = u[i];
    }
  
  //	f = Force(t, x, v);
  f = (*Phase::Force)(t, x, v);
  
  dudt[3] = f[0];
  dudt[4] = f[1];
  dudt[5] = f[2];
  
  dudt[6] = f.dot(v);
}


void SNewton(double t, Eigen::VectorXd& xx, Eigen::VectorXd& vv,
	     Eigen::VectorXd& dudt)
{
  static Eigen::Vector3d f, x, v;
  
  for (int i=0; i<3; i++)
    {
      x[i] = xx[i];
      v[i] = vv[i];
    }
  
  //	f = Force(t, x, v);
  f = (*Phase::Force)(t, x, v);
  
  for (int i=0; i<3; i++) dudt[i] = f[i];
  
  dudt[6] = f.dot(v);
}


Phase Phase::integrate_to(double tfin)
{
  static Phase newpoint;
  double dt;
  
  newpoint = *this;
  
  dt = Sign(tfin-t) * 0.01 * get_Xscale()/get_Vscale();
  
  do 
    {
      if (fabs(tfin-newpoint.t) <= fabs(dt)) dt = tfin-newpoint.t;
      if (fabs(dt) < TINY) break;
      newpoint = newpoint.advance(dt);
    } while ((newpoint.t < tfin)==(t < tfin) && 
	     sqrt(newpoint.Position().dot(newpoint.Position())) < 1.0e6 * xscale );
  
  return newpoint;
}

Phase Phase::integrate_to_symplectic(double tfin, double scale)
{
  static Phase newpoint;
  double dt;
  
  newpoint = *this;
  
  dt = Sign(tfin-t) * scale * get_Xscale()/get_Vscale();
  
  do 
    {
      if (fabs(tfin-newpoint.t) <= fabs(dt)) dt = tfin-newpoint.t;
      if (fabs(dt) < TINY) break;
      newpoint = newpoint.advance_symplectic(dt);
    } while ((newpoint.t < tfin)==(t < tfin) && 
	     sqrt(newpoint.Position().dot(newpoint.Position())) < 1.0e6 * xscale );
  
  return newpoint;
}


Phase Phase::advance(double& dt)
{
  static Phase newpoint;
  double tt, hnew, hold;
  
  Eigen::VectorXd u(8), dudt(8), uscale(8);
  
  for (int i=0; i<3; i++)
    {
      u[i]   = x[i];
      u[i+3] = v[i];
    }
  u[6] = work;
  tt = t;
  
  Newton(tt, u, dudt);
  
  /* scale for phase_space position */
  
  for (int i=0; i<6; i++) 
    uscale[i] = fabs(u[i]) + fabs(dt*dudt[i]) + TINY;
  
  
  // scale for work done
  
  uscale[6] = fabs(u[6]) + fabs(dt*dudt[6]) + 1.0;
  
  integrator(u, dudt, 7, tt, dt, 
	     tolerance, uscale, hold, hnew, Newton);
  
  dt = hnew;
  
  
  for (int i=0; i<3; i++)
    {
      newpoint.x[i] = u[i];
      newpoint.v[i] = u[i+3];
    }
  newpoint.work = u[6];
  newpoint.t = tt;
  newpoint.m = m;
  
  return newpoint;
}


Phase Phase::advance_symplectic(double& dt)
{
  static Phase newpoint;

  Eigen::VectorXd x0(4), v0(4);
  Eigen::VectorXd x1(4), v1(4);
  
  for (int i=0; i<3; i++) {
    x0[i] = x[i];
    v0[i] = v[i];
  }
  
  
  mapping(x0, v0, x1, v1, t, dt, SNewton);
  
  for (int i=0; i<3; i++) {
    newpoint.x[i] = x1[i];
    newpoint.v[i] = v1[i];
  }
  newpoint.t = t + dt;
  newpoint.m = m;
  
  return newpoint;
}



Phase Phase::rotate_frame(double omega)
{
  static Phase newpoint;
  Eigen::Vector3d W;
  
  W.setZero();
  W[2] = omega;
  
  newpoint.x = x;
  newpoint.v = v - W.cross(x);
  newpoint.m = m;
  newpoint.t = t;
  newpoint.work = work;
  
  
  return newpoint;
}


Phase Phase::rotate_view(double theta)
{
  static Phase newpoint;
  Eigen::Matrix3d R;
  
  
  R = z_Rotation(theta);	
  
  newpoint.x = R * x;
  newpoint.v = R * v;
  newpoint.m = m;
  newpoint.t = t;
  newpoint.work = work;
  
  return newpoint;
}


Eigen::Matrix3d z_Rotation(double theta)
{
  double c, s;
  Eigen::Matrix3d R;
  R.setZero();
  
  c = cos(theta);
  s = sin(theta);
  
  R(0, 0) = c;
  R(0, 1) = -s;
  R(1, 0) = s;
  R(1, 1) = c;
  R(2, 2) = 1.0;
  
  return R;
}

void Phase::binwrite(std::ostream& fp)
{
  binwrite_vector3(x, fp);
  binwrite_vector3(v, fp);
  fp.write((char *)&m, sizeof(double));
  fp.write((char *)&t, sizeof(double));
  fp.write((char *)&work, sizeof(double));
}

Eigen::VectorXd binread_vectorX(std::istream& in)
{
  int sz;
  in.read((char *)&sz, sizeof(int));
  
  Eigen::VectorXd x(sz);
  in.read((char *)x.data(), sz*sizeof(double));
  
  return x;
}

void binwrite_vectorX(Eigen::VectorXd& x, std::ostream& out)
{
  int sz = x.size();
  out.write((const char *)&sz, sizeof(int));
  
  out.write((const char *)x.data(), sz*sizeof(double));
}

Eigen::Vector3d binread_vector3(std::istream& in)
{
  int sz;
  in.read((char *)&sz, sizeof(int));
  
  Eigen::Vector3d x;
  in.read((char *)x.data(), 3*sizeof(double));
  
  return x;
}

void binwrite_vector3(Eigen::Vector3d& x, std::ostream& out)
{
  int sz = 3;
  out.write((const char *)&sz, sizeof(int));
  
  out.write((const char *)x.data(), sz*sizeof(double));
}


void Phase::binread(istream& fp)
{
  x = binread_vector3(fp);
  v = binread_vector3(fp);
  fp.read((char *)&m, sizeof(double));
  fp.read((char *)&t, sizeof(double));
  fp.read((char *)&work, sizeof(double));
}


void Phase::set_Integrator(char *s)
{
  if (!strcmp(s, "rk4"))
    {
      integrator = rkqc;
      strcpy(integrator_name, s);
      return;
    }
  
  if (!strcmp(s, "bs"))
    {
      integrator = bsstep;
      strcpy(integrator_name, s);
      return;
    }
  
  cerr << "did not recognise integrator: " << s << "\n";
  cerr << "known integrators are:\n";
  cerr << "bs = Bulirsch-Stoer\nrk4 = Runge-Kutta 4th\n";
  exit(0);
}

char *Phase::get_Integrator(void)
{
  return integrator_name;
}
