#pragma implementation

#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <numerical.h>
#include <Vector.h>
#include <phase.h>

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

//	Three_Vector Force(double, Three_Vector &, Three_Vector &);

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

Three_Vector Phase::Default_Force(double, Three_Vector &, Three_Vector &)
{
	cerr << "no potential registered" << endl;
	cerr << "you must call Phase::Register_Potential()" << endl;
	exit(0);
}

double Phase::Default_Potential(double, Three_Vector &)
{
	cerr << "no potential registered" << endl;
	cerr << "you must call Phase::Register_Potential()" << endl;
	exit(0);
}

double Phase::Jacobi(void)
{
	double centrifugal, w;

	w = Frame.omega;
	centrifugal = 0.5*m*w*w*(x[1]*x[1] + x[2]*x[2]);

	return Energy() - centrifugal;
}


double Phase::Accuracy(Phase &initial)
{
	return (m*work - 0.5*m*(v*v - initial.v*initial.v))/
		(0.5*m*initial.v*initial.v + TINY);
}



void Phase::print(ostream& out)
{
  out << t << " "
      << m << " "
      << x[1] << " "
      << x[2] << " "
      << x[3] << " "
      << v[1] << " "
      << v[2] << " "
      << v[3] << " "
      << '\n';
}

void Newton(double t, double *u, double *dudt)
{
	static Three_Vector f, x, v;
	int i;

	for (i=1; i<=3; i++)
	{
		dudt[i] = u[i+3];
		v[i] = u[i+3];
		x[i] = u[i];
	}

	//	f = Force(t, x, v);
	f = (*Phase::Force)(t, x, v);

	dudt[4] = f[1];
	dudt[5] = f[2];
	dudt[6] = f[3];

	dudt[7] = f * v;
}


void SNewton(double t, double *xx, double *vv, double *dudt)
{
	static Three_Vector f, x, v;
	int i;

	for (i=1; i<=3; i++)
	{
		x[i] = xx[i];
		v[i] = vv[i];
	}

	//	f = Force(t, x, v);
	f = (*Phase::Force)(t, x, v);

	dudt[1] = f[1];
	dudt[2] = f[2];
	dudt[3] = f[3];

	dudt[7] = f * v;
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
		newpoint = newpoint.advance(&dt);
	} while ((newpoint.t < tfin)==(t < tfin) && 
		sqrt(SQR(newpoint.Position())) < 1.0e6 * xscale );

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
		newpoint = newpoint.advance_symplectic(&dt);
	} while ((newpoint.t < tfin)==(t < tfin) && 
		sqrt(SQR(newpoint.Position())) < 1.0e6 * xscale );

	return newpoint;
}


Phase Phase::advance(double *dt)
{
	static Phase newpoint;
	double u[8], dudt[8], uscale[8];
	double tt, hnew, hold;
	int i;

	for (i=1; i<=3; i++)
	{
		u[i] = x[i];
		u[i+3] = v[i];
	}
	u[7] = work;
	tt = t;

	Newton(tt, u, dudt);

	/* scale for phase_space position */

	for (i=1; i<=6; i++) 
		uscale[i] = fabs(u[i]) + fabs(*dt*dudt[i]) + TINY;
	

	/* scale for work done */

	uscale[7] = fabs(u[7]) + fabs(*dt*dudt[7]) + 1.0;

	integrator(u, dudt, 7, &tt, *dt, 
		tolerance, uscale, &hold, &hnew, Newton);

	*dt = hnew;


	for (i=1; i<=3; i++)
	{
		newpoint.x[i] = u[i];
		newpoint.v[i] = u[i+3];
	}
	newpoint.work = u[7];
	newpoint.t = tt;
	newpoint.m = m;

	return newpoint;
}


Phase Phase::advance_symplectic(double *dt)
{
	static Phase newpoint;
	double x0[4], v0[4];
	double x1[4], v1[4];
	int i;

	for (i=1; i<=3; i++) {
	  x0[i] = x[i];
	  v0[i] = v[i];
	}


	mapping(x0, v0, x1, v1, t, *dt, SNewton);

	for (i=1; i<=3; i++) {
	  newpoint.x[i] = x1[i];
	  newpoint.v[i] = v1[i];
	}
	newpoint.t = t + *dt;
	newpoint.m = m;

	return newpoint;
}



Phase Phase::rotate_frame(double omega)
{
	static Phase newpoint;
	Three_Vector W;
	W.zero();
	W[3] = omega;

	newpoint.x = x;
	newpoint.v = v - Cross(W,x);
	newpoint.m = m;
	newpoint.t = t;
	newpoint.work = work;


	return newpoint;
}


Phase Phase::rotate_view(double theta)
{
	static Phase newpoint;
	Matrix R(1, 3, 1, 3);


	R = z_Rotation(theta);	

	newpoint.x = (Three_Vector) (R * (Vector) x);	
	newpoint.v = (Three_Vector) (R * (Vector) v);
	newpoint.m = m;
	newpoint.t = t;
	newpoint.work = work;

	return newpoint;
}


Matrix z_Rotation(double theta)
{
	double c, s;
	Matrix R(1, 3, 1, 3);
	R.zero();

	c = cos(theta);
	s = sin(theta);

	R[1][1] = c;
	R[1][2] = -s;
	R[2][1] = s;
	R[2][2] = c;
	R[3][3] = 1.0;

	return R;
}

void Phase::binwrite(ostream& fp)
{
	x.binwrite(fp);
	v.binwrite(fp);
	fp.write((char *)&m, sizeof(double));
	fp.write((char *)&t, sizeof(double));
	fp.write((char *)&work, sizeof(double));
}

void Phase::binread(istream& fp)
{
	x.binread(fp);
	v.binread(fp);
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







	
	



	





















