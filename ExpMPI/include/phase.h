class Three_Vector;
class Phase;
class Ensemble;

typedef Three_Vector (*force_func_ptr)(double, Three_Vector &, Three_Vector &);
typedef double (*pot_func_ptr)(double, Three_Vector &);


class Phase
{
protected:
	Three_Vector x;			/* position vector */
	Three_Vector v;			/* velocity vector */
	double m;			/* mass */
	double t;			/* current time */
	double work;			/* total work done so far */

	static double vscale;		/* typical length scale */
	static double xscale;		/* typical speed scale */
	static double tolerance;	/* desired integration accuracy */
	static int potential_defined;	/* has a potential been set? */
	static ode_integrator integrator;	/* type of integrator */
	static symp_integrator mapping;	/* type of sympletic integrator */
	static char integrator_name[8];	/* name of ODE integrator */
	static char symp_integrator_name[8];	/* name of SYMP integrator */
public:
	/* constructors */

	Phase(void);
	Phase(const Phase &);


	/* specification of model */

	static void Register_Potential(pot_func_ptr p, force_func_ptr f)
		{potential_defined=1; Potential = p; Force = f;}
	static Three_Vector Default_Force(double, Three_Vector &, Three_Vector &);
	static double Default_Potential(double, Three_Vector &);
	static pot_func_ptr Potential;
	static force_func_ptr Force;


	
	/* assignment operator */

	Phase &operator=(const Phase &);



	/* access to protected data */

	Three_Vector &Position(void)	{return x;}
	Three_Vector &Velocity(void)	{return v;}
	double &Position(int i)	{return x[i];}
	double &Velocity(int i)	{return v[i];}
	double &Time(void)	{return t;}
	double &Work(void) 	{return work;}
	double &Mass(void)	{return m;}



	/* set and get the scale dimensions and tolerance */

	void set_Tolerance(double tol)	{tolerance = tol;}
	void set_Xscale(double xx) 	{xscale = xx;}
	void set_Vscale(double vv) 	{vscale = vv;}

	double get_Tolerance(void) 	{return tolerance;}
	double get_Xscale(void) 	{return xscale;}
	double get_Vscale(void) 	{return vscale;}


	/* set or retreive the ode integrator type */

	void set_Integrator(char *);
	char *get_Integrator(void);




	/* print a phase-space point */

	void print(ostream& out);


	/* rotate frame or point about z axis */

	Phase rotate_frame(double);
	Phase rotate_view(double);



	/* various functions of the coordinates */

	double Speed(void) 	{return sqrt(v*v);}
	double Energy(void)	{return 0.5*m*v*v + m*(*Potential)(t, x);}
	double Jacobi(void);
	double Accuracy(Phase &);
	Three_Vector Angular_Momentum(void)	{return m*Cross(x, v);}



	/* numerical integration */

	Phase integrate_to(double);
	Phase advance(double *);
	Phase x_return_map(void);

	/* symplectic mapping */

	Phase integrate_to_symplectic(double, double);
	Phase advance_symplectic(double *);

	/* predictor-corrector stuff */

	void set_PC(double, int);
	Phase PC_advance(double);


	void binread(istream& fp);
	void binwrite(ostream& fp);

	friend class Ensemble;
};


void Newton(double, double *, double *);
double v_circ(double);
double vertical_kappa(double);
double epicyclic_kappa(double);
double Oort_A(double);
double Oort_B(double);

class Frame_Rotation
{
public:
	double omega, omegasq, omega2;
	void set_omega(double w) 	
	{
		omega = w;
		omegasq = omega*omega;
		omega2 = 2.0*omega;
	}
	void set_corotation(double rc)
	{
		if (rc <= 0.0) omega=0.0;
		else omega = v_circ(rc)/rc;
		omegasq = omega*omega;
		omega2 = 2.0*omega;
	}
};


Matrix z_Rotation(double);
typedef int (*freezing_function)(Phase &);



class Ensemble
{
protected:
	int Nstars;
	Phase *stars;
	double t;

public:
	Ensemble(void);
	~Ensemble(void);
	Ensemble(int);
	Ensemble(Ensemble &);

	void setsize(int);
	int get_Nstars(void) 	{return Nstars;}

	Ensemble &operator=(Ensemble &);

	Phase &operator[](int);

	Ensemble integrate_to(double);
	Ensemble Step_Positions(double);
	Ensemble Step_Velocities(double, Three_Vector *);
	Ensemble integrate_to(double, freezing_function);

	void write_log(ostream&);
	void write_snapshot(ostream&);
	void read_snapshot(istream&);
	void read_tidesfile(istream&);
	void write_orbits(ostream&);
	void read_orbits(istream&);

	void settime(double);
	void set_masses(double);
	void scale_masses(double);
	void scale_positions(double);
	void scale_speeds(double);
	void rotate_frame(double);
	void rotate_view(double);

	void translate(Three_Vector &);

	double total_Energy(void);
	double total_Kinetic(void);
	double total_Potential(void);
	double Virial(void);
	Three_Vector CM_Position(void);
	Three_Vector CM_Velocity(void);
	Three_Vector total_Angular_Momentum(void);

	Matrix Moment_Tensor(double);
	Matrix Inertia_Tensor(double);
	Three_Vector Principal_Axes(double, Matrix &);
	Three_Vector Solid_Body_Frequency(double);

	void binread(istream& fp);
	void binwrite(ostream& fp);

};



double inverse_interp_zero(double, double, double, double,
	double, double );

void forward_interp(double, double, double, double,
	double, double, double, double *, double *);















	
		

		
