




#include "expand.h"



/*
	reduce masses of all particles according to the
	function
		
				  	             -t/t_mloss
		m(t) = m(0) (1 - frac_mloss (1.0 - e	         ))

	The initial masses are stored in the array initial_mass[],
	declared globally.


To add this code, it was necessary to modify the following
files:
	parse.c		(reading parameters frac_mloss and t_mloss)
	init.c		(store initial mass of all particles)
	step.c		(call mass_loss() at the beginning of each step)
	expand.h	(declarations for *initial_mass, frac_mloss, and
			t_mloss)



	KL 7/3/92
*/



void mass_loss(void)
{
	int i;
	double f;


	f = 1.0 - frac_mloss * (1.0 - exp(-tnow/t_mloss));

	for (i=1; i<=nbodies; i++)
	{
		mass[i] = f * initial_mass[i];
	}
}
		











	


