#include "itc_model.h"
#include "itc_sim.h"
#include "energies.h"

struct mWorkspace model;
struct sWorkspace sim;

int setup(int size, int cyclic )
{
	model.size		= size;
	model.cyclic	= cyclic;

	if( setupModelWorkspace(&model) > 0 )
		return 0;
	if( setupSimWorkspace(&sim,model) > 0 )
		return 0;

	return 0;
}

int calc( int n, double temp, double* P, double* L, double* Q, double *params )
{
	int status = 0;

	/*
		n		: Number of points (size of P,L, and Q arrays)
		temp	: Experimental temperature
		P		: Array of (total) protein concentrations
		L		: Array of (total) ligand concentrations
		Q		: Array to hold (total) enthalpic heats at each set of conditions
		params	: Array of model-specific parameters
	*/

	model.temp = temp;
	assignEnergies(&model, &sim, params);

	for(int i=0; i<n; i++)
	{
		model.Ptot = P[i];
		model.Ltot = L[i];

		status = setFree( &model );
		if( status != 0 )
			return status;

		Q[i] = getQ( sim, model );
	}

    return status;
}

int close()
{
	freeSimWorkspace(sim,model);
	freeModelWorkspace(model);

	return 0;
}