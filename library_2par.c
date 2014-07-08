//
//  main.c
//  itc_sim
//
//  Created by Elihu Ihms on 2/18/14.
//  Copyright (c) 2014 Elihu Ihms. All rights reserved.
//

#include "itc_model.h"
#include "itc_sim.h"
#include "energies_2par.h"

int calc( double* P, double* L, double* Q, int n, int size, int cyclic, double temp, double *params )
{
	struct mWorkspace model;
	struct sWorkspace sim;

	model.size		= size;
	model.cyclic	= cyclic;
	model.temp		= temp;

	setupModelWorkspace(&model);
	setupSimWorkspace(&sim,model);

	assignEnergies(&model, &sim, params[0], params[1]);

	for(int i=0; i<n; i++)
	{
		model.Ptot = P[i];
		model.Ltot = L[i];
		setFree( &model );

		Q[i] = getQ( sim, model );
	}

	freeSimWorkspace(sim,model);
	freeModelWorkspace(model);

    return 0;
}