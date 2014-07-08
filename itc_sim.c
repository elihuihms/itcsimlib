/*
 *  itc_sim.c
 *  itcsim
 *
 *  Created by Elihu Ihms on 2/24/14.
 *  Copyright 2014 The Ohio State University. All rights reserved.
 *
 */

#include "itc_sim.h"
#include <stdlib.h>
#include <math.h>

int setupSimWorkspace( struct sWorkspace *sim, struct mWorkspace model )
{
	sim->model = &model;
	sim->enthalpies = (double *)malloc(pow(2, model.size) * sizeof(double) );
	if( sim->enthalpies == NULL )
		return 1;

	for(int i=0; i<pow(2, model.size); i++)
		sim->enthalpies[i]=0;
	
	return 0;
}

int freeSimWorkspace( struct sWorkspace sim, struct mWorkspace model )
{
	free( sim.enthalpies );
	return 0;
}

double getQ( struct sWorkspace sim, struct mWorkspace model )
{
	double Q=0;
	for(int i=0; i<pow(2,model.size); i++)
		Q += sim.enthalpies[i] * model.probs[i];

	return Q;
}