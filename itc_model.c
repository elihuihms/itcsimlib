/*
 *  stat_thermo.c
 *  itcsim
 *
 *  Created by Elihu Ihms on 2/24/14.
 *  Copyright 2014 The Ohio State University. All rights reserved.
 *
 */

#include "itc_model.h"

#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

int setupModelWorkspace( struct mWorkspace *w )
{
	w->energies = (double *)malloc(pow(2, w->size) * sizeof(double) );
	if( w->energies == NULL )
		return 1;
	w->probs = (double *)malloc(pow(2, w->size) * sizeof(double) );
	if( w->probs == NULL )
		return 1;
	w->bound = (int *)malloc(pow(2, w->size) * sizeof(int) );
	if( w->bound == NULL )
		return 1;

	w->configs = (int **)malloc(pow(2, w->size) * sizeof(int*) );
	if( w->configs == NULL )
		return 1;

	for(int i=0; i<pow(2, w->size); i++){
		w->configs[i] = (int *)malloc( w->size * sizeof(int) );
		if( w->configs[i] == NULL )
			return 1;
	}

	int counter;
	for(int i=0; i<pow(2, w->size); i++)
	{
		w->energies[i]=0;
		w->bound[i]=0;

		/* use the binary representation of the number to describe the config */
		counter=i;
		for (int j = w->size -1; j >= 0; --j){
			w->configs[i][j] = (counter & 1);
			counter >>= 1;
		}

		/* determine the number of occupied sites in each configuration */
		for(int j=0; j<w->size; j++)
			if( w->configs[i][j] > 0 )
				w->bound[i]++;
	}

	/* initialize the root finder */
	w->fsolver_s = gsl_root_fsolver_alloc( gsl_root_fsolver_brent );
	w->fsolver_F.function = &getFree;

	return 0;
}

int freeModelWorkspace( struct mWorkspace w )
{
	for(int i=0; i<w.size; i++){
		free( w.configs[i] );
	}
	free( w.configs );
	free( w.energies );
	free( w.probs );
	free( w.bound );

	gsl_root_fsolver_free( w.fsolver_s );

	return 0;
}

void setProbabilities( struct mWorkspace *w )
{
	double	R = 8.3144621; // J/(K*mol)
	double	sum = 0;

	for(int i=0; i<pow(2,w->size); i++)
	{
		w->probs[i] = exp( (-1 * w->energies[i]) / ( R * w->temp ) ) * pow(w->Lfree,w->bound[i]);
		sum += w->probs[i];
	}

	for(int i=0; i<pow(2,w->size); i++)
		w->probs[i] /= sum;

	return;
}

double getOccupation( struct mWorkspace w, int bound )
{
	double	ret=0;
	for(int i=0; i<pow(2,w.size); i++)
		if(w.bound[i]==bound)
			ret += w.probs[i];

	return ret;
}

double getFree( double Lfree, void *params )
{
	struct mWorkspace *w = (struct mWorkspace *)params;

	w->Lfree = Lfree;
	setProbabilities( w ); /* get the bound fraction of each protein state */

	double bound = 0; /* concentration of sites in bound state */
	for(int i=0; i<pow(2,w->size); i++)
		bound += w->probs[i] * w->Ptot * w->bound[i];

	return w->Ltot - (bound + w->Lfree);
}

int setFree( struct mWorkspace *w )
{
	int status;
	int iter = 0, max_iter = 100;

	w->fsolver_F.params = w;
	gsl_root_fsolver_set( w->fsolver_s, &w->fsolver_F, 0, w->Ltot );

	double r,a,b;
	do
	{
		iter++;
		status = gsl_root_fsolver_iterate( w->fsolver_s );
		r = gsl_root_fsolver_root( w->fsolver_s );
		a = gsl_root_fsolver_x_lower( w->fsolver_s );
		b = gsl_root_fsolver_x_upper( w->fsolver_s );

		status = gsl_root_test_interval( a, b, 0.0, 0.0 );

	}while( status == GSL_CONTINUE && iter < max_iter );

	return status;
}

double getNbar( struct mWorkspace w )
{
	double bound=0; /* concentration of sites in bound state */
	for(int i=0; i<pow(2,w.size); i++)
		bound += w.probs[i] * w.Ptot * w.bound[i];

	return bound/(w.Ptot*w.size);
}

void printConfig( struct mWorkspace w, int index )
{
	for(int i=0; i<w.size; i++)
		printf("%i",w.configs[index][i]);
	printf("\n");
	return;
}