#ifndef ITC_MODEL_H
#define ITC_MODEL_H

#include <gsl/gsl_roots.h>

struct mWorkspace
{
	int		size;
	double	temp;
	double	Ptot;
	double	Ltot;
	double	Pfree;
	double	Lfree;
	int		cyclic;
	int**	configs;
	int*	bound;
	double* energies;
	double*	probs;
	gsl_root_fsolver*		fsolver_s;
	gsl_function			fsolver_F;
};

int setupModelWorkspace( struct mWorkspace *w );

int freeModelWorkspace( struct mWorkspace w );

void setProbabilities( struct mWorkspace *w );

double getOccupation( struct mWorkspace w, int bound );

double getFree( double Lfree, void *params );

int setFree( struct mWorkspace *w );

double getNbar( struct mWorkspace w );

void printConfig( struct mWorkspace w, int index );

#endif