#ifndef _itc_sim_h_
#define _itc_sim_h_

#include "itc_model.h"

struct sWorkspace
{
	struct	mWorkspace *model;
	double*	enthalpies;
} sWorkspace;

int setupSimWorkspace( struct sWorkspace *sim, struct mWorkspace model );

int freeSimWorkspace( struct sWorkspace sim, struct mWorkspace model );

double getQ( struct sWorkspace sim, struct mWorkspace model );

#endif
