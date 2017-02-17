#ifndef ITC_SIM_H
#define ITC_SIM_H

#include "itc_model.h"

struct sWorkspace
{
	struct	mWorkspace *model;
	double*	enthalpies;
};

int setupSimWorkspace( struct sWorkspace *sim, struct mWorkspace model );

int freeSimWorkspace( struct sWorkspace sim, struct mWorkspace model );

double getQ( struct sWorkspace sim, struct mWorkspace model );

#endif
