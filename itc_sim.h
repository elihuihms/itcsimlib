/*
 *  itc_sim.h
 *  itcsim
 *
 *  Created by Elihu Ihms on 2/24/14.
 *  Copyright 2014 The Ohio State University. All rights reserved.
 *
 */

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
