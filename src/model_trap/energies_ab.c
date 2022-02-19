/*
	configuration energies/enthalpies as described in Elihu's thesis
*/

#include "energies.h"
#include "itc_model.h"
#include "itc_sim.h"
#include <math.h>

int permute( int i, int n );
int permute( int i, int n )
{
	if( i < 0 )
		return i + n;
	else if( i >= n )
		return i % n;
	else
		return i;
}

void assignEnergies( struct mWorkspace *w, struct sWorkspace *sim, double *params )
{
	/*

	Configurations:
	 1				(dGA if i == even)
	 1				(dGB if i == odd)
	010				(+0)
	011 or 110		(+dGC)
	111				(+dGC +dGC)

	param[0]		dGA
	param[1]		dGB
	param[2]		dGC
	param[3]		dHA
	param[4]		dHB
	param[5]		dHC
	*/

	/* assign energies for each configuration */
	for(int i=0; i<pow(2, w->size); i++)
	{
		w->energies[i]=0;
		sim->enthalpies[i]=0;

		for(int j=0; j<w->size; j++)
		{
			if( w->configs[i][j] > 0 )
			{
				/*
					Set the intrinsic energy of the binding site.
				*/
				if (i % 2 == 0) {
					/* Even (A) sites. */
					w->energies[i]		+= params[0];
					sim->enthalpies[i]	+= params[3];
				} else {
					/* Odd (B) sites. */
					w->energies[i]		+= params[1];
					sim->enthalpies[i]	+= params[4];
				}
				
				/*
					Add any coupling energies.
				*/
				if( w->configs[i][ permute(j+1, w->size) ] > 0 ) {
					/* The right-hand site is occupied. */
					w->energies[i]		+= params[2];
					sim->enthalpies[i]	+= params[5];
				}
				/*
					Coupling energies are not exclusive!
				*/
				if( w->configs[i][ permute(j-1, w->size) ] > 0 ) {
					/* The left-hand site is occupied. */
					w->energies[i]		+= params[2];
					sim->enthalpies[i]	+= params[5];
				}
			}
		}
	}

	return;
}
