/*
	configuration energies/enthalpies as described in Ian's thesis
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

void assignEnergies(struct mWorkspace *w, struct sWorkspace *sim, double *params)
{
	/*
	param[0]		dG0
	param[1]		dG1
	param[2]		dG2
	param[3]		dH0
	param[4]		dH1
	param[5]		dH2
	*/

	/* assign energies for each configuration */
	for(int i=0; i<pow(2, w->size); i++)
	{
		w->energies[i]=0;
		sim->enthalpies[i]=0;

		/* start w/ presence of ligand dG */
		for(int j=0; j<w->size; j++)
		{
			if( w->configs[i][j] > 0 )
			{

				if( w->configs[i][ permute(j+1,w->size) ] > 0 && w->configs[i][ permute(j-1,w->size) ] > 0 )
				{
					w->energies[i]		+= params[2];
					sim->enthalpies[i]	+= params[5];
				}
				else if( w->configs[i][ permute(j+1,w->size) ] > 0 || w->configs[i][ permute(j-1,w->size) ] > 0 )
				{
					w->energies[i]		+= params[1];
					sim->enthalpies[i]	+= params[4];
				}
				else
				{
					w->energies[i]		+= params[0];
					sim->enthalpies[i]	+= params[3];
				}
			}
		}
	}
	return;
}