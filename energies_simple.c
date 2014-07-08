/*
	configuration energies/enthalpies as described in Ian's thesis
*/

#include "itc_model.h"
#include "itc_sim.h"

void assignEnergies(
	struct mWorkspace *w, struct sWorkspace *sim,
	double Gbind,
	double Hbind )
{
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
				w->energies[i]		+= Gbind;
				sim->enthalpies[i]	+= Hbind;
			}
		}
	}
	return;
}