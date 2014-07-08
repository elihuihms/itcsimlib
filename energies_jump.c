/*

*/

#include "itc_model.h"
#include "itc_sim.h"

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

void assignEnergies(
	struct mWorkspace *w, struct sWorkspace *sim,
	double Gbind, double Gnear, double Gfar,
	double Hnone, double Hnear, double Hfar )
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
				w->energies[i]			+= Gbind; /* presence of ligand dG */
				sim->enthalpies[i]		+= Hnone; /* binding to site enthalpy */

				if( w->configs[i][ permute(j+1, w->size) ] > 0 ){
					w->energies[i]		+= Gnear;
					sim->enthalpies[i]	+= Hnear;
				}

				if( w->configs[i][ permute(j-1, w->size) ] > 0 ){
					w->energies[i]		+= Gnear;
					sim->enthalpies[i]	+= Hnear;
				}

				if( w->configs[i][ permute(j+2, w->size) ] > 0 ){
					w->energies[i]		+= Gfar;
					sim->enthalpies[i]	+= Hfar;
				}

				if( w->configs[i][ permute(j-2, w->size) ] > 0 ){
					w->energies[i]		+= Gfar;
					sim->enthalpies[i]	+= Hfar;
				}
			}
		}
	}

	return;
}