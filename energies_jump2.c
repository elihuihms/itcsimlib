/*

*/

#include "itc_model.h"
#include "itc_sim.h"

int permute( int i, int n );
int permute( int i, int n )
{
	if( i < 0 )
		return i + (i % n);
	else if( i >= n )
		return i % n;
	else
		return i;
}

void assignEnergies(
	struct mWorkspace *w, struct sWorkspace *sim,
	double Gbind, double Ga, double Gb, double Gc,
	double Hnone, double Ha, double Hb, double Hc )
{
	double tmp1, tmp2;

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
				tmp1 = Gbind; /* presence of ligand dG */
				tmp2 = Hnone; /* binding to site enthalpy */

				if( w->configs[i][ permute(j+1, w->size) ] > 0 ){
					tmp1	*= Ga;
					tmp2	*= Ha;
				}

				if( w->configs[i][ permute(j-1, w->size) ] > 0 ){
					tmp1	*= Ga;
					tmp2	*= Ha;
				}

				if( w->configs[i][ permute(j+2, w->size) ] > 0 ){
					tmp1	*= Gb;
					tmp2	*= Hb;
				}

				if( w->configs[i][ permute(j-2, w->size) ] > 0 ){
					tmp1	*= Gb;
					tmp2	*= Hb;
				}

				if( w->configs[i][ permute(j+3, w->size) ] > 0 ){
					tmp1	*= Gc;
					tmp2	*= Hc;
				}

				if( w->configs[i][ permute(j-3, w->size) ] > 0 ){
					tmp1	*= Gc;
					tmp2	*= Hc;
				}

				w->energies[i] += tmp1;
				sim->enthalpies[i] += tmp2;
			}
		}
	}
	return;
}