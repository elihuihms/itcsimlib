/*
	Combination of SK and ian's models of Trp-TRAP binding
*/

#include "itc_model.h"
#include "itc_sim.h"

void assignEnergies( struct mWorkspace *w, struct sWorkspace *sim,
					double Gbind, double Ga, double Gb, double Gc,
					double Hnone, double Ha, double Hb, double Hc )
{
	double tmp1, tmp2;

	/* assign energies for each configuration */
	for(int i=0; i<pow(2, w->size); i++)
	{
		w->energies[i]=0;
		sim->enthalpies[i]=0;

		for(int j=1; j<w->size -1; j++)
		{
			if( w->configs[i][j] > 0 )
			{
				tmp1 = Gbind; /* presence of ligand dG */
				tmp2 = Hnone; /* binding to site enthalpy */

				if( (w->configs[i][j+1]>0) && (w->configs[i][j-1]>0) ) /* both neighboring sites are occupied */
				{
					tmp1 *= Gc;
					tmp2 *= Hc;
				}
				else
				{
					if( w->configs[i][j+1] > 0 ){
						tmp1 *= Gb; /* coupling to an occupied site */
						tmp2 *= Hb;
					}else{
						tmp1 *= Ga; /* coupling to an unoccupied site */
						tmp2 *= Ha;
					}

					if( w->configs[i][j-1] > 0 ){
						tmp1 *= Gb;
						tmp2 *= Hb;
					}else{
						tmp1 *= Ga;
						tmp2 *= Ha;
					}
				}

				w->energies[i] += tmp1;
				sim->enthalpies[i] += tmp2;
			}
		}

		if( w->cyclic > 0 )
		{
			if( w->configs[i][0] > 0 )
			{
				tmp1 = Gbind;
				tmp2 = Hnone;

				if( (w->configs[i][1]>0) && (w->configs[i][w->size -1]>0) )
				{
					tmp1 *= Gc;
					tmp2 *= Hc;
				}
				else
				{
					if( w->configs[i][1] > 0 ){
						tmp1 *= Gb;
						tmp2 *= Hb;
					}else{
						tmp1 *= Ga;
						tmp2 *= Ha;
					}

					if( w->configs[i][w->size -1] > 0 ){
						tmp1 *= Gb;
						tmp2 *= Hb;
					}else{
						tmp1 *= Ga;
						tmp2 *= Ha;
					}
				}

				w->energies[i] += tmp1;
				sim->enthalpies[i] += tmp2;
			}
			if( w->configs[i][w->size -1] > 0 )
			{
				tmp1 = Gbind;
				tmp2 = Hnone;

				if( (w->configs[i][0]) && (w->configs[i][w->size -2]>0) )
				{
					tmp1 *= Gc;
					tmp2 *= Hc;
				}
				else
				{
					if( w->configs[i][0] ){
						tmp1 *= Gb;
						tmp2 *= Hb;
					}else{
						tmp1 *= Ga;
						tmp2 *= Ha;
					}

					if( w->configs[i][w->size -2] > 0 ){
						tmp1 *= Gb;
						tmp2 *= Hb;
					}else{
						tmp1 *= Ga;
						tmp2 *= Ha;
					}
				}

				w->energies[i] += tmp1;
				sim->enthalpies[i] += tmp2;
			}
		}
	}

	return;
}

