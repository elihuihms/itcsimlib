/*
	Combination of SK and ian's models of Trp-TRAP binding
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

void assignEnergies( struct mWorkspace *w, struct sWorkspace *sim,
					double Ga, double Gb, double Gc, double Gd,
					double Ha, double Hb, double Hc, double Hd )
{
	double a1, a2;

	a1 = Gb/Ga; /* equivalent to SK's alpha1 */
	a2 = Gc/Ga; /* equivalent to SK's alpha2 */

	/* assign energies for each configuration */
	for(int i=0; i<pow(2, w->size); i++)
	{
		w->energies[i]=0;
		sim->enthalpies[i]=0;

		for(int j=0; j<w->size; j++)
		{
			if( w->configs[i][j] > 0 )
			{
				if( (w->configs[i][ permute(j+1,w->size) ]>0) && (w->configs[i][ permute(j-1,w->size) ]>0) ) /* both neighboring sites are occupied */
				{
					w->energies[i]		+= Gd;
					sim->enthalpies[i]	+= Hd;
				}
				else
				{
					if( (w->configs[i][ permute(j+1,w->size) ]>0) || (w->configs[i][ permute(j-1,w->size) ]>0) ){
						w->energies[i]		+= Gb; /* coupling to an occupied site */
						sim->enthalpies[i]	+= Hb;
					}else{
						w->energies[i]		+= Ga; /* coupling to an unoccupied site */
						sim->enthalpies[i]	+= Ha;
					}
				}
			}
		}
	}

	return;
}

