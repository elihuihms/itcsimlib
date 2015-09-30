/*
	configuration energies/enthalpies as described in Elihu's thesis
*/

#include "itc_model.h"
#include "itc_sim.h"
#include "energies.h"

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
	 1				(K)
	010				(a)
	011 or 110		(b)
	111				(c)

	param[0]		dG0
	param[1]		dGa
	param[2]		dGb
	param[3]		dGc
	param[4]		dH0
	param[5]		dHa
	param[6]		dHb
	param[7]		dHc
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
				w->energies[i]			+= params[0]; /*  1  */
				sim->enthalpies[i]		+= params[4];

				if( (w->configs[i][ permute(j-1,w->size) ] == 0) && (w->configs[i][ permute(j+1,w->size) ] == 0) ){
					w->energies[i]		+= params[1]; /* 010 */
					sim->enthalpies[i]	+= params[5];
				}else if( (w->configs[i][ permute(j-1,w->size) ] == 1) && (w->configs[i][ permute(j+1,w->size) ] == 0) ){
					w->energies[i]		+= params[2]; /* 110 */
					sim->enthalpies[i]	+= params[6];
				}else if( (w->configs[i][ permute(j-1,w->size) ] == 0) && (w->configs[i][ permute(j+1,w->size) ] == 1) ){
					w->energies[i]		+= params[2]; /* 011 */
					sim->enthalpies[i]	+= params[6];
				}else if( (w->configs[i][ permute(j-1,w->size) ] == 1) && (w->configs[i][ permute(j+1,w->size) ] == 1) ){
					w->energies[i]		+= params[3]; /* 111 */
					sim->enthalpies[i]	+= params[7];
				}
			}
		}
	}

	return;
}
