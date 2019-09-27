#
# This script optimizes the heats of dilution for each experiment given an existing global fit.
# Note that this script will first attempt to use the compiled C model (for speed), but will fall back to the native python model if they were not built during installation.
#

from itcsimlib import ITCSim,ITCFit
from itcsimlib.model_trap import IK
from itcsimlib.model_ising import NonAdditive
from itcsimlib.utilities import *

from scipy import optimize

if __name__ == "__main__": # multiprocessing execution guards (only necessary for Windows, and if sim threads > 1)

	sim = ITCSim(T0=273.15+40, units='kcal', verbose=False, threads=6)
	sim.add_experiment_file( 'data/20C-01-TRAPstk-TrpA.txt', skip=[0])
	sim.add_experiment_file( 'data/35C-01-TRAPstk-TrpA.txt', skip=[0])
	sim.add_experiment_file( 'data/40C-01-TRAPstk-TrpA.txt', skip=[0,1])
	sim.add_experiment_file( 'data/45C-01-TRAPstk-TrpA.txt', skip=[0])
	sim.add_experiment_file( 'data/55C-01-TRAPstk-TrpB.txt', skip=[0])
	sim.add_experiment_file( 'data/65C-01-TRAPstk-TrpB.txt', skip=[0,85,91])
	
	try: # attempt to use the much-faster compiled model
		sim.set_model( IK() )
	except (ImportError, OSError):
		sim.set_model( NonAdditive(nsites=11, circular=1, lattice_name="TRAP", ligand_name="Trp") )

	sim.set_model_params(dGX=-7.70, dGY=-8.74, dGZ=-8.85, dHX=-13.5, dHY=-19.7, dHZ=-18.2, dCpX=0.10, dCpY=-0.39, dCpZ=-0.35)	
	sim.run()
	sim.make_plots(hardcopy=True,hardcopyprefix="dilution_pre_",hardcopytype='pdf')
	
	def optimize_dilution( x ):
		sim.remove_all_experiments()
		sim.add_experiment_file( 'data/20C-01-TRAPstk-TrpA.txt', skip=[0],	Q_dil=x[0] )
		sim.add_experiment_file( 'data/35C-01-TRAPstk-TrpA.txt', skip=[0],	Q_dil=x[1] )
		sim.add_experiment_file( 'data/40C-01-TRAPstk-TrpA.txt', skip=[0,1],	Q_dil=x[2] )
		sim.add_experiment_file( 'data/45C-01-TRAPstk-TrpA.txt', skip=[0],	Q_dil=x[3] )
		sim.add_experiment_file( 'data/55C-01-TRAPstk-TrpB.txt', skip=[0],	Q_dil=x[4] )
		sim.add_experiment_file( 'data/65C-01-TRAPstk-TrpB.txt', skip=[0,85,91],	Q_dil=x[5] )			
		chisq = sim.run()

		print(x,chisq)
		return chisq

	print(optimize.fmin_powell(func=optimize_dilution, x0=[0.0]*6))

	sim.make_plots(hardcopy=True,hardcopyprefix="dilution_post_",hardcopytype='pdf')
	sim.done()