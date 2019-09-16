#!/usr/bin/env python

#
# This script compares the critical bootstrap and chi-square approaches for obtaining confidence boundaries across all of the heat capacity terms.
# Note that this script will first attempt to use the compiled C model (for speed), but will fall back to the native python model if they were not built during installation.
#

from itcsimlib.itc_experiment import ITCExperimentBase

old_init = ITCExperimentBase.__init__ 
def new_init(self, T, V0, injections, dQ, Cell, Syringe, skip=[], ddQ=[], Q_dil=0.0, cellRef=None, syringeRef=None, title=None, units='J'):
	old_init(self, T, V0, injections, dQ, Cell, Syringe, skip, ddQ, Q_dil, cellRef, syringeRef, title, units)
	
	# overwrite heat of dilutions with old-style itcsimlib values
	for i in range(self.npoints):
		if i==0:
			self.dQ_dil[i] = self.Concentrations[i][self.syringeRef]
		else:
			self.dQ_dil[i] = self.Concentrations[i][self.syringeRef] -self.Concentrations[i-1][self.syringeRef]
		self.dQ_dil[i] *= -1.0*self.V0*self.Q_dil*self.Syringe[self.syringeRef]*self.injections[i]

ITCExperimentBase.__init__ = new_init

from itcsimlib import ITCSim,ITCFit
from itcsimlib.model_trap import IK
from itcsimlib.model_ising import NonAdditive
from itcsimlib.utilities import *
       
if __name__ == "__main__": # multiprocessing execution guards (only necessary for Windows, and if sim threads > 1)
	sim = ITCSim(T0=273.15+40, units='kcal', verbose=False, threads=6)
	
	try: # attempt to use the much-faster compiled model
		sim.set_model( IK() )
	except (ImportError, OSError):
		sim.set_model( NonAdditive(nsites=11, circular=1, lattice_name="TRAP", ligand_name="Trp") )

	sim.set_model_params(dGX=-7.70, dGY=-8.74, dGZ=-8.85, dHX=-13.5, dHY=-19.7, dHZ=-18.2, dCpX=0.10, dCpY=-0.39, dCpZ=-0.35)
	sim.add_experiment_file( 'data/20C-01-TRAPstk-TrpA.txt', skip=[0], Q_dil=4.17)
	sim.add_experiment_file( 'data/35C-01-TRAPstk-TrpA.txt', skip=[0], Q_dil=-6.10 )
	sim.add_experiment_file( 'data/40C-01-TRAPstk-TrpA.txt', skip=[0,1], Q_dil=-2.79 )
	sim.add_experiment_file( 'data/45C-01-TRAPstk-TrpA.txt', skip=[0], Q_dil=-3.15 )
	sim.add_experiment_file( 'data/55C-01-TRAPstk-TrpB.txt', skip=[0], Q_dil=0.37 )
	sim.add_experiment_file( 'data/65C-01-TRAPstk-TrpB.txt', skip=[0,85,91], Q_dil=-3.56 )

	fit = ITCFit( sim, verbose=True, method='powell' )

	""" Uncomment for additional optimization
	print("Initial chisquare: %.2f" % sim.run())
	sim.make_plots(hardcopy=True,hardcopyprefix="estimate_pre_",hardcopytype='pdf')	
	opt,chisq = fit.optimize( params=sim.get_model_params() )
	print("Optimized chisquare: %.2f" % chisq)
	sim.set_model_params( **opt )
	sim.run()
	sim.make_plots(hardcopy=True,hardcopyprefix="estimate_post_",hardcopytype='pdf')
	write_params_to_file( 'estimate_optimized.txt', opt )
	"""

	est_bootstrap = fit.estimate( params=['dCpX', 'dCpY', 'dCpZ'], method='bootstrap', bootstraps=100 )
	est_critchisq = fit.estimate( params=['dCpX', 'dCpY', 'dCpZ'], method='sigma', stdevs=1 )
	
	for param in ['dCpX', 'dCpY', 'dCpZ']:
		bounds_bootstrap = sorted( est_bootstrap[ param ] )
		bounds_critchisq = est_critchisq[ param ].sorted()
		write_data_to_file( 'estimate_optimized.txt', "Bootstrap %s = %.3f (%.3f to %.3f)\n"%(param, sim.get_model_params()[ param ], bounds_bootstrap[0], bounds_bootstrap[1]) )
		write_data_to_file( 'estimate_optimized.txt', "CriticalC %s = %.3f (%.3f to %.3f)\n"%(param, sim.get_model_params()[ param ], bounds_critchisq[0], bounds_critchisq[1]) )
	
	sim.done()