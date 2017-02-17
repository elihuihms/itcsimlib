#!/usr/bin/env python

#
# This script uses the critical chi-square approach to generate upper and lower confidence boundaries for the 1998 Saroff and Kiefer TRAP+Tryptophan binding model parameters described by .
#

from itcsimlib.itc_experiment import ITCExperimentBase

old_init = ITCExperimentBase.__init__ 
def new_init(self, T, V0, injections, dQ, Cell, Syringe, skip=[], ddQ=[], Q_dil=0.0, cellRef=None, syringeRef=None, title=None, units='J'):
	old_init(self, T, V0, injections, dQ, Cell, Syringe, skip, ddQ, Q_dil, cellRef, syringeRef, title, units)
	
	# overwrite heat of dilutions with old-style itcsimlib values
	for i in xrange(self.npoints):
		if i==0:
			self.dQ_dil[i] = self.Concentrations[i][self.syringeRef]
		else:
			self.dQ_dil[i] = self.Concentrations[i][self.syringeRef] -self.Concentrations[i-1][self.syringeRef]
		self.dQ_dil[i] *= -1.0*self.V0*self.Q_dil*self.Syringe[self.syringeRef]*self.injections[i]

ITCExperimentBase.__init__ = new_init

from itcsimlib import ITCSim,ITCFit
from itcsimlib.model_trap import SKa
from itcsimlib.utilities import *

sim = ITCSim(T0=273.15+40,verbose=True,threads=6)
sim.set_model( SKa() )
sim.set_model_params(
	dG0=-32157.324179654042,
	dGb=-1915.3395538189311,
	dH0=-88297.98259332242,
	dHb=6385.4580531436668,
	dCp0=-1663.5157839533472,
	dCpb=96.05928066671521
	)
	
sim.add_experiment_file( 'data/20C-01-TRAPstk-TrpA.txt', skip=[0],	Q_dil=8.524E+4  )
sim.add_experiment_file( 'data/35C-01-TRAPstk-TrpA.txt', skip=[0],	Q_dil=2.174E+5 )
sim.add_experiment_file( 'data/40C-01-TRAPstk-TrpA.txt', skip=[0,1],Q_dil=6.521E+4 )
sim.add_experiment_file( 'data/45C-01-TRAPstk-TrpA.txt', skip=[0],	Q_dil=3.057E+5 )
sim.add_experiment_file( 'data/55C-01-TRAPstk-TrpB.txt', skip=[0],	Q_dil=1.237E+5 )
sim.add_experiment_file( 'data/65C-01-TRAPstk-TrpB.txt', skip=[0,85,91],	Q_dil=3.673E+5 )

print sim.run()
sim.make_plots(hardcopy=True,hardcopyprefix="pre_",hardcopytype='pdf')

fit = ITCFit( sim, verbose=True )
opt,chisq = fit.optimize( params=sim.get_model_params() )
sim.set_model_params( **opt )

sim.make_plots(hardcopy=True,hardcopyprefix="post_",hardcopytype='pdf')

write_params_to_file( 'fit.txt', opt )

for p in sim.get_model_params():
	ret = fit.estimate_sigma( params=[p], method='bisect', stdevs=1 )
	write_data_to_file( 'est_%s.log'%p, str(ret) )

sim.done()