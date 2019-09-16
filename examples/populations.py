#!/usr/bin/env python

#
# This script uses an Ising model to dump the best-fit populations/abundancies of different TRAP+Tryptophan stoichiometries at each experimental temperature.
# Note that a similar effect can be obtained by generating a mass spec experiment and simulating the populations
#

from itcsimlib import ITCSim,ITCFit
from itcsimlib.model_ising import Ising
from itcsimlib.thermo	import *

# monkeypatch Ising class
def populations_Q(self,T0,T,concentrations):
	"""Return the enthalpy of the system at each of the specified concentrations """
	
	# set the free energies (and enthalpic energies if necessary) of each configuration
	self.set_energies(T0,T)
	self.precision=1E-12
	
	f = open("populations_nonadditive_%i.txt"%T,'w')
	
	# calculate the enthalpy at each set of conditions
	Q = [0.0]*len(concentrations)
	for i,c in enumerate(concentrations):
		# set the weights (probabilities) of each lattice configuration
		self.set_probabilities(c['Lattice'],c['Ligand'],T)
		
		f.write("%i\t%f\t"%(i,c['Ligand']/c['Lattice']))
		for j in range(self.nsites+1):
			f.write("%f\t"%(sum([self.weights[k] for k in range(self.nconfigs) if self.bound[k]==j])))
		f.write("\n")
		
		# enthalpy is sum of all weighted enthalpies of the lattices
		Q[i] = sum( [self.weights[j] * self.enthalpies[j] for j in range(self.nconfigs)] )

	f.close()

	return Q
		
Ising.Q = populations_Q

# import NonAdditive model, which will inheirt our monkeypatched Ising model and modified Q function.
from itcsimlib.model_ising import NonAdditive

if __name__ == "__main__": # multiprocessing execution guards (only necessary for Windows, and if sim threads > 1)
	sim = ITCSim(T0=273.15+40, units='kcal', verbose=True, threads=6)
	sim.set_model( NonAdditive(nsites=11) )
	sim.set_model_params(dGX=-7.70, dGY=-8.74, dGZ=-8.85, dHX=-13.5, dHY=-19.7, dHZ=-18.2, dCpX=0.10, dCpY=-0.39, dCpZ=-0.35)

	for f in ['data/20C-01-TRAPstk-TrpA.txt', 'data/35C-01-TRAPstk-TrpA.txt', 'data/40C-01-TRAPstk-TrpA.txt', 'data/45C-01-TRAPstk-TrpA.txt', 'data/55C-01-TRAPstk-TrpB.txt', 'data/65C-01-TRAPstk-TrpB.txt']:
		sim.remove_all_experiments()
		sim.add_experiment_file(f)
		for e in sim.experiments:
			e.change_component_name('TRAP','Lattice')
			e.change_component_name('Trp','Ligand')

		sim.run()
		sim.make_plots(hardcopy=True,hardcopytype='pdf')

	sim.done()
