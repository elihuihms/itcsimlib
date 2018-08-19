#!/usr/bin/env python

#
# This script uses an Ising model to dump the best-fit populations/abundancies of different TRAP+Tryptophan stoichiometries at each experimental temperature.
# Note that a similar effect can be obtained by generating a mass spec experiment and simulating the populations
#

from itcsimlib import ITCSim,ITCFit
from itcsimlib.model_ising import Ising
from itcsimlib.thermo	import *

# monkeypatch Ising class
def Q(self,T0,T,concentrations):
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
		for j in xrange(self.nsites+1):
			f.write("%f\t"%(sum([self.weights[k] for k in xrange(self.nconfigs) if self.bound[k]==j])))
		f.write("\n")
		
		# enthalpy is sum of all weighted enthalpies of the lattices
		Q[i] = sum( [self.weights[j] * self.enthalpies[j] for j in xrange(self.nconfigs)] )

	f.close()

	return Q
		
Ising.Q = Q

from itcsimlib.model_ising import NonAdditive

if __name__ == "__main__": # multiprocessing execution guards (only necessary for Windows, and if sim threads > 1)
	sim = ITCSim(T0=273.15+40,verbose=True,threads=6)
	sim.set_model( NonAdditive(nsites=11) )

	dG0,dGoe,dGoo,dH0,dHoe,dHoo,dCp0,dCpoe,dCpoo = -3.14694E+04,-4.98796E+03,-5.51985E+03,-4.89957E+04,-3.23101E+04,-2.58665E+04,4.13250E+02,-1.98629E+03,-1.83856E+03
	sim.set_model_params(
		dG0,dG0+dGoe,dG0+dGoo,
		dH0,dH0+dHoe,dH0+dHoo,
		dCp0,dCp0+dCpoe,dCp0+dCpoo)

	files =(
		( 'data/20C-01-TRAPstk-TrpA.txt', [0],	3.17924838 ),
		( 'data/35C-01-TRAPstk-TrpA.txt', [0],	-8.49621114 ),
		( 'data/40C-01-TRAPstk-TrpA.txt', [0,1],-5.60920276 ),
		( 'data/45C-01-TRAPstk-TrpA.txt', [0],	-6.31553711 ),
		( 'data/55C-01-TRAPstk-TrpB.txt', [0],	-1.92624856 ),
		( 'data/65C-01-TRAPstk-TrpB.txt', [0,85,91],	-5.54462468 ))

	for f,skip,Q_dil in files:
		sim.remove_all_experiments()

		sim.add_experiment_file( f, skip=skip, Q_dil=Q_dil )
		for e in sim.experiments:
			e.change_component_name('TRAP','Lattice')
			e.change_component_name('Trp','Ligand')

		print sim.run()
		sim.make_plots(hardcopy=True,hardcopytype='pdf')

	sim.done()
