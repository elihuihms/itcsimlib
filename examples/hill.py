#!/usr/bin/env python

#
# This script writes the free ligand and lattice saturation (theta) to a file at each titration point.
#

from itcsimlib import ITCSim
from itcsimlib.model_ising import *

class Hill_Printer_Model(NonAdditive):

	def Q(self,T0,T,concentrations):
		
		# set the free energies (and enthalpic energies if necessary) of each configuration
		self.set_energies(T0,T)
		
		fh = open("hill_%i.txt"%T,'a')

		Q = [0.0]*len(concentrations)
		for i,c in enumerate(concentrations):
			# set the weights (probabilities) of each lattice configuration
			freeL = self.set_probabilities(c['Lattice'],c['Ligand'],T)
			
			# get the partial saturation
			theta = 0.0
			for j in xrange(self.nconfigs):
				theta += self.bound[j] * self.weights[j]
				
			fh.write("%i	%0.5E	%f\n"%(i,freeL,theta))
			
			# enthalpy is sum of all weighted enthalpies of the lattices
			Q[i] = sum( [self.weights[j] * self.enthalpies[j] for j in xrange(self.nconfigs)] )

		fh.close()

		return Q

model = Hill_Printer_Model(nsites=11,circular=1)
sim = ITCSim(T0=273.15+40,verbose=True,threads=1)

injection_vols = [0.1]*100
injection_vols.extend( [1.0]*200 )
injection_vols.extend( [10]*10 )

sim.add_experiment_synthetic(
	T=273.15+40,
	V0=1416.6,
	injections=injection_vols,
	Cell={"Lattice":1E-6},
	Syringe={"Ligand":100E-6},
	noise=1.0,
	title='simulated')

sim.set_model( model )
sim.set_model_params(-31237.753,-36510.231,-36983.080,-47678.518,-81769.489,-75924.637,405.754,-1606.945,-1494.137)
sim.run()
sim.make_plots(hardcopy=True,hardcopyprefix="hill_")

sim.done()

