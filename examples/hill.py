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
		fh.write("#	TotalLattice	TotalLigand 	FreeLigand	Saturation\n")
		Q = [0.0]*len(concentrations)
		for i,c in enumerate(concentrations):
			# set the weights (probabilities) of each lattice configuration
			freeL = self.set_probabilities(c['Lattice'],c['Ligand'],T)
			
			# get the partial saturation
			theta = 0.0
			for j in range(self.nconfigs):
				theta += self.bound[j] * self.weights[j]
				
			fh.write("%i	%0.5E	%0.5E	%0.5E	%f\n"%(i,c['Lattice'],c['Ligand'],freeL,theta))
			
			# enthalpy is sum of all weighted enthalpies of the lattices
			Q[i] = sum( [self.weights[j] * self.enthalpies[j] for j in range(self.nconfigs)] )

		fh.close()

		return Q

model = Hill_Printer_Model(nsites=11, circular=1)
sim = ITCSim(T0=273.15+40, units='kcal', verbose=True, threads=1) # threads=1, as if multiple experiments are being simulated simultaneously, results from both will be written to the files out of order

injection_vols = [0.1]*100
injection_vols.extend( [1.0]*100 )

sim.add_experiment_synthetic(
	T=273.15+40,
	V0=1416.6,
	injections=injection_vols,
	Cell={"Lattice":1E-6},
	Syringe={"Ligand":200E-6},
	title='simulated'
)

sim.set_model( model )
sim.set_model_params(dGX=-7.70, dGY=-8.74, dGZ=-8.85, dHX=-13.5, dHY=-19.7, dHZ=-18.2, dCpX=0.10, dCpY=-0.39, dCpZ=-0.35)
sim.run()
sim.make_plots(hardcopy=True,hardcopyprefix="hill_")

sim.done()

