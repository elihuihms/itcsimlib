#!/usr/bin/env python

#
# This script provides the number of sites that are empty, or have 0, 1, or 2 neighbors in a circular lattice.
#

from itcsimlib import ITCSim
from itcsimlib.thermo import *
from itcsimlib.model_ising import NonAdditive

class NonAdditiveRecorder(NonAdditive):
	def __init__(self,log=None,nsites=3,circular=1,**kwargs):
		NonAdditive.__init__(self,nsites,circular,**kwargs)
		
		self.log = log
		self.log_counter = 1
		self.neighbors = None

	def get_neighbors(self):
		neighbors = [ [0,0,0,0] for i in range(self.nconfigs) ]
		
		for i in range(self.nconfigs):			
			for j in range(self.nsites):
				
				if self.get_site_occupancy(i,j): # is site occupied?
					if self.get_site_occupancy(i,j+1): # is the next neighboring site occupied?
						if self.get_site_occupancy(i,j-1): # is previous neighboring site occupied?
							neighbors[i][3]+=1
						else:
							neighbors[i][2]+=1
							
					elif self.get_site_occupancy(i,j-1):
						neighbors[i][2]+=1
					else:
						neighbors[i][1]+=1
				else:
					neighbors[i][0]+=1
			
			# normalize
			for j in range(4):
				neighbors[i][j] *= self.weights[i]
				
		return neighbors
		
	def set_probabilities(self,totalP,totalL,T):
	
		# set the weights of each configuration
		NonAdditive.set_probabilities(self,totalP,totalL,T)
		if self.log == None:
			return
		
		# get the configuration-weighted number of sites either unoccupied, with no neighbors, one, or two neighboring sites occupied
		# must be called after set_probabilities()
		neighbors = self.get_neighbors()

		probabilities = [0.0,0.0,0.0,0.0]
		for i in range(self.nconfigs):
			for j in range(4):
				probabilities[j] += (neighbors[i][j] / self.nsites)

		# write probability information to logfile			
		handle = open(self.log,'a+')
		if self.log_counter == 1:
			handle.write("#inj	ratio	0	1	2	3\n")
		handle.write( "%i\t"%(self.log_counter) )
		handle.write( "%f\t"%(totalL / totalP) )
		handle.write( "%f\t"%(probabilities[0]) )
		handle.write( "%f\t"%(probabilities[1]) )
		handle.write( "%f\t"%(probabilities[2]) )
		handle.write( "%f\t"%(probabilities[3]) )
		handle.write( "\n" )
		handle.close()
		
		self.log_counter+=1

for name in ['20C-01-TRAPstk-TrpA','35C-01-TRAPstk-TrpA','40C-01-TRAPstk-TrpA','45C-01-TRAPstk-TrpA','55C-01-TRAPstk-TrpB','65C-01-TRAPstk-TrpB']:
	sim = ITCSim(T0=273.15+40, units='kcal', verbose=True, threads=1) # threads=1, as if multiple experiments are being simulated simultaneously, results from both will be written to the files out of order
	sim.set_model( NonAdditiveRecorder(log="log_%s.txt"%(name), nsites=11, circular=1, lattice_name="TRAP", ligand_name="Trp") )
	sim.set_model_params(dGX=-7.70, dGY=-8.74, dGZ=-8.85, dHX=-13.5, dHY=-19.7, dHZ=-18.2, dCpX=0.10, dCpY=-0.39, dCpZ=-0.35)
	sim.add_experiment_file( "./data/%s.txt"%(name) )
	sim.run()
	sim.make_plots(hardcopy=True)
	sim.done()
