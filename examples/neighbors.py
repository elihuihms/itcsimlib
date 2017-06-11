#!/usr/bin/env python

#
# This script provides the number of sites that are empty, or have 0, 1, or 2 neighbors in a circular lattice.
#

from itcsimlib import ITCSim
from itcsimlib.thermo import *
from itcsimlib.model_ising import NonAdditive

class NonAdditiveRecorder(NonAdditive):
	def __init__(self,nsites=3,circular=1,log=None):
		NonAdditive.__init__(self,nsites,circular)
		
		self.log = log
		self.log_counter = 1
		self.neighbors = None

	def get_neighbors(self):
		neighbors = [ [0,0,0,0] for i in xrange(self.nconfigs) ]
		
		for i in xrange(self.nconfigs):			
			for j in xrange(self.nsites):
				
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
			for j in xrange(4):
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
		for i in xrange(self.nconfigs):
			for j in xrange(4):
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

names = ['20C-01-TRAPstk-TrpA','35C-01-TRAPstk-TrpA','40C-01-TRAPstk-TrpA','45C-01-TRAPstk-TrpA','55C-01-TRAPstk-TrpB','65C-01-TRAPstk-TrpB']
skip = [ [0],[0],[0,1],[0],[0],[0,85,91] ]

for i in xrange(1):

	sim = ITCSim(T0=273.15+40,verbose=True,threads=1)
	sim.set_model( NonAdditiveRecorder(nsites=11,circular=1,log="log_%s.txt"%(names[i])) )
	sim.set_model_params(-31237.753,-36510.231,-36983.080,-47678.518,-81769.489,-75924.637,405.754,-1606.945,-1494.137)

	sim.add_experiment_file( "./data/%s.txt"%(names[i]),	skip=skip[i] )
	sim.get_experiment(0).change_component_name('TRAP','Lattice')
	sim.get_experiment(0).change_component_name('Trp','Ligand')

	sim.run()
	sim.make_plots(hardcopy=True)
	sim.done()
