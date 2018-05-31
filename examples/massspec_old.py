#!/usr/bin/env python

#
# This script demonstrates how an itcsimlib experiment can be extended in order to fit mass spectrometric data
#

import os

from itcsimlib import *
from itcsimlib.model_ising import *
from numpy import *

_MATPLOTLIB_BACKEND = None #None for default

class MSExperiment():
	def __init__(self, file):
		self.title 			= os.path.splitext(os.path.basename(file))[0]
		self.chisq			= None

		fh = open(file)
		self.T = float(fh.readline().split()[2])
		PConc = map(float,fh.readline().split()[2:])
		LConc = map(float,fh.readline().split()[2:])
		fh.close()
		
		self.npoints = len(PConc)
		assert len(LConc) == self.npoints
		
		data = genfromtxt( file, usecols=xrange(1,self.npoints+1) )

		assert len(data)/2 == self.npoints # first half is measured values, second half is sigmas
		self.npops = len(data[0])

		self.Concentrations,self.PopIntens,self.PopErrors,self.PopFits = [],[],[],[]
		for i in xrange(self.npoints):
			self.Concentrations.append({})
			self.Concentrations[i]['Lattice'] = PConc[i] / 1.0E6
			self.Concentrations[i]['Ligand'] = LConc[i] / 1.0E6
			self.PopIntens.append( data[i][:] )
			self.PopErrors.append( data[i+self.npoints][:] )
			self.PopFits.append( [0.0]*self.npops )
	
	def make_plot(self,hardcopy=False,hardcopydir='.',hardcopyprefix='', hardcopytype='png'):
		try:
			if _MATPLOTLIB_BACKEND != None:
				import matplotlib
				matplotlib.use(_MATPLOTLIB_BACKEND)
			import matplotlib.pyplot as pyplot
		except:
			pyplot = None

		if pyplot == None: return
		if hardcopy: fig = pyplot.figure()

		pyplot.clf()
		pyplot.title(self.title)
		pyplot.ylabel("Experimental - Fit Abundance")
		pyplot.xlabel("Mass Populations")
		
		ax1 = pyplot.subplot(3,1,1)
		ax2 = pyplot.subplot(3,1,2)
		ax3 = pyplot.subplot(3,1,3)

		xax_positions,xax_labels = [],[]
		width,space,left = 0.25,0.5,0.0
		for i in xrange(self.npops):
			bars1 = ax1.bar( [left + (j*width) for j in xrange(self.npoints)], [self.PopIntens[j][i] for j in xrange(self.npoints)], width=width, edgecolor='r' )
			bars2 = ax2.bar( [left + (j*width) for j in xrange(self.npoints)], [self.PopFits[j][i] for j in xrange(self.npoints)], width=width, edgecolor='r' )
			bars3 = ax3.bar( [left + (j*width) for j in xrange(self.npoints)], [self.PopIntens[j][i] - self.PopFits[j][i] for j in xrange(self.npoints)], width=width, edgecolor='r' )
			xax_positions.append( left + (self.npoints*width)/2.0 )
			xax_labels.append("%i"%i)
			left += (j*width)+space
	
			for j in xrange(self.npoints):
				#color = float(j) / self.npoints	
				#bars1[j].set_color( (0, 0, color) )	
				#bars2[j].set_color( (0, color, 0) )	
				#bars3[j].set_color( (color, 0, 0) )
				color = 1.0 - (float(j) / self.npoints)
				bars1[j].set_color( (color, color, 1) )	
				bars2[j].set_color( (color, 1, color) )	
				bars3[j].set_color( (1, color, color) )

				bars1[j].set_edgecolor( (0,0,0) )	
				bars2[j].set_edgecolor( (0,0,0) )	
				bars3[j].set_edgecolor( (0,0,0) )	
		
		ax1.set_ylabel("Experimental")
		ax1.set_xticks( xax_positions )
		ax1.set_xticklabels( xax_labels)
		ax2.set_ylabel("Fit")
		ax2.set_xticks( xax_positions )
		ax2.set_xticklabels( xax_labels)
		ax3.set_ylabel("Residuals")
		ax3.set_xticks( xax_positions )
		ax3.set_xticklabels( xax_labels)
			
		pyplot.draw()

		if hardcopy:
			fig.savefig( os.path.join(hardcopydir,"%s%s.%s"%(hardcopyprefix,self.title,hardcopytype)), bbox_inches='tight')
			pyplot.close(fig)
		else:
			pyplot.show()
			
	def export_to_file(self,path):
		fh = open(path, 'w')
		fh.write("# Experimental data\n")
		for i in xrange(self.npoints):
			fh.write("%i	%s\n"%(i,"\t".join(["%f"%f for f in self.PopIntens[i]])))
		fh.write("# Fit data\n")
		for i in xrange(self.npoints):
			fh.write("%i	%s\n"%(i,"\t".join(["%f"%f for f in self.PopFits[i]])))
		fh.close()
			
	def get_chisq(self, data, writeback=False):
		self.chisq = 0.0
		for i in xrange(self.npoints):
			for j in xrange(self.npops):
				self.chisq += (data[i][j] - self.PopIntens[i][j])**2 / self.PopErrors[i][j]**2
		
		if writeback:
			self.PopFits = data
		
		self.chisq = self.chisq / (self.npoints * self.npops)

		return self.chisq
		
class NonAdditiveMS(Ising):
	def __init__(self,nsites=11,circular=1):
		Ising.__init__(self,nsites,circular)
		self.units = 'J'
		self.add_parameter( 'dGX',	'dG',	description='Free energy change upon binding to a site flanked by two unoccupied' )
		self.add_parameter( 'dGY',	'dG',	description='Free energy change upon binding to a site flanked by one occupied' )
		self.add_parameter( 'dGZ',	'dG',	description='Free energy change upon binding to a site flanked by two occupied' )

	def set_energies(self,T0,T):
		for i in xrange(self.nconfigs):
			self.gibbs[i] = 0.0
			
			for j in xrange(self.nsites):
				if self.get_site_occupancy(i,j): # is site occupied?
					
					if self.get_site_occupancy(i,j+1): # is the next neighboring site occupied?
						if self.get_site_occupancy(i,j-1): # is previous neighboring site occupied?
							self.gibbs[i]+= self.params['dGZ']
							
						elif self.circular:
							self.gibbs[i]+= self.params['dGY']
							
					elif self.get_site_occupancy(i,j-1):
						self.gibbs[i]+= self.params['dGY']
						
					elif self.circular:
						self.gibbs[i]+= self.params['dGX']

	def Q(self,T0,T,concentrations):
		pops = []
		print T
		self.set_energies(T0,T)
		for i,c in enumerate(concentrations):
			pops.append( [0.0]*(self.nsites+1) )
			
			self.set_probabilities(c['Lattice'],c['Ligand'],T)
			for j in xrange(self.nconfigs):
				pops[i][self.bound[j]] += self.weights[j]

		return pops

if __name__ == "__main__":

	model = NonAdditiveMS()
	model.set_params(-27000,-27000,-30000)
	
	sim = ITCSim(verbose=True)

	sim.set_model(model)

	sim.add_experiment( MSExperiment('data/TRAP_populations_EDDA_old.txt') )	
	
	sim.run()

	sim.make_plots(hardcopy=True,hardcopytype='png')

	sim.experiments[0].export_to_file("massspec_old.fit")

	sim.done()