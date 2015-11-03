import numpy
import numpy.random

from os		import path
from uuid	import uuid4

from thermo	import *
from utilities	import savitzky_golay

_MATPLOTLIB_BACKEND = None #None for default

class ITCExperimentBase:
	"""A class that encapsulates a specific ITC experiment.

	Attributes:
		title (string): A descriptive title for the experiment.
		T (float): The experimental temperature (in Kelvin).
		V0 (float): The volume of the calorimeter's cell, in microliters.
		injections (list of floats): The amount (in microliters) injected at each injection point.
		dQ_exp (list of floats): The experimental (observed) enthalpies at each injection, in Joules.
		dQ_fit (list of floats): The predicted enthalpies (if any) at each injection, in Joules.
		Cell (dict of floats): The starting concentrations of components in the cell, in Molar units.
		Syringe (dict of floats): The concentrations of components in the syringe solution
		Concentrations (list of dicts): The concentrations of each component in the cell at each titration point, in Molar units.
		skip (list of ints): Titration points to exclude during fitting.
		ddQ (list of floats): List of estimated errors in each titration point enthalpy.
		spline_pts (int): Number of points to use in the Savitsky-Golay filter used to estimate errors in experimental enthalpies.
		spline_order (int): Order of the SG filter.
		Q_dil (float): Heat of dilution for syringe solution.
		chisq (float): If a fit has been generated, the reduced chi-squared goodness-of-fit value.
	"""

	def __init__(self, T, V0, injections, dQ, Cell, Syringe, skip=[], ddQ=[], Q_dil=0.0, cellRef=None, syringeRef=None, title=None, units='J'):
		"""Constructor for the ITCExperiment object.

		Args:
			T (float): The experimental temperature (in Kelvin).
			V0 (float): The volume of the calorimeter's cell, in microliters.
			injections (list of floats): The volumes injected at each titration point.
			Cell (dict of floats): The starting concentration of the named components in the cell.
			Syringe (dict of floats): The concentration of the named components in the syringe.
			skip (list of ints): Titration points to exclude during fitting.
			ddQ (list of floats): List of estimated errors in each titration point enthalpy.
			Q_dil (float): Heat of dilution for syringe solution.
			cellRef (string): The cell reference component.
			syringeRef (string): The syringe reference component.
			spline_pts (int): Number of points to use in the Savitsky-Golay filter used to estimate errors in experimental enthalpies.
			spline_order (int): Order of the SG filter.
			title (string): A descriptive title for the experiment.
		"""
		
		self.T			= T
		self.V0			= V0
		self.injections	= injections
		self.npoints	= len(self.injections)
		self.Cell		= Cell
		self.Syringe	= Syringe
		self.skip		= skip
		self.Q_dil		= Q_dil
		self.units		= units

		if title:
			self.title		= title
		else:
			self.title 		= uuid4()
	
		# reference components
		if cellRef == None:
			self.cellRef	= self.Cell.keys()[0]
		else:
			assert cellRef in Cell.keys()
			self.cellRef	= cellRef
		if syringeRef == None:
			self.syringeRef	= self.Syringe.keys()[0]
		else:
			assert syringeRef in Syringe.keys()
			self.syringeRef	= syringeRef

		# initialize the list of dicts used to track the concentrations of the components (in the cell) at each titration point
		self.Concentrations = []
		for i in xrange(self.npoints):
			self.Concentrations.append({})
			for s in Cell:
				self.Concentrations[-1][s] = 0.0
			for s in Syringe:
				self.Concentrations[-1][s] = 0.0

		# calculate the total ligand and macromolecule concentration in the cell at each titration point
		# this uses the dilution formula described in Microcal's data processing in Origin manual
		for i in xrange(self.npoints):
			dV = sum(self.injections[0:i])
			
			# these are += because it's possible that a component could be in both the syringe and cell solutions
			for s in Syringe:
				self.Concentrations[i][s] += (Syringe[s]*self.injections[i]/V0) + ((Syringe[s]*dV/V0) * (1.0/(1.0+(dV/(2.0*V0)))))
			for s in Cell:
				self.Concentrations[i][s] += Cell[s] * ( (1-(dV/(2.0*V0))) / (1.0+(dV/(2.0*V0))) )

		# convert raw data (in calories) to joules per mol of injectant
		assert len(dQ) == self.npoints
		self.dQ_exp = numpy.array([J_from_cal(dQ[i]) for i in xrange(self.npoints)],dtype='d')
		self.dQ_fit	= None
		self.ddQ = [None]*self.npoints
		self.chisq	= None

	def make_plot(self,hardcopy=False,hardcopydir='.',hardcopyprefix='', hardcopytype='png'):
		"""Generate a plot of the experimental data, and the fit if present.

		Arguments:
			hardcopy (boolean): Display the fit to the screen, or write it to a file?
			hardcopydir (string): The directory to write the hardcopy to.
			hardcopyprefix (string): A prefix to append to the experiment's title when generating the output plot filename.
			hardcopytype (string): The output format for the plot (availability is dependent upon what backends matplotlib was compiled with).

		Returns:
			None

		"""
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
		pyplot.ylabel("%s/mol of %s"%(self.units,self.syringeRef))
		pyplot.xlabel("%s / %s"%(self.syringeRef,self.cellRef))

		tmpx = [ self.Concentrations[i][self.syringeRef]/self.Concentrations[i][self.cellRef] for i in xrange(self.npoints) if i not in self.skip ]
		tmpy = [ convert_from_J(self.units,self.dQ_exp[i]) for i in xrange(self.npoints) if i not in self.skip ]
		tmpd = [ convert_from_J(self.units,self.ddQ[i]) for i in xrange(self.npoints) if i not in self.skip ]

		pyplot.errorbar(tmpx,tmpy,yerr=tmpd,c='#000000',fmt='s')

		if getattr(self,'spline',False):
			tmpz = [ convert_from_J(self.units,self.spline[i]) for i in xrange(self.npoints) if i not in self.skip ]
			pyplot.plot(tmpx,tmpz,c='g')

		if getattr(self,'dQ_fit',False):
			tmpy = [ convert_from_J(self.units,self.dQ_fit[i]) for i in xrange(self.npoints) if i not in self.skip ]
			pyplot.plot(tmpx,tmpy,c='r')

		if len(self.skip) > 0:
			tmpx = [ self.Concentrations[i][self.syringeRef]/self.Concentrations[i][self.cellRef] for i in self.skip ]
			tmpy = [ convert_from_J(self.units,self.dQ_exp[i]) for i in self.skip ]
			pyplot.errorbar(tmpx,tmpy,yerr=0,c='g',fmt='s')

		pyplot.draw()
		if hardcopy:
			fig.savefig( path.join(hardcopydir,"%s%s.%s"%(hardcopyprefix,self.title,hardcopytype)), bbox_inches='tight')
			pyplot.close(fig)
		else:
			pyplot.show()
			
	def export_data(self,path,units=None):
		"""Export the attributes of the experiment to a file.

		Args:
			path (string): The filesystem path to write the output to.
			units (string): The units to use in the export, must be either "J" for Joules (the default), or "cal" for calories.

		Returns:
			None
		"""
		
		if units:
			assert units in ('cal','kcal','J')
		else:
			units = self.units

		from datetime import datetime

		h = open(path,'w')
		h.write("# itcsim export of %s\n#\n"%(self.title))
		h.write("# Date	%s\n"%(datetime.ctime(datetime.today())))
		h.write("# units %s\n"%(units))
		h.write("# T	%.5f\n"%(self.T))
		h.write("# V0	%.5f\n"%(self.V0))
		for s in self.Cell:
			h.write("# Cell	%s	%.5f\n"%(s,self.Cell[s])) 
		for s in self.Syringe:
			h.write("# Syringe	%s	%.5f\n"%(s,self.Syringe[s])) 
		h.write("# Q_dil	%.5E\n"%(self.Q_dil))
		h.write("# skip: %s\n"%(",".join(map(str,self.skip))))
		h.write("#\n# Ivol	dQ_exp	ddQ_exp dQ_fit	spline	skipped	%s	%s"%("\t".join(self.Cell.keys()),"\t".join(self.Syringe.keys())))
		
		if self.spline == None:
			spline = [0 for i in xrange(self.npoints)]
		else:
			spline = self.spline[:]

		for i in xrange(self.npoints):
			cell	= ["%.5f"%(self.Concentrations[i][s]) for s in self.Cell]
			syringe	= ["%.5f"%(self.Concentrations[i][s]) for s in self.Syringe]
			h.write("%i	%.5f	%.5E	%.5E	%.5E	%.5f	%.5E	%.5E	%.5E	%.5E	%s	%s	%s\n"%(
				i,
				self.injections[i],
				self.Concentrations[i][self.cellRef],
				self.Concentrations[i][self.syringeRef],
				convert_from_J(units,self.dQ_exp[i]),
				convert_from_J(units,self.ddQ[i]),
				convert_from_J(units,self.dQ_fit[i]),
				convert_from_J(units,spline[i]),
				(i in self.skip),
				"\t".join(cell),
				"\t".join(syringe)
			))

		h.close()
		
	def get_chisq(self, Q, writeback=False):
		"""Calculate the goodness-of-fit between the provided data and the experimental data.

		Args:
			Q (list of floats): The predicted total heat at each injection point.

		Returns:
			(float): The goodness of the fit, as a reduced chi-square.
		"""

		dV = 0.0
		for i in xrange(self.npoints):
			dV += self.injections[i]
			Q[i] += self.Q_dil*self.injections[i]*(1.0 - 1.0/(1.0+self.V0/dV)) # add heat of dilution
			Q[i] *= self.V0 * self.Concentrations[i][self.cellRef] # normalize total heat content to macromolecule concentration in the cell volume
		
		# obtain the change in cell heat b/t each titration point
		dQ = [0.0]*self.npoints
		for i in xrange(self.npoints):
			if i==0:
				dQ[i] = Q[i] + ( (self.injections[i]/self.V0)*(Q[i]/2.0) )
			else:
				dQ[i] = Q[i] + ( (self.injections[i]/self.V0)*((Q[i]+Q[i-1])/2.0) ) - Q[i-1]
				
		if writeback:
			self.dQ_fit = dQ[:]

		# calculate reduced chi-square value using error estimate
		self.chisq = sum([(self.dQ_exp[i] - dQ[i])**2 / self.ddQ[i]**2 for i in xrange(self.npoints) if i not in self.skip])
		self.chisq /= (self.npoints -len(self.skip)) # reduced chi-square corrected for n degrees of freedom
		
		return self.chisq
		
class ITCExperiment(ITCExperimentBase):
	"""
	spline (list of floats): The smoothed spline (if any) to the experimental enthalpies, in Joules.
	"""
	
	def __init__(self, spline_pts=7, spline_order=1, *args, **kwargs ):
		ITCExperimentBase.__init__(self,*args,**kwargs)
	
		if ddQ != []:
			assert len(ddQ) == self.npoints
			self.ddQ = numpy.array(ddQ[:],dtype='d')
			self.spline = None
		else: # estimate errors by fitting spline
			counter,tmp = 0,[0]*self.npoints
			for i in xrange(self.npoints):
				if i not in self.skip:
					tmp[i] = counter
					counter += 1

			spl = savitzky_golay([self.dQ_exp[i] for i in xrange(self.npoints) if i not in self.skip], spline_pts, spline_order )
			err = numpy.std([self.dQ_exp[i] - spl[tmp[i]] for i in xrange(self.npoints) if i not in self.skip])
			self.ddQ = [err]*self.npoints
			self.spline = [spl[tmp[i]] for i in xrange(self.npoints)]

class ITCExperimentSynthetic(ITCExperimentBase):
	def __init__(self, injections, noise=0.0, *args, **kwargs):
		ITCExperimentBase.__init__(self,injections=injections,dQ=[0.0]*len(injections),*args,**kwargs)
		
		self.noise = noise
		self.initialized = False
		
	# monkeypatch chisq to updated dQ if uninitialized
	def get_chisq(self, *args, **kwargs):
	
		if self.noise == 0 and self.initialized: # no noise is nonsensical to fit
			return 0.0
			
		self.ddQ = [1.0]*self.npoints	
		ret = ITCExperimentBase.get_chisq(self, *args, **kwargs)

		if not self.initialized:
			if self.noise == 0:
				self.dQ_exp = self.dQ_fit[:]
			else:
				self.dQ_exp = [numpy.random.normal(self.dQ_fit[i],self.noise) for i in xrange(self.npoints)]
				
			self.dQ_fit, self.initialized = None, True
			ret = ITCExperimentBase.get_chisq(self, *args, **kwargs)
			
		if self.noise == 0:
			self.ddQ = [0.0]*self.npoints
			
		return ret






