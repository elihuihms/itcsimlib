import scipy
import uuid
import os

from thermo	import *
from thermo import _UNITS
from utilities	import savitzky_golay

class ITCExperimentBase:
	"""A class that encapsulates a specific ITC experiment.

	Attributes
	----------
	chisq : float
		If a fit has been generated, the reduced chi-squared goodness-of-fit value.
	T : float
		The experimental temperature (in Kelvin).
	Concentrations : list of dicts
		The concentrations of the components (in the cell) at each titration point
			
	Notes
	-----
		Energies are stored in Joules internally, and then converted (if necessary) to the desired unit when retrieved by the user.
	"""

	def __init__(self, T, V0, injections, dQ, Cell, Syringe, skip=[], dQ_err=[], Q_dil=0.0, cellRef=None, syringeRef=None, title=None, units='J'):
		"""Constructor for the ITCExperiment object.

		Arguments
		---------
		T : float
			The experimental temperature (in Kelvin).
		V0 : float
			The volume of the calorimeter's cell, in microliters.
		injections : list of floats
			The volumes injected at each titration point.
		dQ : list of floats 
			The evolved heat at each injection (in calories, unnormalized)
		Cell : dict of floats
			The starting concentration of the named components in the cell.
		Syringe : dict of floats
			The concentration of the named components in the syringe.
		skip : list of ints
			Titration points to exclude during fitting.
		dQ_err : list of floats
			List of estimated error in each titration point enthalpy.
		Q_dil : float
			Heat of dilution for syringe solution.
		cellRef : string
			The cell reference component.
		syringeRef : string
			The syringe reference component.
		title : string
			A descriptive title for the experiment.
		units : string
			The units of dQ heats.
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
		else: # assign a pseudounique title
			self.title 		= uuid.uuid4()
	
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
		self.dDQ_conc = [0.0]*self.npoints # fractional concentration of syringe solution in cell to use for dilution calculations
		for i in xrange(self.npoints):
			dV = sum(self.injections[0:i])
			
			self.dDQ_conc[i] = (1.0*self.injections[i]/V0) + ((1.0*dV/V0) * (1.0/(1.0+(dV/(2.0*V0)))))
			
			# these are += because it's possible that a component could be in both the syringe and cell solutions
			for s in self.Syringe:
				self.Concentrations[i][s] += (self.Syringe[s]*self.injections[i]/V0) + ((self.Syringe[s]*dV/V0) * (1.0/(1.0+(dV/(2.0*V0)))))
			for s in self.Cell:
				self.Concentrations[i][s] += self.Cell[s] * ( (1-(dV/(2.0*V0))) / (1.0+(dV/(2.0*V0))) )

		# heat of dilution will be proportional to the difference in concentration between syringe solution and cell solutions
		self.dQ_dil = [0.0]*self.npoints
		for i in xrange(self.npoints):
			if i == 0:
				self.dQ_dil[i] = (1.0 -self.dDQ_conc[i]) * self.Q_dil
			else:
				self.dQ_dil[i] = (1.0 -self.dDQ_conc[i] -self.dDQ_conc[i-1]) * self.Q_dil

		assert len(dQ) == self.npoints

		# convert raw data (in calories) to joules (note that this is not normalized per mol of injectant!)
		self.dQ_exp = scipy.array([J_from_cal(dQ[i]) for i in xrange(self.npoints)],dtype='d')
		self.dQ_fit	= None
		self.spline = None

		if dQ_err == []:
			self.dQ_err = None
		else:
			assert len(dQ_err) == self.npoints
			self.dQ_err = scipy.array([J_from_cal(dQ_err[i]) for i in xrange(self.npoints)],dtype='d')

		self.chisq	= None
		self.initialized = False # will be set to True by implementors of this base class
	
	def __str__(self):
		"""Stringify the experiment to be suitable for display to the user
		
		Parameters
		----------
		None
		
		Returns
		-------
		string
			String containing the description and other experimental info.
		"""
		
		ret = "Title: %s\n"%(self.title)
		if (self.chisq != None):
			ret+= "Chisq: %f\n"%(self.chisq)
		ret+= "Temperature: %fK\n"%(self.T)
		ret+= "%i Injections (%i skipped)\n"%(self.npoints,len(self.skip))
		ret+= "Dilution enthalpy: %.3E\n"%(self.Q_dil)
		ret+= "Cell components:\n"
		for s in self.Cell:
			ret+="\t%s (%.3E M)"%(s,self.Concentrations[0][s])
			if s == self.Cell.keys()[0]:
				ret+=" (reference)"
			ret+="\n"
		ret+= "Syringe components:\n"
		for s in self.Syringe:
			ret+="\t%s (%.3E M)"%(s,self.Concentrations[0][s])
			if s == self.Syringe.keys()[0]:
				ret+=" (reference)"
			ret+="\n"
		return ret

	def change_component_name(self,old_name,new_name):
		"""Update the name of a component in the experiment to the specified one.
		
		Arguments
		---------
		old_name : string
			The name of the component in the experiment object to update
		new_name : string
			The new name of the component to use
			
		Returns
		-------
		None
		
		Notes
		-----
			Because syringe and cell components are referred to by name, problems can arise where the component is called by one name in the model, and another in the experiment.
			This method allows the user to update an experimental component by name to match the one used in the model, e.g. from the specific "Tryptophan" specified in the experiment file to the generic "Ligand" used in the model.
		"""

		assert old_name in self.Concentrations[0]
		for i in xrange(self.npoints):
			self.Concentrations[i][new_name] = self.Concentrations[i][old_name]
			del self.Concentrations[i][old_name]
		
		if old_name in self.Cell:
			self.Cell[new_name] = self.Cell[old_name]
			del self.Cell[old_name]
			
		if old_name in self.Syringe:
			self.Syringe[new_name] = self.Syringe[old_name]
			del self.Syringe[old_name]
			
		if self.cellRef == old_name:
			self.cellRef = new_name
		if self.syringeRef == old_name:
			self.syringeRef = new_name
		
	def make_plot(self,hardcopy=False,hardcopydir='.',hardcopyprefix='', hardcopytype='png'):
		"""Generate a plot of the experimental data, and the fit if present.

		Arguments
		----------
		hardcopy : boolean
			Display the fit to the screen, or write it to a file?
		hardcopydir : string
			The directory to write the hardcopy to.
		hardcopyprefix : string
			A prefix to append to the experiment's title when generating the output plot filename.
		hardcopytype : string
			The output format for the plot (availability is dependent upon what backends matplotlib was compiled with).

		Returns
		-------
		None
		"""
		if not self.initialized:
			raise Exception("No data to plot. If experiment is synthetic, call sim.run() first.")
		
		from __init__ import MATPLOTLIB_BACKEND
		try:
			matplotlib.get_backend()
		except:
			if MATPLOTLIB_BACKEND != None:
				import matplotlib
				matplotlib.use(MATPLOTLIB_BACKEND)
		
		import matplotlib.pyplot as pyplot

		if hardcopy: fig = pyplot.figure()

		pyplot.clf()
		pyplot.title(self.title)
		pyplot.ylabel("%s/mol of %s"%(self.units,self.syringeRef))
		pyplot.xlabel("%s / %s"%(self.syringeRef,self.cellRef))

		# For convention, normalize the heat evolved as per mol of injected reference ligand
		tmpx = [ self.Concentrations[i][self.syringeRef]/self.Concentrations[i][self.cellRef] for i in xrange(self.npoints) if i not in self.skip ]
		tmpy = [ convert_from_J(self.units,self.dQ_exp[i])/self.Syringe[self.syringeRef]/self.injections[i] for i in xrange(self.npoints) if i not in self.skip ]
		
		if self.dQ_err is not None:
			tmpd = [ convert_from_J(self.units,self.dQ_err[i]) for i in xrange(self.npoints) if i not in self.skip ]
			pyplot.errorbar(tmpx,tmpy,yerr=tmpd,c='#000000',fmt='s')

		if self.spline is not None:
			tmpz = [ convert_from_J(self.units,self.spline[i])/self.Syringe[self.syringeRef]/self.injections[i] for i in xrange(self.npoints) if i not in self.skip ]
			pyplot.plot(tmpx,tmpz,c='g')

		if self.dQ_fit is not None:
			tmpy = [ convert_from_J(self.units,self.dQ_fit[i])/self.Syringe[self.syringeRef]/self.injections[i] for i in xrange(self.npoints) if i not in self.skip ]
			pyplot.plot(tmpx,tmpy,c='r')

		if len(self.skip) > 0:
			tmpx = [ self.Concentrations[i][self.syringeRef]/self.Concentrations[i][self.cellRef] for i in self.skip ]
			tmpy = [ convert_from_J(self.units,self.dQ_exp[i])/self.Syringe[self.syringeRef]/self.injections[i] for i in self.skip ]
			pyplot.errorbar(tmpx,tmpy,yerr=0,c='g',fmt='s')
			
		if self.Q_dil != 0:
			tmpx = [ self.Concentrations[i][self.syringeRef]/self.Concentrations[i][self.cellRef] for i in xrange(self.npoints) if i not in self.skip ]
			tmpy = [ convert_from_J(self.units,self.dQ_dil[i])/self.Syringe[self.syringeRef]/self.injections[i] for i in xrange(self.npoints) if i not in self.skip ]
			pyplot.plot(tmpx,tmpy,c='b')

		pyplot.draw()
		if hardcopy:
			fig.savefig( os.path.join(hardcopydir,"%s%s.%s"%(hardcopyprefix,self.title,hardcopytype)), bbox_inches='tight')
			pyplot.close(fig)
		else:
			pyplot.show()
			
	def export_to_file(self,path,units='cal',full=False):
		"""Export the attributes of the experiment to a file.

		Arguments
		---------
		path : string
			The filesystem path to write the output to.
		units : string
			The units to use in the export, must be either "J" for Joules (the default), or "cal" for calories.

		Returns
		-------
		None
			
		Notes
		-----
			Most instruments (for historical reasons) report binding enthalpies in calories. Thus, this is the default behavior for this method.
		"""
				
		from datetime import datetime

		assert units in _UNITS

		h = open(path,'w')
		h.write("# itcsim export of %s\n#\n"%(self.title))
		h.write("# Date	%s\n"%(datetime.ctime(datetime.today())))
		h.write("# units %s\n"%(units))
		h.write("# T	%.5f\n"%(self.T))
		h.write("# V0	%.5f\n"%(self.V0))
		for s in self.Cell:
			h.write("# Cell	%s	%.5E\n"%(s,self.Cell[s])) 
		for s in self.Syringe:
			h.write("# Syringe	%s	%.5E\n"%(s,self.Syringe[s])) 
		h.write("# Q_dil	%.5E\n"%(self.Q_dil))
		h.write("# skip %s\n"%(",".join(map(str,self.skip))))
		
		if full:
			h.write("#\n# Ivol	dQ_exp	dQ_err_exp dQ_fit	dQ_spline	skipped	%s	%s\n"%("\t".join(self.Cell.keys()),"\t".join(self.Syringe.keys())))
		else:
			h.write("#\n# Ivol	dQ_exp\n")
		
		if self.dQ_fit is None:
			fit = [0.0]*self.npoints
		else:
			fit = self.dQ_fit

		if self.spline is None:
			spline = [0.0]*self.npoints
		else:
			spline = self.spline

		if self.dQ_err is None:
			dQ_err = [0.0]*self.npoints
		else:
			dQ_err = self.dQ_err

		for i in xrange(self.npoints):
			cell	= ["%.5E"%(self.Concentrations[i][s]) for s in self.Cell]
			syringe	= ["%.5E"%(self.Concentrations[i][s]) for s in self.Syringe]
			if full:
				h.write("%.5f	%.5E	%.5E	%.5E	%.5E	%i	%s	%s\n"%(
					self.injections[i],
					convert_from_J(units,self.dQ_exp[i]),
					convert_from_J(units,dQ_err[i]),
					convert_from_J(units,dQ_fit[i]),
					convert_from_J(units,spline[i]),
					(i in self.skip),
					"\t".join(cell),
					"\t".join(syringe)
				))
			else:
				h.write("%.5f	%.5f\n"%(self.injections[i],convert_from_J('cal',self.dQ_exp[i])) )

		h.close()
		
	def get_chisq(self, Q, writeback=False):
		"""Calculate the goodness-of-fit between the provided data and the experimental data.

		Arguments
		---------
		Q : list of floats
			The predicted total heat at each injection point.
		writeback : boolean
			Update the experiment's simulated heat attribute (dQ_fit) with the provided Qs, as well as the chisq value?
			
		Returns
		-------
		float
			The goodness of the fit, as a reduced chi-square or as a sum-of-squares if a experimental error in dQ was not provided.
		"""
		
		dV = 0.0
		for i in xrange(self.npoints):
			dV += self.injections[i]
			Q[i] *= self.V0 * self.Concentrations[i][self.cellRef] # normalize total heat content to macromolecule concentration in the cell volume

		# obtain the change in cell heat between each titration point
		dQ = [0.0]*self.npoints
		for i in xrange(self.npoints):
			if i==0:
				dQ[i] = Q[i] + ( (self.injections[i]/self.V0)*(Q[i]/2.0) )
			else:
				dQ[i] = Q[i] + ( (self.injections[i]/self.V0)*((Q[i]+Q[i-1])/2.0) ) - Q[i-1]
			
			dQ[i] += self.dQ_dil[i] # add heat of dilution

		if self.dQ_err is None:
			dQ_err = [1.0]*self.npoints
		else:
			dQ_err = self.dQ_err

		# calculate reduced chi-square / SOS value using error estimate
		chisq = sum([(self.dQ_exp[i] - dQ[i])**2 / dQ_err[i]**2 for i in xrange(self.npoints) if i not in self.skip])
		chisq /= (self.npoints -len(self.skip)) # reduced chi-square corrected for n degrees of freedom
			
		if writeback:
			self.dQ_fit = dQ[:]
			self.chisq = chisq

		return chisq
		
class ITCExperiment(ITCExperimentBase):
	"""An object that encapsulates an empirical (i.e. "real") ITC experiment.
	
	Attributes
	----------
	spline : list of floats
		The smoothed spline (if any) to the experimental enthalpies, in Joules.
	"""
	
	def __init__(self, spline_pts=7, spline_order=1, *args, **kwargs ):
		"""The constructor for the empirical ITCExperiment class
	
		Arguments
		---------
		spline_pts : int
			Number of points to use in the Savitsky-Golay filter used to estimate errors in experimental enthalpies.
		spline_order : int
			Order of the Savitsky-Golay filter.
		*args
			Positional arguments to pass along to ITCExperimentBase.
		**kwargs
			Keyword arguments to pass along to ITCExperimentBase.
		"""

		ITCExperimentBase.__init__(self,*args,**kwargs)
	
		if self.dQ_err is None: # estimate errors by fitting spline
			counter,tmp = 0,[0]*self.npoints
			for i in xrange(self.npoints):
				if i not in self.skip:
					tmp[i] = counter
					counter += 1

			spl = savitzky_golay([self.dQ_exp[i] for i in xrange(self.npoints) if i not in self.skip], spline_pts, spline_order )
			err = scipy.std([self.dQ_exp[i] - spl[tmp[i]] for i in xrange(self.npoints) if i not in self.skip])
			self.dQ_err = [err]*self.npoints
			self.spline = [spl[tmp[i]] for i in xrange(self.npoints)]
		
		self.initialized = True

class ITCExperimentSynthetic(ITCExperimentBase):
	"""An object that encapsulates a synthetic (i.e. simulated) ITC experiment.
	
	Attributes
	----------
	noise : float
		The standard deviation (in percent) of normally-distributed noise to add to the pure simulated total evolved heat.
	initialized: boolean
		A flag that is set once the class's dQ_fit attribute has been set (see notes).
		
	Notes
	-----
		The first time the get_chisq() method of this class is called, the "experimental" dQ attribute is populated with the per-injection heats simulated by the model (+ the amount of specified noise).
		Subsequent calls will return the actual reduced chi-square difference between the experiment's calculated dQ and the model predictions, as expected.

		If noise is set to 0 (or None), then the get_chisq() method of this class will always return 1.0, as without realistic noise, the reduced chisquare is meaningless. 
	"""

	def __init__(self, injections, noise=None, *args, **kwargs):
		"""The constructor for the ITCExperimentSynthetic class
		
		Arguments
		---------
		injections : list of floats
			The injection volumes (in uL) for each point in the titration.
		noise : float
			The standard deviation (in percent) of normally-distributed noise to add to the pure simulated total evolved heat.
		*args
			Positional arguments to pass along to ITCExperimentBase
		**kwargs
			Keyword arguments to pass along to ITCExperimentBase
	
		"""
		
		ITCExperimentBase.__init__(self,injections=injections,dQ=[0.0]*len(injections),*args,**kwargs)
		
		if noise == None:
			self.noise = None	
		else:
			self.noise = convert_to_J(self.units,noise)*self.Syringe[self.syringeRef]
		self.initialized = False
		
	def get_chisq(self, Q, writeback):
		"""Monkeypatches the parent ITCExperimentBase's get_chisq() class method to update the class's dQ_exp attribute the first time the method is called.
		
		Arguments
		---------
		Q : list of floats
			The predicted total heat at each injection point.
		writeback : boolean
			Update the experiment's simulated heat attribute (dQ_fit) with the provided Qs?

		Returns
		-------
		float
			The goodness of the fit, as a reduced chi-square, unless the class's noise attribute is zero or none in which case 1.0 is returned (see class notes).
		"""
		if not self.initialized:
			
			if self.noise:
				self.dQ_err = [self.noise]*self.npoints
				ret = ITCExperimentBase.get_chisq(self, Q[:], writeback=True)
				self.dQ_exp = [scipy.random.normal(self.dQ_fit[i],self.noise) for i in xrange(self.npoints)]
			elif self.noise == 0 or self.noise == None:
				ITCExperimentBase.get_chisq(self, Q[:], writeback=True)
				self.chisq = 1.0
				self.dQ_exp = self.dQ_fit[:]
				ret = 1.0

			self.initialized = True			
			return ret

		if self.noise:
			return ITCExperimentBase.get_chisq(self, Q[:], writeback)
		else:
			ITCExperimentBase.get_chisq(self, Q[:], writeback)
			return 1.0