import math
import scipy
import scipy.optimize

from itc_model	import ITCModel
from thermo		import _R
from thermo		import *

try:
	_tmp = OrderedDict()
except:
	from ordered_dict	import OrderedDict

class Ising(ITCModel):
	"""An ITC model based on ligand binding to an Ising lattice.
	
	Attributes
	----------
	nconfigs : int
		The total number of possible configurations of the lattice.
	configs : list of list of ints
		The occupation state (1=bound, 0=unbound) of each site in each lattice configuration.
	bound : list of ints
		The number of ligands bound to each lattice configuration.
	weights : list of floats
		The absolute (normalized) probability of each configuration.
	gibbs : list of floats
		The free energy of each configuration.
	enthalpies: list of floats
		The enthalpy of each configuration.
	precision: float
		The precision in ligand concentration required for convergence during set_probabilities()
	parameter_symbols : dict of sympy symbols
		Convenience container for the sympy symbols corresponding to the model parameter names, ultimately used to construct the symbolic configuration expressions.
	config_expressions : list of sympy expressions
		The symbolic expression of the configuration's free energy.
	"""
		
	def __init__(self,nsites=3,circular=True):
		"""The constructor for the base Ising binding model.
		
		Arguments
		---------
		nsites : int
			The number of potential binding sites in the lattice.
		circular : boolean
			Is the lattice circular?
		"""
		ITCModel.__init__(self)
		
		self.nsites,self.circular = nsites,circular
		
		self.nconfigs	= 2**self.nsites
		self.configs	= [ [int(s) for s in ("{0:0%ib}"%(self.nsites)).format(i)] for i in xrange(self.nconfigs) ]
		self.bound		= [ c.count(1) for c in self.configs ] # number of bound sites
		self.weights	= [0.0]*self.nconfigs # probability of each config
		self.gibbs		= [0.0]*self.nconfigs # free energy of each config
		self.enthalpies	= [0.0]*self.nconfigs # enthalpic energy of each config
		self.precision	= 1E-9 # the precision in ligand concentration required for convergence during set_probabilities()
		
		self.parameter_symbols = {} # model parameters used during partition function generation
		self.config_expressions = [ 0 for i in xrange(self.nconfigs) ] # expressions of configuration free energies using the parameter symbols
		
		if self.circular:		
			self.add_component('Lattice',description='A circular lattice with %i binding sites'%(self.nsites))
		else:
			self.add_component('Lattice',description='A linear lattice with %i binding sites'%(self.nsites))
		self.add_component('Ligand',description='A lattice-binding ligand')

	def add_parameter(self, name, type, **kwargs):
		"""Wrapper for the typical ITC model add_parameter, with the added tweak that a sympy symbol is created for eventually generating the model's partition function."""
		ITCModel.add_parameter(self, name, type, **kwargs)
		
		try:
			import sympy
		except:
			return	
		self.parameter_symbols[name] = sympy.symbols(name)

	def get_site_occupancy(self,config,site):
		"""Return whether a given site is occupied or not, correctly accounting for circular lattices and lattice indicies <0 or >n.
		
		Arguments
		---------
		config : int
			The lattice configuration index
		site : int
			The lattice site index
		
		Returns
		-------
		boolean or None
			Whether the site is occupied (bound), or None if the lattice is nonlinear and the site index is meaningless
		"""
		if site < 0:
			if not self.circular:
				return None
			return self.configs[config][site+self.nsites] == 1
		if site >= self.nsites:
			if not self.circular:
				return None
			return self.configs[config][site%self.nsites] == 1
		return self.configs[config][site] == 1

	def set_probabilities(self,totalP,totalL,T):
		"""Set the normalized weights (probabilities) of each configuration at the specified free energies and component concentrations.

		Arguments
		---------
		totalP : float
			The total concentration of binding lattices.
		totalL : float
			The total concentration of ligands.
		T : float
			The experimental temperature.
			
		Returns
		-------
		float
			The free concentration of ligand.
		"""

		def _freeL_dev(freeL): # 
			# Return the deviation between predicted and actual free ligand concentration
			
			# set the probability of each configuration at the free ligand concentration
			self.weights = [math.exp( (-1.0 * self.gibbs[i]) / ( _R * T ) ) * freeL**self.bound[i] for i in xrange(self.nconfigs) ]
			total = sum(self.weights)
			self.weights = [ self.weights[i] / total for i in xrange(self.nconfigs) ]
			
			# concentration of sites in bound state
			bound = sum( [totalP * self.weights[i] * self.bound[i] for i in xrange(self.nconfigs)] )

			return totalL - (freeL + bound)
		
		# find where the deviation between actual and test free ligand is zero. Use zero free and total ligand as our bracketing guesses
		freeL = scipy.optimize.brentq( _freeL_dev, 0.0, totalL, xtol=self.precision, disp=True )
		
		# make sure before we quit that we set the weights at the correct free ligand conc
		_freeL_dev( freeL )
		
		return freeL
		
	def Q(self,T0,T,concentrations):
		"""Return the enthalpy of the system at each of the specified component concentrations.
		
		Arguments
		---------
		T0 : float
			The reference temperature of the simulation.
		T : float
			The temperature of the experiment to simulate.
		concentrations : list of dicts
			The concentrations of each component at each titration point.
		
		Returns
		-------
		list of floats
			The total enthalpy of the system at each injection point.
		"""
		
		# set the free energies (and enthalpic energies if necessary) of each configuration
		self.set_energies(T0,T)

		# calculate the enthalpy at each set of conditions
		Q = [0.0]*len(concentrations)
		for i,c in enumerate(concentrations):
			# set the weights (probabilities) of each lattice configuration
			self.set_probabilities(c['Lattice'],c['Ligand'],T)
			
			# enthalpy is sum of all weighted enthalpies of the lattices
			Q[i] = sum( [self.weights[j] * self.enthalpies[j] for j in xrange(self.nconfigs)] )

		return Q
			
	def get_partition_function(self,substitute_Ks=True,full_simplify=True):
		"""Return the partition function of the binding model as a sympy expression.
		
		Arguments
		---------
		substitute_Ks : boolean
			Improve condensation of the partition function by converting free energies to binding constants, where possible.
		full_simplify : boolean
			Perform a final reduction of terms for the compiled partition function?
		
		Returns
		-------
		sympy object
			The symbolic partition function for the model.
		"""
		import sympy

		self.set_energies(273.15,273.15) # ensure that this is run at least once to populate config_terms

		L,R,T = sympy.symbols("L R T") # ligand, gas constant, temp
				
		config_expressions = [ None for i in xrange(self.nconfigs) ] # convert all configuration free energies to effective K(a)s
		for i in xrange(self.nconfigs):
			config_expressions[i] = sympy.exp(self.config_expressions[i] / (R*T))

		bound_expressions = [ 0 for i in xrange(self.nsites+1) ] # sum configuration Ks at each stoichiometry
		for i in xrange(self.nconfigs):
			bound_expressions[self.bound[i]] += config_expressions[i]
				
		ret = 0
		for i in xrange(self.nsites+1):
			if substitute_Ks:
				bound_expressions[i] = sympy.expand(bound_expressions[i]) # expand first in order to effectively combine terms later	
				for p in self.params: # replace simple intrinsic or multiplicative binding factors
					bound_expressions[i] = bound_expressions[i].subs( sympy.exp(self.parameter_symbols[p] / (R*T)), sympy.symbols("K_%s"%p) )

			bound_expressions[i] = sympy.simplify(bound_expressions[i]) # simplify each stoichiometric expression

			bound_expressions[i] = bound_expressions[i] * (L**i) # don't forget ligand concentration
			ret = ret + bound_expressions[i] # build full partition function
		
		if full_simplify:
			ret = sympy.simplify(ret)
		
		return ret

	def set_energies(self,T0,T):
		"""Set the free and enthalpic energy of each lattice configuration using the current parameter values. This method is a stub that child classes should replace.
		
		Arguments
		---------
		T0 : float
			The reference temperature of the simulation.
		T : float
			The current temperature of the system.
		
		Returns
		-------
		None
		"""
		self.reset_partition_function()
		
		for i in xrange(self.nconfigs):
			self.gibbs[i], self.enthalpies = 0.0, 0.0
			
		raise NotImplementedError("Valid ITC Ising models should implement this!")

class FullAdditive(Ising):
	"""An Ising-type model, in which ligands bind to either a linear or circular lattice. Coupling can occur to both unoccupied and occupied lattice points."""

	def __init__(self,nsites=3,circular=1):
		Ising.__init__(self,nsites,circular)
		
		self.add_parameter( 'dG0',	'dG',	description='Intrinsic free energy change upon binding to a site.' )
		self.add_parameter( 'dGa',	'dG',	description='Additional free energy change upon binding to a site flanked by an unoccupied site' )
		self.add_parameter( 'dGb',	'dG',	description='Additional free energy change upon binding to a site flanked by an occupied site' )
		self.add_parameter( 'dH0',	'dH',	description='Intrinsic enthalpy change upon binding to a site.' )
		self.add_parameter( 'dHa',	'dH',	description='Additional enthalpy change upon binding to a site flanked by an unoccupied site' )
		self.add_parameter( 'dHb',	'dH',	description='Additional enthalpy change upon binding to a site flanked by an occupied site' )
		self.add_parameter( 'dCp0',	'dCp',	description='Intrinsic heat capacity change upon binding to a site.' )
		self.add_parameter( 'dCpa',	'dCp',	description='Additional heat capacity change upon binding to a site flanked by an unoccupied site' )
		self.add_parameter( 'dCpb',	'dCp',	description='Additional heat capacity change upon binding to a site flanked by an occupied site' )

	def set_energies(self,T0,T):
		dG0 = dG_vant_Hoff( self.params['dG0'], self.params['dH0'], self.params['dCp0'], T, T0 )
		dGa = dG_vant_Hoff( self.params['dGa'], self.params['dHa'], self.params['dCpa'], T, T0 )
		dGb = dG_vant_Hoff( self.params['dGb'], self.params['dHb'], self.params['dCpb'], T, T0 )
		dH0 = dH_vant_Hoff( self.params['dH0'], self.params['dCp0'], T, T0 )
		dHa = dH_vant_Hoff( self.params['dHa'], self.params['dCpa'], T, T0 )
		dHb = dH_vant_Hoff( self.params['dHb'], self.params['dCpb'], T, T0 )
		
		for i in xrange(self.nconfigs):
			self.gibbs[i],self.enthalpies[i] = 0.0,0.0
			self.config_expressions[i] = 0
			
			for j in xrange(self.nsites):

				if self.get_site_occupancy(i,j): # is site occupied?
					self.gibbs[i]+=dG0
					self.enthalpies[i]+=dH0
					self.config_expressions[i] += self.parameter_symbols['dG0']
					
					if self.get_site_occupancy(i,j+1): # is the next neighboring site occupied?
						self.gibbs[i]+=dGb
						self.enthalpies[i]+=dGb
						self.config_expressions[i] += self.parameter_symbols['dGb']

					elif self.circular:
						self.gibbs[i]+=dGa
						self.enthalpies[i]+=dGa
						self.config_expressions[i] += self.parameter_symbols['dGa']
					
					if self.get_site_occupancy(i,j-1): # is previous neighboring site occupied?
						pass # Note: this avoids double counting, and thus is implemented as in Saroff & Kiefer

					elif self.circular:
						self.gibbs[i]+=dGa
						self.enthalpies[i]+=dGa
						self.config_expressions[i] += self.parameter_symbols['dGa']
		return

class HalfAdditive(Ising):
	"""An Ising-type model, in which ligands bind to either a linear or circular lattice. Coupling only occurs between occupied lattice points."""

	def __init__(self,nsites=3,circular=1):
		Ising.__init__(self,nsites,circular)
		
		self.add_parameter( 'dG0',	'dG',	description='Intrinsic free energy change upon binding to a site.' )
		self.add_parameter( 'dGb',	'dG',	description='Additional free energy change upon binding to a site flanked by an occupied site' )
		self.add_parameter( 'dH0',	'dH',	description='Intrinsic enthalpy change upon binding to a site.' )
		self.add_parameter( 'dHb',	'dH',	description='Additional enthalpy change upon binding to a site flanked by an occupied site' )
		self.add_parameter( 'dCp0',	'dCp',	description='Additional heat capacity change upon binding to a site flanked by an unoccupied site' )
		self.add_parameter( 'dCpb',	'dCp',	description='Additional heat capacity change upon binding to a site flanked by an occupied site' )

	def set_energies(self,T0,T):
		dG0 = dG_vant_Hoff( self.params['dG0'], self.params['dH0'], self.params['dCp0'], T, T0 )
		dGb = dG_vant_Hoff( self.params['dGb'], self.params['dHb'], self.params['dCpb'], T, T0 )
		dH0 = dH_vant_Hoff( self.params['dH0'], self.params['dCp0'], T, T0 )
		dHb = dH_vant_Hoff( self.params['dHb'], self.params['dCpb'], T, T0 )
		
		for i in xrange(self.nconfigs):
			self.gibbs[i],self.enthalpies[i] = 0.0,0.0
			self.config_expressions[i] = 0

			for j in xrange(self.nsites):

				if self.get_site_occupancy(i,j): # is site occupied?
					self.gibbs[i]+=dG0
					self.enthalpies[i]+=dH0
					self.config_expressions[i] += self.parameter_symbols['dG0']
					
					if self.get_site_occupancy(i,j+1): # is the next neighboring site occupied?
						self.gibbs[i]+=dGb
						self.enthalpies[i]+=dGb
						self.config_expressions[i] += self.parameter_symbols['dGb']
		return

class NonAdditive(Ising):
	"""An Ising-type model, in which ligands bind to either a linear or circular lattice. Binding energy depends upon whether zero, one, or both neighboring sites are occupied."""

	def __init__(self,nsites=3,circular=1):
		Ising.__init__(self,nsites,circular)
		
		self.add_parameter( 'dGX',	'dG',	description='Free energy change upon binding to a site flanked by two unoccupied' )
		self.add_parameter( 'dGY',	'dG',	description='Free energy change upon binding to a site flanked by one occupied' )
		self.add_parameter( 'dGZ',	'dG',	description='Free energy change upon binding to a site flanked by two occupied' )
		self.add_parameter( 'dHX',	'dH',	description='Enthalpy change upon binding to a site flanked by two unoccupied' )
		self.add_parameter( 'dHY',	'dH',	description='Enthalpy change upon binding to a site flanked by one occupied' )
		self.add_parameter( 'dHZ',	'dH',	description='Enthalpy change upon binding to a site flanked by two occupied' )
		self.add_parameter( 'dCpX',	'dCp',	description='Change in heat capacity upon binding to a site flanked by two unoccupied' )
		self.add_parameter( 'dCpY',	'dCp',	description='Change in heat capacity upon binding to a site flanked by one occupied' )
		self.add_parameter( 'dCpZ',	'dCp',	description='Change in heat capacity upon binding to a site flanked by two occupied' )

	def set_energies(self,T0,T):
		dGX = dG_vant_Hoff( self.params['dGX'], self.params['dHX'], self.params['dCpX'], T, T0 )
		dGY = dG_vant_Hoff( self.params['dGY'], self.params['dHY'], self.params['dCpY'], T, T0 )
		dGZ = dG_vant_Hoff( self.params['dGZ'], self.params['dHZ'], self.params['dCpZ'], T, T0 )
		dHX = dH_vant_Hoff( self.params['dHX'], self.params['dCpX'], T, T0 )
		dHY = dH_vant_Hoff( self.params['dHY'], self.params['dCpY'], T, T0 )
		dHZ = dH_vant_Hoff( self.params['dHZ'], self.params['dCpZ'], T, T0 )
		
		for i in xrange(self.nconfigs):
			self.gibbs[i],self.enthalpies[i] = 0.0,0.0
			self.config_expressions[i] = 0
			
			for j in xrange(self.nsites):
				if self.get_site_occupancy(i,j): # is site occupied?
					
					if self.get_site_occupancy(i,j+1): # is the next neighboring site occupied?
						if self.get_site_occupancy(i,j-1): # is previous neighboring site occupied?
							self.gibbs[i]+= dGZ
							self.enthalpies[i]+= dHZ
							self.config_expressions[i] += self.parameter_symbols['dGZ']
							
						elif self.circular:
							self.gibbs[i]+= dGY
							self.enthalpies[i]+= dHY
							self.config_expressions[i] += self.parameter_symbols['dGY']
							
					elif self.get_site_occupancy(i,j-1):
						self.gibbs[i]+= dGY
						self.enthalpies[i]+= dHY
						self.config_expressions[i] += self.parameter_symbols['dGY']
						
					elif self.circular:
						self.gibbs[i]+= dGX
						self.enthalpies[i]+= dHX
						self.config_expressions[i] += self.parameter_symbols['dGX']
						
		return
