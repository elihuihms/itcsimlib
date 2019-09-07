"""Helper classes and functions that permit use of DRAKON to build statistical thermodynamics models.


"""

from .thermo import *
from .model_ising import Ising


class DRAKONIsingModel(Ising):
	"""A wrapper/simplification interface for an Ising model for use in DRAKON flow diagram-based models.
	
	Attributes
	----------
	None
	"""
	
	def setup(self):
		raise NotImplementedError("Valid DRAKON models must implement this method")
	
	"""
	def site(self, i, j ):
		raise NotImplementedError("Valid DRAKON models must implement site() or config()")
		
	def configuration(self, *args):
		raise NotImplementedError("Valid DRAKON models must implement site() or config()")
	"""
	
	def __init__(self):
		"""Constructor for the model.
		
		Arguments
		---------
		None
		"""
		
		# convenience functions
		self.occupied = self.get_site_occupancy
		self.neighbor = self.get_site_occupancy
		self.set_parameter = self.set_param
		
		self.setup()

		#assert ("configuration" in dir(self) and "site" in dir(self)) == False
	
	def initialize(self, *args, **kwargs):
		"""Alias for the parent class constructor, used to explicity call init in DRAKON.
		See the Ising.__init__() method for the argument list.
		"""
		
		Ising.__init__(self, *args, **kwargs)
	
	def add_parameter(self, name, type, **kwargs):
		"""Alias for the ITCModel add_parameter() method. This method also adds the parameter to the class attributes.
		
		See the ITCModel.add_parameter() method for the argument list.
		
		Note that it's usually bad practice to pollute the namespace with random attributes.
		Here, it's a calculated risk for more grokable DRAKON diagrams."""
		
		if getattr(self, name, "unassigned") != "unassigned":
			raise Exception("The parameter name \"%s\" is already used in this class."%(args[0]))
		else:
			setattr(self, name, 0.0)
			
		Ising.add_parameter(self, name, type, **kwargs)
			
	def set_param(self, name, value):
		"""Alias for the ITCModel set_parameter() method. Updates the class attribute value as well.
		
		See the ITCModel.add_parameter() method for the argument list.
		"""

		Ising.set_param(self, name, value) # convert units if necessary
		setattr(self, name, self.params[name])
		
	def count_occupied(self, i):
		"""Convenience function for DRAKON models - getter for bound[].
		
		Arguments
		---------
		i : integer
			The index of the configuration.
			
		Returns
		-------
		integer : the number of occupied sites in the configuration.
		"""
		
		return self.bound[i]
		
	def set_energies(self, T0, T):
		"""Set the gibbs and enthalpic energy of each lattice configuration using the DRAKON site() method.
		
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
		self._T0,self._T = T0,T
		
		config_energy_function = getattr(self, "configuration", False)
		
		for i in range(self.nconfigs):
			self.gibbs[i],self.enthalpies[i] = 0.0,0.0
			self.config_expressions[i] = 0
			
			if config_energy_function:
				self.configuration(i)
			else:
				for j in range(self.nsites):
					self.site(i, j)		
		
	def add_dG(self, i, dG, dH=None, dCp=None):
		"""Convenience function for DRAKON models to increment the gibbs free energy of a configuration.
		This function also permits temperature-dependent van't Hoff correction, if dH and dCp are not None.
		Alternatively, an expression consisting of existing model parameters may be provided for each argument, but this will result in slow model evaluation as the expression will need to be eval()'d for every configuration. 
		
		Arguments
		---------
		i : integer
			The index of the configuration
		dG : string
			The name of the free energy change parameter to add to the configuration, the reference free energy change if dH and dCp are provided, or an expression of other existing model parameters.
		dH : string or None
			The name of the  enthalpy change parameter to use in the van't Hoff correction.
		dCp : string or None
			The name of the heat capacity change parameter to use in the van't Hoff correction.
		"""
		
		if dG in self.params:
			self.config_expressions[i] += self.parameter_symbols[dG]
			dG = self.get_param(dG,units="J")
		else: # expression to evaluate
			dG_v,dG_s = dG,dG # dG_values, dG_symbols
			for p in self.params:
				dG_v = dG_v.replace(p,"self.get_param('%s',units='J')"%p)
				dG_s = dG_s.replace(p,"self.parameter_symbols['%s']"%p)
			dG = eval(dG_v)
			self.config_expressions[i] += eval(dG_s)

		if dH in self.params:
			dH = self.get_param(dH,units="J")
		elif dH is None:
			pass
		else:
			for p in self.params:
				dH = dH.replace(p,"self.get_param('%s',units='J')"%p)
			dH = eval(dH)

		if dCp in self.params:
			dCp = self.get_param(dCp,units="J")
		elif dCp is None:
			pass
		else:
			for p in self.params:
				dCp = dCp.replace(p,"self.get_param('%s',units='J')"%p)
			dCp = eval(dCp)
		
		if dH==None or dCp==None:
			self.gibbs[i] += dG
		else:
			self.gibbs[i] += dG_vant_Hoff( dG, dH, dCp, self._T, self._T0 )
		
	def add_dH(self, i, dH, dCp=None):
		"""Convenience function for DRAKON models to increment the enthalpy of a configuration.
		This function permits temperature-dependent van't Hoff correction, if dCp is not None.
		Alternatively, an expression consisting of existing model parameters may be provided, but this may result in slow model execution, as the expression will need to be eval()'d for every configuration. 

		Arguments
		---------
		i : integer
			The index of the configuration			
		dH : string
			The name of the enthalpy change parameter to add to the configuration, the reference enthalpy if dCp is provided, or an expression of other existing model parameters.
		dCp : string or None
			The name of the heat capacity change parameter to use in the van't Hoff correction.
		"""
		
		if dH in self.params:
			dH = self.get_param(dH,units="J")
		else:
			for p in self.params:
				dH = dH.replace(p,"self.get_param('%s',units='J')"%p)
			dH = eval(dH)

		if dCp in self.params:
			dCp = self.get_param(dCp,units="J")
		elif dCp is None:
			pass
		else:
			for p in self.params:
				dCp = dCp.replace(p,"self.get_param('%s',units='J')"%p)
			dCp = eval(dCp)
		
		if dCp==None:
			self.enthalpies[i] += dH
		else:
			self.enthalpies[i] += dH_vant_Hoff( dH, dCp, self._T, self._T0 )
		
