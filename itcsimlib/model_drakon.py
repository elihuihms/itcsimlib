from thermo import *
from model_ising import Ising

class DRAKONIsingModel(Ising):
	"""A wrapper/simplification interface for an Ising model for use in DRAKON flow diagram-based models.
	
	Attributes
	----------
	None
	"""
	
	def setup(self):
		raise NotImplementeError("Valid DRAKON models must implement this method")
		
	def site(self, i, j ):
		raise NotImplementeError("Valid DRAKON models must implement this method")
	
	def __init__(self):
		"""Constructor for the model.
		
		Arguments
		---------
		None
		"""
		
		# convenience functions
		self.occupied = self.get_site_occupancy
		self.neighbor = self.get_site_occupancy
		self.set_param = self.set_parameter
		
		self.setup()
	
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
			
	def set_parameter(self, name, value):
		"""Alias for the ITCModel set_parameter() method. Updates the class attribute value as well.
		
		See the ITCModel.add_parameter() method for the argument list.
		"""

		Ising.set_param(self, name, value) # convert units if necessary
		setattr(self, name, self.params[name])
		
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
		
		for i in xrange(self.nconfigs):
			self.gibbs[i],self.enthalpies[i] = 0.0,0.0
			self.config_expressions[i] = 0
			
			for j in xrange(self.nsites):
				self.site(i, j)
		
	def add_dG(self, i, dG, dH=None, dCp=None):
		"""Convenience function for DRAKON models that allows incrementing the gibbs free energy of a configuration.
		Also permits temperature-dependent van't Hoff correction, if dH and dCp are not None.
		
		Arguments
		---------
		i : integer
			The index of the configuration
		dG : float
			The free energy change to add to the configuration, or the reference free energy change if dH and dCp are provided.
		dH : float or None
			The enthalpy change to use in the van't Hoff correction.
		dCp : float or None
			The heat capacity change to use in the van't Hoff correction.
		"""
		
		if dH==None or dCp==None:
			self.gibbs[i] += dG
		else:
			self.gibbs[i] += dG_vant_Hoff( dG, dH, dCp, self._T, self._T0 )
		
	def add_dH(self, i, dH, dCp=None):
		"""Convenience function for DRAKON models that allows incrementing the enthalpy of a configuration.
		Also permits temperature-dependent van't Hoff correction, if dCp is not None.
		
		Arguments
		---------
		i : integer
			The index of the configuration			
		dH : float
			The enthalpy change to add to the configuration.
		dCp : float or None
			The heat capacity change to use in the van't Hoff correction.
		"""
		
		if dCp==None:
			self.enthalpies[i] += dH
		else:
			self.enthalpies[i] += dH_vant_Hoff( dH, dCp, self._T, self._T0 )
		
	

		