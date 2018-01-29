from thermo import *
from model_ising import Ising

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
		
		self.setup()
	
	def initialize(self, *args, **kwargs):
		"""Alias for the parent class constructor, used to explicity call init in DRAKON.
		See the Ising.__init__() method for the argument list.
		"""
		
		Ising.__init__(self, *args, **kwargs)
	
	def add_parameter(self, name, type=None, **kwargs):
		"""Convenience function for the ITCModel add_parameter() method.
		
		See the ITCModel.add_parameter() method for the argument list.
		"""
			
		if type==None: # a bit cheesy
			type = name[0:2]
			
		Ising.add_parameter(self, name, type, **kwargs)
		
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
		
		config_energy_function = getattr(self, "config", False)
		
		for i in xrange(self.nconfigs):
			self.gibbs[i],self.enthalpies[i] = 0.0,0.0
			self.config_expressions[i] = 0
			
			if config_energy_function:
				self.configuration(i,*self.configs[i])
			else:
				for j in xrange(self.nsites):
					self.site(i, j)				
		
	def add_dG(self, i, dG, dH=None, dCp=None):
		"""Convenience function for DRAKON models to increment the gibbs free energy of a configuration.
		This function also permits temperature-dependent van't Hoff correction, if dH and dCp are not None.
		
		Arguments
		---------
		i : integer
			The index of the configuration
		dG : string
			The name of the free energy change parameter to add to the configuration, or the reference free energy change if dH and dCp are provided.
		dH : string or None
			The name of the  enthalpy change parameter to use in the van't Hoff correction.
		dCp : string or None
			The name of the heat capacity change parameter to use in the van't Hoff correction.
		"""
		
		if dH==None or dCp==None:
			self.gibbs[i] += self.get_param(dG,units="J")
		else:
			self.gibbs[i] += dG_vant_Hoff( self.get_param(dG,units="J"), self.get_param(dH,units="J"), self.get_param(dCp,units="J"), self._T, self._T0 )
		
		self.config_expressions[i] += self.parameter_symbols[dG]
		
	def add_dH(self, i, dH, dCp=None):
		"""Convenience function for DRAKON models to increment the enthalpy of a configuration.
		This function permits temperature-dependent van't Hoff correction, if dCp is not None.
		
		Arguments
		---------
		i : integer
			The index of the configuration			
		dH : string
			The name of the enthalpy change parameter to add to the configuration.
		dCp : string or None
			The name of the heat capacity change parameter to use in the van't Hoff correction.
		"""
		
		if dCp==None:
			self.enthalpies[i] += self.get_param(dH,units="J")
		else:
			self.enthalpies[i] += dH_vant_Hoff( self.get_param(dH,units="J"), self.get_param(dCp,units="J"), self._T, self._T0 )
		
