from thermo import *
from thermo import _UNITS

try:
	_tmp = OrderedDict()
except:
	from ordered_dict	import OrderedDict

class ITCModel():
	"""A base class that serves as the foundation for generation of per-injection binding enthalpies
	
	Attributes
	----------
	units : string
		The units to report binding enthalpies in. Note that internally, all calculations should be performed in SI units (e.g. Joules)
	params : ordered dict of strings
		The parameters of the model
	components : list of strings
		The names of the components that are involved in the binding model (e.g. the lattice/macromolecule and ligand)
	_param_meta : dict of tuples
		A storage variable that contains the current value and meta information for each model parameter
	_component_meta : dict of tuples
		A storage variable that contains information for each model component
	lattice_name : string
		The name of the binding component (for binary systems)
	ligand_name : string
		The name of the bound component (for binary systems)

	Notes
	-----
		Parameter types can be "n" (stoichiometry), "k" (a kinetic constant, not currently used), "dG" (change in free energy), "dH" (change in enthalpy), "dS" (change in entropy), or "dCp" (change in heat capacity)
		By convention, if lattice_name and ligand_name are not provided, they'll be set to the first and second (respectively) components set by add_component()
	"""
	
	def __init__(self, units="J", lattice_name=None, ligand_name=None):
		"""The constructor for the base ITCModel class. Child class constructors should call this parent constructor first and then probably do something else with an argument or two.
		
		Arguments
		---------
		units : string
			The units to use for energetic parameters (J, kcal, cal)
		lattice_name : string
			The name of the binding component (for binary systems)
		ligand_name : string
			The name of the bound component (for binary systems)

		Returns
		-------
		None

		""" 
		self.lattice_name = lattice_name
		self.ligand_name = ligand_name
		self.units = units
		
		self.params = OrderedDict()
		self.components = []
		self._param_meta = {}
		self._component_meta = {}
		
	def add_parameter(self, name, type, bounds=[None,None], default=0.0, linked='', description=''):
		"""Register a model parameter.
		
		Arguments
		---------
		name : string
			The name of the model parameter. Should be descriptive, like "dG1" (e.g. the delta G(ibbs) energy for binding mode 1).
		type : string
			The type, i.e. energy, stoichiometry, etc. See class notes.
		bound : tuple
			The high and low boundaries for the parameter. Depending on the fitting mode, these may or not be enforced.
		default : float
			The default value for the parameter
		linked : string
			(Not used yet)
		description: string
			A description of the model parameter. Only ever reported to the user for their information.
			
		Returns
		-------
		None
		"""
				
		assert name not in self.params.keys()
		assert type in ('n','k','dG','dH','dS','dCp')
		self.params[name] = default
		self._param_meta[name]	= (default,type,bounds,description,linked,type in ('dG','dH','dS','dCp'))
	
	def add_component(self, name, description=''):
		"""Register a model component.
		
		Arguments
		---------
		name : string
			The name of the model component (e.g. "Lattice", "Macromolecule", or "Ligand")
		description : string
			A description of the model component. Only ever reported to the user for their information.
		
		Returns
		-------
		None
		"""

		if self.lattice_name is None:
			self.lattice_name = name
		elif self.ligand_name is None:
			self.ligand_name = name

		self.components.append(name)
		self._component_meta[name] = (description)
	
	# used if we need any special startup or shutdown procedures
	def start(self):
		"""Class method stub for any necessary setup calls the model needs (scratch space, etc.). Actual models will need to overwrite this only if necessary.
		
		Arguments
		---------
		None
		
		Returns
		-------
		None
		"""
		pass
		
	def stop(self):
		"""Class method stub for any necessary shutdown calls the model needs (clearing scratch space, etc.). Actual models will need to overwrite this only if necessary.
		
		Arguments
		---------
		None
		
		Returns
		-------
		None
		"""
		pass
						
	# setters
	def set_units(self,units):
		"""Set the units the model will A) receive parameter values, and B) return dQ heats.
		
		Arguments
		---------
		units : string
			The units to use for model parameters and dQ heats
		
		Returns
		-------
		None
			
		Notes
		-----
			All internal calculations are still to be done in SI units (Joules)!
		"""
		
		assert units in _UNITS
				
		self.units = units		

	def set_param(self, name, value ):
		"""Set the value of the single named parameter.
		
		Arguments
		---------
		name : string
			The name of the model parameter
		value : float
			The value to set the named model parameter too
		
		Notes
		-----
			If the model parameter is an energy, it will be automatically converted to Joules for storage from whatever units are currently specified
		"""
				
		if self._param_meta[name][5]:
			self.params[name] = convert_to_J(self.units,value)
		else:
			self.params[name] = value

	def set_params(self, *args, **kwargs):
		"""Set the values of multiple parameters.
		
		Arguments
		---------
		*args
			Set the values of the model parameters in the order as returned by get_param_names().
		**kwargs
			Set the values of the named model parameters.
			
		Returns
		-------
		None
			
		Notes
		-----
			To prevent confusion, when setting parameter values positionally, all model values must be specified (in the correct order of course).
		"""
		
		for i in xrange(len(args)):
			assert len(args) == len(self.params)
			self.set_param( self.params.keys()[i], args[i] )
		for k,v in kwargs.iteritems():
			self.set_param( k, v )

	# getters
	def get_params(self,units=None):
		"""Return the value of all model parameters.
		
		Arguments
		---------
		units : string
			Use the specified units. If omitted, use the units already specified by the model (if applicable).
		
		Returns
		-------
		OrderedDict
			Parameter name-keyed dict of values.
		"""
		return OrderedDict( (name,self.get_param(name,units)) for name in self.params.keys() )
		
	def get_param(self,name,units=None):
		"""Return the value of the specified model parameter.
		
		Arguments
		---------
		name : string
			The name of the parameter to return.
		units : string
			Use the specified units. If omitted, use the units already specified by the model (if applicable).
		
		Returns
		-------
		float
			The value of the parameter.
		"""
			
		if not self._param_meta[name][5]:
			return self.params[name]
		if units:
			return convert_from_J(units,self.params[name])
		return convert_from_J(self.units,self.params[name])
		
	def get_param_names(self):
		"""Return a list of the parameter names used by the model.
		
		Arguments
		---------
		None
		
		Returns
		-------
		list of strings
			The names of the model parameters.
		"""
		
		return self.params.keys()
		
	def get_component_names(self):
		"""Return a list of the component names present in the model.
		
		Arguments
		---------
		None
		
		Returns
		-------
		list of strings
			The names of the model components.
		"""
		return self.components.keys()

	def get_param_type(self,name):
		"""Return the type of the named model parameter.
		
		Arguments
		---------
		name : string
			Name of the model parameter.
		
		Returns
		-------
		string
			The type of the requested parameter (e.g. "n", "k", "dG", "dH", or "dCp"). See class notes.
		"""
		return self._param_meta[name][1]

	def get_param_bounds(self,name):
		"""Return the boundaries for the named model parameter.
		
		Arguments
		---------
		name : string
			Name of the model parameter.
		
		Returns
		-------
		tuple of floats
			The high and low boundaries (if specified, otherwise they are None) for the requested model parameter.
		"""
		return self._param_meta[name][2]

	def get_param_description(self,name):
		"""Return the description for the named model parameter.
		
		Arguments
		---------
		name : string
			Name of the model parameter.
		
		Returns
		-------
		string
			The description of the model parameter previously set by add_parameter().
		"""
		return self._param_meta[name][3]
									
	def __str__(self):
		"""Stringifies the model and its current state.
		
		Arguments
		---------
		None
		
		Returns
		-------
		string
			A user-interpretable string that describes the parameters of the model, their values, and other info.
		"""
		
		types = {'n':"Stoichiometry",'k':"Rate constant",'dG':"Free energy",'dH':"Enthalpy",'dS':"Entropy",'dCp':"Heat capacity"}
		units = {'n':"sites",'k':"1/s or 1/s/mol",'dG':"%s/mol"%self.units,'dH':"%s/mol"%self.units,'dS':"%s/mol/K"%self.units,'dCp':"%s/mol/K"%self.units}
		ret = "\nModel: %s.%s\n\nDescription:\n%s\n"%(self.__module__,self.__class__.__name__, self.__doc__)
		ret+= "\nComponents:\n"
		ret+= "Index	Name	Description\n"
		for i,c in enumerate(self.components):
			ret+="%i)\t%s\t%s\n"%(i+1,c,self._component_meta[c])
		ret+= "\nParameters:\n"
		ret+= "Index	Param     Type                Value               Description\n"
		for i,k in enumerate(self.params):
			ret+="%i)\t%s%s%s%s\n" % (
				i+1,
				k.ljust(10),
				types[self.get_param_type(k)].ljust(20),
				(("%.3f"%self.get_param(k))+" "+units[self.get_param_type(k)]).ljust(20),
				self._param_meta[k][3]
				)

		return ret

	def Q(self,T0,T,concentrations): # return the enthalpy at each set of component concentrations
		"""Return the total binding heat at each injection predicted by the model and its current parameter values. This is just a stub, and child classes should implement their own calculations of Q.
		
		Arguments
		---------
		T0 : float
			The reference temperature to be used for the model (used in temperature-dependent dG, dH, dCp).
		T : float
			The temperature the titration was performed at.
		concentrations : list of dicts
			The concentration of components at each injection point.
			
		Returns
		-------
		list of floats
			The total heat in the system at each injection point.
		"""
		raise NotImplementedError("Valid ITC models should implement this!")
