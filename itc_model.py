from thermo import *

try:
	_tmp = OrderedDict()
except:
	from ordered_dict	import OrderedDict

class ITCModel():

	def __init__(self):		
		self.params = OrderedDict()
		self.components = []
		self.units = "J"

		self._param_meta, self._component_meta = {},{}
		
	def add_parameter(self, name, type, bounds=[None,None], default=0.0, linked='', description=''):
		assert name not in self.params.keys()
		assert type in ('n','k','dG','dH','dS','dCp')
		self.params[name] = default
		self._param_meta[name]	= (default,type,bounds,description,linked,type in ('dG','dH','dS','dCp'))
	
	def add_component(self, name, description=''):
		self.components.append(name)
		self._component_meta[name] = (description)
	
	# used if we need any special startup or shutdown procedures
	def start(self):
		pass
		
	def stop(self):
		pass
						
	# setters
	def set_units(self,units):
		self.units = units		

	def set_param(self, param, value ):
		if self._param_meta[param][5]:
			self.params[param] = convert_to_J(self.units,value)
		else:
			self.params[param] = value

	def set_params(self, *args, **kwargs):
		for i in xrange(len(args)):
			assert len(args) == len(self.params)
			self.set_param( self.params.keys()[i], args[i] )
		for k,v in kwargs.iteritems():
			self.set_param( k, v )

	# getters
	def get_params(self,units=None):
		return OrderedDict( (name,self.get_param(name,units)) for name in self.params.keys() )
		
	def get_param(self,name,units=None):
		if not self._param_meta[name][5]:
			return self.params[name]
		if units:
			return convert_from_J(units,self.params[name])
		return convert_from_J(self.units,self.params[name])
		
	def get_param_names(self):
		return self.params.keys()
		
	def get_component_names(self):
		return self.components.keys()

	def get_param_type(self,name):
		return self._param_meta[name][1]

	def get_param_bounds(self,name):
		return self._param_meta[name][2]

	def get_param_description(self,name):
		return self._param_meta[name][3]
									
	def __str__(self):
		types = {'n':"Stoichiometry",'k':"Rate constant",'dG':"Free energy",'dH':"Enthalpy",'dS':"Entropy",'dCp':"Heat capacity"}
		units = {'n':"sites",'k':"1/s or 1/s/mol",'dG':"%s/mol"%self.units,'dH':"%s/mol"%self.units,'dS':"%s/mol/K"%self.units,'dCp':"%s/mol/K"%self.units}
		ret = "\nITCModel \"%s\"\nDescription:\n%s\n"%(self.__module__, self.__doc__)
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
		raise NotImplementedError("Valid ITC models should implement this!")
