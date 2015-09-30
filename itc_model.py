try:
	_tmp = OrderedDict()
except:
	from ordered_dict	import OrderedDict

class ITCModel():

	def __init__( self ):
		self.params = OrderedDict()
		self._meta = {}
		self._order = []

	def add_param(self, name, type, bounds=[None,None], default=0.0, linked='', description=''):
		assert name not in self.params.keys()
		assert type in ('n','dG','dH','dS','dCp')
		self.params[name]	= default
		self._meta[name]	= (default,type,bounds,description,linked)
		self._order.append(name)

	# setters
	def set_params(self, *args, **kwargs):
		for i in xrange(len(args)):
			assert len(args) == len(self.params)
			self.params[self._order[i]] = args[i]
		for k,v in kwargs.iteritems():
			self.params[k] = v

	# getters
	def get_param_names(self):
		return self.params.keys()

	def get_param_type(self,name):
		return self._meta[name][1]

	def get_param_bounds(self,name):
		return self._meta[name][2]

	def get_param_description(self,name):
		return self._meta[name][3]

	def get_format(self, T, T_ref):
		return self.params.values()

	def get_worker(self,in_queue,out_queue):
		return None

	def __str__(self):
		types = {'n':"Stoichiometry",'dG':"Free energy",'dH':"Enthalpy",'dS':"Entropy",'dCp':"Heat capacity"}
		units = {'n':"sites",'dG':"J/mol",'dH':"J/mol",'dS':"J/mol/K",'dCp':"J/mol*K"}
		ret = "\nITCModel \"%s\" Description: %s\n"%(self.__module__, self.__doc__)
		ret+= "Param     Type                Value               Description\n"
		for k in self.params:
			ret+= "%s%s%s%s\n" % (
				k.ljust(10),
				types[self.get_param_type(k)].ljust(20),
				(("%.3f"%self.params[k])+" "+units[self.get_param_type(k)]).ljust(20),
				self._meta[k][3]
				)

		return ret
