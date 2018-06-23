import random
import copy
import math

try:
	_tmp = OrderedDict()
except:
	from ordered_dict	import OrderedDict

from itc_model import ITCModel

class GAReducer:
	def __init__(self, fit, nparams, popsize):
		self.fit = fit
		self.sim = fit.get_sim()
		self.model = GAModel(sim.model)
		self.start_params = self.model.parent_model.get_params()

		self.nparams = nparams
		self.popsize = popsize
		
		self.populationA = GAPopulation(self.model, self.nparams, self.popsize)
		self.populationB = GAPopulation(self.model, self.nparams, self.popsize)
		
	def start():
		self.sim.set_model( self.model )
		self.populationA.optimize()
	
	def done():
		# restore the original model and parameter values
		self.sim.set_model( self.model.parent_model )
		self.sim.set_model_params( self.start_params )

class GAModel(ITCModel):
	def __init__(self, model):
		ITCModel.__init__(self)
		self.parent_model = model

	def set_mapping(self, mapping):
		ITCModel.__init__(self)
		self.set_units(self.parent_model.units)
		self._component_meta = self.parent_model._component_meta
		
		self.map = {}
		for i,group in enumerate(mapping):
			type = self.parent_model.get_param_type(group[0])
			self.map[ "%s:%i"%(type,i)] = group[i][:]
			self.add_parameter("%s:%i"%(type,i),type)
			self._param_meta[name] = self.parent_model._param_meta[map[i][0]]

	def set_param(self, param, value):
		for p in self.map[param]:
			self.parent_model.set_param(p, value)
			
	def Q(self,T0,T,concentrations):
		return self.parent_model.Q(T0,T,concentrations)
		
class GAPopulation:
	
	def __init__(self, model, nparams, popsize):
		self.model = model
		self.nparams = nparams
		self.popsize = popsize

		self.member_param_map = [ [ [] for j in xrange(self.nparams)] for i in xrange(self.popsize) ]
		self.member_param_val = [ [0.0 for j in xrange(self.nparams)] for i in xrange(self.popsize) ]
		self.member_fitness = [ None for i in xrange(self.popsize) ]
		
		param_names = self.model.parent_model.get_param_names()
		param_types = [self.model.parent_model.get_param_type[p] for p in param_names]
					
		# generate a random assignment of possible parameters to each member of the population
		for member in xrange(self.popsize):
			# randomly assign one output parameter from each type to the initial input parameter array
			# many "output" (model) parameters can be assigned to one input parameter
			assigned_params = 0
			
			for type in self.params: # first make sure that each TYPE of parameter (dG, dH, etc.) is represented at least once in the input parameter map
				self.member_param_map[member][assigned_params].append( params[type].pop() )
				self.member_param_val[member][assigned_params] = self.model.get_param( self.member_param_map[member][assigned_params][0] ) # starting value is the currently assigned model parameter value
				assigned_params += 1
				
			while assigned_params < self.nparams: # randomly assign model parameters to the input parameter map
				type = random.choice(params.keys())
				if len(params[type]) == 0:
					del params[type]
					continue

				self.member_param_map[member][assigned_params].append( params[type].pop() )
				self.member_param_val[member][assigned_params] = self.model.get_param( self.member_param_map[member][assigned_params][0] )
				assigned_params += 1
			
			while len(params.keys()) > 0: # randomly assign any remaining model parameters to the correct type in the input parameter map
				type = random.choice(params.keys())
				if len(params[type]) == 0:
					del params[type]
					continue
					
				index = random.randrange(self.nparams)
				if type == self.model.get_param_type(self.member_param_map[member][index][0]):
					self.member_param_map[member][index].append( params[type].pop() )
	"""		for i in xrange(self.popsize):
			self.sim.set_graph( self.populationA.get_graph(i) )
			self.sim.set_model_params( **self.populationA.get_params(i) )
			opt,chisq = fit.optimize( params=self.sim.get_model_params().keys() )
"""

	def __getitem__(self, member):
		return self.member_fitness[member]
        
	def __setitem__(self, member, value):
		self.member_fitness[member] = value
	
	def copy(self, population):
		self.member_param_map = copy.deepcopy(population.member_param_map)
		self.member_param_val = copy.deepcopy(population.member_param_val)
		self.member_fitness = copy.deepcopy(population.member_fitness)

	def find_unassigned_params(self, member):
		ret = []
		for param in self.model.get_param_names():
			found = False
			for params in self.member_param_map[member]:
				if param in params:
					found = True
			if not found:
				ret.append(param)
		return ret
			
	def cross( self, member1=None, member2=None ):
		if member1==None:
			member1 = random.randrange(self.popsize)
		if member2==None:
			member2 = random.randrange(self.popsize)
		
		i,j = random.randrange(self.nparams),random.randrange(self.nparams)
		t = self.model.get_param_type(self.member_param_map[member1][i][0])
		while self.model.get_param_type(self.member_param_map[member2][j][0]) != t:
			j = random.randrange(self.nparams)
			
		self.member_param_map[member1][i],self.member_param_map[member1][j] = self.member_param_map[member1][j],self.member_param_map[member1][i]
		
		self.assign_missing_params(member1)
		self.assign_missing_params(member2)
		self.remove_duplicate_params(member1)
		self.remove_duplicate_params(member2)
		
	def mutate( self, member=None ):
		if member==None:
			member1 = random.randrange(self.popsize)
			
		i = random.randrange(self.nparams)
		if len(self.member_param_map[member][i]) > 1:
			p = random.choice(self.member_param_map[member][i])
			t = self.model.get_param_type(p)
			self.member_param_map[member][i].remove(p)
		else:
			return
			
		while True:
			i = random.randrange(self.nparams)
			if t == self.model.get_param_type(self.member_param_map[member][i][0]):
				self.member_param_map[member][i].append(p)
				return
				
	def get_unique_params( self, member ):
		ret = {}
		for index,map in enumerate(self.member_param_map[member]):
			ret[ self.member_param_map[member][index][0] ] = self.member_param_val[member][index]
		return ret

	def get_params( self, member ):
		ret = {}
		for index,map in enumerate(self.member_param_map[member]):
			for p in map:
				ret[p] = self.member_param_val[member][index]
		return ret

	def set_unique_params( self, member, params ):
		for index,map in enumerate(self.member_param_map[member]):
			for p in params:
				if self.member_param_map[member][index][0] == p:
					self.member_param_val[member][index] = params[p]

	def set_params( self, member, params ):
		for p in params:
			self.member_param_val[member] = params[:]
			
	def assign_missing_params(self, member):
		params = random.shuffle(self.model.get_param_names())
		for p in params: # see which model parameters are found in the member's map
			indices = random.shuffle(xrange(self.nparams))
			for i in xrange(self.nparams):
				index = indices.pop()
				if p in self.member_param_map[member][index]:
					params.remove(p)
		for p in params: # add the missing parameters to a member's map that contains parameters of the right type
			indices = random.shuffle(xrange(self.nparams))
			for i in xrange(self.nparams):
				index = indices.pop()
				if model.get_param_type(p) == model.get_param_type(self.model_get_param_map[member][index][0]):
					self.model_get_param_map[member][index].append(p)
		return len(params)

	def remove_duplicate_params(self, member):
		indices = random.shuffle(xrange(self.nparams))
		duplicates,seen_params = 0,[]
		for i in xrange(self.nparams):
			index = indices.pop()
			for p in seen_params:
				if p in self.member_param_map[member][index]:
					self.member_param_map[member][index].remove(p)
					duplicates += 1
			seen_params.extend(self.member_param_map[member][index])
		return duplicates
					
				