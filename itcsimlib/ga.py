import random
import copy
import math
import shelve

try:
	_tmp = OrderedDict()
except:
	from ordered_dict	import OrderedDict

from itc_model import ITCModel

class GAModel(ITCModel):
	def __init__(self, model, mapping):
		ITCModel.__init__(self)
		self.parent_model = model
		self.set_units(self.parent_model.units)
		self._component_meta = self.parent_model._component_meta
		
		self.map = {}
		for type,maps in mapping.iteritems():
			for i,map in enumerate(maps):
				self.map[ "%s:%i"%(type,i) ] = map[0][:]
				self.add_parameter("%s:%i"%(type,i),type,description=" ".join(map[0]))
				self.set_param("%s:%i"%(type,i), map[1])

	def set_param(self, param, value ):
		ITCModel.set_param(self, param, value)
		for p in self.map[param]:
			self.parent_model.set_param(p, value)
			
	def Q(self,T0,T,concentrations):
		try:
			return self.parent_model.Q(T0,T,concentrations)
		except:
			return [0.0]*len(concentrations)
		
class GAPopulation:
	
	def __init__(self, fit, popsize, nparams, params=None, chiTol=1.0, mutFreq=0.1, mutSwap=0.5, crossFreq=0.1, scratchFile=None):
		self.fit = fit
		self.sim = fit.get_sim()
		self.model = self.sim.get_model()
		self.model.precision=1E-20
		self.model_param_values = self.sim.get_model_params()
		self.popsize = popsize
		self.nparams = nparams
		if 	params == None:
			self.variable_params = model.get_param_names()
		else:
			self.variable_params = params
		self.chiTolerance = chiTol
		self.mutationFrequency = mutFreq
		self.mutationSwapPercent = mutSwap
		self.crossFrequency = crossFreq
		self.scratchFile=scratchFile

		self.member_param_map = [ {} for i in xrange(self.popsize*2) ]	
		self.member_chisquare = [ None for i in xrange(self.popsize*2) ]
		self.member_optimized = [ False for i in xrange(self.popsize) ]
		
		if self.scratchFile != None:
			self.scratch = shelve.open( self.scratchFile )
			self.scratch['size'] = self.popsize
		else:
			self.scratch = None
		
		self.randomize()
		
	def __len__(self):
		return self.popsize
	
	def __getitem__(self, index):
		return (self.member_param_map[index],self.member_chisquare[index])
        
	def __setitem__(self, index, value):
		self.member_param_map[index],self.member_chisquare[index] = value
		
	def __del__(self):
		if self.scratch != None:
			self.scratch.close()
			
	def get_model(self):
		return self.model
						
	def evolve( self ):
		for i in xrange(self.popsize):
			
			if random.random() < self.mutationFrequency:
				self.mutate( i )
				self.member_optimized[i] = False			
			#if random.random() < self.crossFrequency:
			#	j = random.randrange(self.popsize)
			#	self.cross( i, j )
			#	self.member_optimized[i] = False
			#	self.member_optimized[j] = False

			if not self.member_optimized[i]:		
				model = GAModel(self.model,self.member_param_map[i])
				self.sim.set_model( model )
				self.fit.set_sim( self.sim )
			
				opt,self.member_chisquare[i] = self.fit.optimize( params=model.get_params() )
				for p,val in opt.iteritems():
					type,number = p.split(":")
					self.member_param_map[i][type][int(number)][1] = val
				self.member_optimized[i] = True
			
			if self.scratch != None:
				self.scratch['pop_%i_map'%i] = self.member_param_map[i]
				self.scratch['pop_%i_chi'%i] = self.member_chisquare[i]
				self.scratch['pop_%i_opt'%i] = self.member_optimized[i]
									
		if not None in self.member_chisquare: # sorting not performed on first generation
			self.sort()
		
		for i in xrange(self.popsize):
			self.member_param_map[i+self.popsize] = self.member_param_map[i]
			self.member_chisquare[i+self.popsize] = self.member_chisquare[i]
			
	def sort( self ):
		tmp = zip(self.member_param_map,self.member_chisquare)
		tmp.sort(key=lambda x: float(x[1])*self.chiTolerance*random.random() ) # fuzzy sort within chi-tolerance value to retard monoculturing
		self.member_param_map = [x[0] for x in tmp]
		self.member_chisquare = [x[1] for x in tmp]
		
	def randomize( self ):
		# generate a dict of all base model parameters keyed by type
		parameters = {}
		for p in self.variable_params:
			type = self.model.get_param_type(p)
			if type in parameters.keys():
				parameters[type].append(p)
			else:
				parameters[type] = [p]
		
		# number of parameters in the base model
		n_base_params = len(self.variable_params)
		
		assert self.nparams < n_base_params
		assert self.nparams > len(parameters)

		# randomly generate parameter trees (mappings) for each population member
		for i in xrange(self.popsize):
			self.member_chisquare[i] = random.random()

			# parameter mappings are also keyed by parameter type for convenience
			self.member_param_map[i] = dict.fromkeys( parameters.keys() )
			unassigned_params = copy.deepcopy(parameters)

			# ensure that at each parameter type is represented at least once in the maps
			j = n_base_params
			for type in unassigned_params:
				self.member_param_map[i][type] = []
				self.member_param_map[i][type].append( [[],0.0] )
				self.member_param_map[i][type][-1][0].append( random.choice(unassigned_params[type]) )
				self.member_param_map[i][type][-1][1] = self.model.get_param(self.member_param_map[i][type][-1][0][0],self.sim.get_units())
				unassigned_params[type].remove(self.member_param_map[i][type][-1][0][0])
				j -= 1
				
				if unassigned_params[type] == []:
					del unassigned_params[type]
					
			# generate the appropriate number of maps for the model
			for k in xrange(self.nparams - (n_base_params-j)):
				type = random.choice(unassigned_params.keys())
				self.member_param_map[i][type].append( [[],0.0] )
				self.member_param_map[i][type][-1][0].append( random.choice(unassigned_params[type]) )
				self.member_param_map[i][type][-1][1] = self.model.get_param(self.member_param_map[i][type][-1][0][0],self.sim.get_units())
				unassigned_params[type].remove(self.member_param_map[i][type][-1][0][0])
				j -= 1
				
				if unassigned_params[type] == []:
					del unassigned_params[type]

			# assign the remaining base model parameters to the existing map groups
			while j > 0:
				for type in random.sample(unassigned_params.keys(),len(unassigned_params.keys())):
					k = random.randrange( len(self.member_param_map[i][type]) )
					self.member_param_map[i][type][k][0].append( random.choice(unassigned_params[type]) )
					self.member_param_map[i][type][k][1]+= self.model.get_param(self.member_param_map[i][type][k][0][-1],self.sim.get_units())
					unassigned_params[type].remove(self.member_param_map[i][type][k][0][-1])
					j -= 1
					
					if unassigned_params[type] == []:
						del unassigned_params[type]
			
			# map parameter values are the average of the parent model parameter values in the map
			for type in self.member_param_map[i]:
				for map in self.member_param_map[i][type]:
					map[1] = map[1]/len(map[0])
	
	def split_map(self, i, iteration=0):
		"""
		Splits a existing, randomly-selected map into two maps.
		This is used to return to a correct number of maps after removal of empty maps.
		"""
		
		if iteration > 100:
			raise Exception("ga:split_map() iteration exceeded 100.\nBusted map:\n%s\n"%(str(self.member_param_map[i])))
		
		# select the map to split
		type = random.choice(self.member_param_map[i].keys())		
		n = len(self.member_param_map[i][type])
		if n == 0: # try again
			self.split_map( i, iteration+1 )

		j = random.randrange(n)
		n = len(self.member_param_map[i][type][j][0])
		if n == 1: # try again
			self.split_map( i, iteration+1 )

		self.member_param_map[i][type].append( [[],self.member_param_map[i][type][j][1]] ) # create the new map and copy the map parameter value

		for k in xrange(int(math.ceil(n/2))): # populate it with randomly-selected elements from the old map
			self.member_param_map[i][type][-1][0].append( random.choice(self.member_param_map[i][type][j][0]) )
			self.member_param_map[i][type][ j][0].remove( self.member_param_map[i][type][-1][0][-1] )
		
	def mutate( self, i ):
		type = random.choice(self.member_param_map[i].keys())
		n = len(self.member_param_map[i][type])
		if n == 1:
			return

		mapA,mapB = random.sample(xrange(n),2)
		elA = random.randrange(len(self.member_param_map[i][type][mapA][0]))
		
		if random.random() < self.mutationSwapPercent:
			elB = random.randrange(len(self.member_param_map[i][type][mapB][0]))
			self.member_param_map[i][type][mapA][0][elA],	\
			self.member_param_map[i][type][mapB][0][elB] =	\
			self.member_param_map[i][type][mapB][0][elB],	\
			self.member_param_map[i][type][mapA][0][elA]
		else:
			self.member_param_map[i][type][mapB][0].append( self.member_param_map[i][type][mapA][0][elA] )
			del self.member_param_map[i][type][mapA][0][elA]
			
			# did we generate an empty map?
			if len(self.member_param_map[i][type][mapA][0]) == 0:
				del self.member_param_map[i][type][mapA]
				self.split_map(i)
			
	def cross( self, iA, iB ):
		type = random.choice(self.member_param_map[iA].keys())
		mapA = random.randrange(len(self.member_param_map[iA][type]))
		mapB = random.randrange(len(self.member_param_map[iB][type]))
		
		tmpA = self.member_param_map[iA][type][mapA][0][:]
		tmpB = self.member_param_map[iB][type][mapB][0][:]
		
		for p in tmpA:
			for map in self.member_param_map[iB][type]:
				if p in map[0]: map[0].remove(p)
		for p in tmpB:
			for map in self.member_param_map[iA][type]:
				if p in map[0]: map[0].remove(p)
		
		self.member_param_map[iA][type][mapA][0] = tmpB
		self.member_param_map[iB][type][mapB][0] = tmpA
		
		# if any maps ended up as zombies, replace them
		for map in self.member_param_map[iA][type]:
			if len(map[0]) == 0:
				self.member_param_map[iA][type].remove(map)
				self.split_map(iA)
		for map in self.member_param_map[iB][type]:
			if len(map[0]) == 0:
				self.member_param_map[iB][type].remove(map)
				self.split_map(iB)