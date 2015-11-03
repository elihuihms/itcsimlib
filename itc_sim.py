from multiprocessing	import cpu_count,Queue

from itc_experiment		import *
from itc_calc			import ITCCalc
from thermo				import *
from utilities			import *

class ITCSim:
	"""
	Contains all the necessary data and methods to collect experiments
	"""

	def __init__(self,T0=298.15,units='J',verbose=False,threads=None):
		self.T0	= T0 # reference temperature
		assert units in thermo._units
		self.units = units
		self.size = 0 # number of experiments
		self.experiments = []
		self.chisq = {}
		self.verbose = verbose

		self.model = None
		self.in_Queue,self.out_Queue = Queue(),Queue()

		if threads==None or threads<1:
			threads = cpu_count()
		self.workers = [None] * threads

	def __str__(self):
		ret = "\nITCSim \"%s\"\n"%(self.__module__)
		ret+= "Current model information:\n"
		ret+= "\n".join( ["\t%s"%s for s in str(self.model).split("\n")] )
		ret+= "\nPer-model chisq values:\n"
		for E in self.get_experiments():
			if E.title in self.chisq:
				ret+= "\t%s	%f\n"%(E.title,self.chisq[E.title])
			else:
				ret+= "\t%s	N/A\n"%(E.title)
		return ret

	# getters
	def get_experiment(self, index=0):
		return self.experiments[index]

	def get_experiments(self):
		return self.experiments

	def get_experiment_by_title(self, title):
		for E in self.experiments:
			if E.title==title:
				return E
		return None

	def get_model(self):
		return self.model

	def get_model_param(self,name):
		return self.model.get_param(name)
	
	def get_model_params(self):
		return self.model.get_params()
		
	def get_chisq(self):
		return sum([self.chisq[t] for t in self.chisq])/self.size

	# setters
	def set_model(self, model):
		if not (None in self.workers):
			self.stop_workers()

		self.model = model
		self.model.set_units(self.units)
		for i in xrange(len(self.workers)):
			self.workers[i] = ITCCalc( self.T0, self.model, self.in_Queue, self.out_Queue )
			self.workers[i].start()

	def set_model_params(self, units=None, *args, **kwargs ):
		if units == None:
			units = self.units
		self.model.set_params( units=units, *args, **kwargs )
		
	def set_model_param(self, param, value, units=None):
		if units == None:
			units = self.units
		self.model.set_param( param, value, units=units)

	###
	def info(self):
		print str(self)

	def done(self):
		self.stop_workers()

	def stop_workers(self):
		# send term signal to workers
		for i in xrange(len(self.workers)):
			self.in_Queue.put( None )

		# make sure they're all shut down
		for i in xrange(len(self.workers)):
			if self.workers[i] != None:
				self.workers[i].join()
			self.workers[i] = None

	def add_experiment( self, experiment ):
		self.experiments.append( experiment )
		self.size +=1
	
	def add_experiment_file( self, file, **kwargs ):
		tmp,data = read_itcsimlib_exp(file)
		# overwrite any file-obtained info with explicit values
		info = tmp.copy()
		info.update(kwargs)
		if len(data) == 2:
			self.add_experiment( ITCExperiment(injections=data[0],dQ=data[1],units=self.units,**info) )
		elif len(data) == 3:
			self.add_experiment( ITCExperiment(injections=data[0],dQ=data[1],ddQ=data[2],units=self.units,**info) )
		return self.get_experiment(self.size-1)
	
	def add_experiment_synthetic( self, injections, *args, **kwargs ):
		self.add_experiment( ITCExperimentSynthetic(injections=injections,units=self.units,*args,**kwargs) )
		return self.experiments[self.size -1]
		
	def remove_experiment( self, experiment ):
		self.experiments.remove(experiment)
		if experiment.title in self.chisq:
			del self.chisq[experiment.title]
		self.size -=1		

	def make_plots(self,indices=None,**kwargs):
		"""
		Generate plots for all experimental datasets
		"""
		for (i,E) in enumerate(self.experiments):
			if(indices==None) or (i in indices):
				E.make_plot(**kwargs)

	def export_data(self,indices=None,**kwargs):
		for (i,E) in enumerate(self.experiments):
			if(indices==None) or (i in indices):
				E.export_data(**kwargs)
				
	def write_params(self, **kwargs):
		write_params_to_file( file, params=self.model.get_params(), units=self.units, **kwargs )
			
	def run( self, experiments=None, writeback=True ):
		"""
		Using the current parameters, generates fits for either the specified experiments or all experimental datasets, returns the normalized/average chi-square goodness of fit
		"""
		
		if experiments:
			for E in experiments:
				self.in_Queue.put( (self.model.get_params(),E) )
		else:
			for E in self.get_experiments():
				self.in_Queue.put( (self.model.get_params(),E) )

		queue_contents = []
		while len(queue_contents) < self.size:
			queue_contents.append( self.out_Queue.get(True) )

		for title,data in queue_contents:
			# in the case of an exception during model execution, title will be None
			if title == None:
				print "\nFatal error during model evalution: %s"%(data)
				self.done()
				return None
			else:
				self.chisq[title] = self.get_experiment_by_title(title).get_chisq(data,writeback)

		return self.get_chisq()


