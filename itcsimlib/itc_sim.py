"""Essential simulator for generating or fitting data with binding models.


"""

import os
import copy
import multiprocessing

from .					import __version__
from .itc_experiment	import *
from .itc_calc			import ITCCalc
from .thermo			import _UNITS
from .utilities			import *


class ITCSim:
	"""This core class evaluates binding models against experimental data.
	
	Attributes
	----------
	chisq : dict of floats
		The average reduced chi squared goodness-of-fit across the experiments in the simulator.
	experiments : list of ITCExperiments
		The list of experiments in the simulation. Do not directly modify this list, use the add and remove experiment class methods.
	model : ITCModel
		The model used by the simulator to generate/fit data.
	"""

	def __init__(self,T0=298.15,units='J',verbose=False,threads=0):
		"""Constructor for the ITCSim class.
		
		Arguments
		---------
		T0 : float
			The reference temperature for the simulation.
		units : string
			The units to use when reporting parameters to the user.
		verbose : boolean
			Write extra information to stdout?
		threads : int
			Number of threads to use when simulating ITC data. Default (0) disables multiprocessing, None uses all available cores.
		"""
		
		self.T0	= T0 # reference temperature
		assert units in _UNITS
		self.units = units
		self.size = 0 # number of experiments
		self.experiments = []
		self.chisq = {}
		self.verbose = verbose

		self.model = None
		self.in_Queue,self.out_Queue = multiprocessing.Queue(),multiprocessing.Queue()

		# Enable/diable multithreading, avoids __name__ guards on Windows
		if threads == 1:
			threads = 0
		elif threads == None:
			threads = multiprocessing.cpu_count()
		self.workers = [None] * threads

	def __str__(self):
		"""Stringify the simulator to be suitable for display to the user.
		
		Parameters
		----------
		None
		
		Returns
		-------
		string
			String containing the current state of the simulator.
		"""

		import sys
		from os.path			import abspath, getmtime
		from datetime			import date

		ret = "################################################################################\n"
		ret+= "{:^80}".format("itcsimlib v.%s"%(__version__))
		ret+= "\n\nExecution Date: %s\n"%(date.today().ctime())
		ret+= "Script: %s\n"%(abspath(sys.argv[0]))
		ret+= "Modification Date: %s\n"%(date.fromtimestamp(getmtime(abspath(sys.argv[0]))).ctime())
		ret+= "################################################################################\n"
		ret+= "Simulation reference temperature: %fK\n"%(self.T0)
		ret+= "Simulation units: %s\n"%(self.units)
		if len(self.chisq)>0:
			ret+= "Simulation chisq: %f\n"%(self.get_chisq())
		ret+= "Simulation model:\n"
		ret+= "\n".join( ["\t%s"%s for s in str(self.model).split("\n")] )
		ret+= "\n################################################################################\n"
		for i,E in enumerate(self.experiments):
			ret+= "Experiment %i:"%(i)
			ret+= "\n".join( ["\t%s"%s for s in str(E).split("\n")] )
			ret+= "\n"
		return ret

	# getters
	def get_experiment_by_title(self, title):
		"""Return the named experiment from simulation.
		
		Arguments
		---------
		title : string
			The title (name) of the experiment to return.
		
		Returns
		-------
		ITCExperiment
			The requested experiment
		"""
		for E in self.experiments:
			if E.title==title:
				return E
		raise KeyError("Experiment \"%s\" not found in simulation."%(title))

	def get_model_param(self,name,units=None):
		"""Return the value of the named simulation model parameter.
		
		Arguments
		---------
		name : string
			The name of the model parameter.
		units : string
			The units to use (if different than the specified model units to use, if applicable).
		
		Returns
		-------
		float
			The requested model parameter value.
		"""
		if units == None:
			units = self.units
		return self.model.get_param(name,units)
	
	def get_model_params(self,units=None):
		"""Return all the parameter values of the current simulation model.
		
		Arguments
		---------
		units : string
			The units to use (if different than the specified model units to use, if applicable).
		
		Returns
		-------
		dict of floats
			The requested model parameter values.
		"""
		if units == None:
			units = self.units
		return self.model.get_params(units)
		
	def get_chisq(self):
		"""Return the total reduced chisquared goodness-of-fit for the simulator.
		
		Arguments
		---------
		None
		
		Returns
		-------
		float
			The (total chisq / # experiments) between the experimental data and the model, if they exist.
		"""
		return sum([self.chisq[E.title] for E in self.experiments]) / self.size

	# setters
	def set_model(self, model):
		"""Set the model for the simulator to use.
		
		Arguments
		---------
		model : ITCModel
			The model to use for the simulator.
		
		Returns
		-------
		None		
		"""
		if not (None in self.workers):
			self.done()

		self.model = model
		self.model.set_units(self.units)

		for i in range(len(self.workers)):
			self.workers[i] = ITCCalc( self.T0, self.model, self.in_Queue, self.out_Queue )
			self.workers[i].start()

	def set_model_params(self, *args, **kwargs ):
		"""Passthrough for the simulator's model set_params()"""
		self.model.set_params( *args, **kwargs )
		
	def set_model_param(self, param, value):
		"""Passthrough for the simulator's model set_param()"""
		self.model.set_param( param, value )
		
	def done(self):
		"""Cleanly shuts down the simulator.
		
		Arguments
		---------
		None
		
		Returns
		-------
		None
		"""

		# send term signal to workers
		for i in range(len(self.workers)):
			self.in_Queue.put( None )

		# make sure they're all shut down
		for i in range(len(self.workers)):
			if self.workers[i] != None:
				self.workers[i].join()
			self.workers[i] = None

	def add_experiment( self, experiment ):
		"""Add a pre-defined experiment to the simulator.
		
		Arguments
		---------
		experiment : ITCExperiment
			The experiment to add.
			
		Returns
		-------
		None
		"""
		self.experiments.append( experiment )
		self.size +=1
		
	def add_experiment_synthetic( self, *args, **kwargs ):
		"""Add a set of synthetic experimental conditions to the simulator.
		
		Arguments
		---------
		*args
			Positional arguments for ITCExperimentSynthetic constructor
		**kwargs
			Keyword arguments for ITCExperimentSynthetic constructor
		
		Returns
		-------
		ITCExperiment
			The experiment created and added to the simulator.
		"""
		self.add_experiment( ITCExperimentSynthetic(units=self.units,*args,**kwargs) )
		return self.experiments[-1]

	def add_experiment_file( self, path, **kwargs ):
		"""Add an experiment present in the provided path to the simulator.
		
		Arguments
		---------
		path : string
			Path to the experiment file.
		**kwargs
			Keyword arguments to the ITCExperiment constructor, can be used to overwrite any file-defined arguments.
		
		Returns
		-------
		ITCExperiment
			The experiment created and added to the simulator.
		"""

		fname,ext = os.path.splitext(path)

		if os.path.isdir( path ):
			experiment = read_nitpic_exp( path, exp_args=kwargs )
		elif ext == ".itcpkl":
			experiment = read_itcsimlib_pkl( path )
		elif ext == ".nitpkl":
			experiment = read_nitpikl_exp( path, exp_args=kwargs )
		elif ext == ".dh" or ext = ".DH":
			experiment = read_origin_exp( path, exp_args=kwargs )
		else:
			experiment = read_itcsimlib_exp( path, exp_args=kwargs )

		if experiment != None:
			experiment.units = self.units
			self.add_experiment( experiment )

		return self.experiments[-1]
			
	def remove_experiment( self, experiment ):
		"""Remove an experiment from the simulator.

		Arguments
		---------
		experiment : ITCExperiment
			The experiment to remove.
		
		Returns
		-------
		None
		"""
		self.experiments.remove(experiment)
		if experiment.title in self.chisq:
			del self.chisq[experiment.title]
		self.size -=1
		
	def remove_all_experiments( self ):
		"""Removes all experiments from the simulator.
		
		Arguments
		---------
		None
		
		Returns
		-------
		None
		"""
		while len(self.experiments) > 0:
			self.remove_experiment( self.experiments[0] )

	def make_plots(self,indices=None,**kwargs):
		"""Call the make_plot() methods of simulator experiments.
		
		Arguments
		---------
		indices : list of ints
			If provided, generate only plots for the specified experimental indices. Otherwise, generate plots for all experiments.
		**kwargs
			Keyword arguments for the ITCExperiment make_plot() method.
		
		Returns
		-------
		None
		"""
		for (i,E) in enumerate(self.experiments):
			if(indices==None) or (i in indices):
				E.make_plot(**kwargs)

	def export_data(self,indices=None,**kwargs):
		"""Call the export_data() methods of simulator experiments.
		
		Arguments
		---------
		indices : list of ints
			If provided, generate exported files for the specified experimental indices. Otherwise, generate exports for all experiments.
		**kwargs
			Keyword arguments for the ITCExperiment export_data() method.
		
		Returns
		-------
		None
		"""
		for (i,E) in enumerate(self.experiments):
			if(indices==None) or (i in indices):
				E.export_data(**kwargs)
				
	def run( self, experiments=None, writeback=True ):
		"""Using the current model parameters, generate fits for either the specified experiments, and return the average reduced chi-squared goodness-of-fit.
		
		Arguments
		---------
		experiments : list of ITCExperiments
			The experiments to run through the simulator. If None, run all experiments in the simulator.
		writeback : boolean
			Update the dQ_fits of the experiments in the simulator?
		
		Returns
		-------
		float
			The average reduced chi-squared goodness-of-fit across the experiments.	
		"""
		if experiments == None:
			experiments = self.experiments

		if len(experiments) == 0:
			print("itc_sim: No experiments to simulate.")
			return None

		# without multiprocessing (avoids requirement for __name__ guards in Windows)
		if len(self.workers) == 0:
			self.model.start()

			for E in experiments:
				data = self.model.Q( self.T0, E.T, E.Concentrations )
				self.chisq[E.title] = E.get_chisq(data,writeback)

			self.model.stop()

		# with multiprocessing
		else:
			for E in experiments:
				self.in_Queue.put( (self.model.get_params(units=self.units),E) )

			queue_contents = []
			while len(queue_contents) < self.size:
				queue_contents.append( self.out_Queue.get(True) )

			for title,data in queue_contents:
				# in the case of an exception during model execution, title will be None
				if title == None:
					print("\nitc_sim: Fatal error during model evalution: %s"%(data))
					self.done()
					return None
				else:
					self.chisq[title] = self.get_experiment_by_title(title).get_chisq(data,writeback)

		return self.get_chisq()

