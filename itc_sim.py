import os
from multiprocessing	import cpu_count,Queue
from scipy				import array,dtype,genfromtxt

from itc_experiment		import ITCExperiment

class ITCSim:
	"""
	Contains all the necessary data and methods to collect experiments
	"""

	def __init__(self,model,T_ref=298.15,verbose=False,threads=None):
		self.model	= model
		self.T_ref	= T_ref # reference temperature
		self.dG		= []
		self.dH		= []
		self.dCp	= []
		self.size	= 0 # number of experiments
		self.experiments = []
		self.verbose = verbose
		self.chisq	= 0.0

		self.in_Queue,self.out_Queue = Queue(),Queue()

		if threads==None or threads<1:
			threads = cpu_count()

		self.workers = [None] * threads
		for i in xrange(threads):
			self.workers[i] = self.model.get_worker( self.in_Queue, self.out_Queue )
			self.workers[i].start()

	def __del__(self):

		# send term signal to workers
		for i in xrange(len(self.workers)):
			self.in_Queue.put( None )

		# make sure they're all shut down
		for i in xrange(len(self.workers)):
			if self.workers[i] != None:
				self.workers[i].join()
			self.workers[i] = None

	def __str__(self):
		ret = "\nITCSim \"%s\"\n"%(self.__module__)
		ret+= "Per-model chisq values:\n"
		for E in self.get_experiments():
			ret+= "%s	%f\n"%(E.title,E.chisq)
		return ret

	def get_experiment(self, index=0):
		return self.experiments[index]

	def get_experiments(self):
		return self.experiments

	def get_experiment_by_title(self, title):
		for E in self.experiments:
			if E.title==title:
				return E
		return None

	def get_experiments_by_temperature(self, temperature):
		ret = []
		for E in self.experiments:
			if E.T == temperature:
				ret.append(E)
		return ret

	def set_params(self, dG=None, dH=None, dCp=None ):
		"""
		Set the various initial parameters for the specified model
		"""

		if len(self.dG) > 0 and dG != None:
			assert( len(self.dG) == len(dG) )
		if len(self.dG) > 0 and dH != None:
			assert( len(self.dG) == len(dH) )
		if len(self.dG) > 0 and dCp != None:
			assert( len(self.dG) == len(dCp) )

		# note, dG & dH values are always at ref temp!
		if dG != None:
			self.dG		=	dG
		if dH != None:
			self.dH		=	dH
		if dCp != None:
			self.dCp	=	dCp

	def add_experiment( self, title, T, V0, M0, L0, DQ, I_vol, reverse=False, skip=[], dil_Q=0.0 ):
		"""
		Add ITC data, along with the necessary experimental and instrument parameters necessary to generate fits
		Calculates the active concentration of protein and ligand at each point
		"""

		for E in self.experiments:
			if E.title == title:
				print "Experiment \"%s\" already present in this simulation" % (title)

		self.experiments.append(
			ITCExperiment(
				title	= title,
				T		= T,
				V0		= V0,
				M0		= M0,
				L0		= L0,
				dQ_exp	= array(DQ,dtype('d')),
				I_vol	= array(I_vol,dtype('d')),
				reverse	= reverse,
				skip	= skip,
				dil_Q	= dil_Q
			)
		)

		self.size +=1

	def load_file( self, path, T, V0, M0, L0, reverse=False, skip=[], dil_Q=0.0, title=None):
		"""
		Read experimental ITC data from a file containing two columns: the observed deltaH and the injection volumes
		"""

		(DQ,I_vol) = genfromtxt( path, unpack=True, usecols=(0,1) )
		if title==None:
			title = os.path.splitext(os.path.basename(path))[0]
		self.add_experiment( title, T, V0, M0, L0, DQ, I_vol, reverse, skip, dil_Q)

	def make_plots(self,indices=None,hardcopy=False,hardcopydir='.',hardcopyprefix='',hardcopytype='png'):
		"""
		Generate plots for all experimental datasets
		"""
		for (i,E) in enumerate(self.experiments):
			if(indices==None) or (i in indices):
				E.show_plot(hardcopy,hardcopydir,hardcopyprefix,hardcopytype)

	def export_data(self,dir='.',indices=None,prefix='export_'):
		for (i,E) in enumerate(self.experiments):
			if(indices==None) or (i in indices):
				E.export_data(prefix+E.title+".txt")

	def make_fits( self, writeback=True ):
		"""
		Using the current model parameters, generates fits for all experimental datasets, returns the normalized chi-square goodness of fit
		"""

		for E in self.get_experiments():
			self.in_Queue.put( (self.model.get_format( E.T, self.T_ref ),E) )

		queue_contents = []
		while len(queue_contents) < self.size:
			queue_contents.append( self.out_Queue.get(True) )

		self.chisq = 0.0
		for title,data in queue_contents:

			# in the case of an exception during model execution, title will be None
			if title == None:
				raise data
			else:
				self.chisq += self.get_experiment_by_title(title).calc_chisq(data,writeback)/self.size

		return self.chisq


