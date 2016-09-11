from random	import choice,random
from scipy 	import optimize,mean,std
from thermo	import *

from utilities		import *
try:
	_tmp = OrderedDict()
except:
	from ordered_dict	import OrderedDict

class ITCFit:
	"""A class for optimizing model parameters to accurately fit experimental data.

	Attributes:
		sim (ITCSim): The ITCSim object containing experimental data to replicate.
		model (ITCModel): The model to use.
		bounds (dict of tuples): A parameter-name keyed dict of low and high bounds to enforce during fitting.
		method (string): The optimization algorithm to use.
		method_args (dict): Arguments to pass to the optimization algorithm.
		verbose (boolean): Print additional information to the console?
		chisq (float): The most recently evaluated goodness-of-fit chisquare.
	"""

	def __init__(self, sim, method='simplex', method_args={}, verbose=False):
		"""Constructor function for the ITCFit object.

		Arguments:
			sim (ITCSim): The ITCSim object that possesses the experimental data to replicate.
			bounds (dict of tuples): A parameter-name keyed dict of low and high bounds to enforce during fitting.
			method (string): The optimization algorithm to use.
			method_args (dict): Arguments to pass to the optimization algorithm.
			verbose (boolean): Print additional information to the console?
		"""

		self.sim = sim
		self.model	= self.sim.get_model()
		self.method	= method
		self.method_args = method_args
		self.verbose = verbose
		self.chisq = 0.0

		# obtain model-defined boundaries to enforce during fitting
		self.bounds = dict( (name,self.model.get_param_bounds(name)) for name in self.model.get_param_names() )
		
	def set_sim(self, sim):
		self.sim = sim
		self.model = self.sim.get_model()
		
	def get_sim(self):
		return self.sim
	
	def add_bounds(self, param, low=None, high=None):
		"""Add low and/or high boundaries for a specified parameter.

		Args:
			param (string): The name of the parameter to create the boundary for
			low (float): The lower boundary for the parameter, or None for no boundary
			high (float): The upper boundary for the parameter, or None for no boundary

		Returns:
			None
		"""

		assert param in self.model.get_param_names()
		self.bounds[param] = (low,high)

	def optimize(self, params=[], callback=None, update_fits=False ):
		"""Optimize the specified parameters.

		Args:
			params (dict of strings): The names of the model parameters to optimize.
			callback (function): A callback function to call at each optimization step, will be passed a vector of the current parameter values.

		Returns:
			(dict,float): A parameter-name keyed dict of optimized values, and the corresponding goodness-of-fit.
		"""

		if self.verbose:
			print "\nTask: Optimizing %s parameters using the %s algorithm\n"%(",".join(params),self.method)
			def _printer(x):
				print "%s (%f)" %(" ".join(map(str,x)),self.sim.chisq)
		else:
			def _printer(x):
				return

		# starting parameter values to restore later
		start_params = self.sim.get_model_params().copy()
		
		# initial param guesses as list
		x0 = [start_params[p] for p in params]
		
		# the target objective function to minimize
		def _target(x,sim):
			for i,p in enumerate(params):
				sim.set_model_param(p, x[i])

			# this is hackish - if we violate a boundary, return a scaled well function based on the last iteration point
			m = self._check_bounds()
			if m > 0:
				return sim.get_chisq() * (1+m)

			return sim.run(writeback=False)
	
		# optimize parameters
		opt = self._fitter( _target, x0, callback )

		ret = OrderedDict( (p,opt[0][i]) for i,p in enumerate(params) )
		self.sim.set_model_params(**ret)

		if update_fits:
			self.sim.run()
		else: # restore initial parameter values
			self.sim.set_model_params(**start_params)
	
		# return the optimized parameters and the chisquare value
		return ret,opt[1]

	def estimate(self, params=[], bootstraps=1000, randomize=0.1, callback=None, logfile=None ):
		"""Generate confidence intervals for optimized parameters.

		Args:
			params (dict of strings): The names of the model parameters to optimize.
			bootstraps (int): The number of synthetic datasets to refit.
			randomize (float): A fraction by which to perturb the starting parameters before optimization.
			callback (function): A callback function to call at each optimization step, will be passed a vector of the current parameter values.
			logfile (string): A file path to write the optimized parameter values for each bootstrap dataset

		Returns:
			(dict of tuples): A parameter-name keyed dict of tuples consisting of the mean and standar deviation of the optimized values
		"""

		if self.verbose:
			print "\nTask: Estimating uncertainty intervals for %s using %i bootstraps\n"%(",".join(params),bootstraps)

		# initialize to make sure we have fit data to use
		self.sim.run()

		# starting parameter values to restore later
		start_params = self.sim.get_model_params().copy()

		# store the original experimental and fit data for restoration later
		dQ_exp,dQ_fit = [],[]
		for E in self.sim.get_experiments():
			assert len(E.dQ_exp) == len(E.dQ_fit)
			dQ_exp.append( E.dQ_exp[:] )
			dQ_fit.append( E.dQ_fit[:] )

		# generates a synthetic dataset from the fit and fit residuals
		def _make_bootstrap(n,exp,fit):
			res = [ exp[i]-fit[i] for i in xrange(n) ]
			return [ f + choice(res) for f in exp ]

		param_values = OrderedDict( (p,[]) for p in params )
		for i in xrange(bootstraps):
			if self.verbose:
				print "Note: Bootstrap %i"%(i)

			# randomize starting point by user-specifiable amount
			for p in params:
				self.sim.set_model_param(p,self.sim.get_model_param(p) * (1+(2*random()-0.5)*randomize))

			# replace the experimental data points with the synthetic data
			for j,E in enumerate(self.sim.get_experiments()):
				E.dQ_exp = _make_bootstrap(E.npoints,dQ_exp[j],dQ_fit[j])

			# re-optimize the parameters using the new, synthetic datasets
			optimized,chisq = self.optimize(params)

			# append the optimized parameter values to our growing array
			for p in optimized:
				param_values[p].append( optimized[p] )

			if callback != None:
				callback(optimized)

			if logfile != None:
				write_params_to_file(logfile,optimized,header=(i==0),post="%.3f"%(chisq))

			# restore initial parameter values for the next iteration
			self.sim.set_model_params(**start_params)

		# restore original experimental data and fit for the experiments
		for i,E in enumerate(self.sim.get_experiments()):
			E.dQ_exp = dQ_exp[i]
			E.dQ_fit = dQ_fit[i]

		# return the mean and standard deviation of the model parameters
		return OrderedDict( (p,(mean(param_values[p]),std(param_values[p]))) for p in params )

	def _check_bounds(self):
		ret = 0
		for k,v in self.sim.get_model_params().iteritems():
			if self.bounds[k][0] != None and v < self.bounds[k][0]:
				ret += (self.bounds[k][0] - v)
				if self.verbose:
					print "Note: Boundary violation for \"%s\" (%f<%f)"%(k,v,self.bounds[k][0])

			elif self.bounds[k][1] != None and v > self.bounds[k][1]:
				ret += (v - self.bounds[k][1])
				if self.verbose:
					print "Note: Boundary violation for \"%s\" (%f>%f)"%(k,v,self.bounds[k][1])

		return ret

	def _fitter(self, func, x0, callback=None):
		# wrapper function for a variety of optimization algorithms
		# note that some of these aren't fully integrated yet
		if self.method == 'simplex':
			ret = optimize.fmin(
				func=func,
				x0=x0,
				args=(self.sim,),
				full_output=True,
				callback=callback,
				**self.method_args
			)
		elif self.method == 'basinhopping':
			kwargs = {'method':'Powell','args':(self.sim,)}
			opt = optimize.basinhopping(
				func=func,
				x0=x0,
				args=(self.sim,),
				minimizer_kwargs=kwargs,
				disp=self.verbose,
				**self.method_args
			)
			ret = opt.x,opt.fun
		elif self.method == 'powell':
			ret = optimize.fmin_powell(
				func=func,
				x0=x0,
				args=(self.sim,),
				disp=self.verbose,
				full_output=True,
				callback=callback,
				**self.method_args
			)
		elif self.method == 'tnc':
			opt = optimize.fmin_tnc(
				func=func,
				x0=x0,
				args=(self.sim,),
				approx_grad=True,
				disp=self.verbose,
				**self.method_args)[0]
			ret = opt,func(x0,*args)
		elif self.method == 'bfgs':
			ret = optimize.fmin_l_bfgs_b(
				func=func,
				x0=x0,
				args=(self.sim,),
				approx_grad=True,
				maxfun=fit.maxfun,
				disp=self.verbose,
				**self.method_args)
		else:
			raise Exception('Unrecognized fitting algorithm')

		return ret