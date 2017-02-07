import random
import scipy
import scipy.optimize

from thermo		import *
from utilities	import *

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
		
	def _noisy_bisect(self, f, a, b, fa, fb, tolerance):
		# From http://www.johndcook.com/blog/2012/06/14/root-finding-with-noisy-functions/
		# Assume a < b, fa = f(a) < 0, fb = f(b) > 0.
		if b-a < tolerance:
			return (a, b)
		mid = 0.5*(a+b)
		fmid = f(mid)
		if fmid < fa or fmid > fb: # Monotonicity violated. Reached resolution of noise.
			return (a, b)
		if fmid < 0:
			a, fa = mid, fmid
		else:
			b, fb = mid, fmid
		return self._noisy_bisect(f, a, b, fa, fb, tolerance)
		
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
			print "\nitc_fit: Optimizing %s parameters using the %s algorithm\n"%(",".join(params),self.method)
			def _printer(x):
				print "itc_fit: %s (%f)" %(" ".join(map(str,x)),self.sim.chisq)
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
		
	def estimate(self, params, method='bootstrap', *args, **kwargs ):
		assert method in ['sigma','bootstrap']
	
		if method == 'sigma':
			return self.estimate_sigma( params, *args, **kwargs )
		elif method == 'bootstrap':
			return self.estimate_bootstrap( params, *args, **kwargs )
	
	def estimate_sigma(self, params=[], opt_params=None, sigma=None, stdevs=1, estimate=0.1, method='bisect', tolerance=0.001 ):
		"""Generate confidence intervals for optimized parameters by finding reduced chi-square +n standard deviations
		
		Notes:
			Requires that the simulation model parameters already be at a (global) minima.
			Using the secant (brent's method) is about four times slower than using a noise-tolerant bisection method.
		
		Args:
			params (list of strings): The names of the model parameters to find intervals for.
			opt_params (dict of strings): The names of the model parameters to optimize. Setting to None ensures that all available model parameters will be optimized.
			stdevs (int) : The number of standard deviations above 
			method (string, "bisect" or "secant") : The method used to find the root of the target function.
			tolerance (float) : The tolerance (in standard deviations) used to find the critical parameter values.

		Returns:
			(dict of tuples): A parameter-name keyed dict of tuples consisting of the mean and standard deviation of the optimized values.
		"""
		
		assert method in ("bisect","secant")
				
		# starting parameter values to restore later
		start_params = self.sim.get_model_params().copy()
		
		# get the remaining model parameters that are to be optimized while the selected param is gridded 
		if opt_params == None:
			model_params = self.sim.get_model_params().copy()
		else:
			model_params = OrderedDict( (p,self.sim.get_model_param(p)) for p in model_params )
		
		if sigma == None: # calculate the expected sigma based on the number of observations
			# Andrae, Rene, Tim Schulze-Hartung, and Peter Melchior. "Dos and don'ts of reduced chi-squared." arXiv preprint arXiv:1012.3754 (2010).
			sigma = stdevs * (scipy.sqrt( 2.0 / sum([ e.npoints - len(e.skip) for e in self.sim.get_experiments() ]) ))
		
		# the critical chisq is the point where the low and high parameter estimates are obtained
		critical_chisq = self.sim.run() + sigma
				
		param_values = {}
		for p in params:
			param_values[p] = [None,None]

			opt_params = model_params.copy()
			try: # remove the parameter to be held constant from the list of those to be optimized
				del opt_params[p]
			except KeyError:
				pass

			def target_function( x ): # return the discrepancy between the critical chisq and the chisq of the best fit when parameter p is fixed to argument x
				self.sim.set_model_param(p,x)
				return critical_chisq - self.optimize(opt_params)[1]

			self.sim.set_model_params(**start_params)
			estimate_counter,param_values[p][0] = estimate,model_params[p] * (1.0-estimate)
			chisq_diff = target_function( param_values[p][0] )
			while chisq_diff > 0: # if necessary, decrease the fixed parameter value until we've exceeded the critical chisq
				estimate_counter = estimate_counter + estimate
				if self.verbose:
					print "itc_fit: Lower guess (%f) for parameter \"%s\" is insufficient (%f from critical chisq). Decreasing parameter value to %f." % (param_values[p][0],p,chisq_diff,model_params[p] * (1.0-estimate_counter))
				param_values[p][0] = model_params[p] * (1.0-estimate_counter)
				chisq_diff = target_function( param_values[p][0] )
			
			# actually find the param value that provides the desired confidence interval
			if method == "secant":
				param_values[p][0] = scipy.optimize.brentq( target_function, model_params[p], param_values[p][0], xtol=tolerance )
			else:
				param_values[p][0] = sum(self._noisy_bisect(target_function, model_params[p], param_values[p][0], sigma, chisq_diff, tolerance))/2.0
			
			self.sim.set_model_params(**start_params)
			estimate_counter,param_values[p][1] = estimate,model_params[p] * (1.0+estimate)
			chisq_diff = target_function( param_values[p][1] )
			while chisq_diff > 0: # if necessary, increase the fixed parameter value until we've exceeded the critical chisq
				estimate_counter = estimate_counter + estimate
				if self.verbose:
					print "itc_fit: Upper guess (%f) for parameter \"%s\" is insufficient (%f from critical chisq). Increasing parameter value to %f." % (param_values[p][1],p,chisq_diff,model_params[p] * (1.0+estimate_counter))
				param_values[p][1] = model_params[p] * (1.0+estimate_counter)
				chisq_diff = target_function( param_values[p][1] )
			
			# actually find the param value that provides the desired confidence interval
			if method == "secant":
				param_values[p][1] = scipy.optimize.brentq( target_function, model_params[p], param_values[p][1], xtol=tolerance )
			else:
				param_values[p][1] = sum(self._noisy_bisect(target_function, model_params[p], param_values[p][1], sigma, chisq_diff, tolerance))/2.0

			# From chapter 2 tutorial data (stdevs=2, method='secant'):
			#{'dG': [-10.727396269988125, -11.070049429908266], 'dCp': [-0.09939296467688882, -0.1524609718306775], 'dH': [-11.493031681092932, -12.022877990152883], 'n': [1.778168415517038, 1.8276370690304506]}
			
		# restore initial parameter values for the next iteration
		self.sim.set_model_params(**start_params)
			
		return param_values
				
	def estimate_bootstrap(self, params=[], bootstraps=1000, randomize=0.1, callback=None, logfile=None ):
		"""Generate confidence intervals for optimized parameters by bootstrapping

		Args:
			params (list of strings): The names of the model parameters to find intervals for.
			bootstraps (int): The number of synthetic datasets to refit.
			randomize (float): A fraction by which to perturb the starting parameters before optimization.
			callback (function): A callback function to call at each optimization step, will be passed a vector of the current parameter values.
			logfile (string): A file path to write the optimized parameter values for each bootstrap dataset

		Returns:
			(dict of tuples): A parameter-name keyed dict of tuples consisting of the mean and standard deviation of the optimized values.
		"""

		if self.verbose:
			print "\nitc_fit: Estimating uncertainty intervals for %s using %i bootstraps\n"%(",".join(params),bootstraps)

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
			return [ f + random.choice(res) for f in exp ]

		param_values = OrderedDict( (p,[]) for p in params )
		for i in xrange(bootstraps):
			if self.verbose:
				print "itc_fit: Bootstrap %i"%(i)

			# randomize starting point by user-specifiable amount
			for p in params:
				self.sim.set_model_param(p,self.sim.get_model_param(p) * (1+(2*random.random()-0.5)*randomize))

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
		return OrderedDict( (p,(scipy.mean(param_values[p]),scipy.std(param_values[p]))) for p in params )

	def _check_bounds(self):
		ret = 0
		for k,v in self.sim.get_model_params().iteritems():
			if self.bounds[k][0] != None and v < self.bounds[k][0]:
				ret += (self.bounds[k][0] - v)
				if self.verbose:
					print "itc_fit: Boundary violation for \"%s\" (%f<%f)"%(k,v,self.bounds[k][0])

			elif self.bounds[k][1] != None and v > self.bounds[k][1]:
				ret += (v - self.bounds[k][1])
				if self.verbose:
					print "itc_fit: Boundary violation for \"%s\" (%f>%f)"%(k,v,self.bounds[k][1])

		return ret

	def _fitter(self, func, x0, callback=None):
		# wrapper function for a variety of optimization algorithms
		# note that some of these aren't fully integrated yet
		if self.method == 'simplex':
			ret = scipy.optimize.fmin(
				func=func,
				x0=x0,
				args=(self.sim,),
				disp=self.verbose,
				full_output=True,
				callback=callback,
				**self.method_args)
		elif self.method == 'powell':
			opt = scipy.optimize.fmin_powell(
				func=func,
				x0=x0,
				args=(self.sim,),
				disp=self.verbose,
				full_output=True,
				callback=callback,
				**self.method_args)
			ret = opt,func(opt[0],self.sim)
		elif self.method == 'tnc':
			opt = scipy.optimize.fmin_tnc(
				func=func,
				x0=x0,
				args=(self.sim,),
				approx_grad=True,
				disp=self.verbose,
				**self.method_args)[0]
			ret = opt,func(opt,self.sim)
		elif self.method == 'bfgs':
			ret = scipy.optimize.fmin_l_bfgs_b(
				func=func,
				x0=x0,
				args=(self.sim,),
				approx_grad=True,
				disp=self.verbose,
				**self.method_args)
		else:
			raise Exception('Unrecognized fitting algorithm')

		return ret