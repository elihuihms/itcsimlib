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

	Attributes
	----------
	model : ITCModel
		The model to use for fitting the experimental data.
	bounds : dict of tuples
		A parameter-name keyed dict of low and high bounds to enforce during fitting (retrieved from the model itself, or explicitly provided by the user).
	chisq : float
		The most recently evaluated goodness-of-fit chisquare.
	"""

	def __init__(self, sim, method='simplex', method_args={}, verbose=False):
		"""Constructor function for the ITCFit object.

		Arguments
		---------
		sim : ITCSim
			The ITCSim object that possesses the experimental data to replicate.
		bounds : dict of tuples
			A parameter-name keyed dict of low and high bounds to enforce during fitting.
		method : string
			The optimization algorithm to use.
		method_args : dict
			Arguments to pass to the optimization algorithm (method specific).
		verbose : boolean
			Print additional information to the console?
		"""

		self.sim = sim
		self.model	= self.sim.model
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
		self.model = self.sim.model
		
	def get_sim(self):
		return self.sim
	
	def add_bounds(self, param, low=None, high=None):
		"""Add low and/or high boundaries for the specified parameter.

		Arguments
		---------
		param : string
			The name of the parameter to create the boundary for
		low : float
			The lower boundary for the parameter, or None for no boundary
		high : float
			The upper boundary for the parameter, or None for no boundary

		Returns
		-------
		None
		
		Notes
		-----
			Implementation of parameter boundaries is fit method dependent and can be unpredictable. Use at your own risk, and sparingly.
		"""

		assert param in self.model.get_param_names()
		self.bounds[param] = (low,high)

	def optimize(self, params=[], callback=None, update_fits=False ):
		"""Optimize the specified parameters.

		Arguments
		---------
		params : dict of strings
			The names of the model parameters to optimize.
		callback : function
			A callback function to call at each optimization step, will be passed a vector of the current parameter values.
		update_fits : boolean
			Update the dQ_fit attributes of the experiments in the simulator object?

		Returns
		_______
		(dict,float)
			A parameter-name keyed dict of optimized values, and the corresponding goodness-of-fit.
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
		"""Wrapper for the two methods of estimating uncertainties in the fitted parameter values
		
		Arguments
		--------
		params : list of strings
			The names of the model parameters to find intervals for.
		method : string
			The method to use for interval estimation (either "sigma" or "bootstrap")
		*args
			Positional arguments to pass to either estimate_sigma or estimate_bootstrap
		**kwargs
			Keyword arguments to pass to either estimate_sigma or estimate_bootstrap
		
		Returns
		-------
		dict of tuples
			A parameter-name keyed dict of tuples consisting of the mean and standard deviation, or high and low values of the provided parameters.
		"""
					
		assert method in ['sigma','bootstrap']
	
		if method == 'sigma':
			return self.estimate_sigma( params, *args, **kwargs )
		else:
			return self.estimate_bootstrap( params, *args, **kwargs )
	
	def estimate_sigma(self, params=[], opt_params=None, sigma=None, stdevs=1, estimate=0.1, rootfinder='bisect', tolerance=0.001 ):
		"""Generate high and low parameter value estimates for optimized parameters according to the provided criterion
				
		Arguments
		---------
		params : list of strings
			The names of the model parameters to find intervals for.
		opt_params : dict of strings
			The names of the model parameters to optimize. Setting to None implies that all available model parameters will be optimized, except the parameter that is being evaluated.
		sigma : float
			The critical chisq value of the overall fit to use for estimating a parameter's expected maximum and minimum (to be used in lieu of the stdevs parameter) 
		stdevs : int
			The number of standard deviations above the best-fit chi-square value to use as the critical cutoff for estimating a parameter's expected maximum and minimum (to be used in lieu of the sigma parameter)
		rootfinder : string
			The algorithm used to find the root of the target function ("bisect" or "secant")
		tolerance : float
			The tolerance (in standard deviations) used to find the critical parameter values. Smaller numbers mean a more accurate estimation of a parameter's estimated maximum and minimum values.

		Returns
		-------
		(dict of tuples)
			A parameter-name keyed dict of tuples consisting of the high and low estimates for each of the parameters in the "params" argument.
		"""
		
		assert rootfinder in ("bisect","secant")
				
		# starting parameter values to restore later
		start_params = self.sim.get_model_params().copy()
		
		# get the remaining model parameters that are to be optimized while the selected param is gridded 
		if opt_params == None:
			model_params = self.sim.get_model_params().copy()
		else:
			model_params = OrderedDict( (p,self.sim.get_model_param(p)) for p in model_params )
		
		if sigma == None: # calculate the expected sigma based on the number of observations
			# Andrae, Rene, Tim Schulze-Hartung, and Peter Melchior. "Dos and don'ts of reduced chi-squared." arXiv preprint arXiv:1012.3754 (2010).
			sigma = stdevs * (scipy.sqrt( 2.0 / sum([ e.npoints - len(e.skip) for e in self.sim.experiments ]) ))
		
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
			if rootfinder == "secant":
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
			if rootfinder == "secant":
				param_values[p][1] = scipy.optimize.brentq( target_function, model_params[p], param_values[p][1], xtol=tolerance )
			else:
				param_values[p][1] = sum(self._noisy_bisect(target_function, model_params[p], param_values[p][1], sigma, chisq_diff, tolerance))/2.0

			# From chapter 2 tutorial data (stdevs=2, rootfinder='secant'):
			#{'dG': [-10.727396269988125, -11.070049429908266], 'dCp': [-0.09939296467688882, -0.1524609718306775], 'dH': [-11.493031681092932, -12.022877990152883], 'n': [1.778168415517038, 1.8276370690304506]}
			
		# restore initial parameter values for the next iteration
		self.sim.set_model_params(**start_params)
			
		return param_values
				
	def estimate_bootstrap(self, params=[], bootstraps=1000, randomize=0.1, callback=None, logfile=None ):
		"""Generate confidence intervals for optimized parameters by bootstrapping (also known as jackknife estimation)

		Arguments
		---------
		params : list of strings
			The names of the model parameters to find intervals for.
		bootstraps : int
			The number of synthetic datasets to refit.
		randomize : float
			A fraction by which to perturb the starting parameters before optimization.
		callback : function
			A callback function to call at each optimization step, will be passed a vector of the current parameter values.
		logfile : string
			A file path to write the optimized parameter values for each bootstrap dataset

		Returns
		-------
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
		for E in self.sim.experiments:
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
			for j,E in enumerate(self.sim.experiments):
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
		for i,E in enumerate(self.sim.experiments):
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