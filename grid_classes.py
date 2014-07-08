from multiprocessing	import Pool
from numpy 				import array,dtype,linspace,logspace
from scipy				import optimize

import grid_functions

def _encode_dG(dG,experiment,sim):
	"""
	Given a set of dG parameters, fit the provided experiment by optimizing dH parameters
	Returns a tuple containing the optimized dH parameters and goodness-of-fit: ([dG0,dG1,...],[dH0,dH1,...],Chisq)
	"""

	def fit_Q(x,dG,sim,experiment):
		return itc.fit(
			array(list(dG)+list(x),dtype('d')),
			sim,
			experiment)

	ret = optimize.fmin_powell(
		fit_Q,
		x0=sim.dH[:],
		full_output=True,
		disp=False,
		args=(dG,sim,experiment))

	return (dG,ret[0],ret[1])

class ITCGrid:

	def __init__(threads=None,xtol=1,ftol=1):
		self.threads = threads # number of processing threads to use

		self.xtol = xtol # convergence tolerance
		self.ftol = ftol

		if(threads!=None):
			self.threads = threads
		else:
			self.threads = cpu_count()

		return

	def grid_Q(self, sim, bounds, savefile):
		"""
		Find optimal dH values for each point on the provided boundaries for the FIRST experimental dataset in the sim
		"""

		if len(bounds) != len(sim.dG):
			print "Error:   Boundary size tuple array must match number of dG parameters!"
			return
		for tmp in bounds:
			if len(tmp) != 4:
				print "Error:   Each boundary must be a 4-component tuple: (low,high,steps,logstep)"
				return

		ret = _fit_Q(bounds,sim)

		f = open(savefile,'w')
		for e in ret:
			for f in e:
				f.write("%.5E\t" % (f))
			f.write("\n")
		f.close()
		return


	def fit_Q(dG_bounds,sim):
		"""
		Provided a list of tuples (start,end,steps,logstep) for each dG model parameter, optimizes dH parameters for each coordinate on the grid for the first experiment in the given sim
		Returns the optimized dH parameters and goodness-of-fit as a tuple at each grid coordinate as: ([dG0,dG1,...],[dH0,dH1,...],Chisq)
		"""

		# build a list of the values to test for each parameter
		param_pts = []
		for (low,high,steps,logstep) in dG_bounds:
			if logstep:
				param_pts.append( logspace(low,high,steps) )
			else:
				param_pts.append( linspace(low,high,steps) )

		def get_point( points, iteration ):
			"""
			Given a multidimensional grid of points, provides a unique combination for a given iteration
			"""
			gridpt = []
			for i in range(len(param_pts)):
				gridpt.append( param_pts[i][ iteration % len(param_pts[i]) ] )
				iteration /= len(param_pts[i])
			return gridpt

		# how many points are on the grid?
		n_points = 1
		for i in range(len(param_pts)):
			n_points *= len(param_pts[i])

		# set up the multiprocessing pool
		pool = Pool(processes=sim.threads)
		handles,ret = [],[]

		# optimize dH values for each point the parameter grid
		for i in range(n_points):
			T = sim.experiments.keys()[0] # temperatures
			handles.append( pool.apply_async(
				func = _encode_dG,
				args = ( get_point(param_pts,i), sim.experiments[T][0], sim )
			))

		pool.close()

		grid = [h.get() for h in handles]
		return grid

