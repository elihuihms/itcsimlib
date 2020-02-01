"""Provides discrete sample parameter values selection during data fitting.


"""

import numpy


class ITCGrid:
	"""A class for either discretely spacing parameters either for different starting conditions or holding them fixed during optimization.

	Attributes
	----------
	fit : ITCFit
		The fitter used for optimization.
	sim : ITCSim
		The simulator (and associated model) used for fitting.
	callback : function
		A function to be called with the current parameter vector after optimization at each grid point.
	verbose : boolean
		Whether or not to print additional information to the console.
	"""

	def __init__(self, fit, start=0, end=None, callback=None, verbose=False ):
		"""The constructor function for the ITCGrid object.

		Arguments
		---------
		fit : ITCFit
			The fitter used for optimization.
		start : integer
			The point index on the grid to start at.
		end : integer
			The point index on the grid to end at.
		callback : function
			A function to be called with the current parameter vector after optimization at each grid point.
		verbose : boolean
			Whether or not to print additional information to the console.
		"""

		self.fit	= fit
		self.sim	= self.fit.sim
		self.verbose	= verbose
		self.callback	= callback

		self._start_index	= start
		self._end_index	= end

		self._grid_size	= 1
		self._grid_pts	= []
		self._grid_order	= []
		self._grid_results	= []

	def add_axis(self, param, start, stop, steps, logspace=False):
		"""Add a parameter discretization axis to the the grid

		Arguments
		---------
		param : string
			The name of the model parameter.
		start : float
			The starting value of the model parameter.
		stop : float
			The ending value of the model parameter.
		steps : integer
			The number of steps to insert between the start and stop.
		logspace : boolean
			Space the steps logarithmically?

		Returns
		-------
		None
		"""

		assert param in self.sim.get_model_params().keys()
		
		if logspace:
			self._grid_pts.append( numpy.logspace(start,stop,steps) )
		else:
			self._grid_pts.append( numpy.linspace(start,stop,steps) )
		self._grid_size *= len(self._grid_pts[-1])
		self._grid_order.append(param)

	def define_axis(self, param, points):
		"""Add a parameter discretization axis to the the grid using a set of points.

		Arguments
		---------
		param : string
			The name of the model parameter.
		points : list of floats
			The points at which to sample the parameter

		Returns
		-------
		None
		"""

		assert param in self.sim.get_model_params().keys()
		
		self._grid_pts.append(points)
		self._grid_size *= len(self._grid_pts[-1])
		self._grid_order.append(param)

	def get_axis_names(self):
		"""Returns the names of the parameters the grid is being evaluated over
		
		Arguments
		---------
		None
		
		Returns
		-------
		list of strings
			The parameter names that constitute the axes of the grid
		
		"""
		return self._grid_order

	def _get_point(self, index):
		gridpt = []
		for i in range(len(self._grid_pts)):
			gridpt.append( self._grid_pts[i][ int(index % len(self._grid_pts[i])) ] )
			index /= len(self._grid_pts[i])

		return gridpt

	def optimize(self, params=[], **kwargs ):
		"""Optimize the model at each point on the grid defined by the parameter axes

		Arguments
		---------
		params : list of strings
			The names of the parameters to optimize.
		**kwargs
			Keyword arguments to pass to the ITCFit optimizer at each grid point

		Returns
		-------
		(list of tuples)
			A list of tuples, where each tuple consists of the grid point and the resulting optimized model parameters
		"""

		assert self._grid_size > 1

		if self.verbose:
			print("\nitc_grid: Optimizing %s parameters over %i grid points of %s parameters\n"%(",".join(params),self._grid_size,",".join(self._grid_order)))

		# archive original model params
		start_params = self.sim.get_model_params().copy()

		# initialize storage
		self._grid_results = [[]]*self._grid_size

		# determine endpoints
		if self._end_index == None:
			self._end_index = self._grid_size

		for i in range(self._start_index, self._end_index):
			point = self._get_point(i)

			# reset starting model params
			self.sim.set_model_params( **start_params )

			for (j,p) in enumerate(self._grid_order):
				self.sim.set_model_param(p,point[j])

			try:
				self._grid_results[i] = (point,)+self.fit.optimize( params=params, **kwargs )
			except Exception as e:
				print("itc_grid: Error, caught exception at grid point index %i: %s"%(i,str(e)))
				continue

			if self.callback != None:
				self.callback( *self._grid_results[i] )

		# restore original model params
		self.sim.set_model_params( **start_params )

		return self._grid_results
		