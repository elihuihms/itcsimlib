from multiprocessing	import Pool
from numpy 				import array,dtype
from scipy				import optimize

import itc_lib
from thermo_functions	import *

def fit_H(x,fit,sim):
	"""
	Uses a pool of workers to find optimal dH parameters for a given set of dG parameters using the provided model for each experiment
	"""
	pool = Pool(processes=fit.threads)

	handles = []
	for T in sim.experiments.keys(): # temperatures
		for E in sim.experiments[T]:
			handles.append( pool.apply_async(
				func = itc_lib.sse,
				args = (array(list(E.dG)+list(x),dtype('d')),sim,E)
			))

	pool.close()
	pool.join()

	return sum( [h.get() for h in handles] )

def fit_GH(x,fit,sim):
	"""
	Uses a pool of workers to find optimal dG and dH parameters for the provided model for each experiment
	"""
	pool = Pool(processes=fit.threads)

	handles = []
	for T in sim.experiments.keys(): # temperatures
		for E in sim.experiments[T]:
			handles.append( pool.apply_async(
				func = itc_lib.sse,
				args = (array(x,dtype('d')),sim,E)
			))

	pool.close()
	pool.join()

	return sum( [h.get() for h in handles] )

def fit_GHC(x,fit,sim):
	"""
	Uses  pool of workers to find optimal dG, dH, and dCp parameters for the provided model that best fit each experiment in the sim
	"""

	n = len(x)/3
	params = [0.0]*(n*3)
	pool = Pool(processes=fit.threads)

	handles = []
	for T in sim.experiments.keys(): # temperatures
		for E in sim.experiments[T]:

			#print "\nT:%.3f (%s)" % (T,E.title)
			for i in xrange(n):
				dG  = x[i]
				dH  = x[n+i]
				dCp = x[(2*n)+i]

				params[i]		= dG_vant_Hoff( dG, dH, dCp, T, sim.T_ref )
				params[i+n]		= dH_vant_Hoff( dH, dCp, T, sim.T_ref )
				params[i+(2*n)]	= dCp

			handles.append( pool.apply_async(
				func = itc_lib.sse,
				args = (array(params,dtype('d')),sim,E)
			))

	pool.close()
	pool.join()

	return sum( [h.get() for h in handles] )
