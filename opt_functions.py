import sys

from multiprocessing	import Pool
from numpy 				import array,dtype
from scipy				import optimize

import itc_lib

def _optimize_dH(dG,dH,experiment,sim):
	"""
	Fit the provided experiment by optimizing dH parameters at a given set of dGs
	Returns a tuple containing the optimized dH parameters and goodness-of-fit: (experiment title,[dH0,dH1,...],Chisq)
	"""

	def sse(x):
		return itc_lib.sse( array(list(dG)+list(x),dtype('d')), sim, experiment)

	ret = optimize.fmin_powell(
		func=sse,
		x0=dH[:],
		full_output=True,
		disp=False)

	sys.stdout.write("#")
	sys.stdout.flush()

	return (experiment.title,ret[0],ret[1])

def fit_G_opt_H(fit,sim):
	"""
	Finds optimal dG for all datasets, and optimizes dH for each individually.
	Returns a dict of experiment titles and optimized dH parameters
	"""

	# initialize experiment-title-keyed dict of dH parameters
	dH = {}
	for T in sim.experiments.keys():
		for E in sim.experiments[T]:
			dH[ E.title ] = list(E.dH)

	def _update_dH(x,dH,sim):
		pool = Pool(processes=fit.threads)

		sys.stdout.write("dG: ")
		for f in x:
			sys.stdout.write("%.5E " % f)
		sys.stdout.write("\ndH optimization: ")
		sys.stdout.flush()

		handles = []
		for T in sim.experiments.keys():
			for E in sim.experiments[T]:
				handles.append( pool.apply_async(
					func = _optimize_dH,
					args = (x,dH[ E.title ],E,sim)
				))

		pool.close()
		pool.join()

		sys.stdout.write("\n")
		sys.stdout.flush()

		sum = 0
		for tmp in [h.get() for h in handles]:
			sum += tmp[2]
			dH[ tmp[0] ] = tmp[1]

		return sum

	tmp = optimize.fmin(
		_update_dH,
		x0=sim.dG,
		args=(dH,sim),
		xtol=fit.xtol,
		ftol=fit.ftol,
		disp=True)

	return (tmp,dH)
