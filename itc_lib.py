import os

from ctypes	import cdll,c_double
from numpy	import zeros,dtype

from thermo_functions	import dQ_calc

def sse(params,sim,experiment,writeback=False):
	"""
	Given a set of parameters and an experiment, returns the reduced SSE between the experimental data and resulting fit
	"""

	n = len(experiment.dQ_exp)
	Q_fit = zeros(n,dtype('d'))

	# *NOTE*
	#
	# Loading the dll here incurs a not-insignificant overhead, but ensures
	# that pickling for multiprocess pooling is handled correctly
	#
	# *NOTE*
	itcsimlib = cdll.LoadLibrary(
			os.path.join(
				os.path.dirname( os.path.realpath(__file__) ),
				sim.lib))

	itcsimlib.calc(
		experiment.M_conc.ctypes,
		experiment.L_conc.ctypes,
		Q_fit.ctypes,
		n,
		sim.size,
		sim.cyclic,
		c_double( experiment.T ),
		params.ctypes
	)

	# calculate total heat content in cell at each point
	for i in xrange(n):
		Q_fit[i] *= experiment.M_conc[i] * experiment.V0

	# obtain the change in heat b/t each titration point
	dQ_fit = dQ_calc(Q_fit, experiment.V0, experiment.I_vol)

	for i in xrange(n):
		# normalize by injected amount of material per this titration point
		if experiment.reverse:
			dQ_fit[i] /= experiment.M0 * experiment.I_vol[i]
		else:
			dQ_fit[i] /= experiment.L0 * experiment.I_vol[i]

		# incorrect normalization factor (to CHANGE in titration concentration) Ian used in his thesis:
		#if i==0:
		#	dQ_fit[i] /= experiment.V0*(experiment.L_conc[i] - 0.0)
		#else:
		#	dQ_fit[i] /= experiment.V0*(experiment.L_conc[i] - experiment.L_conc[i-1])

		# remove ligand dilution heat
		if i==0:
			dQ_fit[i] -= experiment.V0*(experiment.L_conc[i])*experiment.dil_Q
		else:
			dQ_fit[i] -= experiment.V0*(experiment.L_conc[i]-experiment.L_conc[i-1])*experiment.dil_Q

	if writeback:
		experiment.dQ_fit = dQ_fit

	# obtain SSE
	sum = 0.0
	for i in xrange(n -1):
		if i not in experiment.skip:
			sum += (dQ_fit[i] - experiment.dQ_exp[i])**2

	# normalize SSE by number of points used in comparison
	return sum / (n-len(experiment.skip))