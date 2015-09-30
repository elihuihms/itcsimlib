from multiprocessing	import Process
from scipy				import array,dtype,zeros
from ctypes				import cdll,c_double,c_int

from thermo				import dQ_calc

class ITCCalc(Process):
	"""Worker daemon that uses a provided function to predict titration point enthalpies.
	
	Note:
		Use this version for python-native models.

	Attributes:
		calc (function): A model-defined function used to calculate per-injection enthalpies.
		iQ (Queue): The queue to read incoming model parameters and ITCExperiments from.
		oQ (Queue): The queue to put calculated enthalpies into.
	"""
	
	def __init__(self,calc,in_queue,out_queue):
		"""Initialize the ITCCalc object.
					
		Arguments:
			calc (function): A model-defined function used to calculate per-injection enthalpies.
			iQ (Queue): The queue to read incoming model parameters and ITCExperiments from.
			oQ (Queue): The queue to put calculated enthalpies into.
		"""
		super(ITCCalc, self).__init__()

		self.calc = calc
		self.iQ = in_queue
		self.oQ = out_queue
		self.daemon = True

	def run(self):
		"""Start looking for experiments to process in the input queue.
		
		Arguments:
			None
		
		Returns:
			None
		"""
		
		# pull parameter,experiment tuple from the input queue, without blocking on an empty queue
		for params,experiment in iter(self.iQ.get, None):
			Q_fit = [0.0]*experiment.npoints

			# pass the incoming experiment and provided parameters to our model-specific calculator function
			try:
				self.calc(
					experiment.npoints,
					experiment.T,
					experiment.M_conc,
					experiment.L_conc,
					Q_fit,
					params
				)

			# in the case of an exception, set the title field to warn the calling thread and stuff the whole exception in the queue
			except Exception as exc:
				self.oQ.put( (None,exc) )
			else:
				self._put_Q( Q_fit, experiment )

	def _put_Q(self, Q_fit, experiment):
		# normalize total heat content to macromolecule concentration in the cell volume
		for i in xrange(experiment.npoints):
			Q_fit[i] *= experiment.V0 * experiment.M_conc[i]

		# obtain the change in cell heat b/t each titration point
		dQ_fit = dQ_calc(Q_fit, experiment.V0, experiment.I_vol)

		for i in xrange(experiment.npoints):

			# normalize by injected amount of material per this titration point
			# remove heat of ligand or macromolecule dilution
			if experiment.reverse:
				dQ_fit[i] /= experiment.M0 * experiment.I_vol[i]

				if i==0:
					dQ_fit[i] -= experiment.V0*(experiment.M_conc[i])*experiment.dil_Q
				else:
					dQ_fit[i] -= experiment.V0*(experiment.M_conc[i]-experiment.M_conc[i-1])*experiment.dil_Q
			else:
				dQ_fit[i] /= experiment.L0 * experiment.I_vol[i]

				if i==0:
					dQ_fit[i] -= experiment.V0*(experiment.L_conc[i])*experiment.dil_Q
				else:
					dQ_fit[i] -= experiment.V0*(experiment.L_conc[i]-experiment.L_conc[i-1])*experiment.dil_Q

		self.oQ.put( (experiment.title,dQ_fit) )

class ITCCalcDLL(ITCCalc):
	"""Worker daemon that uses an external cdll to predict titration point enthalpies.
	
	Note:
		Use this version for shared library models.
	
	Attributes:
		calc (function): A model-defined function used to calculate per-injection enthalpies.
		iQ (Queue): The queue to read incoming model parameters and ITCExperiments from.
		oQ (Queue): The queue to put calculated enthalpies into.
	"""

	def __init__(self,path,setup_args,in_queue,out_queue):
		"""Initialize the ITCCalcDLL object.
					
		Arguments:
			path (string): The path to the shared library.
			setup_args (list): List of arguments to pass to the library's "setup" function.
			iQ (Queue): The queue to read incoming model parameters and ITCExperiments from.
			oQ (Queue): The queue to put calculated enthalpies into.
		"""
		super(ITCCalcDLL, self).__init__(None,in_queue,out_queue)

		self.lib = cdll.LoadLibrary(path)
		self.lib.setup(*setup_args)

	def run(self):
		"""Start looking for experiments to process in the input queue.
		
		Arguments:
			None
		
		Returns:
			None
		"""

		for (params,experiment) in iter(self.iQ.get, None):
			Q_fit = zeros(experiment.npoints,dtype('d'))

			status = self.lib.calc(
				c_int( experiment.npoints ),
				c_double( experiment.T ),
				experiment.M_conc.ctypes,
				experiment.L_conc.ctypes,
				Q_fit.ctypes,
				array(params,dtype('d')).ctypes
			)

			# if the libraries' calc function returns nonzero, create an exception and stuff it in the queue
			if status != 0:
				self.oQ.put( (None,Exception("DLL returned a non-zero error code: %i"%(status))) )
			else:
				self._put_Q( Q_fit, experiment )

		self.lib.close()
