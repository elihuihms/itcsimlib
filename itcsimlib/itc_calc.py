"""Base itc calculation routines for handling multiprocessing.


"""

import sys
import traceback
import multiprocessing


class ITCCalc(multiprocessing.Process):
	"""Worker daemon that uses the provided model to predict titration point enthalpies.
	
	Attributes
	----------
	T0 : float
		The reference temperature used in the simulation.
	model : ITCModel
		The model used to generate per-injection enthalpies.
	iQ : Queue
		The queue to read incoming model parameters and ITCExperiments from.
	oQ : Queue
		The queue to submit calculated enthalpies (or errors) to.
	"""
	
	def __init__(self,T0,model,in_queue,out_queue):
		"""Construtor for the ITCCalc object.
					
		Arguments
		---------
		T0 : float
			The reference temperature used in the simulation
		model : ITCModel
			The model used to calculate per-injection enthalpies.
		in_queue : Queue
			The queue to read incoming model parameters and ITCExperiments from.
		out_queue : Queue
			The queue to put calculated enthalpies into.
		"""
		multiprocessing.Process.__init__(self)
		
		self.T0 = T0
		self.model = model
		self.iQ = in_queue
		self.oQ = out_queue
		self.daemon = True
				
	def run(self):
		"""Start processing submitted experiments in the input queue.
		
		Arguments
		---------
		None
		
		Returns
		-------
		None
		"""
		
		# start the model
		try:
			self.model.start()
		except Exception as exc:
			_type, _value, _traceback = sys.exc_info()
			self.oQ.put( (None,traceback.format_exc()) )

		# pull parameter,experiment tuple from the input queue, without blocking on an empty queue
		for params,E in iter(self.iQ.get, None):
			self.model.set_params(**params)

			try: # in the case of an exception, set the title field to warn the calling thread and stuff the whole exception in the queue
				Q = self.model.Q( self.T0, E.T, E.Concentrations )
			except Exception as exc:
				_type, _value, _traceback = sys.exc_info()
				self.oQ.put( (None,traceback.format_exc()) )
			else:
				self.oQ.put( (E.title,Q) )

		# done with the model now
		self.model.stop()

		return
