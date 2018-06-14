#!/usr/bin/env python

import unittest
import os
import sys
import numpy as np

try:
	from itcsimlib import *
except ImportError:
	sys.path.append(os.path.abspath(".."))
	from itcsimlib import *

from base import TestITCBase

class TestMSExperiment(TestITCBase):
	def test_experiment_file(self):
		from itcsimlib.mass_spec import MSExperiment
		E = MSExperiment('./data/massspec_1.txt')
		self.assertIsNotNone( E )
				
	def test_experiment_plot(self):
		from itcsimlib.mass_spec import MSExperiment
		E = MSExperiment('./data/massspec_1.txt')
		E.make_plot(hardcopy=True, hardcopydir=self.test_dir, hardcopytype="png")
		self.assertTrue( os.path.isfile(os.path.join(self.test_dir,"massspec_1.png")) )
		
	def test_experiment_chisq(self):
		from itcsimlib.itc_experiment import ITCExperiment
		from itcsimlib.mass_spec import MSExperiment
		E = MSExperiment('./data/massspec_1.txt')
		Q = np.zeros(shape=(E.npoints,E.npops))
		self.assertEqual(round(E.get_chisq(Q),1), 292.8)

class TextMSModel(TestITCBase):
	def setUp(self):
		TestITCBase.setUp(self)
		self.sim = ITCSim(T0=298.15,units="kcal",verbose=True,threads=1)
		
	def tearDown(self):
		TestITCBase.tearDown(self)
		self.sim.done()

	def test_convert_model(self):
		from itcsimlib.model_ising import NonAdditive
		from itcsimlib.mass_spec import convert_model

		MSNonAdditive = convert_model(NonAdditive)
		M = MSNonAdditive(nsites=11,circular=1)
		M.set_params(dGX=-27000,dGY=-27000,dGZ=-30000)

		self.assertIsNotNone(M)

	def test_run(self):
		from itcsimlib.model_ising import NonAdditive
		from itcsimlib.mass_spec import convert_model,MSExperiment

		MSNonAdditive = convert_model(NonAdditive)
		M = MSNonAdditive(nsites=11,circular=1)
		M.set_params(dGX=-27000,dGY=-27000,dGZ=-30000)

		self.sim.remove_all_experiments()
		self.sim.add_experiment( MSExperiment('./data/massspec_1.txt') )
		self.sim.set_model( M )
		
		self.assertEqual(round(self.sim.run(),1), 43.9)
	
if __name__ == '__main__':
	unittest.main()