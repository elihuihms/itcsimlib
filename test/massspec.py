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

from base import TestITCBase,get_test_data

class TestMSExperiment(TestITCBase):
	def test_experiment_file(self):
		from itcsimlib.mass_spec import MSExperiment
		E = MSExperiment(get_test_data('massspec_1.txt'))
		self.assertIsNotNone( E )
				
	def test_experiment_plot(self):
		from itcsimlib.mass_spec import MSExperiment
		E = MSExperiment(get_test_data('massspec_1.txt'))
		E.make_plot(hardcopy=True, hardcopydir=self.test_dir, hardcopytype="png")
		self.assertTrue( os.path.isfile(os.path.join(self.test_dir,"massspec_1.png")) )

	def test_experiment_export(self):
		from itcsimlib.mass_spec import MSExperiment
		E = MSExperiment(get_test_data('massspec_1.txt'))
		E.export_to_file(path=self.getFilePath(new=True))
		self.assertTrue( os.path.isfile(self.getFilePath()) )

	def test_experiment_chisq_fixedsigma(self):
		from itcsimlib.itc_experiment import ITCExperiment
		from itcsimlib.mass_spec import MSExperiment
		E = MSExperiment(get_test_data('massspec_1.txt'),sigma=0.05)
		Q = np.zeros(shape=(int(E.npoints/E.npops),E.npops))
		self.assertEqual(round(E.get_chisq(Q),1), 292.8)

	def test_experiment_chisq_filesigma(self):
		from itcsimlib.itc_experiment import ITCExperiment
		from itcsimlib.mass_spec import MSExperiment
		E = MSExperiment(get_test_data('massspec_2.txt'),sigma=None)
		Q = np.zeros(shape=(int(E.npoints/E.npops),E.npops))
		self.assertEqual(round(E.get_chisq(Q),1), 282.0)
		self.assertEqual(round(E.sigma,2), 0.01)

class TestMSModelBase(TestITCBase):
	def setUp(self):
		from itcsimlib.model_ising import NonAdditive
		from itcsimlib.mass_spec import MSModel
		TestITCBase.setUp(self)

		self.sim = ITCSim(T0=298.15,verbose=True,threads=1)
		self.model = MSModel(NonAdditive(nsites=11,circular=1))
		self.model.set_params(dGX=-27000,dGY=-27000,dGZ=-30000)
		self.sim.set_model(self.model)
				
	def tearDown(self):
		TestITCBase.tearDown(self)
		self.sim.done()

class TestExpMSModel(TestMSModelBase):
	def test_converted_model(self):
		from itcsimlib.mass_spec import MSExperiment

		self.sim.remove_all_experiments()
		self.sim.add_experiment( MSExperiment(get_test_data('massspec_1.txt')) )
		self.assertEqual(round(self.sim.run(),1), 43.9)

		self.sim.set_model_params(dGX=-27000,dGY=-25000,dGZ=-30000)
		self.assertEqual(round(self.sim.run(),1), 135.5)

class MSExperimentSynthetic(TestMSModelBase):
	def test_synthetic_experiment(self):
		from itcsimlib.mass_spec import MSExperimentSynthetic
		E = MSExperimentSynthetic(lattice_concs=[1E-6]*20, ligand_concs=[1E-6*i for i in range(20)])

		self.sim.remove_all_experiments()
		self.sim.add_experiment( E )
		self.assertTrue(0 < self.sim.run() < 2)

		self.sim.set_model_params(dGX=-30000,dGY=-30000,dGZ=-30000)
		self.assertTrue(self.sim.run() > 2)

	def test_population_plots(self):
		from itcsimlib.mass_spec import MSExperimentSynthetic
		E = MSExperimentSynthetic(lattice_concs=[1E-6]*20, ligand_concs=[1E-6*i for i in range(20)], title="massspec_popplot")

		self.sim.remove_all_experiments()
		self.sim.add_experiment( E )
		self.sim.run()

		self.sim.experiments[0].make_population_plot(dataset='fit', hardcopy=True, hardcopydir=self.test_dir, hardcopyprefix="fit_", hardcopytype="png")
		self.assertTrue( os.path.isfile(os.path.join(self.test_dir,"fit_massspec_popplot.png")) )
		self.sim.experiments[0].make_population_plot(dataset='experimental', hardcopy=True, hardcopydir=self.test_dir, hardcopyprefix="exp_", hardcopytype="png")
		self.assertTrue( os.path.isfile(os.path.join(self.test_dir,"exp_massspec_popplot.png")) )
		self.sim.experiments[0].make_population_plot(dataset='residuals', hardcopy=True, hardcopydir=self.test_dir, hardcopyprefix="res_", hardcopytype="png")
		self.assertTrue( os.path.isfile(os.path.join(self.test_dir,"res_massspec_popplot.png")) )
	
if __name__ == '__main__':
	unittest.main()