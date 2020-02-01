#!/usr/bin/env python

import unittest
import os
import random
import sys
import shutil
import tempfile
import uuid

try:
	from itcsimlib import *
except ImportError:
	sys.path.append(os.path.abspath(".."))
	from itcsimlib import *

from itcsimlib.utilities import *


def get_test_data(filename):
	path = os.path.dirname(os.path.abspath(__file__)) # ../itcsimlib/test
	path = os.path.join(path,"data") # ../itcsimlib/test/data
	return os.path.join(path,filename)

class TestITCBase(unittest.TestCase):
	def setUp(self):
		self.test_dir = tempfile.mkdtemp()
		
	def tearDown(self,clean=True):
		if clean:
			shutil.rmtree(self.test_dir)
		else:
			print("Temporary file path: %s"%self.test_dir)
		
	def getFilePath(self,new=False):
		if new:
			self.test_file = os.path.join(self.test_dir,str(uuid.uuid4()))
		return self.test_file
	
	def getFileContents(self,path):
		with open(path, 'r') as content_file:
			return content_file.read()

class TestITCExperiment(TestITCBase):

	def test_synthetic_experiment(self):
		from itcsimlib.itc_experiment import ITCExperimentSynthetic
		self.assertIsNotNone( ITCExperimentSynthetic(
				T=298.15,
				V0=1416.6,
				injections=[5.0]*50,
				Cell={"Macromolecule":1E-6},
				Syringe={"Ligand":30E-6},
				noise=1.0,
				title='Test_Experiment_1'))
				
	def test_experiment_export(self):
		from itcsimlib.itc_experiment import ITCExperiment
		E = read_itcsimlib_exp(get_test_data('base_1.txt'))
		E.export_to_file( self.getFilePath(True) )
		self.assertTrue( os.path.isfile(self.getFilePath()) )

	def test_experiment_plot(self):
		from itcsimlib.itc_experiment import ITCExperiment
		E = read_itcsimlib_exp(get_test_data('base_1.txt'), exp_args={"skip":[0,1,2,10]})
		E.dQ_err = [5E3*random.random() for dQ in E.spline]
		E.dQ_fit = E.spline
		E.make_plot(hardcopy=True, hardcopydir=self.test_dir, hardcopytype="png")
		self.assertTrue( os.path.isfile(os.path.join(self.test_dir,"base_1.png")) )
		
	def test_experiment_chisq(self):
		from itcsimlib.itc_experiment import ITCExperiment
		E = read_itcsimlib_exp(get_test_data('base_1.txt'))
		Q = [0.0]*E.injections
		self.assertEqual( round(E.get_chisq(Q),1), 691.5 )

class TestITCSIM(TestITCBase):
	def setUp(self):
		TestITCBase.setUp(self)
		self.sim = ITCSim(T0=298.15,units="kcal",verbose=True,threads=1)
		
	def tearDown(self):
		TestITCBase.tearDown(self)
		self.sim.done()

class TestSIMThreading(TestITCSIM):

	def test_run_singlethread(self):
		from itcsimlib.model_independent import OneMode
		
		self.sim.remove_all_experiments()
		self.sim.add_experiment_synthetic(
				T=298.15,
				V0=1416.6,
				injections=[5.0]*50,
				noise=1,
				Cell={"Macromolecule":1E-6},
				Syringe={"Ligand":30E-6},
				title='Test_Experiment_1')
		self.sim.set_model( OneMode() )
		self.sim.set_model_params(n=1.805,dG=-10.94,dH=-11.75,dCp=0.0)
		self.sim.run()
		self.sim.set_model_params(n=2,dG=-11,dH=-12,dCp=0.0)
		self.sim.run()

		self.assertTrue( self.sim.run() > 1 )

	def test_run_multithread(self):
		from itcsimlib.model_independent import OneMode

		multi = ITCSim(T0=298.15,units="kcal",verbose=True,threads=2)
		multi.add_experiment_synthetic(
				T=298.15,
				V0=1416.6,
				injections=[5.0]*50,
				noise=1,
				Cell={"Macromolecule":1E-6},
				Syringe={"Ligand":30E-6},
				title='Test_Experiment_1')
		multi.add_experiment_synthetic(
				T=298.15,
				V0=1416.6,
				injections=[5.0]*50,
				noise=0.05,
				Cell={"Macromolecule":2E-6},
				Syringe={"Ligand":60E-6},
				title='Test_Experiment_2')
		multi.set_model( OneMode() )
		multi.set_model_params(n=1.805,dG=-10.94,dH=-11.75,dCp=0.0)
		multi.run()
		multi.set_model_params(n=2,dG=-11,dH=-12,dCp=0.0)

		self.assertTrue( multi.run() > 1 )

		multi.done()

class TestITCFit(TestITCSIM):
	def setUp(self):
		TestITCSIM.setUp(self)
		from itcsimlib.model_independent import OneMode

		self.sim.add_experiment_file(get_test_data('base_1.txt'))
		self.sim.add_experiment_file(get_test_data('base_2.txt'))
		self.sim.set_model( OneMode() )
		self.sim.set_model_params(n=1.80546,dG=-10.9522,dH=-12.4281,dCp=0.0)
	
	def test_fit_optimize_simplex(self):
		fit = ITCFit( self.sim, method='simplex', method_args={"maxiter":1} )
		self.assertEqual( round(fit.optimize(params=['n','dG','dH'])[1],3), 2.695 )
	
	def test_fit_optimize_powell(self):
		fit = ITCFit( self.sim, method='powell', method_args={"maxiter":1} )
		self.assertEqual( round(fit.optimize(params=['n','dG','dH'])[1],3), 2.695 )

	def test_fit_optimize_tnc(self):
		fit = ITCFit( self.sim, method='tnc', method_args={"maxfun":1} )
		self.assertEqual( round(fit.optimize(params=['n','dG','dH'])[1],3), 2.695 )
	
	def test_fit_optimize_bfgs(self):
		fit = ITCFit( self.sim, method='bfgs', method_args={"maxiter":1} )
		self.assertEqual( round(fit.optimize(params=['n','dG','dH'])[1],3), 2.695 )

	def test_fit_estimate_bootstrap_(self):
		fit = ITCFit( self.sim, method='simplex', method_args={"maxiter":1} )
		self.assertIsNotNone( fit.estimate(params=['n','dG','dH'], method='bootstrap', bootstraps=5) )

	def test_fit_estimate_sigma_bisect(self):
		fit = ITCFit( self.sim, method='simplex', method_args={"maxiter":1} )
		self.assertIsNotNone( fit.estimate(params=['n','dG','dH'], method='sigma', rootfinder='bisect', stdevs=2) )

	def test_fit_estimate_sigma_secant(self):
		fit = ITCFit( self.sim, method='simplex', method_args={"maxiter":1} )
		self.assertIsNotNone( fit.estimate(params=['n','dG','dH'], method='sigma', rootfinder='secant', stdevs=2) )
		
class TestITCGrid(TestITCSIM):
	def setUp(self):
		TestITCSIM.setUp(self)
		from itcsimlib.model_independent import OneMode

		self.sim.add_experiment_file(get_test_data('base_1.txt'))
		self.sim.set_model( OneMode() )
		self.sim.set_model_params(n=1.805,dG=-10.94,dH=-11.75,dCp=0.0)		
		self.fit = ITCFit( self.sim, method='simplex', method_args={"maxiter":1} )

	def test_2d_grid(self):
		grid = ITCGrid( self.fit )
		grid.add_axis(param='dG',start=-12,stop=-10,steps=2)
		grid.add_axis(param='dH',start=-13,stop=-11,steps=2)
		self.assertIsNotNone( grid.optimize(params="n") )
		
	def test_defined_grid(self):
		grid = ITCGrid( self.fit )
		grid.define_axis(param='dG',points=(-11.5,-10.94,-10.5))
		grid.define_axis(param='dH',points=(-12.25,-11.75,-11.50))
		self.assertIsNotNone( grid.optimize(params="n") )
	
if __name__ == '__main__':
	unittest.main()