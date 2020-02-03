
from itcsimlib.model_independent import OneMode
from itcsimlib.itc_fit import *

from .itc_sim import TestITCSIM


class TestITCFit(TestITCSIM):
	def setUp(self):
		super().setUp()

		self.sim.add_experiment_file(self.get_data('base_1.txt'))
		self.sim.add_experiment_file(self.get_data('base_2.txt'))
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
		