import os
import pytc

from itcsimlib.model_independent import OneMode
from itcsimlib.pytc_integration import ISL2pytcExperiment

from . import TestITCBase

class TestPytcIntegration(TestITCBase,):

	def test_fit(self):
		model = OneMode(units="kcal")
		model.set_params(n=1.0, dG=-8.0, dH=-4.0)

		a = ISL2pytcExperiment(self.get_data('utilities_3.dh'), model, model_params=["n", "dG", "dH"])
#		b = ISL2pytcExperiment(self.get_data('utilities_3.dh', model, model_params=["n", "dG", "dH"])
#		c = ISL2pytcExperiment(self.get_data('utilities_3.dh', model, model_params=["n", "dG", "dH"])

		#fitter = pytc.fitters.BayesianFitter
		fitter = pytc.fitters.MLFitter

		g = pytc.GlobalFit()
		g.add_experiment(a)
#		g.add_experiment(b)
#		g.add_experiment(c)
		g.fit()

		fig, axes = g.plot()

		path = os.path.join(self.test_dir,"pytc_test.png")
		fig.savefig(path, bbox_inches='tight')
		self.assertTrue( os.path.isfile(path) )
		