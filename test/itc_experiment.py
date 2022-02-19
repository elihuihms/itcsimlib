import os
import random

from itcsimlib.utilities import *
from itcsimlib.itc_experiment import ITCExperiment, ITCExperimentSynthetic

from . import TestITCBase


class TestITCExperiment(TestITCBase):

	def test_synthetic_experiment(self):
		self.assertIsNotNone( ITCExperimentSynthetic(
				T=298.15,
				V0=1416.6,
				injections=[5.0]*50,
				Cell={"Macromolecule":1E-6},
				Syringe={"Ligand":30E-6},
				noise=1.0,
				title='Test_Experiment_1'))
				
	def test_experiment_export(self):
		E = read_itcsimlib_exp(self.get_data('base_1.txt'))
		E.export_to_file( self.getFilePath(True) )
		self.assertTrue( os.path.isfile(self.getFilePath()) )

	def test_experiment_plot(self):
		E = read_itcsimlib_exp(self.get_data('base_1.txt'), exp_args={"skip":[0,1,2,10]})
		E.dQ_err = [5E3*random.random() for dQ in E.spline]
		E.dQ_fit = E.spline
		E.make_plot(hardcopy=True, hardcopydir=self.test_dir, hardcopytype="png")
		self.assertTrue( os.path.isfile(os.path.join(self.test_dir,"base_1.png")) )
		
	def test_experiment_chisq(self):
		E = read_itcsimlib_exp(self.get_data('base_1.txt'))
		Q = [0.0]*E.injections
		self.assertEqual( round(E.get_chisq(Q),1), 691.5 )

