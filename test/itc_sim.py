from itcsimlib import ITCSim
from itcsimlib.model_independent import OneMode

from . import TestITCBase


class TestITCSIM(TestITCBase):
	def setUp(self):
		super().setUp()
		self.sim = ITCSim(T0=298.15,units="kcal",verbose=True,threads=1)
		
	def tearDown(self):
		super().tearDown()
		self.sim.done()

class TestSIMThreading(TestITCSIM):

	def test_run_singlethread(self):
		
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
