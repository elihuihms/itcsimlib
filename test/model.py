#!/usr/bin/env python

import unittest
import os
import sys

try:
	from itcsimlib import *
except ImportError:
	sys.path.append(os.path.abspath(".."))
	from itcsimlib import *
from itcsimlib.model_independent import *
from itcsimlib.model_ising import *

from base import TestITCSIM

class TestModel(TestITCSIM):
	def reset_simulation(self,cell="Lattice",syringe="Ligand"):
		self.sim.remove_all_experiments()
		self.sim.add_experiment_synthetic(
				T=298.15,
				V0=1416.6,
				injections=[5.0]*50,
				noise=0.1,
				Cell={cell:1E-6},
				Syringe={syringe:30E-6},
				title='Test_Experiment_1')
		self.sim.add_experiment_synthetic(
				T=298.15+10,
				V0=1416.6,
				injections=[5.0]*50,
				noise=0.1,
				Cell={cell:1E-6},
				Syringe={syringe:30E-6},
				title='Test_Experiment_2')

class TestIndependentModels(TestModel):
	
	def test_onemode_model(self):
		self.reset_simulation(cell="Macromolecule")
		self.sim.set_model( OneMode() )
		self.sim.set_model_params(n = 0.805, dG=-10.94, dH=-11.75)
		self.sim.run() # first run populates the experimental data of a synthetic experiment
		self.sim.set_model_params(n = 1.805, dG=-10.94, dH=-11.75)
		self.assertTrue( self.sim.run() > 1.0 )
	
	def test_twomode_model(self):
		self.reset_simulation(cell="Macromolecule")
		model = NModes(modes=2)
		model.precision = 1E-11
		self.sim.set_model( model )
		self.sim.set_model_params(
			n1 = 0.805, n2 = 1.000,
			dG1=-10.94, dG2=-10.94,
			dH1=-11.75, dH2=-11.75)
		self.sim.run()
		self.sim.set_model_params(
			n1 = 1.805, n2 = 1.000,
			dG1=-10.94, dG2=-10.94,
			dH1=-11.75, dH2=-11.75)
		self.assertTrue( self.sim.run() > 1.0 )
	
class TestIsingModels(TestModel):
	def test_halfadditive_model(self):
		self.reset_simulation()
		self.sim.set_model( HalfAdditive(nsites=5,circular=1) )
		self.sim.set_model_params(
			dG0 = -10, dGb = -1,
			dH0 = -12, dHb = -2,
			dCp0= -1, dCpb =-0.5)
		self.sim.run()
		self.sim.set_model_params(
			dG0 = -10, dGb = -1,
			dH0 = -12, dHb = -2,
			dCp0= 0.0, dCpb=0.0)
		self.assertTrue( self.sim.run() > 1.0 )
		
	def test_fulladditive_model(self):
		self.reset_simulation()
		self.sim.set_model( FullAdditive(nsites=5,circular=1) )
		self.sim.set_model_params(
			dG0 = -10, dGa = 1, dGb = -1,
			dH0 = -12, dHa = 2, dHb = -2,
			dCp0= -1, dCpa =0.5,dCpb=-0.5)
		self.sim.run()
		self.sim.set_model_params(
			dG0 = -10, dGa = 1, dGb = -1,
			dH0 = -12, dHa = 2, dHb = -2,
			dCp0= 0.0, dCpa=0.0,dCpb=0.0)
		self.assertTrue( self.sim.run() > 1.0 )
		
	def test_nonadditive_model(self):
		self.reset_simulation()
		self.sim.set_model( NonAdditive(nsites=5,circular=1) )
		self.sim.set_model_params(
			dGX = -9,  dGY =-10, dGZ = -11,
			dHX = -10, dHY =-12, dHZ = -14,
			dCpX= -1, dCpY =-2, dCpZ = -3)
		self.sim.run()
		self.sim.set_model_params(
			dGX = -8, dGY = -9, dGZ = -10,
			dHX = -8, dHY =-10, dHZ = -12,
			dCpX= -0, dCpY=-1, dCpZ = -2)
		self.assertTrue( self.sim.run() > 1.0 )
		
if __name__ == '__main__':
	unittest.main()