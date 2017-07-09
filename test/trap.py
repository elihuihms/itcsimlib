#!/usr/bin/env python

import unittest
import sys

try:
	from itcsimlib import *
except ImportError:
	sys.path.append(os.path.abspath(".."))
	from itcsimlib import *
from itcsimlib.utilities import *
from itcsimlib.model_trap import *

from model import TestModel

class TestTRAPModels(TestModel):
	def test_SK_model(self):
		self.reset_simulation(cell="TRAP",syringe="Trp")
		self.sim.set_model( SK() )
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
	
	def test_SKa_model(self):
		self.reset_simulation(cell="TRAP",syringe="Trp")
		self.sim.set_model( SKa() )
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
	
	def test_IK_model(self):
		self.reset_simulation(cell="TRAP",syringe="Trp")
		self.sim.set_model( IK() )
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
		
	def test_IKi_model(self):
		self.reset_simulation(cell="TRAP",syringe="Trp")
		self.sim.set_model( IKi() )
		self.sim.set_model_params(
			dG0 = -10, dGoe = 1, dGoo = -1,
			dH0 = -12, dHoe = 2, dHoo = -2,
			dCp0= -1, dCpoe =0.5,dCpoo=-0.5)
		self.sim.run()
		self.sim.set_model_params(
			dG0 = -10, dGoe = 1, dGoo = -1,
			dH0 = -12, dHoe = 2, dHoo = -2,
			dCp0= 0.0, dCpoe=0.0,dCpoo=0.0)
		self.assertTrue( self.sim.run() > 1.0 )

if __name__ == '__main__':
	unittest.main()