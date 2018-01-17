#!/usr/bin/env python

import unittest
import os
import sys

try:
	from itcsimlib import *
except ImportError:
	sys.path.append(os.path.abspath(".."))
	from itcsimlib import *
from itcsimlib.model_drakon import *

from model import TestModel

class DRAKONModel(DRAKONIsingModel):
    def setup(self):
        self.initialize(nsites=11,circular=True)
        self.add_parameter("dG_0",type="dG")
        self.add_parameter("dG_oe",type="dG")
        self.add_parameter("dG_oo",type="dG")
        self.add_parameter("dH_0",type="dH")
        self.add_parameter("dH_oe",type="dH")
        self.add_parameter("dH_oo",type="dH")
        self.add_parameter("dCp_0",type="dCp")
        self.add_parameter("dCp_oe",type="dCp")
        self.add_parameter("dCp_oo",type="dCp")
    def site(self, i, j):
        if self.occupied(i,j) == True:
            self.add_dG(i, self.dG_0, dH=self.dH_0, dCp=self.dCp_0 )
            self.add_dH(i, self.dH_0, dCp=self.dCp_0 )
            if self.occupied(i,j-1) == True:
                if self.occupied(i,j+1) == True:
                    self.add_dG(i, self.dG_oo, dH=self.dH_oo, dCp=self.dCp_oo )
                    self.add_dH(i, self.dH_oo, dCp=self.dCp_oo )
                else:
                    self.add_dG(i, self.dG_oe, dH=self.dH_oe, dCp=self.dCp_oe )
                    self.add_dH(i, self.dH_oe, dCp=self.dCp_oe )
                    if self.occupied(i,j+1) == True:
                        self.add_dG(i, self.dG_oe, dH=self.dH_oe, dCp=self.dCp_oe )
                        self.add_dH(i, self.dH_oe, dCp=self.dCp_oe )
            else:
                if self.occupied(i,j+1) == True:
                    self.add_dG(i, self.dG_oe, dH=self.dH_oe, dCp=self.dCp_oe )
                    self.add_dH(i, self.dH_oe, dCp=self.dCp_oe )
	
class TestDRAKONModel(TestModel):
	def test_drakon_model(self):
		self.reset_simulation()
		self.sim.set_model( DRAKONModel() )
		self.sim.set_model_params(
			dG_0 = -10, dG_oe = -1, dG_oo = -1.5,
			dH_0 = -12, dH_oe = -2, dH_oo = -2.5,
			dCp_0= -1, dCp_oe=-0.5, dCp_oo=-0.75)
		self.sim.run()
		self.sim.set_model_params(
			dG_0 = -10, dG_oe = -1, dG_oo = -1.5,
			dH_0 = -12, dH_oe = -2, dH_oo = -2.5,
			dCp_0= 0.0, dCp_oe=0.0, dCp_oo=0.0)
		self.assertTrue( self.sim.run() > 1.0 )
		
if __name__ == '__main__':
	unittest.main()