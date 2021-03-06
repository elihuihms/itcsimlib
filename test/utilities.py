import unittest
import os
import sys
import pickle

try:
	from itcsimlib import *
except ImportError:
	sys.path.append(os.path.abspath(".."))
	from itcsimlib import *

from base import TestITCBase,get_test_data
from itcsimlib.utilities import *

class TestITCUtilities(TestITCBase):
	def test_rw_params_to_file(self):
		params = {'n':1.80546,'dG':-10.9522,'dH':-12.4281,'dCp':0.0}
		write_params_to_file( self.getFilePath(True), params )
		self.assertEqual( params, read_params_from_file( self.getFilePath() ) )

class TestITCConvertors(TestITCBase):
	def test_rw_experiment_file(self):
		E1 = read_itcsimlib_exp(get_test_data("utilities_1.txt"))
		write_itcsimlib_pickle( self.getFilePath(True), E1 )
		E2 = read_itcsimlib_pickle( self.getFilePath() )
		self.assertEqual( E1.npoints, E2.npoints ) 

	def test_read_nitpic(self):
		E = read_nitpic_exp(get_test_data("utilities_2.nitpkl"))
		self.assertEqual( round(sum(E.dQ_exp),1), -831.3 )

if __name__ == '__main__':
	unittest.main()