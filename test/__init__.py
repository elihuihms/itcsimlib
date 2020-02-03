import os
import sys
import tempfile
import unittest
import shutil
import uuid

sys.path.insert(1,os.path.abspath("."))


class TestITCBase(unittest.TestCase):
	def setUp(self):
		self.test_dir = tempfile.mkdtemp()

	def get_data(self, filename):
		path = os.path.dirname(os.path.abspath(__file__)) # ../itcsimlib/test
		path = os.path.join(path,"data") # ../itcsimlib/test/data
		return os.path.join(path,filename)
			
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

from .itc_experiment import *
from .itc_fit import *
from .itc_grid import *
from .itc_model import *
from .itc_sim import *
from .mass_spec import *
#from .model_drakon import *
from .model_trap import *
from .pytc_integration import *
from .utilities import *