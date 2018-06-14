#!/usr/bin/env python

import sys
import os

try:
	import itcsimlib
except ImportError:
	sys.path.append(os.path.abspath(".."))

from base import *
from model import *
from trap import *
from drakon import *
from mass_spec import *

unittest.main()