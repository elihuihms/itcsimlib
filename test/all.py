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

unittest.main()