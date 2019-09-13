import os
import sys

if sys.version_info < (3, 0, 0):
	sys.stderr.write("Error: itcsimlib requires Python 3.0 or newer.\n")
	sys.exit(-1)

from distutils.core import setup, Extension
from distutils.version import LooseVersion
from distutils.command import build

import scipy
if LooseVersion(scipy.__version__) < LooseVersion("0.11"):
	print("Error: itcsimlib requires scipy 0.11 or higher.")

import matplotlib
if LooseVersion(matplotlib.__version__) < LooseVersion("1.3"):
	print("Error: itcsimlib requires matplotlib 1.3 or higher.")

import sympy
if LooseVersion(sympy.__version__) < LooseVersion("1.3"):
	print("Error: itcsimlib requires matplotlib 1.3 or higher.")

import pyx
if LooseVersion(pyx.__version__) < LooseVersion("1.3"):
	print("Error: itcsimlib requires matplotlib 1.3 or higher.")

import itcsimlib

model_sources = ["src/model_trap/itc_model.c","src/model_trap/itc_sim.c","src/model_trap/itc_calc.c"]
model_libraries = ['m','gsl','gslcblas']
model_sk = Extension("itcsimlib.model_trap_sk",
	libraries = model_libraries, extra_compile_args=['-std=c99'],
	sources=["src/model_trap/energies_sk.c"]+model_sources,
	)		
model_ik = Extension("itcsimlib.model_trap_ik",
	libraries = model_libraries, extra_compile_args=['-std=c99'],
	sources=["src/model_trap/energies_ik.c"]+model_sources,
	)		
model_nn = Extension("itcsimlib.model_trap_nn",
	libraries = model_libraries, extra_compile_args=['-std=c99'],
	sources=["src/model_trap/energies_nn.c"]+model_sources,
	)		

class check_c_build(build.build):
	user_options=build.build.user_options + [("build-c-models",None,"Compile TRAP+Tryptophan models written in C")]

	def initialize_options(self, *args, **kwargs):
		self.build_c_models = None
		build.build.initialize_options(self, *args, **kwargs)

	def run(self, *args, **kwargs):
		if self.build_c_models :
			self.distribution.ext_modules = [model_sk,model_ik,model_nn]
		build.build.run(self, *args, **kwargs)
		
setup(
	name				= 'itcsimlib',
	version				= itcsimlib.__version__,
	license				= itcsimlib.__license__,
	author				= itcsimlib.__author__,
	author_email		= itcsimlib.__author_email__,
	description			= u'itcsimlib: isothermal titration calorimetry simulation using statistical thermodynamics',
	long_description	= open('README.txt').read(),
	url					= u'https://github.com/elihuihms/itcsimlib',
	download_url		= u'https://github.com/elihuihms/itcsimlib/archive/master.zip',
	platforms			= 'any',
	classifiers			= [
		'Topic :: Scientific/Engineering',
		'Topic :: Scientific/Engineering :: Chemistry',
		'Topic :: Scientific/Engineering :: Physics',
		''
		'Development Status :: 5 - Stable',
		'Intended Audience :: Science/Research',
		'License :: OSI Approved :: GNU General Public License v2 (GPLv2)'
	],
	packages=['itcsimlib'],
	ext_modules = [],
	cmdclass={
		'build':check_c_build,
	},
)