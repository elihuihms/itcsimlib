#!/usr/bin/env python

import os
import sys
import setuptools
import distutils.ccompiler as C
import itcsimlib

if sys.version_info < (2, 6, 0):
	sys.stderr.write("itcsimlib requires Python 2.6 or newer.\n")
	sys.exit(-1)

setuptools.setup(
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
	install_requires	= [
		'matplotlib>=1.3',
		'scipy>=0.11',
	],
	classifiers			= [
		'Topic :: Scientific/Engineering',
		'Development Status :: 3 - Alpha',
		'Intended Audience :: Science/Research',
		'License :: OSI Approved :: GNU General Public License v2 (GPLv2)'
	],
	packages=['itcsimlib'],
)

# Compile C-based models used in 2017 Bioinformatics paper
compiler = C.new_compiler(verbose=1)
compiler.set_libraries(["m","gslcblas","gsl"])
compiler_arguments = ["-std=c99","-O2"]

if sys.platform.startswith('linux'):
	compiler_arguments.append("-fPIC")

objects = compiler.compile(['itcsimlib/model_trap/itc_sim.c','itcsimlib/model_trap/itc_model.c'],
	extra_preargs=compiler_arguments,
	output_dir="./build")

compiler.set_link_objects(objects)

for s in ['energies_ik','energies_sk','energies_nn']:
	model_object = compiler.compile(["itcsimlib/model_trap/%s.c"%s],
		extra_preargs=compiler_arguments,
		output_dir="./build")

	compiler.link(s, model_object, "%s.so"%s, output_dir="./build/lib/itcsimlib")
