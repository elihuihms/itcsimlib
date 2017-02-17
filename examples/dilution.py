#!/usr/bin/env python

#
# This script adds a parameter for each dataset (in addition to the model parameters used to fit all datasets) to account the apparent heat of dilution.
# This is very computationally expensive, and will likely lead to problems if your starting model parameter values are poor.
#

from itcsimlib import ITCSim,ITCFit
from itcsimlib.model_trap import *
from itcsimlib.utilities import *

from scipy import optimize

sim = ITCSim(T0=273.15+40,verbose=True,threads=6)
sim.set_model( IKi() )

def optimize_dilution( x ):

	for e in sim.get_experiments():
		sim.remove_experiment( e )
		
	sim.add_experiment_file( 'data/20C-01-TRAPstk-TrpA.txt', skip=[0],	Q_dil=x[0] )
	sim.add_experiment_file( 'data/35C-01-TRAPstk-TrpA.txt', skip=[0],	Q_dil=x[1] )
	sim.add_experiment_file( 'data/40C-01-TRAPstk-TrpA.txt', skip=[0,1],	Q_dil=x[2] )
	sim.add_experiment_file( 'data/45C-01-TRAPstk-TrpA.txt', skip=[0],	Q_dil=x[3] )
	sim.add_experiment_file( 'data/55C-01-TRAPstk-TrpB.txt', skip=[0],	Q_dil=x[4] )
	sim.add_experiment_file( 'data/65C-01-TRAPstk-TrpB.txt', skip=[0,85,91],	Q_dil=x[5] )
		
	sim.set_model_params(*x[6:])
	chisq = sim.run()

	print x,chisq
	return chisq

# The first six are heats of dilution per injected microliter of syringe solution, the following nine parameters are the model parameters
x0 = [3.55846570e+00,-8.36032988e+00,-5.55249866e+00,-6.12875969e+00,-1.87734018e+00,-5.54906632e+00,-3.14784532e+04,-4.98085506e+03,-5.51168479e+03,-4.90956383e+04,-3.22566429e+04,-2.58227800e+04,4.12174495e+02,-1.98322174e+03,-1.83544551e+03]

sim.make_plots(hardcopy=True,hardcopyprefix="dilution_pre_",hardcopytype='pdf')

print optimize.fmin( func=optimize_dilution, x0=x0, full_output=True )

sim.make_plots(hardcopy=True,hardcopyprefix="dilution_post_",hardcopytype='pdf')
sim.done()