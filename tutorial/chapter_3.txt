#
# Chapter 3 : Using The Grid
#

from itcsimlib import *
from itcsimlib.model_independent import *
from itcsimlib.utilities import *

# In the previous chapter we used itcsimlib to identify a combination of model parameters that optimized a fit to three sets of experimental data.
#
# Just as important as finding best-fit parameter values is to understand the confidence intervals for those values. Itcsimlib can do this through bootstrapped error analysis.

# As usual, set up the simulator, load the experimental datasets we used in chapter two, and initialize the ITCFit object:
sim = ITCSim(T0=298.15,units="kcal",verbose=True)
sim.set_model( OneMode() )
sim.add_experiment_file('chapter_2_example_288K.txt', skip=[0])
sim.add_experiment_file('chapter_2_example_298K.txt', skip=[0])
sim.add_experiment_file('chapter_2_example_308K.txt', skip=[0])

fit = ITCFit( sim, verbose=True )

# Next, retrieve and apply the optimized parameters we generated in the previous chapter:
tmp = read_params_from_file( "chapter_2_parameters.txt", row=1 )
sim.set_model_params( **tmp )

# We use the "estimate" method of itcsimlib's ITCFit object to evaluate the bootstrapped datasets.
# For the sake of time, we'll estimate the variance in just the dCp parameter using 20 bootstraps.
# Note that in reality, we'd want to run at least several hundred, and we'd want to let all of the model parameters float.
# We'll also save a record of the bootstrap fits to a file, using the estimate() method's logfile argument:
params = fit.estimate( params=['dCp'], bootstraps=20, logfile="bootstraps_chapter_3.txt" )

# The params dictionary object contains the average value and standard deviation of the dCP parameter, which we can either save to a file or print to the console for our information:
print "dCp Mean: %f kcal +/- %f kcal"%params['dCp']

# Sometimes it's useful to systematically vary one or more parameters, such as sampling different starting conditions for fitting, or to get a better understanding of how model parameter are correlated.
# Let's see what happens when we systematically vary the stoichiometry of our model from 0.5 to 2.5 ligands per macromolecule.

# Itcsimlib makes this easy through the ITCGrid object, and to which we pass our fit object that we want to use:
grid = ITCGrid( fit, verbose=True )

# Add an axis for our grid. In this case, our grid is only one dimension (different stoichiometry values, so it only has one axis).
# Also, let's evaluate our dCp range over 11 linearly-sampled steps:
grid.add_axis( param='n', start=0.5, stop=2.5, steps=11, logspace=False )

# Now, we optimize our parameters at each grid point.
# First, let's use our range of stochiometry values simply as starting points (i.e. not held to the grid value). We'll call grid.optimize() like this:
results = grid.optimize( params=['n','dG','dH','dCp'] )

# Iterate over the results provide by the grid, and write them to a file:
for gridpt,values,chisq in results:
	write_params_to_file( "grid_1_chapter_3.txt", values, pre=gridpt[0], post=chisq, header=False )

# You should notice from this grid search that regardless of the starting stoichiometry, fits always converge to around 1.8.

# Now, let's hold the stoichiometry fixed at each pre-set grid point, and only optimize the other parameters.
# This is done by simply leaving out the dCp parameter from the optimization parameters:
results = grid.optimize( params=['dG','dH','dCp'] )

for gridpt,values,chisq in results:
	write_params_to_file( "grid_2_chapter_3.txt", values, pre=gridpt[0], post=chisq, header=False )

# You will notice that the chi-square values for stoichiometries far from 1.8 are high, but become much smaller as the fixed stoichiometry approaches 1.8

# That concludes this chapter, in chapter 4 we will cover developing and using custom models.

sim.done()