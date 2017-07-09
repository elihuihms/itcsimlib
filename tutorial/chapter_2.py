#
# Chapter 2 : Model Fitting
#

from itcsimlib import *

sim = ITCSim(T0=298.15,units="kcal",verbose=True)

# Model fitting is when we use a theoretical model to fit experimental data, in the hope that:
# 1) Our model can faithfully replicate the experimental data (i.e a "good" fit)
# 2) The model parameters that provide the aforementioned good fit tell us something informative

# Itcsimlib can help you satisfy only the first problem. The second, and much more philosophical one, is an "exercise left to the reader".

# To import experimental data, itcsimlib expects a fairly specific format, and it's worth looking at the experimental data file "chapter_2_experiment.txt"
# Unlike Python scripts, itcsimlib does not entirely ignore the lines in an experimental data file that start with a pound (#) sign, as it uses these lines to store some valuable information.
# Let's look at the first several lines of chapter_2_example.txt in the tutorial "data" folder:

# units kcal
# T	298.15000
# V0	1416.60000
# Cell	Macromolecule	1.00000E-06
# Syringe	Ligand	2.00000E-05
# Q_dil	0.00000E+00

# T, as you may gather, is the temperature the experiment took place at.
# V0, as we saw in Chapter 1, is the cell volume.
# Lines that start with "Syringe" inform itcsimlib what model components were present in the syringe solution. Here, we see that the Ligand component was present at 20 micromolar.
# Lines that start with "Cell", on the other hand, inform itcsimlib what the (starting) concentration of the named component was in the cell. Over the course of the titration, this will obviously be diluted a bit as we titrate in syringe components.
# Note that it's possible to have the same component in both the syringe and the cell, potentially even at different concentrations. However, components that are not specifically acted upon by the model will be ignored.

sim.add_experiment_file('chapter_2_example_298K.txt', skip=[0])

# NOTE! Any header fields in the experiment file that start with the pound sign followed by a space are automatically passed on by the add_experiment_file() function, so we don't have to manually provide arguments. If we do, i.e. add_experiment_file("file.txt", T=298.15), we can override the experimental temperature specified in the file header.
# Here, the only extra argument we are providing is the "skip" argument, which tells itcsimlib to ignore specific injection points, in this case just the first one.

# To fit this data, let's use a very simple model that describes a single binding mode:
from itcsimlib.model_independent import OneMode

sim.set_model( OneMode() )

# Because this model has only one binding modes it has only four parameters:
# Stoichiometry: n, i.e. the number of sites that bind with this mode.
# Affinity: dG, i.e. the free energy change upon binding.
# Enthalpy: dH, i.e. the enthalpy change upon binding.
# Heat capacity: dCp, i.e. the heat capacity change upon binding, which affects both dG and dH.

# Using a number of algorithms of your choosing, itcsimlib can optimize any (or all) of these parameters to values that minimize the discrepancy between the specified model experimental data.
# If you want to make fancier models, e.g. reducing the number of parameters by assuming that both modes have the same change in heat capacity, see tutorial Chapter 4.

# Numerical optimization of model parameters IS NOT MAGIC. If you have bad starting conditions, then numerical optimization can either fail outright, or you'll converge to a local minimum. Even if you start with a decent set of starting conditions, you may still converge to a local minimum.
# There's simply not enough space in this tutorial to describe the numerous potential pitfalls presented by numerical optimization/fitting, but grid fitting, in which parameters are systematically varied is described in tutorial chapter ?.

# Set the starting conditions to a reasonable guess of 2 binding sites with ~1uM Kd and -10kcal/mol enthalpy change. Since we're only analyzing a single temperature, our change in heat capacity can be zero.
sim.set_model_params(n=2, dG=-10, dH=-10, dCp=0)

# It's always a good idea to make an initial starting condition fit plot. Don't forget to run the simulator first!
sim.run()
sim.make_plots(hardcopy=True,hardcopyprefix="initial_")

# The actual parameter optimization, i.e. fitting, is carried out by an ITCFit object. We instantiate one and tell it what simulation to optimize:
fit = ITCFit( sim, verbose=True )

# We tell the fitter's optimize method what parameters we'd like to optimize.
# The optimize method will return the optimized parameters and the final chi-square goodness-of-fit. The latter we can use later, e.g. by writing to a file, or using as the starting conditions for a different fit.
optimized_params,chisq = fit.optimize( params=['n','dG','dH'])

# Now, let's add some more experimental data acquired at different experimental temperatures.
sim.add_experiment_file('chapter_2_example_288K.txt', skip=[0])
sim.add_experiment_file('chapter_2_example_308K.txt', skip=[0])

# Update our model with the parameters we optimized with just the fit to the 298K data:
sim.set_model_params( **optimized_params )

# Now that we have multiple temperatures, we'll likely need to take into account the potential effects a change in heat capacity can have on the affinity and enthalpy of binding.
# We can see how much of an effect ignoring a change in heat capacity has by simulating the model using the parameters optimized at 298K at the other experimental temperatures.
sim.run()
sim.make_plots(hardcopy=True,hardcopyprefix="dCp=0_")

# Those look pretty bad, the binding enthalpy especially is significantly overestimated at low temps, and underestimated at high temperatures.
# Note that the green line in the fit plots is the spline itcsimlib uses to estimate noise in the experimental data.

# We can use the two additional experimental temperatures, however, to determine dCp:
optimized_dCp,chisq = fit.optimize( params=['dCp'])

# With the right change in heat capacity, fits are much improved.
sim.set_model_params( **optimized_dCp )
sim.run()
sim.make_plots(hardcopy=True,hardcopyprefix="optimized_")

# We can also write the parameters to a file in order to have a permanent copy:
from utilities			import write_params_to_file
write_params_to_file( "chapter_2_parameters.txt", params=sim.get_model_params() )

# Don't forget to call this at the end of each script:
sim.done()

# Congrats. You've reached the end of Chapter 2. In the next chapter we'll cover error estimation and gridding model parameters.
