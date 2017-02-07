#
# Chapter 1 : The Basics
#

# Every line that starts with a pound sign (#) is a comment, and are not interpreted by itcsimlib.
# When writing your own analysis/simulation scripts, it's a good idea to explain your reasoning at each step with comments, both for the benefit of yourself later on and for those attempting to understand your work.

# Each itcsim analysis script starts, logically, by importing itcsimlib. If the "itcsimlib" directory is in the same folder as your python script, or already added to your environment's $PATH, this is pretty simple:
from itcsimlib import *

# That was pretty easy. Now let's actually create an ITC simulator by instantiating an ITCSim object:
sim = ITCSim(T0=298.15,units="kcal",verbose=True)

# The arguments necessary to create the simulator should be intuitively straightforward.
# Because many thermodynamic parameters vary with temperature, you'll need to specify a reference temperature (the "T0" argument), which becomes especially important when you're working with data acquired at multiple experimental temperatures.
# The units argument is also intuitive, and can be "J" or "kJ" for Joules or kilojoules, "cal", or "kcal" for calories and kilocalories.
# Lastly, the verbose argument tells the simulator whether or not to print extra information when you run the script. You probably want that info, so here the argument is set to True

# Next up, let's add a theoretical ITC experiment:
sim.add_experiment_synthetic(
	T=298.15,
	V0=1416.6,
	injections=[5.0]*50,
	Cell={"Macromolecule":1E-6},
	Syringe={"Ligand":30E-6},
	noise=1.0,
	title='Test_Experiment_1')

# This theoretical (or in itcsimlib-speak, "synthetic") experiment takes place at 25C (298.15K), has a starting cell volume of 1416.6uL, and is comprised of fifty 5uL injections.
# The cell contains our binding macromolecule at a starting concentration of 1uMol concentration, while our syringe contains our ligand at 30uMol. Here, we've chosen to add 1kcal/mol simulated noise. We've also given the experiment the catchy name "Test_Experiment_1".

# To model the binding curve, you'll need specify a binding model. Let's use an independent two-mode description of binding. We'll need to import it first:
from itcsimlib.model_independent import NModes

# We make a new instance of this model, and tell it we'd like to have two independent modes via the modes=2 argument.
# We then pass this new model to our simulator's set_model() method, which can be done in conjunction with the previous command in Python in one line like this...
sim.set_model( NModes(modes=2) )

# ...or in two separate steps on different lines, like this:
my_model = NModes(modes=2)
sim.set_model( my_model )

# Models have parameters that dictate their behavior. If we want to see the details of the model, all we have to do is treat it as a string, in this case we just want to print the model info to our console when we run this script.
print sim.get_model()

# Once our model is set, we'll need to set it's parameters:
sim.set_model_params(n1=1, n2=3, dG1=-11, dG2=-10, dH1=-15, dH2=+10)
# Remember, we specified that our units are in kcal back when we created the simulator.
# With these arguments, we've said that one binding mode occurs with a stoichiometry of 1:1 ligand to macromolecule, and with a free energy change of -11kcal/mol, and an exothermic change in enthalpy of -15kcal/mol
# The second binding mode occurs with a stoichiometry of 3 ligands per molecule, a slightly weaker -10 kca/mol free energy change, and an endothermic change in enthalpy of +10 kcal/mol.

# Now that we've updated the model parameters, let's verify that they've been set.
print sim.get_model()

# Great, but just setting the parameters doesn't really do anything. Now we need to run the simulator, which will populate our synthetic experiment. itcsimlib makes this pretty easy:
sim.run()

# Hopefully, our simulation goes off without a hitch. If it does, all of the experiments that we've assigned will now contain theoretical binding data (in this case, we have just the one experiment). Let's take a look at the simulated curves:
sim.make_plots(hardcopy=True)

# By setting hardcopy to True, itcsimlib will write plots to a file instead of just temporarily displaying them in a plot window.

# When we're all done, we need to tell the simulator so it can stop any processing threads and clean up.
sim.done()

# Congratulations. You've reached the end of Chapter 1. In the next chapter we'll cover fitting experimental data with a model.
