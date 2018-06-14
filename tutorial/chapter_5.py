#
# Chapter 5 : Using the GUI
#

from itcsimlib import *
from itcsimlib.manipulator import Manipulator

# In this last chapter, we'll utilize the itcsimlib Manipulator object, which allows the user to modify model parameters in real-time to see their effect on fits.

sim = ITCSim(T0=298.15,units="kcal",verbose=True)

# Let's add four experiments at different concentrations and temperatures:
sim.add_experiment_synthetic(
	T=298.15,
	V0=1416.6,
	injections=[5.0]*50,
	Cell={"Macromolecule":0.5E-6},
	Syringe={"Ligand":10E-6},
	noise=0.5,
	title='Test_Experiment_1')

sim.add_experiment_synthetic(
	T=298.15,
	V0=1416.6,
	injections=[5.0]*50,
	Cell={"Macromolecule":1E-6},
	Syringe={"Ligand":20E-6},
	noise=0.5,
	title='Test_Experiment_2')
	
sim.add_experiment_synthetic(
	T=278.15,
	V0=1416.6,
	injections=[5.0]*50,
	Cell={"Macromolecule":0.5E-6},
	Syringe={"Ligand":10E-6},
	noise=0.5,
	title='Test_Experiment_3')

sim.add_experiment_synthetic(
	T=278.15,
	V0=1416.6,
	injections=[5.0]*50,
	Cell={"Macromolecule":1E-6},
	Syringe={"Ligand":20E-6},
	noise=0.5,
	title='Test_Experiment_4')

# Let's use a simple two-site model, and give it some starting conditions:
from itcsimlib.model_independent import NModes
sim.set_model( NModes(modes=2) )
sim.set_model_params( n1=1,n2=1,dG1=-10,dG2=-11,dH1=-10,dH2=-10,dCp1=-3,dCp2=-3 )
sim.run()

# To change model parameters in real-time, simply instantiate a Manipulator object:
Manipulator( sim )

# When you run this script, move the sliders around to see how the different parameter values affect the fits.
# If you move a slider to the end of the scale, it will automatically rescale and permit a greater range.
# You can also type precise desired values into the entry field near each parameter name.

# This concludes the tutorial. If you have suggestions or comments, please contact the developer at mail@elihuihms.com!
sim.done()
