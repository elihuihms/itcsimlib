#
# Chapter 4 : Custom Models
#

from itcsimlib import *

# In this chapter, we'll cover one of the most important parts of itcsimlib: how to use it to develop and test your own models.

# All itcsimlib models must be derived from, or "inherit", either the ITCModel class directly, or a class that is itself a descendent of an ITCModel.

# In this example, we will develop a very simple statistical thermodynamics model in which ligands bind to a macromolecule that has five binding sites. The first three binding events will occur with one free energy and enthalpy change, while the binding of the last two ligands to the macromolecule will have a different free energy and enthalpy change. This could mimic either a positively or negatively-cooperative "switch" type behavior, depending on whether the free energy change for the last two binding events is more or less favorable than the first three binding events.

# We will use a built-in statistical thermodynamics Ising lattice class available in itcsimlib as the basis for our model. This Ising model primitive already inherits the ITCModel class and we import it like so:
from itcsimlib.model_ising import Ising

# This Ising class is not actually a valid model...yet. The critical element of any statistical thermodynamics model is the assignment of energies to the different potential states the model adopts. The Ising class already enumerates all possible configurations, we just need to set the energies for each.

# Our model inherits the Ising class through the following Python convention:
class MyModel(Ising):
	"""This is a simple statistical thermodynamics model describing binding to a five-site macromolecule. The first three ligands bind with one characteristic free energy and enthalpy, while the fourth and fifth bind with a different affinity and enthalpy."""
	# itcsimlib automatically converts any ITCModel's class docstring (seen above) into the description for the model. This is not required by itcsimlib, but is a good documentation practice.
	
	def __init__(self):
		# We use Python's somewhat awkwardly named "__init__" method for initializing the model, and for all itcsimlib models, it must in turn first initialize the parent class it is descended fromt:
		Ising.__init__(self,nsites=5,circular=True)
		# Here, we tell the Ising class the number of sites in our macromolecule, and whether or not the lattice representation of the macromolecule is circular, and from there the Ising class automatically enumerates all possible configurations.
		# In this example, we are not actually using the lattice, but simply the number of ligands already bound to the five-point lattice.
		
		# Next, we need to register the parameters that our model requires. Because we have two distinct binding affinities and enthalpies, our model will have four parameters. If our model needed to be apply across multiple temperatures, we would also need at least one additional parameter describing the change in heat capacity for the system at each binding process.
		
		self.add_parameter( 'dG1',	type='dG',	description='Free energy change upon binding the first, second, or third ligands.' )
		self.add_parameter( 'dH1',	type='dH',	description='Enthalpy change upon binding the first, second, or third ligands.' )
		self.add_parameter( 'dG2',	type='dG',	description='Free energy change upon binding the fourth or fifth ligands.' )
		self.add_parameter( 'dH2',	type='dH',	description='Enthalpy change upon binding the fourth or fifth ligands.' )
		
		# The ITCModel's add_parameter() method must have at least the parameter's name and type (dG, dH, dCp, etc.) as arguments.
		# As with the model's description, an argument's description is only used when you examine the model's info, as demonstrated in tutorial chapter 1. However, describing your model's parameters is quite useful for future reference by you or others.
		
		# That's it for the model initialization. If our model had multiple binding components, we'd need to specify those during initialization as well, but that level of complexity goes beyond a simple tutorial.
		return None
	
	def set_energies(self,T0,T):
		# The Ising parent model expects you to define a set_energies() method, that after execution will populate the model's "self.gibbs[]" and "self.enthalpies[]" attributes with the correct energies for each configuration.
		# The Ising model generates several useful class attributes that our set_energies() method can use. For our five-site lattice, there are 2^5 or 32 potentially energetically-distinct configurations of ligands bound to the lattice. This number of configs is determined by the Ising class during initialization and stored in the "self.nconfigs" attribute.
		
		# During model execution, the parameter values set by either the user or the fitting routines are made available through the model's "self.params[]" attribute dictionary.
		# Also, when your model is evaluated during simulation/fitting, the simulations reference temperature (T0), as well as the experimental temperature (T) are provided as arguments. This is to permit any temperature dependence of your model's parameters.
		
		for i in range(self.nconfigs): # Iterate over each of the configurations
			
			# Another useful model attribute set by the parent Ising model is the "bound[]" list, which is the number of ligands bound to the indexed configuration.
			
			if self.bound[i] == 0: # Reference state. If no ligands bound, free energy and enthalpy are zero.
				self.gibbs[i] = 0.0
				self.enthalpies[i] = 0.0
			
			elif self.bound[i] < 4: # If less than three ligands bound, then the change in the free energy and enthalpic changes for each binding event are the same
				self.gibbs[i] = self.params['dG1'] * self.bound[i]
				self.enthalpies[i] = self.params['dG1'] * self.bound[i]
				
			elif self.bound[i] == 4: # If four ligands are bound, the first three still have the same free energy change as before, but the fourth will bind with a different free  and enthalpic energy changes
				self.gibbs[i] = (self.params['dG1'] * 3) + self.params['dG2'] 
				self.enthalpies[i] = (self.params['dH1'] * 3) + self.params['dH2']
			
			elif self.bound[i] == 5: # And for five bound ligands, the fifth binds with the same free energy and enthalpy changes as the fourth.
				self.gibbs[i] = (self.params['dG1'] * 3) + (self.params['dG2'] * 2)
				self.enthalpies[i] = (self.params['dH1'] * 3) + (self.params['dH2'] * 2)
		
		# That's it! The Ising model class will take care of everything else.
		return None

sim = ITCSim(T0=298.15,units="kcal",verbose=True)

sim.add_experiment_synthetic(
	T=298.15,
	V0=1416.6,
	injections=[5.0]*50,
	Cell={"Lattice":1E-6},
	Syringe={"Ligand":50E-6},
	noise=0.0,
	title='MyModel_Test')

# Assign our custom model
sim.set_model( MyModel() )

# Print the full information for our model, just for fun.
print sim.model

# Set the values for the parameters we defined in the model. Let's first model positive cooperativity by making the free energy change for the fourth and fifth binding steps more favorable:
sim.set_model_params(dG1=-10, dH1=-10, dG2=-12, dH2=-5)

sim.run()

sim.make_plots(hardcopy=True,hardcopyprefix="Positive_")

# Now we will model negative cooperativity by making the last two binding steps less favorable:
sim.set_model_params(dG1=-10, dH1=-10, dG2=-9, dH2=-5)

sim.run()

# Let's make a second plot, and compare the titration curves for our model when it's exhibiting positive vs. negative cooperativity.
sim.make_plots(hardcopy=True,hardcopyprefix="Negative_")

# All finished. In the next chapter, we'll look at how you can get a better handle on your model by changing parameters in real-time and seeing the effects.
sim.done()
