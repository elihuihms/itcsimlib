import sympy

from itcsimlib import *

sim = ITCSim(T0=298.15,units="kcal",verbose=True)

sim.add_experiment_synthetic(
	T=298.15,
	V0=1416.6,
	injections=[5.0]*50,
	Cell={"Lattice":1E-6},
	Syringe={"Ligand":40E-6},
	title='drakon_004')

from koshland_1966 import Model

model = Model()
model.precision = 1E-15

sim.set_model( model )

print sim.model
print sympy.latex(sim.model.get_partition_function())

sim.set_model_params(-10,0,0,-1,0,0) #dG_st, dG_AB, dG_BB, #dH_st, dH_AB, dH_BB
sim.run()

sim.set_model_params(-10,-1,0,-1,0,0)
sim.run()

#sim.make_plots(hardcopy=True)
sim.make_plots()

sim.done()