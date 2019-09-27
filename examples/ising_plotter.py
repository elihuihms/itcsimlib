import sympy
sympy.init_printing()

from itcsimlib.model_ising import *

circle_model = NonAdditive(nsites=8, circular=1, units='kcal')
circle_model.set_params(dGX=-7,dGY=-8,dGZ=-9)

print("Circular lattice:")
print(sympy.latex(circle_model.get_partition_function()))
circle_model.draw_lattices( "ising_8site_circle")

print("Linear lattice:")
linear_model = NonAdditive(nsites=8, circular=0, units='kcal')
linear_model.set_params(dGX=-7,dGY=-8,dGZ=-9)

print(sympy.latex(linear_model.get_partition_function()))
linear_model.draw_lattices( "ising_8site_linear")
