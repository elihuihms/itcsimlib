from itcsimlib.model_independent import OneMode
from itcsimlib.itc_fit import ITCFit
from itcsimlib.itc_grid import ITCGrid

from .itc_sim import TestITCSIM


class TestITCGrid(TestITCSIM):
	def setUp(self):
		super().setUp()

		self.sim.add_experiment_file(self.get_data('base_1.txt'))
		self.sim.set_model( OneMode() )
		self.sim.set_model_params(n=1.805,dG=-10.94,dH=-11.75,dCp=0.0)		
		self.fit = ITCFit( self.sim, method='simplex', method_args={"maxiter":1} )

	def test_2d_grid(self):
		grid = ITCGrid( self.fit )
		grid.add_axis(param='dG',start=-12,stop=-10,steps=2)
		grid.add_axis(param='dH',start=-13,stop=-11,steps=2)
		self.assertIsNotNone( grid.optimize(params="n") )
		
	def test_defined_grid(self):
		grid = ITCGrid( self.fit )
		grid.define_axis(param='dG',points=(-11.5,-10.94,-10.5))
		grid.define_axis(param='dH',points=(-12.25,-11.75,-11.50))
		self.assertIsNotNone( grid.optimize(params="n") )
	