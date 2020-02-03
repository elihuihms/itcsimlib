import pytc

from .itc_sim import ITCSim


class ISL2pytcExperiment(pytc.experiments.base.BaseITCExperiment):

	class ISL2pytcModel(pytc.indiv_models.base.ITCModel):
		
		def __init__(self,**kwargs):
			self._isl_simulator = kwargs.pop("simulator")
			self._isl_model_params = kwargs.pop("model_params")
			
			assert(len(self._isl_simulator.experiments) == 1)
			self._isl_experiment = self._isl_simulator.experiments[0]
			
			super().__init__(**kwargs)

		@property
		def dQ(self):
			self._isl_simulator.set_model_params( **{p:self._params[p].value for p in self._isl_model_params} )
			self._isl_simulator.run(writeback=True)
			return self._isl_experiment.dQ_fit

		def param_definition():
			"""Param definitions are built on the fly using _initialize_param below."""
			pass

		def _initialize_param(self):
			super()._initialize_param(
				param_names=self._isl_model_params[:],
				param_guesses=[self._isl_simulator.get_model_param(p) for p in self._isl_model_params])

		def update_values(self,param_values):
			super().update_values(param_values)

	def __init__(self,path,model,model_params=None,reverse=False,experiment_kwargs={}):
		tmp_sim = ITCSim()
		tmp_exp = tmp_sim.add_experiment_file(path,**experiment_kwargs)
		tmp_sim.done()

		self._isl_simulator = ITCSim(
			T0 = tmp_exp.T,
			units = model.units.replace("/mol",""),
			threads = 1
		)
		self._isl_simulator.set_model(model)
		self._isl_experiment = self._isl_simulator.add_experiment_file(path,**experiment_kwargs)

		if model_params is None:
			model_params = list(self._isl_simulator.get_model_params().keys())

		if reverse:
			self._isl_experiment.change_component_name( self._isl_experiment.syringeRef, model.lattice_name )	
			self._isl_experiment.change_component_name( self._isl_experiment.cellRef, model.ligand_name )	
		else:
			self._isl_experiment.change_component_name( self._isl_experiment.cellRef, model.lattice_name )	
			self._isl_experiment.change_component_name( self._isl_experiment.syringeRef, model.ligand_name )	

		if self._isl_experiment.skip:
			shot_start = min(self._isl_experiment.skip)
		else:
			shot_start = 0

		pytc_experiment_kwargs = {
			"dh_file":		str(self._isl_experiment.title),
			"model": 		self.ISL2pytcModel,
			"shot_start":	shot_start,
			"units":		"%s/mol"%(self._isl_simulator.units),
			"uncertainty":	0.0,
			"simulator":	self._isl_simulator,
			"model_params":	model_params,
			"S_syringe":	self._isl_experiment.Syringe.get(self._isl_experiment.cellRef,0.0),
			"T_cell":		self._isl_experiment.Cell.get(self._isl_experiment.syringeRef,0.0),
		}

		super().__init__(**pytc_experiment_kwargs)

	@property
	def isl_simulator(self):
		return self._isl_simulator

	@property
	def isl_model(self):
		return self._isl_simulator.model

	@property
	def isl_experiment(self):
		return self._isl_experiment
	
	def _read_heats_file(self):
		self.temperature = self._isl_experiment.T
		self.stationary_cell_conc = self._isl_experiment.Cell[self._isl_experiment.cellRef]
		self.titrant_syringe_conc = self._isl_experiment.Syringe[self._isl_experiment.syringeRef]
		self.cell_volume = self._isl_experiment.V0

		self._shots = self._isl_experiment.injections
		self._heats = self._isl_experiment.dQ_exp
		if self._isl_experiment.dQ_err is None:
			self._heats_stdev = [0.0] * self._isl_experiment.npoints
		else:
			self._heats_stdev = self._isl_experiment.dQ_err