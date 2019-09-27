:orphan:

.. _itcsimlib-tutorial:

Tutorials
---------

Tutorials for itcsimlib are provided as five Jupyter notebooks (https://jupyter.org/) and cover a range of different topics:

`Chapter 1 <https://github.com/elihuihms/itcsimlib/blob/python3/tutorial/chapter_1.ipynb>`_
	Provides a basic framework for an itcsimlib analysis script, and demonstrates how to use the simulation, model, and theoretical experiment classes (:py:class:`~itcsimlib.itc_sim.ITCSim`, :py:class:`~itcsimlib.itc_model.ITCModel`, :py:class:`~itcsimlib.itc_experiment.ITCExperimentSynthetic`, respectively).

`Chapter 2 <https://github.com/elihuihms/itcsimlib/blob/python3/tutorial/chapter_2.ipynb>`_
	Demonstrates how to read in empirical data to create experiments using the :py:class:`~itcsimlib.itc_experiment.ITCExperiment` class, use an indpendent-modes model such as :py:class:`~itcsimlib.model_independent.NModes` to fit such data, and how to obtain fit parameter confidence intervals.

`Chapter 3 <https://github.com/elihuihms/itcsimlib/blob/python3/tutorial/chapter_3.ipynb>`_
	Describes how to use the :py:class:`~itcsimlib.itc_grid.ITCGrid` class to either restrain model parameters to discretely sampled values during fitting, or to sample starting conditions prior to optimization using multidimensional grids.

`Chapter 4 <https://github.com/elihuihms/itcsimlib/blob/python3/tutorial/chapter_4.ipynb>`_
	Covers the creation of Ising-type binding models based on the :py:class:`~itcsimlib.model_ising.Ising` model class, using the arrangement of occupied sites in the lattice to determine binding free energy and enthalpies.

`Chapter 5 <https://github.com/elihuihms/itcsimlib/blob/python3/tutorial/chapter_5.ipynb>`_
	Demonstrates the use of the :py:class:`~itcsimlib.manipulator.Manipulator` object, which enables the user to change model parameters in real time and see the effect on model fits.

If a given topic of interest isn't described here, please take a look in the :ref:`itcsimlib-examples` section, where more advanced topics and additional tools provided by itcsimlib are demonstrated.

As always, if you need help or have constructive criticism, please drop the author a line at mail@elihuihms.com.
