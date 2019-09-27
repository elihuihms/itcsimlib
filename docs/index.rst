itcsimlib
=========

itcsimlib uses statistical thermodynamics (ST) models to simulate and fit binding data, in particular those obtained from isothermal titration calorimetry (ITC). Because statistical thermodynamics models can calculate the prevalence of individual lattice+ligand configurations, itcsimlib can be readily adapted to fit mass spectrometry data as well.

itcsimlib also provides a number of tools for understanding / analyzing Ising ST models that may be of general interest to those studing complex binding processes. These include automated generation of partition functions and figures depicting lattice + ligand configuration with their degeneracies.

Getting Started
===============

* :ref:`itcsimlib-installation` instructions for required dependencies and optional compilation. 

* :ref:`itcsimlib-tutorial` that demonstrate basic itcsimlib usage and script development.

* :ref:`itcsimlib-examples` for additional functionality, tools, and advanced topics.

Reference
===================

itcsimlib doesn't possess a GUI. You will need to write scripts in Python that make use of the itcsimlib classes described below. 

Essentials
----------

.. autosummary::
   :toctree: modules

   itcsimlib.itc_sim
   itcsimlib.itc_experiment
   itcsimlib.model_independent
   
Fitting
-------

.. autosummary::
   :toctree: modules

   itcsimlib.itc_fit
   itcsimlib.itc_grid
   itcsimlib.manipulator

Model Development
-----------------

.. autosummary::
   :toctree: modules

   itcsimlib.itc_model
   itcsimlib.model_ising
   itcsimlib.model_drakon
   itcsimlib.model_trap
   itcsimlib.mass_spec

Nuts and Bolts
--------------

.. autosummary::
   :toctree: modules

   itcsimlib.thermo
   itcsimlib.utilities
   itcsimlib.itc_calc

Indices
==================

* :ref:`genindex`
