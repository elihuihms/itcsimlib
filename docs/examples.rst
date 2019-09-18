:orphan:

.. _itcsimlib-examples:

Example Scripts
---------------

These examples provide a number of real-life applications of itcsimlib or show off specialized usage that may not be covered in the :ref:`itcsimlib-tutorial`. The author can provide no assurances regarding their quality or suitability.

Note that the `examples/data` directory contains the six experimental TRAP+Tryptophan titrations used in the author’s 2017 Biophysical Journal paper, and provide a good example of the file format itcsimlib uses.

`dilution.py <https://github.com/elihuihms/itcsimlib/blob/master/examples/dilution.py>`_
	Provides an example of a brute force method to estimate heat-of-dilutions in the case of poor baseline control. Should be used as a last resort, and not recommended for poorly-constrained models.

`drakon <https://github.com/elihuihms/itcsimlib/blob/master/examples/drakon/Readme.ipynb>`_
	Demonstrates how the graphical programming language DRAKON (https://drakonhub.com/) may be used to construct Ising-type binding models for use in itcsimlib.

`estimate.py <https://github.com/elihuihms/itcsimlib/blob/master/examples/estimate.py>`_
	Compares the two different methods for estimating fit parameter uncertainties in itcsimlib; critical chi-square or bootstrapping (jackknife).

`hill.py <https://github.com/elihuihms/itcsimlib/blob/master/examples/hill.py>`_
	Generates a hill-like binding isotherm by expanding a model to write the free ligand and lattice saturation (theta) at each titration point.

`ising_plotter.py <https://github.com/elihuihms/itcsimlib/blob/master/examples/ising_plotter.py>`_
	Shows off itcsimlib's ability to make figures enumerating the different lattice + ligand configurations for linear and circular lattices, as well as their energetic degeneracy at different model parameters.

`mass_spec.ipynb <https://github.com/elihuihms/itcsimlib/blob/master/examples/mass_spec.ipynb>`_
	Extends an Ising model to simulate and fit mass spectrometry data, as well as generate useful population plots. 

`neighbors.py <https://github.com/elihuihms/itcsimlib/blob/master/examples/neighbors.py>`_
	Provides an example of how a model may be expanded in order to provide statistics about the number of sites that are empty, or have 0, 1, or 2 neighbors in the system lattices.

`populations.py <https://github.com/elihuihms/itcsimlib/blob/master/examples/populations.py>`_
	Extends a model to provide data for comparing the best-fit populations/abundancies of different lattice+ligand stoichiometries at different experimental temperatures.

If you'd like to contribute your own example of itcsimlib usage, please contact the author at mail@elihuihms.com.
