
## itcsimlib : A statistical thermodynamics isothermal titration calorimeter simulator

# Introduction

itcsimlib is python module that focuses (although not exclusively) on the use of statistical thermodynamics models to simulate, fit, and otherwise interpret isothermal titration calorimetry (ITC) data. Because statistical thermodynamics models individually simulate the prevalence of lattice+ligand configurations, itcsimlib can also be readily extended to fit mass spectrometry data.

Note: itcsimlib doesn't possess a GUI. You will need to write scripts in Python that make use of itcsimlib classes. Users of XPLOR and other programmatic analysis tools will find this quite familiar.

If you're not familiar with Python specifically, don't worry! The provided scripts (in the tutorial and examples directories) should give you a good introduction. However, if you have absolutely no programming background at all, you may want to try out one of the many introduction to Python tutorials available on the web.

Although itcsimlib comes with several binding models that may fit your data well (including both Ising and non-Ising based ones), to really leverage itcsimlib to learn the most about your particular system, you'll probably want to write your own models. Itcsimlib will save you a lot of time dealing with all of the non-exciting and non-modeling aspects of model evaluation.

If at any point you need help, have constructive criticism, or wish to contribute some of your own ideas or models to itcsimlib, please contact the author at mail@elihuihms.com.

# Requirements

Itcsimlib requires the scipy (v0.11 or above) and matplotlib modules (v1.3 or above). For automated generation of symbolic partition functions, you'll also need the sympy module. You can install all of these piecemeal, or better yet, get a complete Python environment, including most of the commonly-used modules pre-packaged together by installing the Anaconda Python stack, available at https://www.continuum.io/downloads.

# Installing itcsimlib

You have two options:

1. Use the included distutils setup.py, i.e. "python setup.py install".

2. Move the itcsimlib directory (the one that contains __init__.py and all the other python sources) to a directory along with your experimental data. You can edit your shell's PYTHONPATH environmental variable to point to a directory that contains itcsimlib, thus avoiding having to maintain multiple copies.

Note that if you want to compile the optional TRAP+Tryptophan binding models that are written in C, you'll either want to run the setup script with the "--build-c-models" flag, or use the traditional configure/make scripts (see "model_trap" in the itcsimlib directory). Keep in mind that compiling these models will additionally require the GNU scientific library: https://www.gnu.org/software/gsl/.

# Acknowledging itcsimlib

It is my sincere hope that itcsimlib assists your research efforts. If it does, please consider citing the following paper in your manuscript: Ihms, Elihu C. et al. "Mechanistic models fit to variable temperature calorimetric data provide insights into cooperativity.” Biophysical Journal (2017)
