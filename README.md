
## itcsimlib : A statistical thermodynamics isothermal titration calorimeter simulator

# Introduction

itcsimlib is python module that uses statistical thermodynamics models to simulate and fit binding data, in particular those obtained from isothermal titration calorimetry (ITC). Because statistical thermodynamics models can calculate the prevalence of individual lattice+ligand configurations, itcsimlib can be readily adapted to fit mass spectrometry data.

Note: itcsimlib doesn't possess a GUI. You will need to write scripts in Python that make use of itcsimlib classes, in a manner that users of [XPLOR-NIH](https://nmr.cit.nih.gov/xplor-nih/) and other programmatic analysis tools will find quite familiar. 

Full instructions for getting started via tutorials, API documentation, and advanced examples can be found at https://itcsimlib.readthedocs.io/

Although itcsimlib comes with several binding models that may fit your data, you'll eventually want to write your own. You may want to try out [itcsimlib-blockly](https://github.com/elihuihms/itcsimlib-blockly), which uses Google's [Blockly](https://developers.google.com/blockly/) visual programming language to build itcsimlib statistical thermodynamics models.

If at any point you need help, have constructive criticism, or wish to contribute some of your own ideas or models to itcsimlib, please contact the author at mail@elihuihms.com.

# Requirements

See `requirements.txt` for all required dependencies. Note that Tkinter is necessary only if you want to use the Manipulator class to visualize the effects of changing model parameters in realtime. Ultimately, the easiest and quickest route to get started is probably to install the Anaconda Python stack available at https://www.continuum.io/downloads.

# Installing itcsimlib

```
cd itcsimlib
pip install -r requirements.txt
python setup.py install
```

If you want to compile the optional TRAP + Tryptophan models that are written in C, you'll either want to run the setup script with the `--build-c-models` flag, or use the traditional configure/make scripts (see `src/model_trap`). Keep in mind that these models require the GNU scientific library: https://www.gnu.org/software/gsl/.

# Acknowledging itcsimlib

It is my sincere hope that itcsimlib assists your research efforts. If it does, please consider citing the following paper in your manuscript: Ihms, Elihu C. et al. "Mechanistic models fit to variable temperature calorimetric data provide insights into cooperativity.” Biophysical Journal (2017) https://doi.org/10.1016/j.bpj.2017.02.031
