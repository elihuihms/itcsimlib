#
# Introduction
#

First, it's important to realize that itcsimlib is not an application, it's a collection of routines to model, fit, and otherwise understand isothermal titration calorimetry (ITC) data. As such, to use itcsimlib you will be writing scripts in Python that make use of itcsimlib classes, much in the same way as structural biologists use other programmatic tools like NMRPipe, XPLOR, etc.

If this scares you, there are several programs out there that let you get a decent handle on ITC binding data through a graphical interface. However, if you need to get the absolute most information out of your experimental data (and who doesn't, really?), then itcsimlib is probably what you need. 

If you're not familiar with Python specifically, don't worry! The provided tutorial scripts (in the tutorial folder) should be relatively self explanatory. However, if you have absolutely no programming background at all, you may want to try out one of the many introduction to Python tutorials available on the web before jumping in ab initio.

Keep in mind that to access the true power of itcsimlib, you'll need to eventually write your own models, and that involves writing code. If you discover you're in this camp, then itcsimlib will save you a lot of time dealing with all of the non-exciting and non-modeling aspects of model evaluation and allow you to get right to testing actual hypotheses.

Lastly, if at any point you need help, have constructive criticism, or wish to contribute some of your own ideas and/or code to itcsim, please drop the author a line at elihuihms@elihuihms.com.

#
# Installing itcsimlib
#

The easiest way to get up and running is to just download (or checkout) itcsimlib from https://github.com/elihuihms/itcsimlib, and then place the itcsimlib directory (the one that contains __init__.py and all the other python sources) in a folder along with your experimental data. You can then easily reference both your experimental data and itcsimlib in your analysis scripts as demonstrated in the tutorial chapters.

If you want, you can install itcsimlib somewhere of *your* choice, and then modify your PYTHONPATH environmental variable to point to the installation. If you don't know what an environmental variable is, or don't feel like editing your shell profile config files, just take the easy installation approach.

If you're used to python packages, you'll notice that I have provided a setup.py. However, I wouldn't use it under most circumstances, for the aforementioned reason that the real power of itcsimlib is writing your own models, and "setup.py install" can stick it somewhere relatively ugly like "/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/". An important exception to this would be if you're deploying itcsimlib on a HPC cluster for a model you've already got running.

Finally, because the TRAP/Tryptophan binding models that were the impetus behind this project are written in C, I've given users the option of using the standard distutils C compiler calls ("setup.py build", if you've got a C compiler available), or traditional configure/make scripts (see "model_trap" in the itcsimlib directory). Keep in mind that compiling these models require the excellent Gnu scientific library (GSL, https://www.gnu.org/software/gsl/). Speaking of requirements...

#
# Requirements
#

Basic itcsimlib functionality requires the scipy and matplotlib libraries. You can install these piecemeal, or better yet, get a fresh install of Python and tons of useful Python libraries pre-packaged together by installing the Anaconda Python stack, available at https://www.continuum.io/downloads. The Enthought Python stack also works well: https://www.enthought.com/products/epd/

#
# Acknowledging itcsimlib
#

It is my sincere hope that itcsimlib assists your research efforts. If it does, please cite the following paper in your manuscripts: ???

Also, if you have developed your own models, please consider contacting the maintainer (elihuihms@elihuihms.com) to add your model to the official itcsimlib codebase.