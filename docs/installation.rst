:orphan:

.. _itcsimlib-installation:

Installation
------------
If you do not already have a Python scientific computing environment set up, we recommend you follow `these steps <https://python-for-scientists.readthedocs.io/en/latest/_pages/install_python.html>`_ to set it up. Once configured, 

.. sourcecode:: bash

    git clone https://github.com/elihuihms/itcsimlib.git
    cd itcsimlib
    python3 setup.py develop

If you want to compile the optional TRAP + Tryptophan models that are written in C, you'll either want to run the setup script with the `--build-c-models` flag, or use the traditional configure/make scripts (see `src/model_trap`). Keep in mind that these models require the GNU scientific library: https://www.gnu.org/software/gsl/.

Dependencies
------------
 + `numpy <http://www.numpy.org/>`_
 + `scipy <https://www.scipy.org/>`_
 + `matplotlib <http://matplotlib.org/>`_
 + `sympy <http://www.sympy.org/>`_
 + `pyx <hhttps://pypi.org/project/PyX/>`_
 + `Tkinter <https://docs.python.org/3/library/tk.html>`_ (only for GUI elements)
 + `jupyter <https://jupyter.org/>`_ (useful for scripting and documentation)