.. CactusTool documentation master file, created by
   sphinx-quickstart on Tue Dec  3 13:08:52 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to CactusTool's documentation!
======================================
CactusTool package provides tools to analysis and visualization `CACTUS <http://www.cactuscode.org>`_ simulations. 

Cactus itself cannot visualize the data and requires third-party tools. Currently, the better software includes `Visit <https://wci.llnl.gov/simulation/computer-codes/visit/>`_, `SimulationTools <https://simulationtools.org>`_  and `gnuplot <http://www.gnuplot.info>`_. In addition to these better tools, `PostCactus <https://bitbucket.org/DrWhat/pycactuset>`_ and `PyCactus <https://bitbucket.org/knarrff/pycactus>`_  are also a good choice. Both are pure Python package. Personally like it very much. Python is a powerful tool and there are very rich third-party libraries, such as `re <https://docs.python.org/3/library/re.html>`_, `pandas <https://pandas.pydata.org>`_, `matplotlib <https://matplotlib.org/index.html>`_, and `HDF5 <https://docs.h5py.org/en/stable/>`_. Most gravitational wave data processing are based on python. With these third-party libraries, you can do more amazing things besides visualization. `Zachariah B. Etienne <http://astro.phys.wvu.edu/zetienne/>`_ wrote a outputC module. I think it can be used on Cactus in the future. However, PostCactus and PyCactus are all based on python2, which is not maintained now. Based on a comprehensive consideration of all aspects, I decided to write this library. CactusTool code benefits from PostCactus and PyCactus code. The software framework follows PostCactus. However, a more comprehensive manual is provided. Cactus will be continuously updated and will be more powerful than the former.

CactusTool is not a powerful 3D animation tool. CactusTool generates rich static picture with matplotlib. Although matplotlib can do this, it will make the code more complicated and more inconvenient to use. In my opinion, animation is more convenient to show the simulation and easier to judge whether the simulation evolution process is normal, but it is not helpful for the analysis of the data. This aspect can be made up by `Jupyter notebook <https://jupyter.org>`_. Jupyter Notebook allows you to create and share documents that contain live code, equations, visualizations and narrative text. It will be more convenient to show results. More than this, `ipywidgets <https://ipywidgets.readthedocs.io/en/latest/>`_ allows us to turn Jupyter Notebooks from static documents into interactive dashboards, perfect for exploring and visualizing data, i.e., **provide animation in another way and be interactive.** With `RISE <https://rise.readthedocs.io/en/maint-5.6/>`_, you can instantly turn your jupyter notebook into a slideshow. ipywidgets interactive dashboards will not be part of CactusTool. Everyone's needs are different and it is difficult to unify. This manual will provide some ipywidgets examples and template.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Note/install
   Tutorial/main
   API/CactusTool

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Built by `Yu Liu <https://yuliumt.github.io>`_