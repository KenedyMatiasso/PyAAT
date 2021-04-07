PyAAT
########

.. |python| image:: https://img.shields.io/badge/Made%20with-Python-1f425f.svg
   :target: https://www.python.org/

.. |sphinx| image:: https://img.shields.io/badge/Made%20with-Sphinx-1f425f.svg
   :target: https://www.sphinx-doc.org/

.. |Open Source? Yes!| image:: https://badgen.net/badge/Open%20Source%20%3F/Yes%21/blue?icon=github
   :target: https://github.com/Naereen/badges/


.. |PyPI status| image:: https://img.shields.io/pypi/status/ansicolortags.svg
   :target: https://pypi.python.org/pypi/ansicolortags/

.. |license| image:: https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square
   :target: https://github.com/KenedyMatiasso/PyAAT/blob/main/LICENSE

.. |docs| image:: https://img.shields.io/badge/docs-latest-brightgreen.svg?style=flat-square
   :target: https://pypi.org/project/PyAAT/#description
   
.. |mybinder| image:: https://mybinder.org/badge.svg
   :target: https://mybinder.org/v2/gh/KenedyMatiasso/PyAAT/3e355dbb77a3db0cfaef2a0a941bc9f79cfb32ed


|python| |sphinx| |Open Source? Yes!| |PyPI status| |license| |mybinder| |docs|

:Name: Python Aerospace Analysis Toolbox- PyAAT
:Author: Kenedy Matiasso Portella
:Website: https://github.com/KenedyMatiasso/PyAAT
:Version: 0.0.dev5

Python Aerospace Analysis Toolbox (PyAAT) is a open-Source python-based toolbox for modeling, analysis and simulation of aerospace systems.

About PyAAT
**********************
Python Aerospace Analysis Toolbox is a python package aimed to support aerospace engineers along the design process.
It includes a series of tools and models for modeling, analysis and simulation of aeronautic and space systems, it includes:

Atmosphere
===========
* International Standard Atmosphere (ISA)
* Committee on Extension to the Standard Atmosphere (COESA 1976)
* Simplified Model
* Check out an example using `atmosphere models <https://hub.gke2.mybinder.org/user/kenedymatiasso-pyaat-ldhiu3mo/notebooks/examples/atmosphere_example.ipynb>`_.

Gravity
========
* Vertical Constant
* Newton Gravity
* High order gravity model (represented by spherical harmonics)
* Check out an example using `gravity models <https://hub.gke2.mybinder.org/user/kenedymatiasso-pyaat-ldhiu3mo/notebooks/examples/gravity_example.ipynb>`_.

Wind
=====
* Not included yet

Aerodynamics
=============
* Based on stability derivatives

Propulsion
===========
* Turbo-propeller
* Turbo-fan
* Turbo-jet
* Check out an example using  `propulsion model <https://hub.gke2.mybinder.org/user/kenedymatiasso-pyaat-ldhiu3mo/notebooks/examples/propulsion_example.ipynb>`_.


Body dynamics (rigid body only)
=================================
* 6 Degrees of fredom


Instalation
**********************

The most recent stable version of PyAAT can be installed directly from pip repository using comand:

$ ``pip install pyaat``

If you want to install the most recent functional version (no guarantee of stability), you can install directly from our github repository using:

$  ``pip install git+https://github.com/KenedyMatiasso/PyAAT``
