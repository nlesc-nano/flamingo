.. image:: https://github.com/nlesc-nano/flamingo/workflows/build%20with%20conda/badge.svg
   :target: https://github.com/nlesc-nano/flamingo/actions
.. image:: https://readthedocs.org/projects/flamingo-molecular-properties-calculator/badge/?version=latest
   :target: https://flamingo-molecular-properties-calculator.readthedocs.io/en/latest/?badge=latest
.. image:: https://codecov.io/gh/nlesc-nano/flamingo/branch/master/graph/badge.svg?token=N6CU1B82X0
   :target: https://codecov.io/gh/nlesc-nano/flamingo
.. image:: https://zenodo.org/badge/300545275.svg
   :target: https://zenodo.org/badge/latestdoi/300545275

########
flamingo
########

Compute and filter molecular properties. See `documentation <https://flamingo-molecular-properties-calculator.readthedocs.io/en/latest/>`_.

Installation
============

- Download miniconda for python3: miniconda_.

- Install according to: installConda_.

- Create a new virtual environment using the following commands:

  - ``conda create -n flamingo``

- Activate the new virtual environment

  - ``conda activate flamingo``

To exit the virtual environment type  ``conda deactivate``.


.. _dependecies:

Dependencies installation
-------------------------

- Type in your terminal:

  ``conda activate flamingo``

Using the conda environment the following packages should be installed:


- install RDKit_ and H5PY_:

  - `conda install -y -q -c conda-forge  h5py rdkit`

.. _installation:

Package installation
--------------------
Finally install the package:

- Install **flamingo** using pip:
  - ``pip install nlesc-flamingo``

Now you are ready to use *flamingo*.


Contributing
************

If you want to contribute to the development of flamingo,
have a look at the `contribution guidelines <CONTRIBUTING.rst>`_.

License
*******

Copyright (c) 2020-2021, Netherlands eScience Center

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.



Credits
*******

This package was created with `Cookiecutter <https://github.com/audreyr/cookiecutter>`_ and the `NLeSC/python-template <https://github.com/NLeSC/python-template>`_.

.. _installConda: https://conda.io/projects/conda/en/latest/user-guide/install/index.html
.. _RDKit: https://www.rdkit.org
.. _H5PY: https://www.h5py.org/
.. _miniconda: https://docs.conda.io/en/latest/miniconda.html
