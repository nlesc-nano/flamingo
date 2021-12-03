##########
Change Log
##########

All notable changes to this project will be documented in this file.
This project adheres to `Semantic Versioning <http://semver.org/>`_.


0.3.0 [03/12/2021]
******************
* Release flamingo on pypi (#64)
* Run all the filters in parallel (#38)
* Filter molecules with a single functional group (#43)
* Add interface to cosmo-rs (#5)


0.2.1 [14/01/2021]
******************
Change
-----
* Use all the available CPU to compute bulkiness with CAT by calling the `imap_unordered Pool's method <https://docs.python.org/3/library/multiprocessing.html#multiprocessing.pool.Pool.imap_unordered>`_.
* Remove the `batch_size` input parameter and fix it to 1000.


0.2.0 [03/11/2020]
******************

Added
-----
* Workflow to compute properties using `CAT <https://github.com/nlesc-nano/CAT>`_


0.1.0 [14/10/2020]
******************

Added
-----
* Move `swan functionality to compute and filter properties<https://github.com/nlesc-nano/swan/issues/44>`_ to **flamingo**
