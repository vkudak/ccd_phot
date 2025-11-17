Installation
============

.. _installation:

How to install the scripts
---------------------------

You can just download the ``ccd_phot`` from GitHub page or clone whole repository.
You need only Python 3 (3.10 or higher) to be installed in your system.

.. note::
    To download use this link https://github.com/vkudak/ccd_phot/archive/refs/heads/master.zip

To clone the project:

.. code:: console

    $ git clone https://github.com/vkudak/ccd_phot.git



In such case you can always update project and have latest code with fixes and new features.

To update project:

.. code:: console

   $ git pull

Than navigate to you ``CCD_PHOT`` folder and install all requirements


It is highly recommended to create new virtual environment
and install all dependencies in new ``venv``.

.. code:: console

   $ cd ccd_phot
   $ python3 -m venv venv
   $ pip install -r requirements.txt

Or install with system python

.. code:: console

   $ cd ccd_phot
   $ pip install -r requirements.txt



After that read this section :ref:`start_scripts`
