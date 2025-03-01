.. ccd_phot documentation master file, created by
   sphinx-quickstart on Sat Mar  1 21:58:28 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

CCD_PHOT DOCUMENTATION
======================

Python scripts to process photometry obtained on CCD QHY-174M

Installation
============

.. _installation:

How to install the scripts
---------------------------

You can just download the ``ccd_phot`` from GitHub page or clone whole repository

.. code:: console

    $ git clone https://github.com/vkudak/ccd_phot.git

In such case you can always update project using

.. code:: console

   $ git pull

Than you must install all requirements

.. code:: console

   $ pip install -r requirements.txt


.. note::

   It is highly recommended to create new virtual environment
   and install all dependencies in new ``VENV``.

   .. code:: python

       python3 -m venv <myenvpath>
       pip install -r requirements.txt

   Another option -- add ``ccd_phot`` directory path to system environment,
   so you can use scripts right from the command line if all requirements are satisfied.


Use for photometry of RSO
=========================

Prepare FITS files
-------------------
Your RSO photometry FITS files should be in some directory.
They should be named in order of creation, otherwise the script will lose the target.

Avoid names like ``file_1.fits``, ``file_10.fits``, etc.
In this case use ``fits_sort`` option in ``config_sat.ini`` file, see :ref:`sat-config-label`.

Use: ``file_00001.fits``, ``file_00002.fits`` and so on.

FITS files do not need to contain WCS solution, ``sat_phot.py`` script do not use it at all.

Starting from ver 1.0.1 ``sat_phot.py`` script can try to find object automatically on the first FITS file.
If script fails than you need to add additional keywords to the first FITS file header
``OBJX: <float>`` & ``OBJY: <float>`` to tell the script where the target initially is situated on the frame.

.. _sat-config-label:

Prepare sat_config file
-----------------------

Example of ``sat_config.ini`` file can be found in project root directory.
Photometry script will try to find ``config_sat.ini`` in working DIR or some other INI file that will be read as config file.
Config file has few sections (see below).

Most important section is STD section.
Here you should set the system zero-point & coefficient of extinction. Can be obtained in the process of system calibration
(see :ref:`sys-calibration-label`).
Filter (band) and all other parameters are optional.

You can sort your FITS files according to some KEY from FITS headers, if there is some filename issues.

.. code-block::

   [NAME]
   NORAD = 43657
   COSPAR = 18082A
   NAME = COSMOS-2528

   [TLE]
   tle_file = path_to_tle_file

   [STD]
   A = 14.5                # System zero-point
   k = 0.1                 # coefficient of extinction in selected filter
   gate = 10               # Window to search the object (in pixels)
   Filter = R              # Filter name
   Time_format = UT        # UT/JD optional parameter, default value = UT
   min_real_mag = 15       # Min mag that can be obtained on photometry system
   max_center_err = 2.0    # Maximum allowed error in center definition
   fits_sort = DATE-OBS    # Header field to sort files by (without quotes).
                           # Default sorted by FITS name, without this option at all or set None.

   [APERTURE]
   ap = 8
   an_in = 12
   an_out = 16
   saturated = False       # reject saturated data. Default is False

   [PLOT]
   plot_errors = False     # Plot errors with ccd_plot.py script
   auto_plot = True        # Plot graphs when result is ready
   python = python3        # python3, python38, or something else
   plot_command = /path_to_file/ccd_plot.py # abspath to ccd_plot.py

   [SITE]
   Name = OBS_NAME
   lat = <float>
   lon = <float>
   h = <float>

Usage
-----
There is two option of ``sat_phot.py`` to start processing FITS files:

   Start script from its location
       PATH parameter where FITS files are located, must be given

       .. code-block:: console

         python3 sat_phot.py path\to\fits\files

   Start script from FITS files folder
       In this case your path to ``sat_phot.py`` should be added to ENVIRONMENT variables. PATH parameter should be ``.``.

       .. code-block:: console

         sat_phot.py .


In general ``sat_phot.py`` has few
options:

   **[-c or - - config ]**    Specify config file [Optional]

   **[-s or - - start_file]** Specify the name of the FITS file to start from [Optional]. Could be useful if photometry script fails in the middle of the file list. This can happen if the object is too faint. In such case you can use this option to start from required FITS. But first rename result file or it will be overwritten.

   **[path]**               Path to FITS files [**REQUIRED**]



.. _sys-calibration-label:

Use for photometry system calibration
=====================================

Select standards
----------------
Observe standards fields in your filters 30 FITS files of each field will be enough.
You should observe different standard fields on different zenith angles (or same field on different zenith angles).
Preferred are `Landold Equatorial Standard <https://www.eso.org/sci/observing/tools/standards/Landolt.html>`_
but script works with `Gaia EDR3 and Johnson-Kron-Cousins standards Catalogue <https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A%2BA/664/A109>`_
that contain [Ref]_ other standard fields such as M67 for example.


.. [Ref] Pancino, E., “The Gaia EDR3 view of Johnson-Kron-Cousins standard stars: the curated Landolt and Stetson collections”, *Astronomy and Astrophysics*, vol. 664, Art. no. A109, EDP, 2022. doi:10.1051/0004-6361/202243939.



Prepare config
--------------

Additional scripts
==================

ccd_plot.py

time2jd

jd2time

star_phot_flux.py






.. toctree::
   :maxdepth: 2
   :caption: Contents:

