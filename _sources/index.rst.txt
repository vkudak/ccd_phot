.. ccd_phot documentation master file, created by
   sphinx-quickstart on Sat Mar  1 21:58:28 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


CCD_PHOT DOCUMENTATION
======================

This is not a package, its a collection of scripts !

Python scripts to process photometry of RSO objects obtained on CCD QHY-174M-GPS (tested on ZWO analog ASI 174) and SharpCap Software.

Tested Hardware & Software
--------------------------
* QHY 174M-GPS
* ZWO ASI 174 MM
* SharpCap (for work with QHY GPS module)

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


Photometry of RSO
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
   lat = 48.000
   lon = 22.000
   h = 230.000

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

Photometry system calibration
=====================================

Select standards
----------------
Observe and grab standards fields in your filters. 30 FITS files of each field in one filter will be enough.
You should observe different standard fields on different zenith angles (or same field on different zenith angles) to
obtain good fit of coefficient of extinction K.

Preferred are `Landold Equatorial Standard <https://www.eso.org/sci/observing/tools/standards/Landolt.html>`_
but script works with `Gaia EDR3 and Johnson-Kron-Cousins standards Catalogue <https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A%2BA/664/A109>`_
that contain other standard fields such as M67 for example [Ref]_.

All FITS files should be organized in folders by filter. For example you will have `V` folder with all FITS files
captured in V band and `B` folder where all FITS fileas are captured in B band.

.. [Ref] Pancino, E., “The Gaia EDR3 view of Johnson-Kron-Cousins standard stars: the curated Landolt and Stetson collections”,
         *Astronomy and Astrophysics*, vol. 664, Art. no. A109, EDP, 2022. `DOI <https://doi.org/10.1051/0004-6361/202243939>`_.



Prepare config
--------------

Example of config file is provided with scripts

.. code-block::

   [Stars_Stand]
   K = 0.35           # coefficient of extinction if Band. Can be None (will be calculated) or value
   Filter = V         # Band
   snr = 1.0          # min star SNR
   max_m_fit = 14.0   # max V_MAG for process on FITS files
   max_m_calc = 13.0  # max MAG for calculations
   r_max_val = 0.25   # delete stars with R2 more that this
   dark_frame = path\to\dark\frame  # Path to dark frame
   dark_stable = 0    # If no dark frame than you can subtract some value
   FoV = 1.5          # Field of View in degrees (radius)

   [APERTURE]
   r_ap = 5
   an_in = 10
   an_out = 14

   [astrometry.net]
   scale_lower = 2.2         # lower limit arcsec/pix
   scale_upper = 2.6         # upper limit arcsec/pix
   api_key = xxxxxxxxxxx     # Astrometry.net api key

   [SITE]
   Name = OBS_NAME
   lat = 48.000
   lon = 22.000
   h = 231.000

Config file should be named `config_stars.ini` and created in directory with FITS files
(in each dir, if you have many dirs for different filters).

Run calibration
---------------

Similar to `sat_phot.py` system calibration scripts can be run in same ways.
The only difference here is that we have two scripts they are stored in `sys_calibrate` folder:

   `stars_process.py`  path/to/fits/files
      Run to process all FITS files. It will produce `ref_stars_cat.json` and logs.
      Only parameter that must be given is path to folder with FITS files.

   `stars_calc_params.py` path/to/ref_stars_cat.json
      Calculate system parameters and build the graphs.
      System parameters will appear in stdout and in log file.
      As the parameter set the path to `ref_stars_cat.json` file.

If you wish to recalculate sys parameters (Z_x, C_x, K_x) with other parameters,
for example change max MAG of selected stars, you can do it without processing FITS files all data are already stored
in `ref_stars_cat.json` file.

References
----------

Calibration process is based on papers:

1. Hu, S. M., Han, S. H., Guo, D. F., & Du, J. J. (2014).
The photometric system of the One-meter Telescope at Weihai Observatory of Shandong University.
Research in Astronomy and Astrophysics, 14(6), 719.
`Link <https://iopscience.iop.org/article/10.1088/1674-4527/14/6/010/meta>`__

2. Huang, F., Li, J. Z., Wang, X. F., Shang, R. C., Zhang, T. M., Hu, J. Y., ... & Jiang, X. J. (2012).
The photometric system of the Tsinghua-NAOC 80-cm telescope at NAOC Xinglong Observatory.
Research in Astronomy and Astrophysics, 12(11), 1585.
`Link <https://iopscience.iop.org/article/10.1088/1674-4527/12/11/012/meta>`__

3. Kinoshita, D., Chen, C. W., Lin, H. C., Lin, Z. Y., Huang, K. Y., Chang, Y. S., & Chen, W. P. (2005).
Characteristics and performance of the CCD photometric system at Lulin observatory.
Chinese Journal of Astronomy and Astrophysics, 5(3), 315.
`Link <https://iopscience.iop.org/article/10.1088/1009-9271/5/3/011/meta>`__


Additional scripts
==================

ccd_plot
-----------

   `ccd_plot.py`  path/to/res_file
      Plot the LC of RSO created by `sat_phot.py` script. All parameters are set in config_sat.ini


res_time2jd
-----------
   `res_time2jd.py`  path/to/res_file
      Convert time format in the LC result file created by `sat_phot.py` script.
      Convert from UT to JD.

res_jd2time
-----------
   `res_jd2time.py`  path/to/res_file
      Convert time format in the LC result file created by `sat_phot.py` script.
      Convert from JD to UT.

star_phot_flux
--------------
   `star_phot_flux.py`  path/to/fits/files
      Act similat to `sat_phot.py` but produce flux for selected star.
      Start must be marked as target in first FITS file header (KEYs `OBJX` & `OBJY`).


Citing CCD_PHOT
===============

If you use CCD_PHOT for a project that leads to a publication, whether directly or as a dependency of another package,
please include the following acknowledgment:

.. code-block::

   This research made use of CCD_PHOT scripts used for
   photometry of RSO and photometry system calibration (Kudak V. <YEAR>).

BibTex files for all Photutils versions can be found at `https://doi.org/10.5281/zenodo.11066198 <https://doi.org/10.5281/zenodo.11066198>`_

For example CCD_PHOT v1.0.1 should be cite Kudak V. 2024 with BibTex entry (https://zenodo.org/records/11066198/export/bibtex)

.. code-block::

   @software{viktor_kudak_2024_11066198,
     author       = {Viktor Kudak},
     title        = {vkudak/ccd\_phot: 1.0.1},
     month        = apr,
     year         = 2024,
     publisher    = {Zenodo},
     version      = {1.0.1},
     doi          = {10.5281/zenodo.11066198},
     url          = {https://doi.org/10.5281/zenodo.11066198},
   }


Related works
--------------

1. Kudak, V., & Perig, V. (2022). QHY-174M-GPS camera as the device for photometry of artificial satellites.
*Artificial Satellites. Journal of Planetary Geodesy*, 57(1), 47-57. `Link <https://doi.org/10.2478/arsa-2022-0003>`__

.. toctree::
   :maxdepth: 3
   :caption: Contents:


