.. ccd_phot documentation master file, created by
   sphinx-quickstart on Sat Mar  1 21:58:28 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


CCD_PHOT DOCUMENTATION
======================

This is not a package, its a collection of scripts !

Python scripts to process photometry of RSO objects obtained on CCD QHY-174M-GPS (tested on ZWO analog ASI 174) and SharpCap Software.

.. toctree::
   :maxdepth: 3
   :caption: Contents:

   installation
   start_scripts
   rso_photometry
   sys_calibration
   additional_scripts


Tested Hardware & Software
--------------------------
* QHY 174M-GPS
* ZWO ASI 174 MM
* SharpCap (for work with QHY GPS module)


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
*Artificial Satellites. Journal of Planetary Geodesy*, 57(1), 47-57. `DOI <https://doi.org/10.2478/arsa-2022-0003>`__

