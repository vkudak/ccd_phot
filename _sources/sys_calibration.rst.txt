
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
that contain other standard fields such as M67 for example [1]_.

Calibration process is based on papers [2]_, [3]_, [4]_.

All FITS files should be organized in folders by filter. For example you will have `V` folder with all FITS files
captured in V band and `B` folder where all FITS fileas are captured in B band.


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

* `stars_process.py`  path/to/fits/files

  Run to process all FITS files. It will produce `ref_stars_cat.json` and logs. Only parameter that must be given is path to folder with FITS files.

* `stars_calc_params.py` path/to/ref_stars_cat.json

  Calculate system parameters and build the graphs.
  System parameters will appear in stdout and in log file.
  As the parameter set the path to `ref_stars_cat.json` file.

If you wish to recalculate sys parameters (Z_x, C_x, K_x) with other parameters,
for example change max MAG of selected stars, you can do it without processing FITS files all data are already stored
in `ref_stars_cat.json` file.

References
----------

Calibration process is based on papers:

.. [1] Pancino, E., “The Gaia EDR3 view of Johnson-Kron-Cousins standard stars: the curated Landolt and Stetson collections”,
         *Astronomy and Astrophysics*, vol. 664, Art. no. A109, EDP, 2022. `DOI <https://doi.org/10.1051/0004-6361/202243939>`_.


.. [2]  Hu, S. M., Han, S. H., Guo, D. F., & Du, J. J. (2014).
    The photometric system of the One-meter Telescope at Weihai Observatory of Shandong University.
    Research in Astronomy and Astrophysics, 14(6), 719.
    `DOI <https://iopscience.iop.org/article/10.1088/1674-4527/14/6/010/meta>`__

.. [3]  Huang, F., Li, J. Z., Wang, X. F., Shang, R. C., Zhang, T. M., Hu, J. Y., ... & Jiang, X. J. (2012).
    The photometric system of the Tsinghua-NAOC 80-cm telescope at NAOC Xinglong Observatory.
    Research in Astronomy and Astrophysics, 12(11), 1585.
    `DOI <https://iopscience.iop.org/article/10.1088/1674-4527/12/11/012/meta>`__

.. [4]  Kinoshita, D., Chen, C. W., Lin, H. C., Lin, Z. Y., Huang, K. Y., Chang, Y. S., & Chen, W. P. (2005).
    Characteristics and performance of the CCD photometric system at Lulin observatory.
    Chinese Journal of Astronomy and Astrophysics, 5(3), 315.
    `DOI <https://iopscience.iop.org/article/10.1088/1009-9271/5/3/011/meta>`__
