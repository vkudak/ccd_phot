Additional scripts
==================

ccd_plot
-----------
`ccd_plot.py`  path/to/res_file
   Plot the LC of RSO created by `sat_phot.py` script. All parameters are set in config_sat.ini


res_time2jd
-----------
`res_time2jd.py`  path/to/res_file or dir with res_files, can read a files by mask (example: *.ph*)
   Convert time format in the LC result file created by `sat_phot.py` script.
   Convert from UT to JD.

res_jd2time
-----------
`res_jd2time.py`  path/to/res_file or dir with res_files, can read a files by mask (example: *.ph*)
   Convert time format in the LC result file created by `sat_phot.py` script.
   Convert from JD to UT.

star_phot_flux
--------------
`star_phot_flux.py`  path/to/fits/files
   Act similat to `sat_phot.py` but produce flux for selected star.
   Start must be marked as target in first FITS file header (KEYs `OBJX` & `OBJY`).


res2eossa
---------
`res2eossa.py`  path/to/phX/files or dir with phX files
   Converter to EOSSA v3.1.1
   documentation `here <https://amostech.com/TechnicalPapers/2021/EOSSA_File_Specification_v3.1.1_release4_DistroA.pdf>`_
