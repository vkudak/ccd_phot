Photometry of RSO
=========================

Prepare FITS files
-------------------
Your RSO photometry FITS files should be in some directory.
They should be named in order of creation, otherwise the script will lose the target
In this case use ``fits_sort`` option in ``config_sat.ini`` file to sort the FITS files, see :ref:`sat-config-label`.

Avoid names like ``file_1.fits``, ``file_10.fits``, etc. Preferred name patters is ``{date_folder}/{norad_folder}/{norad}_{HHMMSS of start}_{count}.fits``

Or simply use: ``file_00001.fits``, ``file_00002.fits`` and so on.

FITS files do not need to contain WCS solution, ``sat_phot.py`` script do not use it at all.

Starting from ver 1.0.1 ``sat_phot.py`` script can try to find object automatically on the first FITS file.
If automatic script fails than it look for add additional keywords in the first FITS file header
``OBJX: <float>`` & ``OBJY: <float>`` or for config file (in section ``[OBJ_POS]``) to tell the script where the target initially is situated on the frame.

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

.. literalinclude:: ../../config_sat_example.ini
   :language: ini


In general ``sat_phot.py`` has few options:

* **[path]**               Path to FITS files [**REQUIRED**]
* **[-c or - - config ]**    Specify config file [Optional]
* **[-s or - - start_file]** Specify the name of the FITS file to start from [Optional]. Could be useful if photometry script fails in the middle of the file list. This can happen if the object is too faint. In such case you can use this option to start from required FITS. But first rename result file or it will be overwritten.



.. admonition:: Warning For Windows users:

    Try to avoid folder and file names with spaces. It can be interpreted as two or more parameters which leeds to incorrect work of the scripts

