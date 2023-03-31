Photometry of artificial satellites with CCD QHY-174M
===================
## ccd_stat_phot

Python scripts to process photometry obtained on CCD QHY-174M

Based on:
- astropy
- photutils
- photometry errors are estimated with wfc3_photometry procedures (https://github.com/spacetelescope/wfc3_photometry)
- astroquery

# Before usage
* You should create a virtual environment and install all requirements.
Run scripts from this ENV in the future.

OR (more like Windows option):

* Add **ccd_phot** directory path to system environment, so you can use scripts right from command line if all requirements all satisfied. 

# Usage
1. Mark satellite on first FITS frame (Add **OBJX** and **OBJY** keywords in header)
2. Prepare `config_sat.ini` in same directory where FITS files are located
3. Start `sat_phot.py` with path to FITS files as parameter (`sat_phot.py . ` OR `sat_phot.py c:\test`)
4. Make resulting plots `ccd_plot.py path\\to\\result_file.phX`
5. Done :)
