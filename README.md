![GitHub Issues or Pull Requests](https://img.shields.io/github/issues/vkudak/ccd_phot)
![GitHub Issues or Pull Requests](https://img.shields.io/github/issues-closed/vkudak/ccd_phot)
![GitHub Release](https://img.shields.io/github/release/vkudak/ccd_phot)

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
* Virtual environment should be created and all requirements installed.
Run scripts from this ENV in the future.

OR (Windows option):

* Add **ccd_phot** directory path to system environment, so you can use scripts right from command line if all requirements all satisfied. 

# Usage
1. Mark satellite on first FITS frame (Add **OBJX** and **OBJY** keywords in header).
    If object is not marked script will try to find it (with **DAOFind** func) as brightest object on the frame.
3. Prepare `config_sat.ini` in same directory where FITS files are located. Config can be namad in other way, but must match the mask `config_satSOMETHING.ini`
3. Start `sat_phot.py` with path to FITS files as parameter (`sat_phot.py . ` OR `sat_phot.py c:\test`)
4. Make resulting plots `ccd_plot.py path\\to\\result_file.phX`
5. Done :)


# Other arguments
usage: `sat_phot.py [-h] [-c CONFIG] [-s_file START_FILE] path`

positional arguments:

  **path**                  `Path to fits files. REQUIRED!`

options:

  **-h, --help**            `show this help message and exit`

  **-c CONFIG, --config CONFIG** `Specify config file`

  **-s_file START_FILE, --start_file START_FILE**
                        `Specify file to start from`
