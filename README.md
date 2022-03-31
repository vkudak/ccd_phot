# ccd_stat_phot
Photometry of artificial satellites with CCD QHY-174M

Python scripts to process photometry obtained on CCD QHY-174M

Based on:
- astropy
- photutils
- photometry errors are estimated with wfc3_photometry procedures (https://github.com/spacetelescope/wfc3_photometry)
- astroquery

# Usage
1. Mark satellite on first Frame
2. prepeare comfig_sat.ini in same directory where fits files are located
3. start sat_phot.py with path to fits files as parameter ("sat_phot.py . " OR "sat_phot.py c:\test" )
4. make result plot "ccd_plot.py path to result_file.txt"
5. Done :)
