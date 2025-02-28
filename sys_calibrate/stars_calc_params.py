#!python3
from astropy.stats import sigma_clipped_stats
from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt
from numpy.ma import masked


import math
# import glob
import sys
import os
import json

from astroquery.astrometry_net import AstrometryNet
from astropy.wcs import WCS
import ephem
import configparser
import warnings

# Додаємо кореневий каталог у PYTHONPATH
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from photometry_with_errors import iraf_style_photometry
from sp_utils import convert_ndarray, read_config_stars, convert_to_numpy, lsqFit, solve_photometric_coefficients

# !!!! http://vizier.u-strasbg.fr/viz-bin/VizieR-4
# Stellar Photometry in Johnson's 11-color system (Ducati, 2002)


'''
Check FILTER and use one of equations:
Uinst = Ustd + Zu + Ku′X + Cu(U − B)std
Binst = Bstd + Zb + Kb′X + Cb(B − V)std
Vinst = Vstd + ZV + Kv′X + Cv(B − V)std
Rinst = Rstd + ZR + Kr′X + Cr(V − R)std
Iinst = Istd + ZI + Ki'X + Ci(V − I)std

Papers to help:
https://www.astro.ncu.edu.tw/~wchen/wp_chen/essay/Kinoshita2005ChJAA5-315%20LOTPI1300B.pdf

https://arxiv.org/pdf/1205.6529.pdf
https://arxiv.org/pdf/1401.4281.pdf
https://adsabs.harvard.edu/pdf/1983PASP...95.1021S
https://www.ias.ac.in/article/fulltext/joaa/012/04/0319-0331

http://dspace.onu.edu.ua:8080/handle/123456789/4660
https://slittlefair.staff.shef.ac.uk/teaching/phy241/lectures/l07/
'''





if len(sys.argv) < 2:
    print("Not enough parameters. Enter path to database file")
    sys.exit()



db_file = sys.argv[1]

# Завантаження JSON-файлу
with open(db_file, "r") as cat_filename:
    database = json.load(cat_filename)
    database = convert_to_numpy(database)


path = os.path.dirname(db_file)
log_file = open(os.path.join(path, 'star_calc.log'), "w")


conf = read_config_stars(os.path.join(path, 'config_stars.ini'), log_file)
if not conf:
    sys.exit()


station = ephem.Observer()
station.lat = conf['site_lat']  #'48.5635505'
station.long = conf['site_lon']  #'22.453751'
station.elevation = conf['site_elev']  #231.1325
band = conf["band"]

for i in range(len(database)):
    # RMS filter inside each star  #############################################################
    dd = True
    while dd:
        database[i]["Flux_mean"] = database[i]["Flux"].mean(axis=0)
        database[i]["Flux_std"] = database[i]["Flux"].std(axis=0)

        database[i]["Flux-mean"] = abs(database[i]["Flux"] - database[i]["Flux_mean"])
        az = []
        for z in range(len(database[i]["Flux"])):
            if database[i]["Flux-mean"][z] > 3 * database[i]["Flux"].std(axis=0):
                az.append(z)

        if len(az) > 0:
            f = np.array(database[i]["Flux"])
            print(database[i]["Name"], "delete..", az)
            log_file.write(database[i]["Name"] + " delete.. " + str(az))
            log_file.write("\n")
            f = np.delete(f, az)
            database[i]["Flux"] = f
            database[i]["Flux_mean"] = database[i]["Flux"].mean(axis=0)
            database[i]["Flux_std"] = database[i]["Flux"].std(axis=0)
            snr_tmp = database[i]["f/b"]
            snr_tmp = np.delete(snr_tmp, az)
            database[i]["f/b"] = snr_tmp
            dd = True
        else:
            dd = False
    #############################################################


m_inst_arr = []
m_std_arr = []
x_arr= []
c_ind_arr = []
star_names_arr = []
c_ind_name = {"V":"B-V", "B":"B-V", "R":"V-R"}

# make arrays
for star in database:
    star['Good'] = "OK; "
    if star["Flux_mean"] < 0:
        star['Good'] += ' Negative Flux;'
    elif star[f"{band}mag"] > conf['max_m_calc']:
        star['Good'] += f" Mag > {conf['max_m_calc']};"
    elif len(star["Flux"]) < 5:
        star['Good'] += ' Got les than 5 Flux values;'
    elif star["f/b"].mean(axis=0) < conf['snr_value']:
        star['Good'] += f" SNR < {conf['snr_value']};"
    else:
        # All good
        m_inst_arr.append(-2.5 * math.log10(star["Flux_mean"]))
        m_std_arr.append(star[f"{band}mag"])
        x_arr.append(star["Mz"])
        c_ind_arr.append(star[c_ind_name[band]])
        star_names_arr.append(star["Name"])


print("Stars left =", len(star_names_arr))
log_file.write("Stars left = %i\n" % len(star_names_arr))
if len(star_names_arr) == 0:
    print('No data to process')
    sys.exit()


if conf["K"] is None:
    Z, K, C, sigma_Z, sigma_K, sigma_C, removed_stars, r2 = (
        solve_photometric_coefficients(m_inst_arr, m_std_arr, x_arr, c_ind_arr, star_names_arr, band,
                                       threshold=conf['r_max_val'],
                                       log_file=log_file,
                                       plot=True,
                                       path=path)
    )

else:
    Z, K, C, sigma_Z, sigma_K, sigma_C, removed_stars, r2 = (
        solve_photometric_coefficients(m_inst_arr, m_std_arr, x_arr, c_ind_arr, star_names_arr, band,
                                       threshold=conf['r_max_val'],
                                       log_file=log_file,
                                       k=conf['K'],
                                       plot=True,
                                       path=path)
    )


print(f"Z_{band} = {Z:6.3f} ± {sigma_Z:6.4f}")
if conf['K'] is None:
    print(f"K_{band} = {K:6.4f} ± {sigma_K:6.4f}")
print(f"C_{band} = {C:6.4f} ± {sigma_C:6.4f}")
print(f"R^2 = {r2:6.4f}")


for i in range(len(database)):
    for name, value in removed_stars:
        if name == database[i]["Name"]:
            database[i]["Good"] = f"Removed with value = {value:6.4f} "


# ### Write stars data
log_file.write("######################--STARS  DATA--###################################\n")
log_file.write(f"Name                  {band}mag     {c_ind_name[band]}     Flux_mean      Flux_std   n_count  (f/b)_maen    Mz        alt      Note\n")
database = sorted(database, key=lambda tstar: tstar["Good"])
for star in database:
    alt = math.asin(1 / star["Mz"])
    alt = math.degrees(alt)
    log_file.write(("{:20s}  {:{width}.{prec}f}   {:{width}.{prec}f}  {:{width2}.{prec2}f}  "+\
                   "{:{width2}.{prec2}f}  {:{width}d}       {:{width}.{prec}f}     {:{width}.{prec}f}     " + \
                   "{:{width}.{prec}f}    {:7s}\n").format(
                   star["Name"].strip(),
                   star[f"{band}mag"],
                   star[c_ind_name[band]],
                   star["Flux_mean"],
                   star["Flux_std"],
                   len(star["Flux"]),
                   star["f/b"].mean(axis=0),
                   star["Mz"],
                   alt,
                   str(star["Good"]),
                   width=5, prec=2, width2=12, prec2=3))
#####################
# '{:{width}.{prec}f}'.format(2.7182, width=5, prec=2)
log_file.write("####################################################################\n")

log_file.write(f"Z_{band} = {Z:6.3f} ± {sigma_Z:6.4f}\n")
if conf['K'] is None:
    log_file.write(f"K_{band} = {K:6.4f} ± {sigma_K:6.4f}\n")
log_file.write(f"C_{band} = {C:6.4f} ± {sigma_C:6.4f}\n")
log_file.write(f"R^2 = {r2:6.4f}")

log_file.close()
