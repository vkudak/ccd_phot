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
from sp_utils import convert_ndarray, read_config_stars, convert_to_numpy, lsqFit

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
# try:
#     db_file = sys.argv[1]
#     database, exp = pickle.load(open(db_file, "rb"))
# except Exception as E:
#     print("Error database file read", db_file)
#     print(repr(E))
#     sys.exit()


path = os.path.dirname(db_file)
log_file = open(os.path.join(path, 'star_calc.log'), "w")


conf = read_config_stars(os.path.join(path, 'config_stars.ini'), log_file)
if not conf:
    sys.exit()

max_m = float(conf['max_m'])

station = ephem.Observer()
station.lat = conf['site_lat']  #'48.5635505'
station.long = conf['site_lon']  #'22.453751'
station.elevation = conf['site_elev']  #231.1325

y_ar = []
x_ar = []
yerr_ar = []
l_ar = []

A_m_list = []
A_mR_list = []
# frf = open("test.txt", "w")
# frfxy = open("test_xy.txt", "w")

for i in range(len(database)):
    # RMS filter#############################################################
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

    # print (database[i]["NOMAD1"], database[i]["Flux_mean"], database[i]["Flux_std"])
    if (database[i]["Flux_mean"] > 0) and \
            (database[i]["Vmag"] < max_m) and \
            (abs(database[i]["V-R"]) < 1.95) and \
            (len(database[i]["Flux"]) > 5):
        m_inst = -2.5 * math.log10(database[i]["Flux_mean"])

        # TODO: make two systems to solve with linReg. One for A and C, another for A and K coefficients

        if conf['c_flag']:
            # yq = database[i]["Rmag"] + m_inst - kr * database[i]["Mz"]
            # yq = database[i]["Vmag"] - m_inst - conf['kr'] * database[i]["Mz"]
            yq = m_inst - database[i]["Vmag"] - conf['kv'] * database[i]["Mz"]
            # yq = database[i]["Vmag"] - m_inst + conf['kv'] * database[i]["Mz"]

            # print(database[i]["Vmag"], m_inst, conf['kv'], database[i]["Mz"])

            # if (yq < 17) and (yq > 14) and (database[i]["f/b"].mean(axis=0) > snr_value):
            if database[i]["f/b"].mean(axis=0) > conf['snr_value']:
                database[i]["yq"] = yq
                y_ar.append(database[i]["yq"])
                x_ar.append(database[i]["B-V"])
                fi = np.array(2.5 * np.log10(database[i]["Flux"]))
                yerr_ar.append(fi.std(axis=0))  # std from  2.5 * log10(Flux)
                l_ar.append(database[i]["Name"])

                # frfxy.write("%17s    %4.2f   %4.2f  %4.2f   %4.2f   %4.2f\n" %
                #             (database[i]["Name"], database[i]["Rmag"], database[i]["V-R"],
                #              database[i]["Flux_mean"], database[i]["V-R"], database[i]["yq"]))
                database[i]["Good"] = ">Good<"
            else:
                database[i]["Good"] = "low s/n"

        else:
            database[i]["A"] = (database[i]["Vmag"] -
                                m_inst -
                                conf['kv'] * database[i]["Mz"] -
                                conf['Cv'] * database[i]["V-R"]
                                )
            if database[i]["f/b"].mean(axis=0) > conf['snr_value']:
                A_m_list.append(database[i]["A"])
                A_mR_list.append(database[i]["Vmag"])
                database[i]["Good"] = True
        # # save to Fedorovich
        # Mz = database[i]["Mz"]
        # el2 = math.degrees(np.arccos(1 / Mz))
        # frf.write("%4.2f  %4.2f  %3.1f  %6i\n" % (database[i]["Rmag"], database[i]["V-R"], el2, int(database[i]["Flux_mean"])))
    else:
        database[i]["Good"] = "Bed"

# frf.close()
# frfxy.close()
if conf['c_flag']:
    print("Stars left =", len(y_ar))
    log_file.write("Stars left = %i\n" % len(y_ar))
else:
    print("Stars left =", len(A_m_list))
    log_file.write("Stars left = %i\n" % len(A_m_list))

if conf['c_flag']:
    y_ar = np.array(y_ar)
    x_ar = np.array(x_ar)

    # am, bm = lreg(x_ar, y_ar) # Klimik formula
    # print(am, bm)

    r_max = 999
    # r_max_val = 0.75
    # print(conf['r_max_val'])
    while (r_max > conf['r_max_val']) and (len(y_ar) > 5):
        # if len(y_ar) > 5:
        c, a, r_max, ind, r2 = lsqFit(y_ar, x_ar)
        c_err = 99
        a_err = 99

        # with errors
        # aa, bb, r2, r_max, ind = fit_lin_reg(x_ar, y_ar, yerr_ar)
        # c, c_err = aa
        # a, a_err = bb
        ####################

        # a2, c2, r_sq = linReg(x_ar, y_ar)
        c_fit = c
        # if c < 0:
        #     c = 1.0 / c
        # print("A=%5.3f  Cr=%5.3f  R^2=%5.3f" % (a2, c2, r_sq))
        # c = c2

        # https://www.annualreviews.org/doi/pdf/10.1146/annurev.astro.41.082801.100251
        lya_v = 544.8 #547
        lya_r = 640.7 #635     lya_eff
        lya_eff = (lya_r * lya_v) / (c * (lya_v - lya_r) + lya_v)
        print("############################LSQ_FIT Results#################################")
        print("A = %2.5f +/- %2.5f, c = %2.5f +/- %2.5f" % (a, a_err, c, c_err))
        print("Lyambda_eff = %5.3f" % lya_eff)
        # log_file.write("A = %3.8f  c = %3.8f\n" % (a, c))

        plt.plot(x_ar, y_ar, "xr")
        p1 = [min(x_ar), max(x_ar)]
        p2 = [a + c_fit * min(x_ar), a + c_fit * max(x_ar)]
        # print(p1)
        # print(p2)
        plt.plot(p1, p2, "k")
        plt.xlabel("B-V")
        plt.ylabel(r'$V_{inst}-V_{st}-K_{v} \cdot M_{z}$')
        # plt.title(fit_file)
        # plt.show()
        plt.savefig(os.path.join(path, "graph_Cr.png"))
        plt.close()

        log_file.write("\n\n")
        log_file.write("####################################################################\n")
        log_file.write("###--------LSQ_FIT Results--(A and Cr all frames)----------------###\n")
        log_file.write("A = %3.5f +/- %3.5f, c = %3.5f +/- %3.5f   R^2=%2.3f\n" % (a, a_err, c, c_err, r2))
        log_file.write("Lyambda_eff = %5.3f nm\n" % lya_eff)
        log_file.write("###--------------------------------------------------------------###\n")
        if r_max > conf['r_max_val']:
            y_ar = np.delete(y_ar, ind)
            x_ar = np.delete(x_ar, ind)
            yerr_ar = np.delete(yerr_ar, ind)
            for i in range(len(database)):
                if database[i]["Name"] == l_ar[ind]:
                    database[i]["Good"] = "Filtered %i r_max= %2.5f" % (ind, r_max)
            log_file.write("## rmax=%5.3f  ind=%i, name = %s\n" % (r_max, ind, l_ar[ind]))
            l_ar = np.delete(l_ar, ind)
        # else:
        #     print("Only %i values. Cand perform LSQ_FIT...skipping frame" % len(y_ar))
        #     log_file.write("Only %i values. Cand perform LSQ_FIT...skipping frame\n" % len(y_ar))
    if r_max > conf['r_max_val']:
        print("Only %i values. Cand perform LSQ_FIT...skipping filtering" % len(y_ar))
        log_file.write("Only %i values. Cand perform LSQ_FIT...skipping filtering\n" % len(y_ar))

    # ### Write stars data
    log_file.write("######################--STARS  DATA--###################################\n")
    log_file.write("Name                  Bmag   Vmag    B-V      Flux_mean      Flux_std   n_count  (f/b)_maen    Mz      Note\n")
    database = sorted(database, key=lambda tstar: tstar["Good"])
    for star in database:
        # log_file.write("%17s  %5.3f  %5.3f  % 2.3f     %8.3f   %5.3f    %i         %2.3f     %2.3f    %s\n" %
        #                (star["Name"], star["Vmag"], star["Rmag"], star["V-R"], star["Flux_mean"], star["Flux_std"],
        #                 len(star["Flux"]), star["f/b"].mean(axis=0), star["Mz"], star["Good"]))

        log_file.write("{:20s}  {:{width}.{prec}f}  {:{width}.{prec}f}   {:{width}.{prec}f}  {:{width2}.{prec2}f}  {:{width2}.{prec2}f}  {:{width}d}       {:{width}.{prec}f}     {:{width}.{prec}f}    {:7s}\n".format(
                       star["Name"].strip(),
                       star["Bmag"],
                       star["Vmag"],
                       star["B-V"],
                       star["Flux_mean"],
                       star["Flux_std"],
                       len(star["Flux"]),
                       star["f/b"].mean(axis=0),
                       star["Mz"],
                       str(star["Good"]),
                       width=5, prec=2, width2=12, prec2=3))
    #####################
    # '{:{width}.{prec}f}'.format(2.7182, width=5, prec=2)

else:
    plt.plot(A_mR_list, A_m_list, "xr")
    plt.xlabel(r'$m_R$')
    plt.ylabel("A")
    plt.title("A value all data")
    # plt.show()
    plt.savefig(os.path.join(path, "graph_A1.png"))
    plt.close()

    print("\nPerforming A rms<%3.5f filter" % conf['rms_val'])
    log_file.write("\nPerforming A rms<%3.5f filter\n" % conf['rms_val'])
    A_m_list, A_mR_list = RMS_del(A_m_list, conf['rms_val'], B=A_mR_list, log_file=log_file)

    plt.plot(A_mR_list, A_m_list, "xr")
    plt.xlabel(r'$m_R$')
    plt.ylabel("A")
    plt.title("A value, filtered %s" % conf['rms_val'])
    # plt.show()
    plt.savefig(os.path.join(path, "graph_A2.png"))
    plt.close()

    log_file.write("NEW A count: %i\n" % len(A_m_list))
    strFormat = len(A_m_list) * '{:5.3f}, '
    formattedList = strFormat.format(*A_m_list)
    log_file.write("Left A values:" + formattedList + "\n")

    A_mean = np.mean(A_m_list, axis=0)
    A_err = np.std(A_m_list, axis=0)

    print("############################--A Results--#################################")
    print("A_mean = %2.5f , A_err = %2.5f " % (A_mean, A_err))

    log_file.write("\n\n")
    log_file.write("####################--A Results--#####################################\n")
    log_file.write("A_mean = %8.5f , A_err =%8.5f\n" % (A_mean, A_err))
    log_file.write("###--------------------------------------------------------------###\n")

log_file.write("####################################################################\n")
log_file.close()
