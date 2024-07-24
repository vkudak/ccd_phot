#!python38
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import configparser

from astropy.time import Time


filename = sys.argv[1]

# cwd = os.getcwd()
wd = os.path.dirname(filename)
# print(wd)


config = configparser.ConfigParser(inline_comment_prefixes="#")
config.read(os.path.join(wd, 'config_sat.ini'))

time_format = config.get('STD', 'Time_format', fallback="UT")
plot_errors = False

try:
    filt = config['STD']['Filter']
    filt = filt.strip("\n")
    filt = filt.strip("\t")
    filt = filt.strip("\r")
    plot_errors = config.getboolean('PLOT', 'plot_errors', fallback=False)
except Exception as e:
    # print(e)
    print("No config file found...\nReading data from header")

    # filt = "None"

    header = {}
    with open(filename) as fres:
        for line in fres:
            if line.startswith("#") and "=" in line:
                hkay, hdata = line[1:].split("=")
                hkay = hkay.strip()
                hdata = hdata.strip()
                hdata = hdata.strip("\n")
                hdata = hdata.strip("\t")
                hdata = hdata.strip("\r")
                header[hkay] = hdata
    if header["Filter"]:
        filt = header["Filter"]
    else:
        filt = "None"
    pass


def read_name(filename):
    cospar = "none"
    norad = "none"
    name = "none"
    dt = "none"

    with open(filename, "r") as myfile:
        head = [next(myfile) for x in range(20)]
    # print(head)
    # print("---------------------------------"
    for line in head:
        line = line.split()
        # print (line)
        if line[1].strip() == "COSPAR":
            cospar = line[3]
        if line[1].strip() == "NORAD":
            norad = line[3]
        if line[1].strip() == "NAME":
            name = " ".join(line[3:])
            name = name.lstrip('0 ')
        if line[1].strip() == "dt":
            dt = line[3]
    return cospar, norad, name, dt


cospar, norad, name, dt = read_name(filename)

dt = float(dt)
# print(cospar, norad, name, dt)

date_time = []

if time_format == "UT":
    flux, f_err, mR, m_err, Az, El = np.genfromtxt(filename, skip_header=True, usecols=(6, 7, 8, 9, 10, 11), unpack=True)
    date, time = np.genfromtxt(filename, unpack=True, skip_header=True, usecols=(0, 1), dtype=None, encoding="utf-8")

    for i in range(0, len(date)):
        # date_time.append(datetime.strptime(date[i].decode('UTF-8') + ' ' + time[i].decode('UTF-8') + "000", "%Y-%m-%d %H:%M:%S.%f"))
        date_time.append(datetime.strptime(date[i] + ' ' + time[i] + "000", "%Y-%m-%d %H:%M:%S.%f"))
else:  # JD
    flux, f_err, mR, m_err, Az, El = np.genfromtxt(filename, skip_header=True, usecols=(5, 6, 7, 8, 9, 10), unpack=True)
    jd = np.genfromtxt(filename, unpack=True, skip_header=True, usecols=(0,), dtype=None, encoding="utf-8")
    jd = [float(x) for x in jd]
    jd = Time(jd, format='jd', scale='utc')

    for i in range(0, len(jd)):
        date_time.append(jd[i].datetime)

    date = [x.date().strftime("%Y-%m-%d") for x in date_time]

# print(date[0], type(date[0]))
# print(date_time[0], type(date_time[0]))

# print(date)

# old style (working!!!!)
# flux, mR, Az, El = np.loadtxt(filename, unpack=True, usecols=(6, 8, 10, 11))
# date, time = np.loadtxt(filename, unpack=True, skiprows=11, usecols=(0, 1), dtype={'names': ('date', 'time'), 'formats': ('S10', 'S12')})

# print (date_time[2])
bands = ["B", "V", "R", "C", "None"]
index = bands.index(filt)
cc = ["blue", "green", "red", "black", "black"]
f_color = cc[index]

## fig im MAG
plt.rcParams['figure.figsize'] = [12, 6]
dm = max(mR) - min(mR)
dm = dm * 0.1
plt.axis([min(date_time), max(date_time), max(mR) + dm , min(mR) - dm])

if plot_errors:
    plt.errorbar(date_time, mR, yerr=m_err, mfc=f_color, fmt='x-', linewidth=0.8, fillstyle="none", markersize=3, capsize=3)
else:
    plt.plot(date_time, mR, "x-", color=f_color, linewidth=0.5, fillstyle="none", markersize=3)


d, t = str(date_time[0]).split(" ")
plt.title("Satellite Name:%s, NORAD:%s, COSPAR:%s \n Date=%s  UT=%s   dt=%2.3f  Filter=%s" % (name, norad, cospar, d, t, dt, filt), pad=6, fontsize=12)
plt.ylabel('m_st')
plt.xlabel('UT')
ax = plt.gca()

# Azimuth axis----------------------------------
ax2 = ax.twiny()
ax2.set_xlim(ax.get_xlim())
numElems = 5
tt_idx = np.round(np.linspace(0, len(date_time) - 1, numElems)).astype(int)
Tt2 = np.array(date_time)
Az2 = np.array(Az)
El2 = np.array(El)

# Az2s = ["%3.2f" % Azt for Azt in Az2]
Az2s = []
Tt2s = []
for kk in range(0, len(Az2)):
    azt = Az2[kk]
    elt = El2[kk]
    t = Tt2[kk]
    Az2s.append("%3.1f; %3.1f"%(azt, elt))
    Tt2s.append(t.strftime("%H:%M:%S"))
Az2s = np.array(Az2s)
ax2.set_xticks(Tt2[tt_idx])  # new_tick_locations
ax2.set_xticklabels(Az2s[tt_idx], fontsize=8)
ax2.set_xlabel(r"Az;h [deg]", fontsize=8)
ax2.tick_params(axis='x', which='major', pad=0)

Tt2s = np.array(Tt2s)
ax.set_xticks(Tt2[tt_idx])  # new_tick_locations
ax.set_xticklabels(Tt2s[tt_idx], fontsize=10)

ax.xaxis.grid()
ax.yaxis.grid()
# ----------------------------------------------------

# plt.show()
Tf = Tt2s[0]
Tf = Tf.split(":")
TF = Tf[0]+Tf[1]+Tf[2]
ymd = date[0].replace("-", "")

figname = os.path.join(wd, norad + "_" + ymd + "_UT" + TF + ".png")
# plt.savefig(norad + "_UT" + TF + ".png")
plt.savefig(figname)
plt.clf()
####


## fig im IMPULS-------------------------------------------------------------------
plt.rcParams['figure.figsize'] = [12, 6]
dm = max(flux) - min(flux)
dm = dm * 0.1
# print(flux[0])
plt.axis([min(date_time), max(date_time), min(flux) - dm, max(flux) + dm])

if plot_errors:
    plt.errorbar(date_time, flux, yerr=f_err, mfc=f_color, fmt='x-', linewidth=0.8, fillstyle="none", markersize=3, capsize=3)
else:
    plt.plot(date_time, flux, "x-", color=f_color, linewidth=0.5, fillstyle="none", markersize=3)

d, t = str(date_time[0]).split(" ")
plt.title("Satellite Name:%s, NORAD:%s, COSPAR:%s \n Date=%s  UT=%s   dt=%2.3f  Filter=%s" % (name, norad, cospar, d, t, dt, filt), pad=6, fontsize=12)
plt.ylabel('Flux')
plt.xlabel('UT')
ax = plt.gca()

# Azimuth axis----------------------------------
ax2 = ax.twiny()
ax2.set_xlim(ax.get_xlim())
numElems = 5
tt_idx = np.round(np.linspace(0, len(date_time) - 1, numElems)).astype(int)
Tt2 = np.array(date_time)
Az2 = np.array(Az)
El2 = np.array(El)

# Az2s = ["%3.2f" % Azt for Azt in Az2]
Az2s = []
Tt2s = []
# flux2s = []
for kk in range(0, len(Az2)):
    azt = Az2[kk]
    elt = El2[kk]
    f = flux[kk]
    t = Tt2[kk]
    Az2s.append("%3.1f; %3.1f" % (azt, elt))
    Tt2s.append(t.strftime("%H:%M:%S"))
    # flux2s.append("%12.6f" % f)
Az2s = np.array(Az2s)
ax2.set_xticks(Tt2[tt_idx])  # new_tick_locations
ax2.set_xticklabels(Az2s[tt_idx], fontsize=8)
ax2.set_xlabel(r"Az;h [deg]", fontsize=8)
ax2.tick_params(axis='x', which='major', pad=0)

Tt2s = np.array(Tt2s)
ax.set_xticks(Tt2[tt_idx])  # new_tick_locations
ax.set_xticklabels(Tt2s[tt_idx], fontsize=10)


# flux2s = np.array(flux2s)
# ff_idx = np.round(np.linspace(0, len(flux) - 1, numElems)).astype(int)
# ax.set_yticks(flux[ff_idx])  # new_tick_locations
# ax.set_yticklabels(flux2s[ff_idx], fontsize=10)

import matplotlib.ticker as ticker
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%8.3f'))

ax.xaxis.grid()
ax.yaxis.grid()
# ----------------------------------------------------

# plt.show()
figname = os.path.join(wd, norad + "_" + ymd + "_UT" + TF + "_Imp.png")
# plt.savefig(norad + "_UT" + TF + "_Imp.png")
plt.savefig(figname)
plt.clf()
####
# plt.savefig('line_plot.jpg', dpi=300, quality=80, optimize=True, progressive=True)  # 50% less size
