import sys
import os
import glob
import numpy as np
from datetime import datetime

from astropy.time import Time


def read_header(fname):
    header_list = []
    with open(fname) as f:
        for line in f:
            if line[0] == "#":
                header_list.append(line)
    # print(header_list)
    return header_list


def read_data(fname):
    date_time = []
    date, time = np.genfromtxt(fname, unpack=True, skip_header=True, usecols=(0, 1), dtype=None, encoding="utf-8")
    x, y, xerr, yerr, flux, flux_err, mag, mag_err, Az, El, Rg = \
        np.genfromtxt(filename, skip_header=True,
                      usecols=(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,),
                      unpack=True)
    fit_file = np.genfromtxt(fname, unpack=True, skip_header=True, usecols=(-1, ), dtype=None, encoding="utf-8")

    for i in range(0, len(date)):
        date_time.append(datetime.strptime(date[i] + ' ' + time[i] + "000", "%Y-%m-%d %H:%M:%S.%f"))

    return date_time, x, y, xerr, yerr, flux, flux_err, mag, mag_err, Az, El, Rg, fit_file


def save_jd_data(filename_jd, header,
                 date_time, x, y, xerr, yerr, flux, flux_err, mag, mag_err, Az, El, Rg, fit_file):
    atime = Time(date_time, format='datetime', scale='utc')
    with open(filename_jd, "w") as f_jd:
        for line in header[:-1]:
            f_jd.write(line)

        f_jd.write(
            "#      JD                     X          Y         Xerr      Yerr             Flux     Flux_err     magR  mag_err     Az(deg)   El(deg)   Rg(Mm)    filename\n")

        for i in range(0, len(date_time)):
            f_jd.write(
                f"{atime[i].jd:0<18}        {x[i]:10.5f} {y[i]:10.5f}  " +
                f"{xerr[i]:8.5f}  {yerr[i]:8.5f}     {'{:13.4f}'.format(flux[i])}  {'{:8.4f}'.format(flux_err[i])}   " +
                f"{mag[i]:6.3f}  {mag_err[i]:6.3f}    {Az[i]:8.3f} {El[i]:8.3f}   {Rg[i]:8.3f}   {fit_file[i]}\n"
            )


if __name__ == "__main__":
    filename = sys.argv[1]

    if os.path.isfile(filename):
        filelist = [filename]
    else:
        # we got path with mask
        filelist = glob.glob(filename)

    
    for filename in filelist:
        wd = os.path.dirname(filename)
        base = os.path.basename(filename)
        fname, ext = os.path.splitext(base)

        header = read_header(filename)
        date_time, x, y, xerr, yerr, flux, flux_err, mag, mag_err, Az, El, Rg, fit_file = read_data(filename)

        f_jd_name = os.path.join(wd, fname + "_jd" + ext)
        # print(f_jd_name)

        save_jd_data(f_jd_name, header, date_time, x, y, xerr, yerr, flux, flux_err, mag, mag_err, Az, El, Rg, fit_file)
