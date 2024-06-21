import sys
import os
import numpy as np

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
    jd = np.genfromtxt(fname, unpack=True, skip_header=True, usecols=(0,), encoding="utf-8")
    x, y, xerr, yerr, flux, flux_err, mag, mag_err, Az, El, Rg = \
        np.genfromtxt(filename, skip_header=True,
                      usecols=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,),
                      unpack=True)
    fit_file = np.genfromtxt(fname, unpack=True, skip_header=True, usecols=(-1, ), dtype=None, encoding="utf-8")

    jd_time = Time(jd, format='jd', scale='utc')

    return jd_time, x, y, xerr, yerr, flux, flux_err, mag, mag_err, Az, El, Rg, fit_file


def save_ut_data(filename_ut, header,
                 jd_time, x, y, xerr, yerr,
                 flux, flux_err, mag, mag_err,
                 Az, El, Rg, fit_file):
    with open(filename_ut, "w") as f_jd:
        for line in header[:-1]:
            f_jd.write(line)

        f_jd.write(
            "#  Date       UT              X          Y         Xerr      Yerr             Flux     Flux_err     magR  mag_err     Az(deg)   El(deg)   Rg(Mm)    filename\n")

        for i in range(0, len(jd_time)):
            date_time = jd_time[i].datetime.strftime("%Y-%m-%dT%H:%M:%S.%f")
            date, time = date_time.split("T")
            f_jd.write(
                f"{date} {time[:12]}   {x[i]:10.5f} {y[i]:10.5f}  " +
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

        # print(fname, fname[-2:])

        if fname[-2:] == "jd":
            fname = fname[:-3]
        # print(fname)
        f_jd_name = os.path.join(wd, fname + "_ut" + ext)
        # print(f_jd_name)

        save_ut_data(f_jd_name, header, date_time, x, y, xerr, yerr, flux, flux_err, mag, mag_err, Az, El, Rg, fit_file)
