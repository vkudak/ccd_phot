#!python3.8
import datetime
from datetime import timedelta

from astropy.stats import sigma_clipped_stats
from astropy.io import fits
import numpy as np
from photutils import CircularAperture, CircularAnnulus
from photutils import aperture_photometry
from photometry_with_errors import *
from sp_utils import *
import math
import glob
import sys
import os
import configparser
import warnings
# import matplotlib
# matplotlib.use('Agg')


def calc_mag(flux, el, rg, A, k, exp, min_mag=15):
    # print (flux)
    if flux < 0:
        mag = min_mag
        return mag
    else:
        m_inst = -2.5 * math.log10(flux/exp)
        mz = 1 / (math.cos(math.pi / 2 - math.radians(el)))  # 1/cos(90-el)
        mr = -5 * math.log10(rg / 1000.0)
        mzr = k * mz

        mag = A + m_inst + mzr + mr
    return mag


if len(sys.argv) < 2:
    print("Not enough parameters. Enter path")
    sys.exit()


warnings.filterwarnings("ignore")

path = sys.argv[1]

# iniList = glob.glob('*config*.ini')

config = configparser.ConfigParser(inline_comment_prefixes="#")

if os.path.isfile(path + '//config_sat.ini'):
    try:
        config.read(path + '//config_sat.ini')

        cospar = config['NAME']['COSPAR']
        norad = config['NAME']['NORAD']
        name = config['NAME']['NAME']

        tle_file = config['TLE']['tle_file']

        A = float(config['STD']['A'])
        k = float(config['STD']['k'])
        gate = int(config['STD']['gate'])
        Filter = config.get('STD','Filter', fallback="F")
        Filter = Filter.strip()
        time_format = config.get('STD', 'Time_format', fallback="UT")
        time_format = time_format.strip()
        min_real_mag = config.getfloat('STD', 'min_real_mag', fallback=15.0)
        max_center_error = config.getfloat('STD', 'max_center_err', fallback=2.0)

        try:
            dark_frame = config['STD']['dark_frame']
        except Exception:
            dark_frame = False

        r_ap = float(config['APERTURE']['ap'])
        an_in = float(config['APERTURE']['an_in'])
        an_out = float(config['APERTURE']['an_out'])

        site_name = config['SITE']['Name']
        site_lat = config['SITE']['lat']
        site_lon = config['SITE']['lon']
        site_elev = float(config['SITE']['h'])

    except Exception as E:
        print("Error in INI file\n", E)
        sys.exit()
else:
    print("Error. Cant find config_sat.ini in " + os.path.join(path, "config_sat.ini"))
    sys.exit()

tle_list = get_tle(tle_file)

list = os.listdir(path)
fl = []
for fn in list:
    f = fn.split('.')
    if f[-1] in ['FIT', 'FITS', 'fit', 'fits']:
        fl.append(fn)
fl.sort()

debug = True

if debug:
    if os.path.isdir(path + "//fig"):
        files = glob.glob(path + "//fig//*.*")
        for f in files:
            os.remove(f)
    else:
        os.makedirs(path + "//fig")

# fl = fl[:10]   # first 10 files  ------  file # -1
# fl = ["Capture_00016.fits"]

# get TIME from first FIT file
header = fits.getheader(path + "//" + fl[0])
date_time = header.get('DATE-OBS')
if date_time == '0001-01-01T00:00:00.0000000':
    date_time = header.get('DATE-END')  # NO GPS time data !!!!!!!!!!!!
ut1 = date_time.split("T")[1]
ut1 = ut1[:2] + ut1[3:5] + ut1[6:8]
ymd = date_time.split("T")[0]
ymd = ymd.replace("-", "")
##############

if norad != "":
    fr = open(path + "//result_" + norad + "_" + ymd + "_UT" + ut1 + ".ph" + Filter, "w")
else:
    fr = open(path + "//result" + "_" + ymd + "_UT" + ut1 + ".ph" + Filter, "w")
# fr.write("     Date              UT                   X                 Y                Xerr          Yerr                 Flux                filename\n")

# from tqdm.auto import tqdm
n_fit = len(fl)
for fit_file in fl:
    perc = fl.index(fit_file)/(n_fit-1) * 100
    print(f"{perc:5.2f}%  filename = {fit_file}", end=" ")
    # hdu = fits.open(path + "//" + fit_file)[0]
    # data = hdu.data
    # header = hdu.header
    data = fits.getdata(path + "//" + fit_file)
    header = fits.getheader(path + "//" + fit_file)

    date_time = header.get('DATE-OBS')
    exp = header.get('EXPTIME')
    exp = float(exp)
    if date_time == '0001-01-01T00:00:00.0000000':
        date_time = header.get('DATE-END')  # NO GPS time data !!!!!!!!!!!!

    # correction DATE_TIME, datetime = datetime + exp / 2
    date_time = datetime.strptime(date_time[:-1], "%Y-%m-%dT%H:%M:%S.%f")
    date_time = date_time + timedelta(seconds=float(exp/2.0))
    # and back to STR
    date_time = date_time.strftime("%Y-%m-%dT%H:%M:%S.%f0")

    width = int(header.get('NAXIS1'))
    height = int(header.get('NAXIS2'))
    t_x, t_y = None, None
    try:
        t_x = float(header.get('OBJX'))
        t_y = float(header.get('OBJY'))
    except Exception:
        pass

    mean, median, std = sigma_clipped_stats(data[20:, :], sigma=3.0)
    # mean, median, std = sigma_clipped_stats(data, sigma=3.0)
    # print (mean, median, std)

    if fit_file == fl[0]:  # make header--------------------------------------------------------------------
        El, Rg, Az, name, nor, cosp, tle_lines = calc_from_tle(site_lat, site_lon, site_elev, tle_list, date_time, cospar, norad, name)
        fr.write("# TLE:\n")
        fr.write("# %s\n" % tle_lines[0])
        fr.write("# %s\n" % tle_lines[1])
        fr.write("# %s\n" % tle_lines[2])

        date, time = date_time.split("T")
        fr.write("# %s %s\n" % (date, time))

        header_last = fits.getheader(path + "//" + fl[-1])
        date_time9 = header.get('DATE-OBS')
        date9, time9 = date_time9.split("T")
        fr.write("# %s %s\n" % (date9, time9))

        # take EXPTIME from the middle of LC
        exp_header = fits.getheader(path + "//" + fl[int(len(fl)/2)])
        fr.write("# dt = %s\n" % exp_header.get('EXPTIME'))
        fr.write("# Filter = %s\n" % Filter)

        fr.write("# COSPAR = %s\n" % cosp)
        fr.write("# NORAD  = %s\n" % nor)
        fr.write("# NAME   = %s\n" % name)

        fr.write("# SITE_NAME   = %s\n" % site_name)
        fr.write("# SITE_LAT   = %s\n" % site_lat)
        fr.write("# SITE_LON   = %s\n" % site_lon)
        fr.write("# SITE_ELEV  = %s\n" % site_elev)

        if time_format == "UT":
            fr.write(f"#  Date       UT              X          Y         Xerr      Yerr             Flux     Flux_err     mag{Filter}  mag_err     Az(deg)   El(deg)   Rg(Mm)    filename\n")
        else:
            fr.write(
                f"#      JD                     X          Y         Xerr      Yerr             Flux     Flux_err     mag{Filter}  mag_err     Az(deg)   El(deg)   Rg(Mm)    filename\n")

    ##################################

    # BEGIN----
    if dark_frame:
        dark_arr = fits.getdata(dark_frame)
        ph_image = substract(data, dark=dark_arr)

        mean2, median2, std2 = sigma_clipped_stats(ph_image[20:, :], sigma=3.0)
        # print (mean2, median2, std2)
        data = substract(ph_image, value=median2)
        # data = ph_image
    else:
        ph_image = data
        data = substract(data, value=median)

    # data = data - median

    # import matplotlib.pyplot as plt
    # from matplotlib.patches import Circle
    # from astropy.visualization import SqrtStretch, LogStretch
    # from astropy.visualization.mpl_normalize import ImageNormalize
    # fig, ax = plt.subplots()
    # # ax.imshow(data, origin='lower')
    # norm = ImageNormalize(stretch=LogStretch())
    # # norm = ImageNormalize(stretch=SqrtStretch())
    # plt.imshow(data[20:,:], cmap='Greys', origin='lower', norm=norm)
    # # ax.scatter(target[0], target[1], s=20, c='red', marker='x')
    # # circle = Circle((target[0], target[1]), 8, facecolor='none', edgecolor='red', linewidth=1, fill=False)
    # # ax.add_patch(circle)
    # plt.show()
    # plt.close('all')
    # sys.exit()

    # gate = 20
    # min_signal = 100

    gate2 = int(gate / 2)

    if t_x is not None:  # target from header
        target = [t_x, height - t_y]
        x0, y0 = int(target[0]), int(target[1])
        error = [0, 0]
        target.append(error[0])
        target.append(error[1])

        try:
            fig_name = path + "//fig//" + fit_file + "_man.png"
            target = fit_m(data, x0, y0, gate=gate2, debug=debug, fig_name=fig_name, centring=True)
        except Exception as E:
            print(E)
            print("Error - curve_fit failed\n")

    else:  # SEARCH target
        target = None
        fig_name = path + "//fig//" + fit_file + ".png"
        try:
            target = fit_m(data, x0, y0, gate=gate, debug=debug, fig_name=fig_name)  # , centring=True)
            # print (target[-1])
            # if target[-1] > min_signal:

            target_original = target

            x0, y0 = int(target[0]), int(target[1])
            err = target[2]
            if (err > 0.9) and (err < 1.5):
                target = fit_m(data, x0, y0, gate=gate2, debug=debug, fig_name=fig_name, centring=True)
            elif err > 1.5:
                target = fit_m(data, x0, y0, gate=gate, debug=debug, fig_name=fig_name)

                x0, y0 = int(target[0]), int(target[1])  # one more time with smaller window
                target = fit_m(data, x0, y0, gate=gate2, debug=debug, fig_name=fig_name, centring=True)

            err = target[2]
            if err > 1:
                target = target_original
                x0, y0 = int(target[0]), int(target[1])
            else:
                x0, y0 = int(target[0]), int(target[1])
            print(f"final = {x0:d},{y0:d}") #x0, y0)

        except Exception as E:
            print(E)
            print("Error - curve_fit failed\n")

    # ------------------------------------ PHOTOMETRY-------------------------------------------------
    # print(target)
    if (target) and (target[2] < max_center_error):  # and (target[-1] > min_signal):
        positions = target[:2]
        aperture = CircularAperture(positions, r=r_ap)
        annulus_aperture = CircularAnnulus(positions, r_in=an_in, r_out=an_out)

        apers = [aperture, annulus_aperture]
        # print (apers)
        # phot_table = aperture_photometry(ph_image, apers)

        # -------------------------------------------------------------
        # bgr_aperture = CircularAperture(positions, r=an_in)
        phot_table = iraf_style_photometry(aperture, annulus_aperture, ph_image, bg_method='mean')
        #-----------------------------------------------------------------------

        # for col in phot_table.colnames:
        #     phot_table[col].info.format = '%.8g'  # for consistent table output
        # # print(phot_table)

        # bkg_mean = phot_table['aperture_sum_1'][0] / annulus_aperture.area
        #
        # # print ("app_sum1=", phot_table['aperture_sum_1'][0])
        # # print ("ann_app_area=", annulus_aperture.area)
        #
        # # print ("bkg_mean = ", bkg_mean)
        # bkg_sum = bkg_mean * aperture.area
        # # print ("bkg_sum = ", bkg_sum)
        #
        # # print("sum and bkgr", phot_table['aperture_sum_0'][0], bkg_sum)
        # final_sum = phot_table['aperture_sum_0'][0] - bkg_sum
        # # print("fin_sum2=", final_sum)
        # phot_table['residual_aperture_sum'] = final_sum
        # # phot_table['residual_aperture_sum'].info.format = '%.8f'  # for consistent table output
        # # print(phot_table['residual_aperture_sum'][0])


        # z = 0
        # if len(phot_table) > 1:
        #     if math.isnan(phot_table['residual_aperture_sum'][z]):
        #         z = 1
        #     for i in range(0, len(phot_table)):
        #         if not math.isnan(phot_table['residual_aperture_sum'][i]):
        #             if phot_table['residual_aperture_sum'][i] > phot_table['residual_aperture_sum'][z]:
        #                 z = i

        date, time = date_time.split("T")
        if len(target) == 5:
            xerr, yerr = target[2], target[3]
        else:
            xerr, yerr = 9, 9

        if (xerr is np.inf) or (yerr is np.inf):
            xerr, yerr = 8, 8


        # flux = phot_table['residual_aperture_sum'][z]

        flux = phot_table['flux'][0]
        # print (flux)
        # print(type(flux))
        flux_err = phot_table['flux_error'][0]
        mag = phot_table['mag']
        mag_err = phot_table['mag_error'][0]

        # print(phot_table['residual_aperture_sum'])
        # print (phot_table)

        El, Rg, Az, name, nor, cosp, tle_lines = calc_from_tle(site_lat, site_lon, site_elev, tle_list, date_time, cospar, norad, name)
        if El < 5:
            print("WARNING! Elevation of satellite < 5 deg. Check settings!")
        mag = calc_mag(flux, El, Rg, A, k, exp, min_mag=min_real_mag)
        # fr.write("%s %s    %8.5f  %8.5f  %8.5f  %8.5f  %12.5f   %s\n" % (date, time[:12], phot_table['xcenter'][z].value, phot_table['ycenter'][z].value, xerr, yerr, flux, fit_file))
        # fr.write("%s %s    %8.5f  %8.5f  %8.5f  %8.5f     %s   %6.3f    %8.3f %8.3f   %8.3f   %s\n" %
            # (date, time[:12], phot_table['xcenter'][z].value, phot_table['ycenter'][z].value, xerr, yerr, '{:13.4f}'.format(flux), mag, Az, El, Rg, fit_file))
        # print(mag)
        # print(mag < min_real_mag, mag > min_real_mag)

        if (mag <= min_real_mag) and (time_format == "UT"):
            fr.write(
                f"{date} {time[:12]}   {phot_table['X'][0]:10.5f} {phot_table['Y'][0]:10.5f}  {xerr:8.5f}  {yerr:8.5f}     {'{:13.4f}'.format(flux)}  {'{:8.4f}'.format(flux_err)}   {mag:6.3f}  {mag_err:6.3f}    {Az:8.3f} {El:8.3f}   {Rg:8.3f}   {fit_file}\n")
        elif (mag <= min_real_mag) and (time_format == "JD"):
            from astropy.time import Time
            date_time_jd = Time(date_time, format='isot', scale='utc')

            fr.write(
                f"{date_time_jd.jd:<23}   {phot_table['X'][0]:10.5f} {phot_table['Y'][0]:10.5f}  {xerr:8.5f}  {yerr:8.5f}     {'{:13.4f}'.format(flux)}  {'{:8.4f}'.format(flux_err)}   {mag:6.3f}  {mag_err:6.3f}    {Az:8.3f} {El:8.3f}   {Rg:8.3f}   {fit_file}\n")
        else:
            print(f"WARNING! mag value < {min_real_mag} mag, skipping this value!")
        # PLOT GENERAL FIT with apperture
        # import matplotlib.pyplot as plt
        # from matplotlib.patches import Circle
        # from astropy.visualization import SqrtStretch, LogStretch
        # from astropy.visualization.mpl_normalize import ImageNormalize
        # fig, ax = plt.subplots()
        # ax.imshow(data, origin='lower')
        # norm = ImageNormalize(stretch=LogStretch())
        # plt.imshow(data, cmap='Greys', origin='lower', norm=norm)
        # # ax.scatter(target[0], target[1], s=20, c='red', marker='x')
        # circle = Circle((target[0], target[1]), 8, facecolor='none', edgecolor='red', linewidth=1, fill=False)
        # ax.add_patch(circle)
        # plt.show()
        # plt.close('all')
