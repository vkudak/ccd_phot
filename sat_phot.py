#!python3
import datetime
from datetime import timedelta

from astropy.stats import sigma_clipped_stats
from astropy.io import fits
from photutils.aperture import CircularAperture, CircularAnnulus, ApertureStats
# from photutils import aperture_photometry
from photometry_with_errors import *
from sp_utils import *
from obj_finder import obj_finder_dao
import glob
import sys
import os
import warnings
import subprocess
import argparse
# import matplotlib
# matplotlib.use('Agg')


parser = argparse.ArgumentParser(description='Photometry of LEO satellites')
# parser.add_argument('-p', '--path', help='Path to fits files', required=True)
parser.add_argument('path', type=str, help='Path to fits files')
parser.add_argument('-c', '--config', help='Specify config file', required=False)
parser.add_argument('-s_file', '--start_file', help='Specify file to start from', required=False)
args = vars(parser.parse_args())

if len(sys.argv) < 2:
    print("Not enough parameters. Enter path")
    sys.exit()


warnings.filterwarnings("ignore")

# path = sys.argv[1]
if args["path"]:
    path = args["path"]
else:
    sys.exit("Not enough parameters. Enter path")

if args["start_file"]:
    start_file = args["start_file"]
else:
    start_file = None

if args["config"]:
    conf = read_config_sat(args["config"])
    conf_name = args["config"]
else:
    print("Search for configuration in working path")
    # search for config files in working dir. Config starts from config_sat_XXXXX.ini
    conf_list = glob.glob(os.path.join(path, 'config_sat*.ini'))
    if len(conf_list) == 0:
        print("No configuration found")
        sys.exit()
    # conf = read_config_sat(os.path.join(path, 'config_sat.ini'))

    # Load first config in list
    conf = read_config_sat(conf_list[0])
    conf_name = conf_list[0]

if not conf:
    print('No config. Exit')
    sys.exit()
else:
    # print(f"Using config '{conf_list[0]}'")
    print(f"Using config '{conf_name}'")

tle_list = get_tle(conf['tle_file'])

fl = get_file_list(path)
# print('Sorting FITS by name...')

if conf['fits_sort'] not in ['name', 'None']:
    print(f'Sorting FITS files by header {conf["fits_sort"]} ...')
    fl = sort_file_list(path, fl, field=conf['fits_sort'])

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

# Start photometry from defined FITS file
if start_file is not None and start_file in fl:
    st_pos = fl.index(start_file)
    fl = fl[st_pos:]

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

if conf['norad'] != "":
    res_path = os.path.join(path, "result_" + conf['norad'] + "_" + ymd + "_UT" + ut1 + ".ph" + conf['Filter'])
else:
    res_path = os.path.join(path, "result" + "_" + ymd + "_UT" + ut1 + ".ph" + conf['Filter'])


fr = open(res_path, "w")
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

    # Fix sec=60
    date_time = fix_datetime(date_time)

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

    # mean, median, std = sigma_clipped_stats(data[20:, :], sigma=3.0)
    # mean, median, std = sigma_clipped_stats(data, sigma=3.0)
    mean, median, std = np.mean(data[5:, :]), np.median(data[5:, :]), np.std(data[5:, :])
    # print(mean, median, std)

    if fit_file == fl[0]:  # make header--------------------------------------------------------------------
        El, Rg, Az, name, nor, cosp, tle_lines, _ = calc_from_tle(conf['site_lat'], conf['site_lon'], conf['site_elev'],
                                                               tle_list,
                                                               date_time,
                                                               conf['cospar'], conf['norad'], conf['name'])
        fr.write("# TLE:\n")
        fr.write("# %s\n" % tle_lines[0])
        fr.write("# %s\n" % tle_lines[1])
        fr.write("# %s\n" % tle_lines[2])

        date, time = date_time.split("T")
        fr.write("# %s %s\n" % (date, time))  # exp corrected start date_time

        header_last = fits.getheader(path + "//" + fl[-1])
        date_time_last = fix_datetime(header_last.get('DATE-OBS'))
        date_time_last = datetime.strptime(date_time_last[:-1], "%Y-%m-%dT%H:%M:%S.%f")
        date_time_last = date_time_last + timedelta(seconds=float(exp / 2.0))
        date_time_last = date_time_last.strftime("%Y-%m-%dT%H:%M:%S.%f")
        date_last, time_last = date_time_last.split("T")
        fr.write("# %s %s\n" % (date_last, time_last))   # exp corrected last date_time

        # take EXPTIME from the middle of LC
        exp_header = fits.getheader(path + "//" + fl[int(len(fl)/2)])
        fr.write("# dt = %s\n" % exp_header.get('EXPTIME'))
        fr.write("# Filter = %s\n" % conf['Filter'])
        fr.write("# Time_format = %s\n" % conf['time_format'])

        fr.write("# COSPAR = %s\n" % cosp)
        fr.write("# NORAD  = %s\n" % nor)
        fr.write("# NAME   = %s\n" % name)

        fr.write("# SITE_NAME   = %s\n" % conf['site_name'])
        fr.write("# SITE_LAT   = %s\n" % conf['site_lat'])
        fr.write("# SITE_LON   = %s\n" % conf['site_lon'])
        fr.write("# SITE_ELEV  = %s\n" % conf['site_elev'])

        if conf['time_format'] == "UT":
            fr.write(f"#  Date       UT              X          Y         Xerr      Yerr             Flux     Flux_err")
            fr.write(f"     mag{conf['Filter']}  mag_err     Az(deg)   El(deg)   Rg(Km)    filename\n")
        else:
            fr.write(f"#      JD                     X          Y         Xerr      Yerr             Flux     Flux_err")
            fr.write(f"     mag{conf['Filter']}  mag_err     Az(deg)   El(deg)   Rg(Km)    filename\n")

        if t_x is None:
            print("\nNo target info in the HEADER of the first FIT file. Trying to search with DAOFind ...")
            t_x, t_y = obj_finder_dao(data[5:, :])
            if t_x and t_y is not None:
                t_y = t_y + 5
                print(f"Object detected at X = {t_x:5.2f} Y = {t_y:5.2f}  OK.  ", end="")
                t_y = height - t_y  # Target_pos coord_system
            else:
                print("No target was found on the image")
                sys.exit()
    ##################################

    # BEGIN----
    if conf['dark_frame']:
        dark_arr = fits.getdata(conf['dark_frame'])
        ph_image = substract(data, dark=dark_arr)

        mean2, median2, std2 = np.mean(ph_image[5:, :]), np.median(ph_image[5:, :]), np.std(ph_image[5:, :])
        data = substract(ph_image, value=median2)
        # # data = ph_image
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

    gate2 = int(conf['gate'] / 2)

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
            target = fit_m(data, x0, y0, gate=conf['gate'], debug=debug, fig_name=fig_name)  # , centring=True)
            # print (target[-1])
            # if target[-1] > min_signal:

            target_original = target

            x0, y0 = int(target[0]), int(target[1])
            err = target[2]
            if (err > 0.9) and (err < 1.5):
                target = fit_m(data, x0, y0, gate=gate2, debug=debug, fig_name=fig_name, centring=True)
            elif err > 1.5:
                target = fit_m(data, x0, y0, gate=conf['gate'], debug=debug, fig_name=fig_name)

                x0, y0 = int(target[0]), int(target[1])  # one more time with smaller window
                target = fit_m(data, x0, y0, gate=gate2, debug=debug, fig_name=fig_name, centring=True)

            err = target[2]
            if err > 1:
                target = target_original
                x0, y0 = int(target[0]), int(target[1])
            else:
                x0, y0 = int(target[0]), int(target[1])
            print(f"final = {x0:d},{y0:d}", end="") #x0, y0)

        except Exception as E:
            print(E)
            print("Error - curve_fit failed\n")

    # ------------------------------------ PHOTOMETRY-------------------------------------------------
    # print(target)
    if (target) and (target[2] < conf['max_center_error']):  # and (target[-1] > min_signal):
        positions = target[:2]
        aperture = CircularAperture(positions, r=conf['r_ap'])
        aper_stats = ApertureStats(data, aperture)
        if aper_stats.max + median == 65535:
            saturated = "*"
            if conf["saturated"]:
                print(f" {saturated}")
            else:
                print(" ", end="")
        else:
            saturated = ""
            print(f" {saturated}")
        annulus_aperture = CircularAnnulus(positions, r_in=conf['an_in'], r_out=conf['an_out'])

        apers = [aperture, annulus_aperture]
        # print (apers)
        # phot_table = aperture_photometry(ph_image, apers)

        # -------------------------------------------------------------
        # bgr_aperture = CircularAperture(positions, r=an_in)
        phot_table = iraf_style_photometry(aperture, annulus_aperture, ph_image, bg_method='mean')
        # -----------------------------------------------------------------------

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

        El, Rg, Az, name, nor, cosp, tle_lines, phase = calc_from_tle(
            conf['site_lat'], conf['site_lon'], conf['site_elev'],
            tle_list, date_time,
            conf['cospar'], conf['norad'], conf['name'])
        if El < 5:
            print("WARNING! Elevation of satellite < 5 deg. Check settings!")
        mag = calc_mag(flux, El, Rg, conf['A'], conf['k'], exp, min_mag=conf['min_real_mag'], phase=phase)
        # fr.write("%s %s    %8.5f  %8.5f  %8.5f  %8.5f  %12.5f   %s\n" % (date, time[:12], phot_table['xcenter'][z].value, phot_table['ycenter'][z].value, xerr, yerr, flux, fit_file))
        # fr.write("%s %s    %8.5f  %8.5f  %8.5f  %8.5f     %s   %6.3f    %8.3f %8.3f   %8.3f   %s\n" %
            # (date, time[:12], phot_table['xcenter'][z].value, phot_table['ycenter'][z].value, xerr, yerr, '{:13.4f}'.format(flux), mag, Az, El, Rg, fit_file))
        # print(mag)
        # print(mag < min_real_mag, mag > min_real_mag)
        if conf['saturated'] == False and (saturated=="*"):
            print(f"WARNING! Aperture is saturated, skipping this value!")
        else:
            if (mag <= conf['min_real_mag']) and (conf['time_format'] == "UT"):
                fr.write(f"{date} {time[:12]}   {phot_table['X'][0]:10.5f} {phot_table['Y'][0]:10.5f}  ")
                fr.write(f"{xerr:8.5f}  {yerr:8.5f}     ")
                fr.write(f"{'{:13.4f}'.format(flux)}  {'{:8.4f}'.format(flux_err)}   {mag:6.3f}  {mag_err:6.3f}    ")
                fr.write(f"{Az:8.3f} {El:8.3f}   {Rg:8.3f}   {fit_file}{saturated}\n")
            elif (mag <= conf['min_real_mag']) and (conf['time_format'] == "JD"):
                from astropy.time import Time
                date_time_jd = Time(date_time, format='isot', scale='utc')

                fr.write(f"{date_time_jd.jd:<23}   ")
                fr.write(f"{phot_table['X'][0]:10.5f} {phot_table['Y'][0]:10.5f}  {xerr:8.5f}  {yerr:8.5f}     ")
                fr.write(f"{'{:13.4f}'.format(flux)}  {'{:8.4f}'.format(flux_err)}   {mag:6.3f}  {mag_err:6.3f}    ")
                fr.write(f"{Az:8.3f} {El:8.3f}   {Rg:8.3f}   {fit_file}{saturated}\n")
            else:
                print(f"WARNING! mag value < {conf['min_real_mag']} mag, skipping this value!")
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
fr.close()
print("Photometry is DONE.")

if conf['auto_plot']:
    print("Plotting Lightcurve...", end="")
    subprocess.run([conf["python"], conf["plot_command"], res_path])
    print("DONE.")
