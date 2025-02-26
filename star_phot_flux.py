#!python3
import datetime
from datetime import timedelta

from astropy.stats import sigma_clipped_stats
from astropy.io import fits
from photutils.aperture import CircularAperture, CircularAnnulus
# from photutils import aperture_photometry
from photometry_with_errors import *
from sp_utils import *
import glob
import sys
import os
import warnings
import subprocess
# import matplotlib
# matplotlib.use('Agg')

'''
Script to get flux for star objects
'''

if len(sys.argv) < 2:
    print("Not enough parameters. Enter path")
    sys.exit()


warnings.filterwarnings("ignore")

path = sys.argv[1]
# iniList = glob.glob('*config*.ini')

conf = read_config_star_flux(os.path.join(path, 'config_star_flux.ini'))

if not conf:
    print('No config. Exit')
    sys.exit()

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

header = fits.getheader(path + "//" + fl[0])
gain = str(header.get('GAIN'))
exp = str(header.get('EXPTIME')).replace(".", "_")

if conf['star_name'] != "":
    res_path = os.path.join(path, "result_" + conf['star_name'] + "_" + ymd + "_UT" + ut1 + "_g"+gain + "_exp"+exp +".ph" + conf['Filter'])
else:
    res_path = os.path.join(path, "result" + "_" + ymd + "_UT" + ut1 + "_g"+gain + "_exp"+exp + ".ph" + conf['Filter'])

fr = open(res_path, "w")
# fr.write("     Date              UT                   X                 Y                Xerr          Yerr                 Flux                filename\n")

# from tqdm.auto import tqdm
n_fit = len(fl)

flux_array, mag_array = [], []

for fit_file in fl:
    perc = fl.index(fit_file)/(n_fit-1) * 100
    print(f"{perc:5.2f}%  filename = {fit_file}", end=" ")
    # hdu = fits.open(path + "//" + fit_file)[0]
    # data = hdu.data
    # header = hdu.header
    data = fits.getdata(path + "//" + fit_file)
    header = fits.getheader(path + "//" + fit_file)

    gain = header.get('GAIN')
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

    mean, median, std = sigma_clipped_stats(data[20:, :], sigma=3.0)
    # mean, median, std = sigma_clipped_stats(data, sigma=3.0)
    # print (mean, median, std)

    if fit_file == fl[0]:  # make header--------------------------------------------------------------------
        # Az, El = get_star_el(star_name=conf['star_name'],
        #                  obs_lat=conf['site_lat'], obs_lon=conf['site_lon'], obs_elev=conf['site_elev'],
        #                  obs_date=date_time)

        fr.write(f"# Star: {conf['star_name']}\n")


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

        fr.write("# SITE_NAME   = %s\n" % conf['site_name'])
        fr.write("# SITE_LAT   = %s\n" % conf['site_lat'])
        fr.write("# SITE_LON   = %s\n" % conf['site_lon'])
        fr.write("# SITE_ELEV  = %s\n" % conf['site_elev'])

        if conf['time_format'] == "UT":
            fr.write(f"#  Date       UT              X          Y         Xerr      Yerr             Flux     Flux_err")
            fr.write(f"     mag{conf['Filter']}  mag_err     Az(deg)   El(deg)     filename\n")
        else:
            fr.write(f"#      JD                     X          Y         Xerr      Yerr             Flux     Flux_err")
            fr.write(f"     mag{conf['Filter']}  mag_err     Az(deg)   El(deg)     filename\n")

    ##################################

    # BEGIN----
    if conf['dark_frame']:
        dark_arr = fits.getdata(conf['dark_frame'])
        ph_image = substract(data, dark=dark_arr)

        mean2, median2, std2 = sigma_clipped_stats(ph_image[20:, :], sigma=3.0)
        # print (mean2, median2, std2)
        data = substract(ph_image, value=median2)
        # data = ph_image
    else:
        ph_image = data
        data = substract(data, value=median)

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
            print(f"final = {x0:d},{y0:d}") #x0, y0)

        except Exception as E:
            print(E)
            print("Error - curve_fit failed\n")

    # ------------------------------------ PHOTOMETRY-------------------------------------------------
    # print(target)
    if (target) and (target[2] < conf['max_center_error']):  # and (target[-1] > min_signal):
        positions = target[:2]
        aperture = CircularAperture(positions, r=conf['r_ap'])
        annulus_aperture = CircularAnnulus(positions, r_in=conf['an_in'], r_out=conf['an_out'])

        apers = [aperture, annulus_aperture]
        # print (apers)
        # phot_table = aperture_photometry(ph_image, apers)

        # -------------------------------------------------------------
        # bgr_aperture = CircularAperture(positions, r=an_in)
        phot_table = iraf_style_photometry(aperture, annulus_aperture, ph_image, bg_method='mean')
        #-----------------------------------------------------------------------

        date, time = date_time.split("T")
        if len(target) == 5:
            xerr, yerr = target[2], target[3]
        else:
            xerr, yerr = 9, 9

        if (xerr is np.inf) or (yerr is np.inf):
            xerr, yerr = 8, 8

        # flux = phot_table['residual_aperture_sum'][z]

        flux = phot_table['flux'][0]

        flux_array.append(flux)

        # print (flux)
        # print(type(flux))
        flux_err = phot_table['flux_error'][0]
        mag = phot_table['mag']
        mag_err = phot_table['mag_error'][0]

        # print(phot_table['residual_aperture_sum'])
        # print (phot_table)

        Az, El = get_star_el(star_name=conf['star_name'],
                         obs_lat=conf['site_lat'], obs_lon=conf['site_lon'], obs_elev=conf['site_elev'],
                         obs_date=date_time)
        if El < 5:
            print("WARNING! Elevation of satellite < 5 deg. Check settings!")
        mag = calc_mag(flux, El, 1000, conf['A'], conf['k'], exp, min_mag=conf['min_real_mag'])
        mag_array.append(mag)

        if (mag <= conf['min_real_mag']) and (conf['time_format'] == "UT"):
            fr.write(f"{date} {time[:12]}   {phot_table['X'][0]:10.5f} {phot_table['Y'][0]:10.5f}  ")
            fr.write(f"{xerr:8.5f}  {yerr:8.5f}     ")
            fr.write(f"{'{:13.4f}'.format(flux)}  {'{:8.4f}'.format(flux_err)}   {mag:6.3f}  {mag_err:6.3f}    ")
            fr.write(f"{Az:8.3f} {El:8.3f}    {fit_file}\n")
        elif (mag <= conf['min_real_mag']) and (conf['time_format'] == "JD"):
            from astropy.time import Time
            date_time_jd = Time(date_time, format='isot', scale='utc')

            fr.write(f"{date_time_jd.jd:<23}   ")
            fr.write(f"{phot_table['X'][0]:10.5f} {phot_table['Y'][0]:10.5f}  {xerr:8.5f}  {yerr:8.5f}     ")
            fr.write(f"{'{:13.4f}'.format(flux)}  {'{:8.4f}'.format(flux_err)}   {mag:6.3f}  {mag_err:6.3f}    ")
            fr.write(f"{Az:8.3f} {El:8.3f}    {fit_file}\n")
        else:
            print(f"WARNING! mag value < {conf['min_real_mag']} mag, skipping this value!")

fr.write("\n")
fr.write("              Mean         StDev\n")
fr.write(f"Flux:      {np.mean(flux_array):8.3f}     {np.std(flux_array):8.3f}\n")
fr.write(f"Mag:        {np.mean(mag_array):8.3f}     {np.std(mag_array):8.3f}\n")

fr.close()
print("Photometry is DONE.")


