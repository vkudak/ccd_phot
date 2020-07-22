#!python38
from astropy.stats import sigma_clipped_stats
from astropy.io import fits
import numpy as np
from photutils import CircularAperture, CircularAnnulus
from photutils import aperture_photometry
from sp_utils import *
import math
import glob
import sys
import os
import configparser


def calc_mag(flux, el, rg, A, k):
    m_inst = -2.5 * math.log10(flux)
    mz = 1 / (math.cos(math.pi / 2 - math.radians(el)))  # 1/cos(90-el)
    mr = -5 * math.log10(rg / 1000.0)
    mzr = k * mz

    mag = A + m_inst + mzr + mr
    return mag


path = sys.argv[1]

config = configparser.ConfigParser()
config.read(path + '//config_sat.ini')
if os.path.isfile(path + '//config_sat.ini'):
    try:
        cospar = config['NAME']['COSPAR']
        norad = config['NAME']['NORAD']
        name = config['NAME']['NAME']

        tle_file = config['TLE']['tle_file']

        A = float(config['STD']['A'])
        k = float(config['STD']['k'])

        r_ap = float(config['APERTURE']['ap'])
        an_in = float(config['APERTURE']['an_in'])
        an_out = float(config['APERTURE']['an_out'])

    except Exception as E:
        print("Error in inin file\n", E)
        sys.exit()
else:
    print ("Error. Cant find config_sat.ini in " + path + '//config_sat.ini')

tle_list = get_tle(tle_file)


list = os.listdir(path)
fl = []
for fn in list:
    f = fn.split('.')
    if f[-1] in ['FIT', 'FITS', 'fit', 'fits']:
        fl.append(l)

fl.sort()

debug = True

if debug:
    if os.path.isdir(path + "//fig"):
        files = glob.glob(path + "//fig//*.*")
        for f in files:
            os.remove(f)
    else:
        os.makedirs(path + "//fig")

fl = fl[:5]   # first 10 files  ------  file # -1


# fl = ["Capture_00016.fits"]

fr = open(path + "//result.txt", "w")
# fr.write("     Date              UT                   X                 Y                Xerr          Yerr                 Flux                filename\n")

for fit_file in fl:
    print ("filename=", fit_file)
    hdu = fits.open(path + "//" + fit_file)[0]
    data = hdu.data
    header = hdu.header
    date_time = header.get('DATE-OBS')
    if date_time == '0001-01-01T00:00:00.0000000':
        date_time = header.get('DATE-END')  # NO GPS time data !!!!!!!!!!!!
    width = int(header.get('NAXIS1'))
    height = int(header.get('NAXIS2'))
    t_x, t_y = None, None
    try:
        t_x = float(header.get('OBJX'))
        t_y = float(header.get('OBJY'))
    except Exception:
        pass

    mean, median, std = sigma_clipped_stats(data, sigma=3.0)

    if fit_file == fl[0]:  # make header--------------------------------------------------------------------
        El, Rg, Az, name, nor, cosp, tle_lines = calc_from_tle(tle_list, date_time, cospar, norad, name)
        fr.write("# TLE:\n")
        fr.write("# %s\n" % tle_lines[0])
        fr.write("# %s\n" % tle_lines[1])
        fr.write("# %s\n" % tle_lines[2])
        fr.write("# NORAD  = %s\n" % nor)
        fr.write("# COSPAR = %s\n" % cosp)
        fr.write("# NAME   = %s\n" % name)

        fr.write("   Date       UT              X          Y         Xerr      Yerr          Flux         mag      filename\n")

    ##################################

    # BEGIN----
    ph_image = data
    data = data - median

    if t_x is not None:  # target from header
        target = [t_x, height - t_y]
        x0, y0 = int(target[0]), int(target[1])
        error = [0, 0]
        target.append(error[0])
        target.append(error[1])

        fig_name = path + "//fig//" + fit_file + "_man.png"
        target = fit_m(data, x0, y0, gate=20, debug=debug, fig_name=fig_name)

    else:  # SEARCH target
        target = None
        fig_name = path + "//fig//" + fit_file + ".png"
        try:
            target = fit_m(data, x0, y0, gate=20, debug=debug, fig_name=fig_name)

            target_original = target

            x0, y0 = int(target[0]), int(target[1])
            err = target[2]
            if (err > 0.9) and (err < 1.5):
                target = fit_m(data, x0, y0, gate=10, debug=debug, fig_name=fig_name, centring=True)
            elif err > 1.5:
                target = fit_m(data, x0, y0, gate=20, debug=debug, fig_name=fig_name)

                x0, y0 = int(target[0]), int(target[1])  # one more time with smaller window
                target = fit_m(data, x0, y0, gate=10, debug=debug, fig_name=fig_name, centring=True)

            err = target[2]
            if err > 1:
                target = target_original
                x0, y0 = int(target[0]), int(target[1])
            else:
                x0, y0 = int(target[0]), int(target[1])
            print ("final=", x0, y0)
        except Exception as E:
            print(E)
            print("Error - curve_fit failed\n")

    # ------------------------------------ PHOTOMETRY-------------------------------------------------
    if (target) and (target[2] < 2):
        positions = target[:2]
        aperture = CircularAperture(positions, r=r_ap)
        annulus_aperture = CircularAnnulus(positions, r_in=an_in, r_out=an_out)

        # from photutils import EllipticalAperture, EllipticalAnnulus
        # aperture = EllipticalAperture(positions, a=4, b=3, theta=-0.25*np.pi)
        # annulus_aperture = EllipticalAnnulus(positions, a_in=12., a_out=14., b_out=8., theta=-0.25*np.pi)

        apers = [aperture, annulus_aperture]
        # print (apers)
        phot_table = aperture_photometry(ph_image, apers)
        for col in phot_table.colnames:
            phot_table[col].info.format = '%.8g'  # for consistent table output
        # print(phot_table)

        bkg_mean = phot_table['aperture_sum_1'] / annulus_aperture.area
        bkg_sum = bkg_mean * aperture.area
        final_sum = phot_table['aperture_sum_0'] - bkg_sum
        phot_table['residual_aperture_sum'] = final_sum
        phot_table['residual_aperture_sum'].info.format = '%.8g'  # for consistent table output

        z = 0
        if len(phot_table) > 1:
            if math.isnan(phot_table['residual_aperture_sum'][z]):
                z = 1
            for i in range(0, len(phot_table)):
                if not math.isnan(phot_table['residual_aperture_sum'][i]):
                    if phot_table['residual_aperture_sum'][i] > phot_table['residual_aperture_sum'][z]:
                        z = i

        date, time = date_time.split("T")
        if len(target) == 4:
            xerr, yerr = target[2], target[3]
        else:
            xerr, yerr = 9, 9

        if (xerr is np.inf) or (yerr is np.inf):
            xerr, yerr = 8, 8

        flux = phot_table['residual_aperture_sum'][z]
        # print (phot_table)

        El, Rg, Az, name, nor, cosp, tle_lines = calc_from_tle(tle_list, date_time, cospar, norad, name)
        mag = calc_mag(flux, El, Rg, A, k)
        # fr.write("%s %s    %8.5f  %8.5f  %8.5f  %8.5f  %12.5f   %s\n" % (date, time[:12], phot_table['xcenter'][z].value, phot_table['ycenter'][z].value, xerr, yerr, flux, fit_file))
        fr.write("%s %s    %8.5f  %8.5f  %8.5f  %8.5f     %s   %8.5f   %s\n" % (date, time[:12], phot_table['xcenter'][z].value, phot_table['ycenter'][z].value, xerr, yerr, '{:010.4f}'.format(flux), mag, fit_file))

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
