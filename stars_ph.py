#!python3.8
from astropy.stats import sigma_clipped_stats
from astropy.io import fits
import numpy as np
from numpy.ma import masked
from photutils import CircularAperture, CircularAnnulus
from photutils import aperture_photometry
from sp_utils import *
import math
# import glob
import sys
import os
# import configparser
# from astroquery.vizier import Vizier
# import astropy.units as u
# import astropy.coordinates as coord
# from astropy.table import Table
from astroquery.astrometry_net import AstrometryNet
from astropy.wcs import WCS
import ephem
import configparser
import warnings


def RMS_del(A, value):
    '''Delete elements of array A until A.RMS>value'''
    A = np.array(A)
    A_del = []
    while A.std(axis=0) > value:
        # rms = A.std(axis=0)
        mean = A.mean(axis=0)
        d = []  # X-mean
        maxx = 0
        for i in range(len(A)):
            d.append(abs(A[i] - mean))
            if d[i] > maxx:
                maxx = d[i]
                imax = i
        A_del.append(A[imax])
        A = np.delete(A, imax)

    strFormat = len(A_del) * '{:5.3f}, '
    formattedList = strFormat.format(*A_del)

    log_file.write("Delete A value(s): " + formattedList + "\n")
    return A


if len(sys.argv) < 2:
    print("Not enouth parameters. Enter path")
    sys.exit()


path = sys.argv[1]

warnings.filterwarnings("ignore")

ploting = False  # plot each frame with appertures

if ploting:
    # for plot
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    # from astropy.visualization import LogStretch
    # from astropy.visualization.mpl_normalize import ImageNormalize

station = ephem.Observer()
station.lat = '48.5635505'
station.long = '22.453751'
station.elevation = 231.1325

config = configparser.ConfigParser(inline_comment_prefixes="#")
config.read(path + '//config_stars.ini')
if os.path.isfile(path + '//config_stars.ini'):
    try:
        Cr = config['Stars_Stand']['C']
        kr = config['Stars_Stand']['K']
        Cr = float(Cr)
        kr = float(kr)
        max_m = config['Stars_Stand']['max_m']
        rms_val = float(config['Stars_Stand']['A_rms'])
        c_flag = config['Stars_Stand'].getboolean('calc_C')

        r_ap = float(config['APERTURE']['r_ap'])
        an_in = float(config['APERTURE']['an_in'])
        an_out = float(config['APERTURE']['an_out'])

    except Exception as E:
        print("Error in inin file\n", E)
        sys.exit()
else:
    print ("Error. Cant find config_stars.ini in " + path + '//config_stars.ini')


# print (path)

list = os.listdir(path)
fl = []
for fn in list:
    f = fn.split('.')
    if f[-1] in ['FIT', 'FITS', 'fit', 'fits']:
        fl.append(fn)

fl.sort()

# print (fl)
# fl = fl[:1]


# #################### BEGIN
ast = AstrometryNet()
# ast.show_allowed_settings()

ast.api_key = "ittzfaqwnrvduhax"


if ploting:
    fig, ax = plt.subplots()


log_file = open(path + '//star_phot.log', "w")

A_general = []
c_general = []

for fit_file in fl:
    print(fit_file)
    log_file.write("####################################################\n")
    log_file.write("filename = " + path + "//" + fit_file + "\n")
    try_again = True
    submission_id = None
    hdu = fits.open(path + "//" + fit_file)
    header = hdu[0].header
    date_time = hdu[0].header.get('DATE-OBS')
    exp = header.get('EXPTIME')
    hdu.close()
    author = None
    try:
        author = header.get('AUTHOR')
        if author == "LKD UZhNU":
            log_file.write('File with WCS\n')
    except Exception:
        pass
    if (author is None) or (author not in ["LKD UZhNU"]):
        while try_again:
            try:
                if not submission_id:
                    print (path + "//" + fit_file)
                    wcs_header = ast.solve_from_image(path + "//" + fit_file, submission_id=submission_id, force_image_upload=True, downsample_factor=2, scale_units="arcsecperpix", scale_type='ul', scale_upper=12.0, scale_lower=8, tweak_order=3)
                else:
                    wcs_header = ast.monitor_submission(submission_id, solve_timeout=120)
            except TimeoutError as e:
                submission_id = e.args[1]
            else:
                # got a result, so terminate
                try_again = False

        if wcs_header:
            # Code to execute when solve succeeds
            print("OK")
            log_file.write("file SOLVED\n")
            with fits.open(path + "//" + fit_file, mode='update') as hdul:
                hdul[0].header = wcs_header
                hdul[0].header.append(('AUTHOR', "LKD UZhNU", 'Solved sucessfuly with astrometry.net'))
                hdul[0].header.append(('DATE-OBS', date_time, "System Clock:Est. Frame Start -OR- GPS:Start Exposure"))
                hdul[0].header.append(('EXPTIME', exp, "EXPOSURE in seconds"))
                hdul.flush()
                hdul.close()

        else:
            # Code to execute when solve fails
            print ("Fail")
            log_file.write("file NOT SOLVED\n")

    # BEGIN star find---------------------------------------------------------------
    # 1. Get RA DEC of image center
    # https://python4astronomers.github.io/astropy/wcs.html
    # 2. Get stars from NOMAD  (done!!!)
    # 3. transform stars RA, DEC -> X,Y. For each star
    # 4. calculate star flux at X,Y
    # 5. make result in flux, Rmag, Mz ??????
    log_file.write('Begin star find procedure...\n')

    header = fits.getheader(path + "//" + fit_file)
    image_tmp = fits.getdata(path + "//" + fit_file)
    mean, median, std = sigma_clipped_stats(image_tmp, sigma=3.0)
    image_tmp = image_tmp - mean

    xc = header["NAXIS1"] / 2.
    yc = header["NAXIS2"] / 2.

    w = WCS(header)
    ra_c, dec_c = w.wcs_pix2world(xc, yc, 1)  # RA DEC of FRAME center
    print("Grab stars from Vizier....")
    log_file.write("Grab stars from Vizier....\n")

    table_res = get_from_NOMAD(ra_c, dec_c, radius="3.0deg", Filter={'Vmag': '<' + max_m})
    table_res = table_res["I/297/out"]  # get NOMAD data

    table_res.remove_columns(["YM", 'r', 'pmRA', 'e_pmRA', 'pmDE', 'e_pmDE', 'Jmag', 'Hmag', 'Kmag', 'R', 'r_Bmag', 'r_Vmag', 'r_Rmag'])
    table_res.sort(["Vmag"])

    len_all = len(table_res)
    log_file.write('Find %i stars\n' % len_all)

    # Write code to delete variable stars from table_res
    print("Deleting variable stars from list....")
    log_file.write('Deleting variable stars...\n')
    table_res = del_var(table_res, Filter={'Vmag': '<' + max_m})
    print ("Stars =", len_all)
    print("Res Stars = Stars - varStars =", len(table_res))
    log_file.write('Stars left -  %i\n' % len(table_res))

    if ploting:
        plt.imshow(image_tmp, cmap='Greys', origin='lower')
    A_list = []
    y_ar, x_ar = [], []
    star_count = 0
    if not c_flag:
        log_file.write("   NOMAD1         Vmag       Rmag         Flux         A        Mz         X           Y\n")
    # "1790-0005788      6.353     5.470    2482736.51558   22.51300  1.31407  894.32825   121.83167"
    for row in table_res:
        if (row["Rmag"] is not None) and (row["Vmag"] is not None):
            ra_s = row["RAJ2000"]
            dec_s = row["DEJ2000"]

            xs, ys = w.wcs_world2pix(ra_s, dec_s, 1)
            # print (xs, ys)
            if (xs > 0) and (ys > 0) and (xs < xc * 2) and (ys < yc * 2):
                star_count = star_count + 1
                # r_ap = 4
                # an_in = 12
                # an_out = 18

                if ploting:
                    fc = Circle((xs, ys), r_ap, facecolor='none', edgecolor='blue', linewidth=1, fill=False)
                    ax.add_patch(fc)

                # fit on xs, ys and measure flux
                try:
                    # targ_star = fit_m(image_tmp, int(xs), int(ys), gate=10, debug=True, fig_name=str(row["Rmag"]) + "_t.png", centring=True)
                    targ_star = fit_m(image_tmp, int(xs), int(ys), gate=5, debug=False, centring=False, silent=True)
                    # Centring is failing !!!!!

                    # targ_star = fit_m(image_tmp, int(xs), int(ys), gate=5, debug=True, fig_name="fig//" + str(row["Rmag"]) + "_t.png", centring=True, silent=True)
                except Exception:
                    print (row["NOMAD1"], "Fail fit Gauss...")
                    log_file.write('%s fail in Gaus fit\n' % row["NOMAD1"])
                    pass

                # Measure aperture flux
                positions = targ_star[:2]
                aperture = CircularAperture(positions, r=r_ap)
                annulus_aperture = CircularAnnulus(positions, r_in=an_in, r_out=an_out)

                apers = [aperture, annulus_aperture]
                phot_table = aperture_photometry(image_tmp, apers)
                for col in phot_table.colnames:
                    phot_table[col].info.format = '%.4g'  # for consistent table output

                bkg_mean = phot_table['aperture_sum_1'] / annulus_aperture.area
                bkg_sum = bkg_mean * aperture.area
                final_sum = phot_table['aperture_sum_0'] - bkg_sum
                phot_table['residual_aperture_sum'] = final_sum
                phot_table['residual_aperture_sum'].info.format = '%.4g'  # for consistent table output

                z = 0
                if len(phot_table) > 1:
                    if math.isnan(phot_table['residual_aperture_sum'][z]):
                        z = 1
                    for i in range(0, len(phot_table)):
                        if not math.isnan(phot_table['residual_aperture_sum'][i]):
                            if phot_table['residual_aperture_sum'][i] > phot_table['residual_aperture_sum'][z]:
                                z = i

                if len(targ_star) == 4:
                    xerr, yerr = targ_star[2], targ_star[3]
                else:
                    xerr, yerr = 9, 9

                if (xerr is np.inf) or (yerr is np.inf):
                    xerr, yerr = 8, 8

                flux = phot_table['residual_aperture_sum'][z]

                # kr = 0.8
                # Cr = 0.005
                vmr = row["Vmag"] - row["Rmag"]

                star = ephem.FixedBody()
                star._ra = ephem.degrees(str(ra_s))
                star._dec = ephem.degrees(str(dec_s))
                station.date = datetime.strptime(date_time[:-1], "%Y-%m-%dT%H:%M:%S.%f")
                star.compute(station)
                el = star.alt  # in radians !!!!!!!!
                Mz = 1 / (math.cos(math.pi / 2 - el))

                if (flux > 0) and (vmr is not masked):
                    if c_flag:
                        m_inst = -2.5 * math.log10(flux)
                        yq = row["Rmag"] - m_inst - kr * Mz
                        y_ar.append(yq)
                        x_ar.append(vmr)
                    else:
                        m_inst = -2.5 * math.log10(flux)
                        A = row["Rmag"] - m_inst - kr * Mz - Cr * vmr  # <------------------ A
                        A_list.append(A)

                        # print ("%8.5f  %8.5f  %10.5f  %8.5f  %8.5f" % (row["Vmag"], math.degrees(el), Mz, ra_s, dec_s))
                        fs = str.format("{0:" ">10.5f}", flux)
                        xx, yy = positions
                        xxs = str.format("{0:" ">8.5f}", xx)
                        yys = str.format("{0:" ">8.5f}", yy)
                        Mzs = str.format("{0:" ">3.5f}", Mz)
                        print ("%s   %8.3f  %8.3f  %15s   %8.5f %8s %10s  %10s" % (row["NOMAD1"], row["Vmag"], row["Rmag"], fs, A, Mzs, xxs, yys))
                        log_file.write("%s   %8.3f  %8.3f  %15s   %8.5f %8s %10s  %10s\n" % (row["NOMAD1"], row["Vmag"], row["Rmag"], fs, A, Mzs, xxs, yys))

                        if ploting:
                            circle = Circle((xx, yy), r_ap, facecolor='none', edgecolor='green', linewidth=1, fill=False)
                            r_in = Circle((xx, yy), an_in, facecolor='none', edgecolor='red', linewidth=1, fill=False)
                            r_out = Circle((xx, yy), an_out, facecolor='none', edgecolor='red', linewidth=1, fill=False)
                            ax.add_patch(circle)
                            ax.add_patch(r_in)
                            ax.add_patch(r_out)
                        # sys.exit()
            # else:
            #     print("not in frame...")

    if c_flag:
        y_ar = np.array(y_ar)
        x_ar = np.array(x_ar)
        c, a, r_max, ind = lsqFit(y_ar, x_ar)

        print ("A and c= ", a, c)

    if not c_flag:
        # log_file.write("Stars in frame = %i\n" % star_count)
        A_list = np.array(A_list)
        mA = np.mean(A_list)
        eA = np.std(A_list)
        print (mA, "+-", eA, "Stars=", len(A_list))
        print ("Filtering, with rms < %3.3f" % rms_val)

        log_file.write("A = %8.5f +/- %8.5f. ###### Stars in frame = %i. A calculated = %i\n" % (mA, eA, star_count, len(A_list)))
        log_file.write("Start A filtering...(err < %3.3f)\n" % rms_val)
        # log_file.write("A count = %i\n" % len(A_list))

        A_list2 = RMS_del(A_list, rms_val)
        mA = np.mean(A_list2)
        eA = np.std(A_list2)
        print (mA, "+-", eA, "Stars=", len(A_list2))
        # log_file.write("A count = %i\n" % len(A_list2))
        log_file.write("A = %8.5f +/- %8.5f. ###### A count = %i\n\n" % (mA, eA, len(A_list2)))

    if ploting:
        plt.show()
        plt.close('all')

    if c_flag:
        A_general.append(a)
        c_general.append(c)
    else:
        A_general.append(mA)


A_general = np.array(A_general)
A_general = RMS_del(A_general, rms_val)

Ag_mean = np.mean(A_general, axis=0)
Ag_err = np.std(A_general, axis=0)

if c_flag:
    # Ag_err = np.std(A_general, axis=0)

    c_general = np.array(c_general)
    c_mean = np.mean(c_general, axis=0)
    c_err = np.std(c_general, axis=0)
    # Ag_err = np.std(A_general, axis=0)[0]


log_file.write("\n\n")
log_file.write("####################################################################\n")
log_file.write("###-------------------A mean for all frames----------------------###\n")
log_file.write("A = %8.5f , sigma =%8.5f\n" % (Ag_mean, Ag_err))
log_file.write("###--------------------------------------------------------------###\n")

if c_flag:
    log_file.write("###----------------C mean for all frames-------------------------###\n")
    log_file.write("C = %8.5f , sigma =%8.5f\n" % (c_mean, c_err))
    log_file.write("###--------------------------------------------------------------###\n")

log_file.write("####################################################################\n")
log_file.close()
