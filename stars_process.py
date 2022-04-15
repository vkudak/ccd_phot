#!python3.8
from astropy.stats import sigma_clipped_stats
from astropy.io import fits
import numpy as np
from numpy.ma import masked
from photutils import CircularAperture, CircularAnnulus
from sp_utils import *
import math
import sys
import os

from astroquery.astrometry_net import AstrometryNet
from astropy.wcs import WCS
import ephem
import configparser
import warnings
from photometry_with_errors import *
import pickle


if len(sys.argv) < 2:
    print("Not enough parameters. Enter path")
    sys.exit()


path = sys.argv[1]

warnings.filterwarnings("ignore")

ploting = False  # plot each frame with appertures


# star_example = {"NOMAD": "1141-0043729", "Vmag": 5, "Rmag": 5.3, "V-R": 0.3, "flux": [1,2,3,4], "flux_err": [1,2,3,4], "flux/bkg": [1,2,3,4], "Mz": [1,2,3,4], "X": [1,2,3,4], "Y": [1,2,3,4]}
# database = ["NOMAD":star_example,]
database = []
a_database = []

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
        kr = config['Stars_Stand']['K']
        kr = float(kr)
        max_m = config['Stars_Stand']['max_m_fit']
        rms_val = float(config['Stars_Stand']['A_rms'])
        c_flag = config['Stars_Stand'].getboolean('calc_C')
        snr_value = config['Stars_Stand'].getfloat('snr')
        if not c_flag:
            Cr = config['Stars_Stand']['C']
            Cr = float(Cr)

        try:
            dark_frame = config['Stars_Stand']['dark_frame']
        except Exception:
            dark_frame = False

        r_ap = float(config['APERTURE']['r_ap'])
        an_in = float(config['APERTURE']['an_in'])
        an_out = float(config['APERTURE']['an_out'])

    except Exception as E:
        print("Error in inin file\n", E)
        sys.exit()
else:
    print("Error. Cant find config_stars.ini in " + path + '//config_stars.ini')


# print (path)

list = os.listdir(path)
fl = []
for fn in list:
    f = fn.split('.')
    if f[-1] in ['FIT', 'FITS', 'fit', 'fits']:
        fl.append(fn)

fl.sort()

# print (fl)
# fl = fl[0:15] + fl[55:70]


# #################### BEGIN
ast = AstrometryNet()
# ast.show_allowed_settings()

ast.api_key = "ittzfaqwnrvduhax"


if ploting:
    fig, ax = plt.subplots()


log_file = open(path + '//star_process.log', "w")

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
                    print(path + "//" + fit_file)
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
            print("Fail")
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

    if dark_frame:
        dark_arr = fits.getdata(dark_frame)
        image_tmp = image_tmp - dark_arr
        ph_image = image_tmp

        minI = np.min(image_tmp)
        if minI < 0:
            print("Warning! Image - Dark has negativ pixels...")
            log_file.write("Warning! Image - Dark has negativ pixels...\n")
    else:
        ph_image = image_tmp
        mean, median, std = sigma_clipped_stats(image_tmp, sigma=3.0)
        image_tmp = image_tmp - mean

    xc = header["NAXIS1"] / 2.
    yc = header["NAXIS2"] / 2.
    exp = header["EXPTIME"]
    exp = float(exp)

    w = WCS(header)
    ra_c, dec_c = w.wcs_pix2world(xc, yc, 1)  # RA DEC of FRAME center
    print("Grab stars from Vizier (Vmag<%s)...." % max_m)
    log_file.write("Grab stars from Vizier (Vmag<%s)....\n" % max_m)

    table_res = get_from_LAND(ra_c, dec_c, radius="3.0deg", Filter={'Vmag': '<' + max_m})
    # print(table_res)
    # sys.exit()
    # table_res = table_res["I/297/out"]  # get NOMAD data
    table_res = table_res[0]

    # table_res.remove_columns(["YM", 'r', 'pmRA', 'e_pmRA', 'pmDE', 'e_pmDE', 'Jmag', 'Hmag', 'Kmag', 'R', 'r_Bmag', 'r_Vmag', 'r_Rmag'])
    table_res.sort(["Vmag"])
    table_res["Rmag"] = table_res["Vmag"] - table_res["V-R"]

    # print(table_res.colnames)
    # sys.exit()

    len_all = len(table_res)
    log_file.write('Find %i stars\n' % len_all)

    # Write code to delete variable stars from table_res
    # print("Deleting variable stars from list....")
    # log_file.write('Deleting variable stars...\n')
    # table_res = del_var(table_res, Filter={'Vmag': '<' + max_m})
    print("Stars =", len_all)
    # print("Res Stars = Stars - varStars =", len(table_res))
    # log_file.write('Stars left -  %i\n' % len(table_res))

    if ploting:
        plt.imshow(image_tmp, cmap='Greys', origin='lower')
    A_list = []
    y_ar, x_ar = [], []
    star_count = 0
    if not c_flag:
        log_file.write("   SimbadName         Vmag       Rmag         Flux         A        Mz         X           Y\n")
    else:
        log_file.write("   SimbadName         Vmag       Rmag      V-R            Flux       Flux_err     bkg          snr      Mz       X         Y\n")
    # "1790-0005788      6.353     5.470    2482736.51558   22.51300  1.31407  894.32825   121.83167"
    for row in table_res:
        if (row["Rmag"] is not None) and (row["Vmag"] is not None):
            ra_s = row["RAJ2000"]
            dec_s = row["DEJ2000"]

            from astropy.coordinates import SkyCoord
            radec = SkyCoord(row["RAJ2000"] + " " + row["DEJ2000"], frame="icrs", unit=(u.hourangle, u.deg))
            # print(radec.ra.deg)
            # sys.exit()

            xs, ys = w.wcs_world2pix(radec.ra.deg, radec.dec.deg, 1)
            # print (xs, ys)
            if (xs > 0) and (ys > 0) and (xs < xc * 2 - 20) and (ys < yc * 2 - 20):   # cut 20 px from edges
                star_count = star_count + 1

                if ploting:
                    fc = Circle((xs, ys), r_ap, facecolor='none', edgecolor='blue', linewidth=1, fill=False)
                    ax.add_patch(fc)

                # fit on xs, ys and measure flux
                try:
                    # targ_star = fit_m(image_tmp, int(xs), int(ys), gate=10, debug=True, fig_name=str(row["Rmag"]) + "_t.png", centring=True)
                    # targ_star = fit_m(image_tmp, int(xs), int(ys), gate=5, debug=False, centring=False, silent=True)
                    figname = "D:\\FTP\\fig.png"
                    targ_star = fit_m(image_tmp, int(xs), int(ys), gate=5, debug=False, fig_name=figname, centring=True, silent=True)

                    positions = targ_star[:2]
                    aperture = CircularAperture(positions, r=r_ap)
                    annulus_aperture = CircularAnnulus(positions, r_in=an_in, r_out=an_out)
                    phot_table = iraf_style_photometry(aperture, annulus_aperture, ph_image)

                    # print(phot_table)
                    # sys.exit()

                    z = 0

                    if len(targ_star) == 4:
                        xerr, yerr = targ_star[2], targ_star[3]
                    else:
                        xerr, yerr = 9, 9

                    if (xerr is np.inf) or (yerr is np.inf):
                        xerr, yerr = 8, 8

                    # flux = phot_table['residual_aperture_sum'][z]
                    # bkg_flux = phot_table['residual_bkg_sum'][z]
                    # snr = flux / bkg_flux

                    flux = phot_table['flux'][z]
                    flux_err = phot_table['flux_error'][z]

                    fb = phot_table['flux_bkg'][z]
                    snr = phot_table['snr'][z]

                    vmr = row["Vmag"] - row["Rmag"]

                    star = ephem.FixedBody()
                    star._ra = ephem.degrees(str(ra_s))
                    star._dec = ephem.degrees(str(dec_s))
                    station.date = datetime.strptime(date_time[:-1], "%Y-%m-%dT%H:%M:%S.%f")
                    star.compute(station)
                    el = star.alt  # in radians !!!!!!!!
                    Mz = 1 / (math.cos(math.pi / 2 - el))

                    if (flux > 0) and (vmr is not masked) and (abs(row["Rmag"] - row["Vmag"]) < 2):
                        if c_flag:
                            m_inst = -2.5 * math.log10(flux/exp)
                            yq = row["Rmag"] - m_inst + kr * Mz
                            if (snr > snr_value):
                                y_ar.append(yq)
                                x_ar.append(vmr)

                            fs = str.format("{0:" ">10.3f}", flux)
                            fes = str.format("{0:" ">8.3f}", flux_err)
                            fbs = str.format("{0:" ">10.3f}", fb)
                            snrs = str.format("{0:" ">10.3f}", snr)
                            xx, yy = positions
                            Mzs = str.format("{0:" ">3.3f}", Mz)
                            xxs = str.format("{0:" ">8.3f}", xx)
                            yys = str.format("{0:" ">8.3f}", yy)
                            if snr < snr_value:  # print "*" on bed star
                                log_file.write("%17s   %8.3f  %8.3f  %8.3f   %15s %10s %12s %5s*  %5s %8s %8s\n" %
                                               (row["SimbadName"], row["Vmag"], row["Rmag"], vmr, fs, fes, fbs, snrs, Mzs, xxs, yys))
                            else:
                                log_file.write("%17s   %8.3f  %8.3f  %8.3f   %15s %10s %12s %5s   %5s %8s %8s\n" %
                                               (row["SimbadName"], row["Vmag"], row["Rmag"], vmr, fs, fes, fbs, snrs, Mzs, xxs, yys))
                        else:
                            m_inst = -2.5 * math.log10(flux/exp)
                            A = row["Rmag"] - m_inst + kr * Mz - Cr * vmr  # <------------------ A
                            A_list.append(A)

                            # print ("%8.5f  %8.5f  %10.5f  %8.5f  %8.5f" % (row["Vmag"], math.degrees(el), Mz, ra_s, dec_s))
                            fs = str.format("{0:" ">10.5f}", flux)
                            xx, yy = positions
                            xxs = str.format("{0:" ">8.5f}", xx)
                            yys = str.format("{0:" ">8.5f}", yy)
                            Mzs = str.format("{0:" ">3.5f}", Mz)
                            # print("%s   %8.3f  %8.3f  %15s   %8.5f %8s %10s  %10s" % (row["NOMAD1"], row["Vmag"], row["Rmag"], fs, A, Mzs, xxs, yys))
                            log_file.write("%17s   %8.3f  %8.3f  %15s   %8.5f %8s %10s  %10s\n" % (row["SimbadName"], row["Vmag"], row["Rmag"], fs, A, Mzs, xxs, yys))

                            if ploting:
                                circle = Circle((xx, yy), r_ap, facecolor='none', edgecolor='green', linewidth=1, fill=False)
                                r_in = Circle((xx, yy), an_in, facecolor='none', edgecolor='red', linewidth=1, fill=False)
                                r_out = Circle((xx, yy), an_out, facecolor='none', edgecolor='red', linewidth=1, fill=False)
                                ax.add_patch(circle)
                                ax.add_patch(r_in)
                                ax.add_patch(r_out)
                            # sys.exit()

                    if snr > snr_value:
                        # check if exist in DB
                        exist = False
                        save_ind = None
                        for ind in range(0, len(database)):
                            if database[ind]["SimbadName"] == (row["SimbadName"]):
                                exist = True
                                save_ind = ind

                        if exist:
                            a_flux = np.array(database[save_ind]["Flux"])
                            a_snr = np.array(database[save_ind]["f/b"])
                            a_xx = np.array(database[save_ind]["X"])
                            a_yy = np.array(database[save_ind]["Y"])
                            # a_A = np.array(database[save_ind]["A"])
                            database[save_ind]["Vmag"] = row["Vmag"]
                            database[save_ind]["Rmag"] = row["Rmag"]
                            database[save_ind]["V-R"] = vmr
                            database[save_ind]["Flux"] = np.append(a_flux, flux)
                            database[save_ind]["f/b"] = np.append(a_snr, snr)
                            database[save_ind]["Mz"] = Mz
                            database[save_ind]["X"] = np.append(a_xx, xx)
                            database[save_ind]["Y"] = np.append(a_yy, yy)
                            # if not c_flag:
                            #     database[save_ind]["A"] = np.append(a_A, A)
                        else:
                            star_e = {
                                "SimbadName": row["SimbadName"],
                                "Vmag": row["Vmag"],
                                "Rmag": row["Rmag"],
                                "V-R": vmr,
                                "Flux": np.array([flux]),
                                "f/b": np.array([snr]),
                                "Mz": Mz,
                                "X": np.array([xx]),
                                "Y": np.array([yy])
                                # "A": np.array([])
                            }
                            # if not c_flag:
                            #     star_e["A"] = np.array([A])
                            database.append(star_e)
                except Exception as e:
                    # except Exception as e:
                    # print(str(e))
                    print(row["SimbadName"], "Fail fit Gauss...")
                    log_file.write('%17s fail in Gauss fit\n' % row["SimbadName"])
                    pass

log_file.write("\n\n")
print("Stars = ", len(database))
log_file.write("Stars total = %i\n" % len(database))


db_filename = os.path.join(path, "res_stars.bin")
with open(db_filename, "wb") as fdb:
    pickle.dump([database, exp], fdb, protocol=pickle.HIGHEST_PROTOCOL)
    # load database
    # database = pickle.load(open("res_stars.bin", "rb"))

print("Database written to file")
log_file.write("Database written to file")




