#!python3
from astropy.stats import sigma_clipped_stats
from astropy.io import fits
import numpy as np
from numpy.ma import masked
from photutils.aperture import CircularAperture, CircularAnnulus
from sp_utils import *
import math
import sys
import os

from astroquery.astrometry_net import AstrometryNet
from astropy.wcs import WCS
import ephem
import warnings
from photometry_with_errors import *
import pickle


if len(sys.argv) < 2:
    print("Not enough parameters. Enter path")
    sys.exit()


path = sys.argv[1]
warnings.filterwarnings("ignore")

ploting = False  # plot each frame with apertures


# star_example = {
#     "NOMAD": "1141-0043729",
#     "Vmag": 5, "Rmag": 5.3,
#     "V-R": 0.3,
#     "flux": [1,2,3,4], "flux_err": [1,2,3,4], "flux/bkg": [1,2,3,4],
#     "Mz": [1,2,3,4],
#     "X": [1,2,3,4], "Y": [1,2,3,4]
# }
# database = ["NOMAD":star_example,]
database = []

log_file = open(path + '//star_process.log', "w")

if ploting:
    # for plot
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    # from astropy.visualization import LogStretch
    # from astropy.visualization.mpl_normalize import ImageNormalize

# config = configparser.ConfigParser(inline_comment_prefixes="#")
# config.read(os.path.join(path,'config_stars.ini'))
# if os.path.isfile(os.path.join(path, 'config_stars.ini')):
#     try:
#         kr = config.getfloat('Stars_Stand', 'K')
#         max_m = config.get('Stars_Stand', 'max_m_fit', fallback="14")
#         rms_val = config.getfloat('Stars_Stand', 'A_rms', fallback=0.05)
#         c_flag = config.getboolean('Stars_Stand', 'calc_C', fallback=True)
#         snr_value = config.getfloat('Stars_Stand', 'snr', fallback=1.2)
#         if not c_flag:
#             Cr = config.getfloat('Stars_Stand', 'C')
#
#         dark_frame = config.get('Stars_Stand', 'dark_frame', fallback=False)
#         dark_stable = config.getfloat('Stars_Stand', 'dark_stable', fallback=0.0)
#
#         r_ap = config.getfloat('APERTURE', 'r_ap')
#         an_in = config.getfloat('APERTURE', 'an_in')
#         an_out = config.getfloat('APERTURE', 'an_out')
#
#         scale_min = config.getfloat('astrometry.net', 'scale_lower', fallback=1)
#         scale_max = config.getfloat('astrometry.net', 'scale_upper', fallback=20)
#         api_key = config.get('astrometry.net', 'api_key', fallback="No key")
#
#         if scale_max == 20 and scale_min == 1:
#             print("No 'astrometry.net' section in INI file. Using default astrometry.net params")
#             log_file.write("No 'astrometry.net' section in INI file. Using default astrometry.net params\n")
#
#         site_name = config.get('SITE', 'Name', fallback="No name")
#         site_lat = config.get('SITE', "lat")
#         site_lon = config.get('SITE', 'lon')
#         site_elev = config.getfloat('SITE', 'h')
#
#     except Exception as E:
#         print("Error in INI file\n", E)
#         sys.exit()
# else:
#     print(f"Error. Cant find config_stars.ini in {os.path.join(path, 'config_stars.ini')}")
#     log_file.write(f"Error. Cant find config_stars.ini in {os.path.join(path, 'config_stars.ini')} \n")
#     sys.exit()

conf = read_config_stars(os.path.join(path, 'config_stars.ini'), log_file)
if not conf:
    sys.exit()

station = ephem.Observer()
station.lat = conf['site_lat']  # '48.5635505'
station.long = conf['site_lon']  # '22.453751'
station.elevation = conf['site_elev']  # 231.1325

fl = get_file_list(path)

if len(fl) == 0:
    print("No FIT files to process. EXIT")
    log_file.write("No FIT files to process. EXIT\n")
    sys.exit()

# print (fl)
# fl = fl[0:15] + fl[55:70]


# #################### BEGIN
ast = AstrometryNet()
# ast.show_allowed_settings()

ast.api_key = conf['api_key']

if ploting:
    fig, ax = plt.subplots()

A_general = []
c_general = []

for fit_file in fl:
    print(fit_file)
    file_with_path = os.path.join(path, fit_file)
    log_file.write("####################################################\n")
    log_file.write("filename = " + file_with_path + "\n")
    try_again = True
    submission_id = None
    hdu = fits.open(file_with_path)  # path + "//" + fit_file)
    header = hdu[0].header
    date_time = hdu[0].header.get('DATE-OBS')
    exp = header.get('EXPTIME')
    hdu.close()
    # author = None
    author = header.get('AUTHOR')

    astrometry_net = False
    # if (header.get('COMMENT') is not None) or (header.get('HISTORY') is not None):
    for as_value in ["Astrometry.net", "astrometry.net"]:
        if word_in_hfield(as_value, header.get('COMMENT')) or word_in_hfield(as_value, header.get('HISTORY')):
            astrometry_net = True

    if not astrometry_net:
        while try_again:
            try:
                if not submission_id:
                    print(file_with_path)
                    wcs_header = ast.solve_from_image(file_with_path,
                                                      submission_id=submission_id,
                                                      force_image_upload=True,
                                                      crpix_center=True,
                                                      downsample_factor=2,
                                                      scale_units="arcsecperpix",
                                                      scale_type='ul',
                                                      scale_upper=conf['scale_max'],
                                                      scale_lower=conf['scale_min'],
                                                      tweak_order=3)
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
            # sys.exit()
            log_file.write("file SOLVED\n")
            with fits.open(file_with_path, mode='update') as hdul:

                hdul[0].header.update(wcs_header)
                hdul.close()

                # hdul[0].header = wcs_header
                # hdul[0].header.append(('AUTHOR', "LKD UZhNU", 'Solved sucessfuly with astrometry.net'))
                # hdul[0].header.append(('DATE-OBS', date_time, "System Clock:Est. Frame Start -OR- GPS:Start Exposure"))
                # hdul[0].header.append(('EXPTIME', exp, "EXPOSURE in seconds"))
                # hdul.flush()
                # hdul.close()

        else:
            # Code to execute when solve fails
            print("Fail solving image")
            log_file.write("FIT file NOT SOLVED\n")

    # BEGIN star find---------------------------------------------------------------
    # 1. Get RA DEC of image center
    # https://python4astronomers.github.io/astropy/wcs.html
    # 2. Get stars from NOMAD  (done!!!)
    # 3. transform stars RA, DEC -> X,Y. For each star
    # 4. calculate star flux at X,Y
    # 5. make result in flux, Rmag, Mz ??????
    log_file.write('Begin star find procedure...\n')

    header = fits.getheader(file_with_path)
    image_tmp = fits.getdata(file_with_path)

    if conf['dark_frame']:
        dark_arr = fits.getdata(conf['dark_frame'])
        image_tmp = image_tmp - dark_arr
        ph_image = image_tmp

        minI = np.min(image_tmp)
        if minI < 0:
            print("Warning! Image - Dark has negative pixels...")
            log_file.write("Warning! Image - Dark has negative pixels...\n")
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
    print("Grab stars from Vizier (Vmag<%s)...." % conf['max_m'])
    log_file.write("Grab stars from Vizier (Vmag<%s)....\n" % conf['max_m'])

    table_res = get_from_Gaia_EDR3_std(ra_c, dec_c, radius="2.0deg", Filter={'Vmag': '<' + conf['max_m']})

    # table_res = get_from_LAND(ra_c, dec_c, radius="2.0deg", Filter={'Vmag': '<' + conf['max_m']})
    # print(table_res)
    # for star in table_res:
    #     print(star)
    # sys.exit()
    # table_res = table_res["I/297/out"]  # get NOMAD data
    table_res = table_res[0]

    # table_res.remove_columns(["YM", 'r', 'pmRA', 'e_pmRA', 'pmDE', 'e_pmDE', 'Jmag', 'Hmag', 'Kmag', 'R', 'r_Bmag', 'r_Vmag', 'r_Rmag'])
    table_res.sort(["Vmag"])
    if "Rmag" not in table_res.colnames:
        print("creating column Rmag")
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

    star_count = 0
    if not conf['c_flag']:
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
            # print (xs, ys, row['SimbadName'])
            if (xs > 0) and (ys > 0) and (xs < xc * 2 - 20) and (ys < yc * 2 - 20):   # cut 20 px from edges
                star_count = star_count + 1

                if ploting:
                    fc = Circle((xs, ys), conf['r_ap'], facecolor='none', edgecolor='blue', linewidth=1, fill=False)
                    ax.add_patch(fc)

                # fit on xs, ys and measure flux
                try:
                    # targ_star = fit_m(image_tmp, int(xs), int(ys), gate=10, debug=True, fig_name=str(row["Rmag"]) + "_t.png", centring=True)
                    # targ_star = fit_m(image_tmp, int(xs), int(ys), gate=5, debug=False, centring=False, silent=True)
                    figname = "D:\\FTP\\fig.png"
                    targ_star = fit_m(image_tmp, int(xs), int(ys), gate=5, debug=False, fig_name=figname, centring=True, silent=True)

                    positions = targ_star[:2]
                    aperture = CircularAperture(positions, r=conf['r_ap'])
                    annulus_aperture = CircularAnnulus(positions, r_in=conf['an_in'], r_out=conf['an_out'])
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

                    vmr = row["V-R"]  #row["Vmag"] - row["Rmag"]

                    star = ephem.FixedBody()
                    star._ra = ephem.hours(str(ra_s))
                    star._dec = ephem.degrees(str(dec_s))
                    station.date = datetime.strptime(date_time[:-1], "%Y-%m-%dT%H:%M:%S.%f")
                    star.compute(station)
                    el = star.alt  # in radians !!!!!!!!
                    Mz = 1 / (math.cos(math.pi / 2 - el))

                    if (flux > 0) and (vmr is not masked) and (abs(row["V-R"]) < 2):
                        if conf['c_flag']:
                            fs = str.format("{0:" ">10.3f}", flux)
                            fes = str.format("{0:" ">8.3f}", flux_err)
                            fbs = str.format("{0:" ">10.3f}", fb)
                            snrs = str.format("{0:" ">10.3f}", snr)
                            xx, yy = positions
                            Mzs = str.format("{0:" ">3.3f}", Mz)
                            xxs = str.format("{0:" ">8.3f}", xx)
                            yys = str.format("{0:" ">8.3f}", yy)
                            if snr < conf['snr_value']:  # print "*" on bed star
                                log_file.write("%17s   %8.3f  %8.3f  %8.3f   %15s %10s %12s %5s*  %5s %8s %8s\n" %
                                               (row["SimbadName"], row["Vmag"], row["Rmag"], vmr, fs, fes, fbs, snrs, Mzs, xxs, yys))
                            else:
                                log_file.write("%17s   %8.3f  %8.3f  %8.3f   %15s %10s %12s %5s   %5s %8s %8s\n" %
                                               (row["SimbadName"], row["Vmag"], row["Rmag"], vmr, fs, fes, fbs, snrs, Mzs, xxs, yys))
                        else:
                            m_inst = -2.5 * math.log10(flux/exp)
                            A = row["Rmag"] - m_inst - (conf['kr'] * Mz) - conf['Cr'] * vmr  # <------------------ A

                            # print ("%8.5f  %8.5f  %10.5f  %8.5f  %8.5f" % (row["Vmag"], math.degrees(el), Mz, ra_s, dec_s))
                            fs = str.format("{0:" ">10.5f}", flux)
                            xx, yy = positions
                            xxs = str.format("{0:" ">8.5f}", xx)
                            yys = str.format("{0:" ">8.5f}", yy)
                            Mzs = str.format("{0:" ">3.5f}", Mz)
                            # print("%s   %8.3f  %8.3f  %15s   %8.5f %8s %10s  %10s" % (row["NOMAD1"], row["Vmag"], row["Rmag"], fs, A, Mzs, xxs, yys))
                            log_file.write("%17s   %8.3f  %8.3f  %15s   %8.5f %8s %10s  %10s\n" % (row["SimbadName"], row["Vmag"], row["Rmag"], fs, A, Mzs, xxs, yys))

                            if ploting:
                                circle = Circle((xx, yy), conf['r_ap'], facecolor='none', edgecolor='green', linewidth=1, fill=False)
                                r_in = Circle((xx, yy), conf['an_in'], facecolor='none', edgecolor='red', linewidth=1, fill=False)
                                r_out = Circle((xx, yy), conf['an_out'], facecolor='none', edgecolor='red', linewidth=1, fill=False)
                                ax.add_patch(circle)
                                ax.add_patch(r_in)
                                ax.add_patch(r_out)
                            # sys.exit()

                    if snr > conf['snr_value']:
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




