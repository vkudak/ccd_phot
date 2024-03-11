import numpy as np
# import astropy.io.fits as fits
import os
import sys
# from astropy.stats import mad_std
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use('Agg')
from matplotlib.patches import Circle
import ephem
from datetime import datetime, timedelta
import math
from astropy import units as u
from astropy.io import fits
from astroquery.vizier import Vizier
import astropy.coordinates as coord
from tqdm import tqdm
from sklearn.linear_model import LinearRegression
import configparser


def substract(image, dark=None, value=None):
    image = np.array(image, dtype=float)

    if dark is not None:
        dark = np.array(dark, dtype=float)
        # print ("substract DARK")
        ph_image = image - dark

    if value is not None:
        # print ("substract VALUE")
        ph_image = image - value

    ph_image[ph_image < 0] = 0
    return ph_image


def linear(x, a, b):
    '''
    Linear function.
    This is just :math:`\\text{linear}(x; a, b) = a x + b`.
    '''
    return a * x + b


def fit_lin_reg(x, y, yerr):
    popt, pconv = curve_fit(linear, x, y, sigma=yerr)
    a, b = popt
    ae, be = np.sqrt(pconv.diagonal())

    residuals = y - linear(x, *popt)
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1 - (ss_res / ss_tot)

    res = []
    for i in range(0, len(x)):
        res.append(abs(y[i] - (a * x[i] + b)))
    res = np.array(res)
    res_max = np.max(res)
    res_min = np.min(res)
    if res_max > abs(res_min):
        ind = np.argmax(res)
    else:
        ind = np.argmin(res)
        res_max = abs(res_min)
    return [a, ae], [b, be], r2, res_max, ind


def linReg(x, y):
    x = np.array(x).reshape((-1, 1))
    y = np.array(y)
    model = LinearRegression()
    model.fit(x, y)
    r_sq = model.score(x, y)
    # print('coefficient of determination:', r_sq)
    # print('intercept:', model.intercept_)
    # print('slope:', model.coef_[0])
    return model.intercept_, model.coef_[0], r_sq


def lsqFit(y, x):
    '''
    y=ax+c
    Return a, c, residual
    '''
    x = np.array(x)
    y = np.array(y)
    A = np.vstack([x, np.ones(len(x))]).T
    res = []
    # linearly generated sequence
    if y != []:
        wb = np.linalg.lstsq(A, y, rcond=None)  # obtaining the parameters
        a, c = wb[0]
        residual_sum = wb[1]
        r2 = 1 - residual_sum / (y.size * y.var())  # = R^2
        # residual = wb[1][0]
    # else:
    #    residual = 999
    for i in range(0, len(x)):
        res.append(abs(y[i] - (a * x[i] + c)))
    res = np.array(res)
    res_max = np.max(res)
    res_min = np.min(res)
    if res_max > abs(res_min):
        ind = np.argmax(res)
    else:
        ind = np.argmin(res)
        res_max = abs(res_min)
    return a, c, res_max, ind, r2


def get_from_NOMAD(RA, DEC, w="0d60m", h="30m", radius=None, Filter={'Rmag': '<13'}):
    center = coord.SkyCoord(ra=RA, dec=DEC, unit=(u.deg, u.deg), frame='icrs')
    Vizier.ROW_LIMIT = 5000
    if radius:
        table = Vizier.query_region(center, width=w, radius=radius, catalog=["NOMAD"], column_filters=Filter)
    else:
        table = Vizier.query_region(center, width=w, height=h, catalog=["NOMAD"], column_filters=Filter)
    return table


def get_from_LAND(RA, DEC, w="0d60m", h="30m", radius=None, Filter={'Rmag': '<13'}):
    center = coord.SkyCoord(ra=RA, dec=DEC, unit=(u.deg, u.deg), frame='icrs')
    Vizier.ROW_LIMIT = 5000
    if radius:
        table = Vizier.query_region(center, width=w, radius=radius, catalog=["Landolt"], column_filters=Filter)
    else:
        table = Vizier.query_region(center, width=w, height=h, catalog=["Landolt"], column_filters=Filter)
    return table


def del_var(table, Filter={'Vmag': '<10'}):
    table_final = table
    Vizier.ROW_LIMIT = 5000
    i = 0
    for row in tqdm(table_final):
        ra = row["RAJ2000"]
        dec = row["DEJ2000"]
        c = coord.SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame='icrs')
        var_table = Vizier.query_region(c, radius='0d0m3s', catalog=["GCVS"],
                                        column_filters=Filter)  # , catalog=["GCVS"]
        var_table_len = len(var_table)
        if var_table_len > 0:
            table_final.remove_row(i)
            # print (table_final.remove_row(i))
            # print ("delete", row)
        i = i + 1
    return table_final


def gaussian(xycoor, x0, y0, sigma, amp):
    '''This Function is the Gaussian Function'''

    x, y = xycoor  # x and y taken from fit function.  Stars at 0, increases by 1, goes to length of axis
    A = 1 / (2 * sigma ** 2)
    eq = amp * np.exp(-A * ((x - x0) ** 2 + (y - y0) ** 2))  # Gaussian
    return eq


def fit(image, mf=None):
    # med = np.median(image)
    # image = image-med
    # image = image[0,0,:,:]

    # max_index = np.where(image >= np.max(image))
    # x0 = max_index[1]  # Middle of X axis
    # y0 = max_index[0]  # Middle of Y axis
    # print(image.shape)

    x0 = int(image.shape[1] / 2)
    y0 = int(image.shape[0] / 2)
    # print (x0, y0)

    x = np.arange(0, image.shape[1], 1)  # Stars at 0, increases by 1, goes to length of axis
    y = np.arange(0, image.shape[0], 1)  # Stars at 0, increases by 1, goes to length of axis
    xx, yy = np.meshgrid(x, y)  # creates a grid to plot the function over

    sigma = np.std(image)  # The standard dev given in the Gaussian
    amp = np.max(image)  # amplitude
    guess = [x0, y0, sigma, amp]  # The initial guess for the gaussian fitting

    low = [0, 0, 0, 0]  # start of data array
    # Upper Bounds x0: length of x axis, y0: length of y axis, st dev: max value in image, amplitude: 2x the max value
    upper = [image.shape[0], image.shape[1], np.max(image), np.max(image) * 2]
    # upper = [image.shape[0], image.shape[1], sigma * 100.0, np.max(image) * 0.5]
    bounds = [low, upper]

    params, pcov = curve_fit(gaussian, (xx.ravel(), yy.ravel()), image.ravel(), p0=guess, bounds=bounds,
                             maxfev=mf)  # optimal fit.  Not sure what pcov is.

    # par, cov, infodict, mesg, ier = optimize.leastsq(residuals, a_guess, args=(x,y),full_output=True)
    # params, pcov = curve_fit(gaussian, (xx.ravel(), yy.ravel()), image.ravel(), p0=guess, bounds=bounds)  #optimal fit.  Not sure what pcov is.
    # print("here")
    return params, pcov


def plotting(image, params, save=False, filename=None, par=None, err=None, tar=None, gate=None):
    try:
        fig, ax = plt.subplots()
        ax.imshow(image, origin='lower')
        ax.scatter(params[0], params[1], s=20, c='red', marker='x')
        circle = Circle((params[0], params[1]), params[2], facecolor='none', edgecolor='red', linewidth=1)

        ax.add_patch(circle)

        if save:
            filename, file_extension = os.path.splitext(filename)
            filename = filename + "_g%i" % gate + file_extension
            # plt.title("parameters= %s \nerrors=%s \ngate=%s" % (par[:2], par[2:], gate), pad=-50, fontsize=10)

            x, y = tar[:2]
            plt.title("parameters= %s \nerrors=%s \ngate=%s  x,y=%2.3f, %2.3f" % (par, err, gate, x, y), pad=-50,
                      fontsize=10)
            plt.savefig(filename)
        else:
            plt.show()
        plt.close('all')
    except Exception as E:
        print(E)


def fit_m(image, x0, y0, gate, debug, fig_name=None, centring=False, silent=False):
    # gate_x = gate_y
    data_fit = image[y0 - gate:y0 + gate, x0 - gate:x0 + gate]
    # plt.imshow(data_fit, cmap='Greys', origin='lower')
    # plt.show()
    if centring:
        if not silent:
            print("centring...", end=" ")
        # ------------------------------------------------------------------ find brighter pixel and center it!!!!
        bx0, by0 = np.unravel_index(data_fit.argmax(), data_fit.shape)
        if not silent:
            # print("bx, by=", bx0, by0)
            print(f"bx, by = {bx0:d},{by0:d}")

        sx = by0 - gate
        sy = bx0 - gate

        data_fit = image[y0 - gate + sy:y0 + gate + sy, x0 - gate + sx:x0 + gate + sx]
        # ----------------------------------------------------------------------------------------
    # print(data_fit.shape)
    # print(data_fit)

    # row_sums = data_fit.sum(axis=1)
    # data_fit = data_fit / row_sums[:, np.newaxis]

    # from sklearn.preprocessing import normalize
    # data_fit = normalize(data_fit, axis=1, norm='l2')
    # print(data_fit)

    # print(np.isfinite(data_fit))
    # import sys
    # sys.exit()

    # data_fit = data_fit[np.isfinite(data_fit)]
    # par = [5,5,0,0]
    # plotting(data_fit, par, save=True, filename=fig_name, par=par, err=None, tar=[5,5], gate=gate)

    # data_fit = data_fit[~np.isnan(data_fit)]
    # data_fit = data_fit *2.5
    # print(data_fit)
    par, pcov = fit(data_fit, mf=20000)

    err = np.sqrt(np.diag(pcov))
    amp = par[-1]
    target = [x0 - gate + par[0], y0 - gate + par[1], err[0], err[1], amp]
    if debug:
        # plotting(data_fit, par, save=True, filename=fig_name, par=target, gate=gate)
        plotting(data_fit, par, save=True, filename=fig_name, par=par, err=err, tar=target, gate=gate)
    return target


def get_tle(tle_path: str) -> list:
    """
    Read file with TLE

    Parameters
    ----------
    :param str tle_path : path to tle file
    """
    tle_list = []
    tle_file = open(tle_path, 'r')
    i = 1
    tle1 = []
    for line in tle_file:
        if line != '':
            tle1.append(line.strip())
            i = i + 1
            if i > 3:
                tle_list.append(tle1)
                tle1 = []
                i = 1
    tle_file.close()
    return tle_list


def calc_from_tle(lat, lon, elev, TLE_list, date_time, COSPAR, NORAD, NAME):
    if COSPAR == '' and NORAD == '' and NAME == '':
        return None
    else:
        try:
            tle = []
            for i in range(0, len(TLE_list)):
                l2 = TLE_list[i][1].split()
                cosp = l2[2]
                nor = l2[1]
                name = TLE_list[i][0]
                if cosp == COSPAR:
                    c = l2[2]
                    no = l2[1][:-1]
                    n = TLE_list[i][0]
                    tle = TLE_list[i]
            if tle == []:
                for i in range(0, len(TLE_list)):
                    l2 = TLE_list[i][1].split()
                    nor = l2[1][:-1]
                    # print 'nor=', nor
                    if nor == NORAD:
                        c = l2[2]
                        no = l2[1][:-1]
                        n = TLE_list[i][0]
                        tle = TLE_list[i]
            if tle == []:
                for i in range(0, len(TLE_list)):
                    l2 = TLE_list[i][1].split()
                    name = TLE_list[i][0]
                    if name == NAME:
                        c = l2[2]
                        no = l2[1][:-1]
                        n = TLE_list[i][0]
                        tle = TLE_list[i]
            if len(tle) > 0:
                # Calculating  El, Rg
                station = ephem.Observer()
                station.lat = lat
                station.long = lon
                station.elevation = elev
                # station.lat = '48.5635505'
                # station.long = '22.453751'
                # station.elevation = 231.1325
                sat = ephem.readtle(tle[0], tle[1], tle[2])
                try:
                    station.date = datetime.strptime(date_time[:-1], "%Y-%m-%dT%H:%M:%S.%f")
                except Exception:
                    print("Error - Wrong date time format...")

                sat.compute(station)
                el = math.degrees(sat.alt)
                rg = sat.range / 1000.0  # km
                az = math.degrees(sat.az)
                n = n.lstrip("0 ")
                return el, rg, az, n, no, c, tle
            else:
                print('\nCan not find TLE for such satellite !!!')
                sys.exit()
        except Exception as ex:
            print("Error occurred. Probably wrong TLE file. \n Error = " + ex.message)


def word_in_hfield(word, h_field):
    """
    Check if word is in header field (h_field)
    Parameters:
        word: search word
        h_field: astropy.io.fits.header._****  Field. (can be list as COMMENT or HISTORY)
    Returns: False or True
    """
    try:
        for c in h_field:
            if word in c:
                return True
            else:
                return False
    except TypeError:
        # return False if h_field does not  exist ('NoneType')
        return False


def read_config_stars(conf_file, log_file):
    """
    Read stars config file
    Parameters
    conf_file: path to config file
    log_file: opened for write text file instance
    Returns
    -------
    res: dictionary with parameters
    """
    config = configparser.ConfigParser(inline_comment_prefixes="#")
    res = {}
    if os.path.isfile(conf_file):
        try:
            config.read(conf_file)
            res['kr'] = config.getfloat('Stars_Stand', 'K')
            res['max_m'] = config.get('Stars_Stand', 'max_m_fit', fallback="14")
            res['rms_val'] = config.getfloat('Stars_Stand', 'A_rms', fallback=0.05)
            res['r_max_val'] = config.getfloat('Stars_Stand', 'r_max_val', fallback=6.25)
            res['c_flag'] = config.getboolean('Stars_Stand', 'calc_C', fallback=True)
            res['snr_value'] = config.getfloat('Stars_Stand', 'snr', fallback=1.2)
            if not res['c_flag']:
                res['Cr'] = config.getfloat('Stars_Stand', 'C')

            res['dark_frame'] = config.get('Stars_Stand', 'dark_frame', fallback=False)
            res['dark_stable'] = config.getfloat('Stars_Stand', 'dark_stable', fallback=0.0)

            res['r_ap'] = config.getfloat('APERTURE', 'r_ap')
            res['an_in'] = config.getfloat('APERTURE', 'an_in')
            res['an_out'] = config.getfloat('APERTURE', 'an_out')

            res['scale_min'] = config.getfloat('astrometry.net', 'scale_lower', fallback=1)
            res['scale_max'] = config.getfloat('astrometry.net', 'scale_upper', fallback=20)
            res['api_key'] = config.get('astrometry.net', 'api_key', fallback="No key")

            if res['scale_max'] == 20 and res['scale_min'] == 1:
                print("No 'astrometry.net' section in INI file. Using default astrometry.net params")
                log_file.write("No 'astrometry.net' section in INI file. Using default astrometry.net params\n")

            res['site_name'] = config.get('SITE', 'Name', fallback="No name")
            res['site_lat'] = config.get('SITE', "lat")
            res['site_lon'] = config.get('SITE', 'lon')
            res['site_elev'] = config.getfloat('SITE', 'h')

            return res

        except Exception as E:
            print("Error in INI file\n", E)
            return False
    else:
        print(f"Error. Cant find config_stars.ini in {conf_file}")
        log_file.write(f"Error. Cant find config_stars.ini in {conf_file} \n")
        return False


def read_config_sat(conf_file):
    """
    Args:
        conf_file (str): Full path to config file

    Returns:
        res (dict, bool): dict of parameters OR False
    """

    config = configparser.ConfigParser(inline_comment_prefixes="#")
    res = {}

    if os.path.isfile(conf_file):
        try:
            config.read(conf_file)

            res['cospar'] = config['NAME']['COSPAR']
            res['norad'] = config['NAME']['NORAD']
            res['name'] = config['NAME']['NAME']

            res['tle_file'] = config.get('TLE', 'tle_file')

            res['A'] = float(config['STD']['A'])
            res['k'] = float(config['STD']['k'])
            res['gate'] = int(config['STD']['gate'])

            band = config.get('STD', 'Filter', fallback="F")
            res['Filter'] = band.strip()

            time_format = config.get('STD', 'Time_format', fallback="UT")
            res['time_format'] = time_format.strip()

            res['min_real_mag'] = config.getfloat('STD', 'min_real_mag', fallback=15.0)
            res['max_center_error'] = config.getfloat('STD', 'max_center_err', fallback=2.0)
            res['fits_sort'] = config.get('STD', 'fits_sort', fallback='None')

            res['dark_frame'] = config.get('STD', 'dark_frame', fallback=False)

            res['r_ap'] = float(config['APERTURE']['ap'])
            res['an_in'] = float(config['APERTURE']['an_in'])
            res['an_out'] = float(config['APERTURE']['an_out'])

            res['site_name'] = config['SITE']['Name']
            res['site_lat'] = config['SITE']['lat']
            res['site_lon'] = config['SITE']['lon']
            res['site_elev'] = float(config['SITE']['h'])

            return res

        except Exception as E:
            print("Error in INI file\n", E)
            return False
    else:
        print("Error. Cant find config_sat.ini in " + conf_file)
        return False


def get_file_list(path):
    """
    Args:
        path (str): directory where to search FITS files
    Returns:
        fl (list): sorted list of files (with paths)
    """
    fl = []
    for file_name in os.listdir(path):
        f = file_name.split('.')
        if f[-1] in ['FIT', 'FITS', 'fit', 'fits']:
            fl.append(file_name)
    fl.sort()
    return fl


def sort_file_list(path, fl, field='DATE-OBS'):
    """
    Sort FITS file list according to field in fits header (default field: 'DATE-OBS')

    Parameters
    ----------
    path: path to files
    fl: file list to sort
    field: parameter to sort by (from header)

    Returns
    -------
    fl: sorted file list
    """
    sort_list = []
    for f in fl:
        header = fits.getheader(os.path.join(path, f))
        tt = header.get(field)
        if tt is not None:
            sort_list.append([f, tt])
        else:
            sort_list.append([f, 0])
    sort_list.sort(key=lambda x: x[1])  # sort by second element in [name, time]
    sort_list = [f[0] for f in sort_list]
    return sort_list


def calc_mag(flux, el, rg, zp, k, exp, min_mag=15):
    """
    Parameters
    ----------
    :param float flux: Flux
    :param float el: elevation in degrees
    :param float rg: Range to satellite
    :param float zp: Zero point of System
    :param float k: coefficient of extinction
    :param float exp: Exposition
    :param min_mag: minimum reachable magnitude

    Returns
    -------
    m (float): standard magnitude

    """
    if flux < 0:
        m = min_mag
        return m
    else:
        m_inst = -2.5 * math.log10(flux / exp)
        mz = 1 / (math.cos(math.pi / 2 - math.radians(el)))  # 1/cos(90-el)
        mr = -5 * math.log10(rg / 1000.0)
        mzr = k * mz

        m = zp + m_inst + mzr + mr
    return m


def fix_datetime(date_time):
    """
    Fixes datetime with wrong seconds (sec = 60)
    TIME format must be '2024-01-10T18:07:60.00000'

    Args:
        :param str date_time: DATE-OBS from header
    Return:
        Corrected datetime as str
    """
    if date_time[17:19] == '60':  # wrong time in sec section
        date_time = date_time[:17] + "00" + date_time[19:]
        date_time = datetime.strptime(date_time[:-1], "%Y-%m-%dT%H:%M:%S.%f")
        date_time = date_time + timedelta(minutes=1)
        date_time = date_time.strftime("%Y-%m-%dT%H:%M:%S.%f")
    return date_time
