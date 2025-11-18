import numpy as np
# import astropy.io.fits as fits
import os
import sys
# from astropy.stats import mad_std
from scipy.optimize import curve_fit
from scipy.stats import pearsonr
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
from sgp4.io import fix_checksum



def convert_ndarray(obj):
    """ Рекурсивно проходить по структурі та конвертує numpy-об'єкти у стандартні типи """
    if isinstance(obj, np.ndarray):
        return obj.tolist()  # Конвертуємо масив у список
    elif isinstance(obj, (np.float32, np.float64)):
        return float(obj)  # Конвертуємо float32/float64 у стандартний float
    elif isinstance(obj, (np.int32, np.int64)):
        return int(obj)  # Конвертуємо int32/int64 у стандартний int
    elif isinstance(obj, list):
        return [convert_ndarray(item) for item in obj]  # Обробляємо список рекурсивно
    elif isinstance(obj, dict):
        return {key: convert_ndarray(value) for key, value in obj.items()}  # Обробляємо словник рекурсивно
    else:
        return obj  # Інші типи залишаємо без змін


def convert_to_numpy(obj):
    """ Рекурсивно перетворює списки назад у numpy-масиви """
    if isinstance(obj, list):
        return np.array(obj) if all(isinstance(i, (int, float)) for i in obj) else [convert_to_numpy(i) for i in obj]
    elif isinstance(obj, dict):
        return {key: convert_to_numpy(value) for key, value in obj.items()}
    else:
        return obj  # Повертаємо без змін, якщо це не список або словник


def solve_photometric_coefficients(m_inst, m_std, x, color_ind, star_names, band,
                                   threshold=0.25,
                                   log_file=None,
                                   k=None,
                                   plot=False,
                                   path=None):
    """
    Calculates the coefficients Z, K, C in some filter,
    using the least squares method with point rejection.

    Parameters:
    m_inst (array): Instrumental Magnitudes.
    M_std (array): Standard Magnitudes (Landolt catalogue or other).
    X (array): Air mass of Star.
    Color_ind (array): Color index of Star (B-V or other).
    Star_names (array): Star names.
    Band(str): Band (Filter).
    Threshold (float): Threshold for rejection of stars (Default = 0.25).
    Log_file (file obj): Log file to write to OR None.
    K (float or None): If None, k is determined, otherwise it is constant.
    Plot (bool): If True, plot the Graph.
    Path (string): Path to save the figure.

    Returns:
    tuple: (Z, K, C, sigma_Z, sigma_K, sigma_C, removed_stars, r_squared, corr_coef)
    """
    m_inst = np.array(m_inst)
    m_std = np.array(m_std)
    x = np.array(x)
    color_ind = np.array(color_ind)
    star_names = np.array(star_names)

    removed_stars = []
    while True:
        # Формуємо матрицю коефіцієнтів A
        if k is None:
            a = np.column_stack((np.ones_like(m_inst), x, color_ind))
        else:
            a = np.column_stack((np.ones_like(m_inst), color_ind))

        # Вектор спостережень b
        if k is None:
            b = m_std - m_inst
        else:

            b = m_std - m_inst - k * x

        # Метод найменших квадратів
        coeffs, residuals, _, _ = np.linalg.lstsq(a, b, rcond=None)

        if k is not None:
            coeffs = np.insert(coeffs, 1, k)

        # Обчислюємо відхилення
        predicted = a @ coeffs[[0, 2]] if k is not None else a @ coeffs
        errors = np.abs(predicted - b)

        # Перевіряємо, чи є точки з великим відхиленням
        max_error_idx = np.argmax(errors)
        max_error_value = errors[max_error_idx]
        if max_error_value < threshold:
            break  # Якщо всі відхилення менші за поріг, виходимо з циклу

        # Зберігаємо інформацію про видалену зірку
        removed_stars.append((star_names[max_error_idx], max_error_value))
        print(f"Star: '{star_names[max_error_idx]}' was removed with error value = {max_error_value:4.3f}")
        if log_file is not None:
            log_file.write(
                f"Star: '{star_names[max_error_idx]}' was removed with error value = {max_error_value:4.3f}\n")

        # Видаляємо точку з найбільшим відхиленням
        m_inst = np.delete(m_inst, max_error_idx)
        m_std = np.delete(m_std, max_error_idx)
        x = np.delete(x, max_error_idx)
        color_ind = np.delete(color_ind, max_error_idx)
        star_names = np.delete(star_names, max_error_idx)

    # Оцінка похибок коефіцієнтів
    residual_variance = residuals / (len(b) - len(coeffs)) if len(residuals) > 0 else np.array([0])
    cov_matrix = np.linalg.inv(a.T @ a) * residual_variance
    errors = np.sqrt(np.diag(cov_matrix))

    if k is not None:
        errors = np.insert(errors, 1, 0)  # Додаємо sigma_K = 0, бо K зафіксований

    # Обчислюємо R²
    ss_total = np.sum((b - np.mean(b)) ** 2)
    ss_residual = np.sum((b - predicted) ** 2)
    r_squared = 1 - (ss_residual / ss_total)

    # Побудова графіка, якщо потрібно
    if plot:
        c = {"B": 'blue', "V": 'green', "R": 'red'}
        my_color = c[f'{band}']
        if path is None:
            path = os.getcwd()

        # Graph 1: Photometric Calibration
        plt.figure(figsize=(12, 8))
        plt.scatter(x, m_std - m_inst, label='Data', color=my_color)
        for i, name in enumerate(star_names):
            plt.annotate(name, (x[i], m_std[i] - m_inst[i]), fontsize=8, alpha=0.7)
        x_range = np.linspace(min(x), max(x), 100)
        fit_line = coeffs[0] + (coeffs[1] * x_range if k is None else k * x_range) + coeffs[2] * np.mean(color_ind)
        plt.plot(x_range, fit_line, color='red', label='Fit')
        plt.xlabel('Airmass')
        plt.ylabel('Magnitude Difference')
        plt.title(f'Photometric Calibration in {band} Filter')
        # plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(path, f"graph_{band}.png"))
        plt.close()

        # Graph 2: Residual Analysis
        plt.figure(figsize=(12, 8))
        residuals = b - predicted
        plt.scatter(predicted, residuals, color=my_color, alpha=0.6)
        plt.axhline(0, color='black', linestyle='--')
        for i, name in enumerate(star_names):
            plt.annotate(name, (predicted[i], residuals[i]), fontsize=8, alpha=0.7)
        plt.xlabel('Fitted Values')
        plt.ylabel('Residuals')
        plt.title('Residual Analysis')
        plt.tight_layout()
        plt.savefig(os.path.join(path, f"residuals_{band}.png"))
        plt.close()

        # Graph 3: Calibrated vs. Catalog Magnitudes
        plt.figure(figsize=(12, 8))
        calibrated_magnitudes = m_inst + coeffs[0] + (coeffs[1] * x if k is None else k * x) + coeffs[2] * color_ind
        plt.scatter(m_std, calibrated_magnitudes, color=my_color, alpha=0.6)
        plt.plot(m_std, m_std, color='black', linestyle='--', label='Ideal Fit')
        plt.xlabel('Catalog Magnitude')
        plt.ylabel('Calibrated Magnitude')
        plt.title('Calibrated vs. Catalog Magnitudes')
        plt.tight_layout()
        plt.savefig(os.path.join(path, f"calibrated_mags_{band}.png"))
        plt.close()

    corr_coef = pearsonr(m_std, calibrated_magnitudes)[0]

    return *coeffs, *errors, removed_stars, r_squared, corr_coef


def RMS_del(A, value, B=None, log_file=None):
    '''Delete elements of array A until A.RMS>value'''
    A = np.array(A)
    if B is not None:
        B = np.array(B)
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
        if B is not None:
            B = np.delete(B, imax)

    strFormat = len(A_del) * '{:5.3f}, '
    formattedList = strFormat.format(*A_del)

    print("N deleted =", len(A_del))
    if log_file:
        log_file.write("N deleted = %i \n" % len(A_del))
        log_file.write("Deleted value(s): " + formattedList + "\n")
    if B is not None:
        return A, B
    else:
        return A


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
    # print('Y=', y)
    if len(y) > 0:
        wb = np.linalg.lstsq(A, y, rcond=None)  # obtaining the parameters
        a, c = wb[0]
        residual_sum = wb[1]
        r2 = 1 - residual_sum / (y.size * y.var())  # = R^2
        # residual = wb[1][0]
    else:
        r2 = 0 #residual = 999

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
    return a, c, res_max, ind, r2[0]


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


def get_from_Gaia_EDR3_std(RA, DEC, w="0d60m", h="30m", radius=None, Filter={'Rmag': '<13'}):
    cat_name = 'J/A+A/664/A109'
    center = coord.SkyCoord(ra=RA, dec=DEC, unit=(u.deg, u.deg), frame='icrs')
    Vizier.ROW_LIMIT = 5000
    if radius:
        table = Vizier.query_region(center, width=w, radius=radius, catalog=[cat_name], column_filters=Filter)
    else:
        table = Vizier.query_region(center, width=w, height=h, catalog=[cat_name], column_filters=Filter)
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
            print(f"bx, by = {bx0:d},{by0:d}", end=" ")

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
                tle[1] = fix_checksum(tle[1])
                tle[2] = fix_checksum(tle[2])

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

                sun = ephem.Sun()
                sun.compute(station)
                my_phase = calc_phase(sun, sat)
                my_phase = math.degrees(my_phase)  # check ?

                return el, rg, az, n, no, c, tle, my_phase
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
            res['K'] = config.getfloat('Stars_Stand', 'K', fallback=None)
            res['max_m_fit'] = config.getfloat('Stars_Stand', 'max_m_fit', fallback="14")
            res['max_m_calc'] = config.getfloat('Stars_Stand', 'max_m_calc', fallback="14")

            res['r_max_val'] = config.getfloat('Stars_Stand', 'r_max_val', fallback=6.25)
            res['snr_value'] = config.getfloat('Stars_Stand', 'snr', fallback=1.2)
            res['band'] = config.get('Stars_Stand', 'filter')

            res['dark_frame'] = config.get('Stars_Stand', 'dark_frame', fallback=False)
            res['dark_stable'] = config.getfloat('Stars_Stand', 'dark_stable', fallback=0.0)
            res['fov'] = config.getfloat('Stars_Stand', 'FoV', fallback=1.5)

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

            res['header_time'] = config.get('STD', 'header_time', fallback="DATE-OBS")
            res['time_moment'] = config.get('STD', 'time_moment', fallback="start")
            if res['time_moment'] not in ['start', 'end', 'middle']:
                res['time_moment'] = 'start'  # fallback if there is a bullshit in config

            res['min_real_mag'] = config.getfloat('STD', 'min_real_mag', fallback=15.0)
            res['max_center_error'] = config.getfloat('STD', 'max_center_err', fallback=2.0)
            res['fits_sort'] = config.get('STD', 'fits_sort', fallback='None')

            res['dark_frame'] = config.get('STD', 'dark_frame', fallback=False)

            res['r_ap'] = float(config['APERTURE']['ap'])
            res['an_in'] = float(config['APERTURE']['an_in'])
            res['an_out'] = float(config['APERTURE']['an_out'])
            res['saturated'] = config.getboolean("APERTURE", "saturated", fallback=True)

            res['site_name'] = config['SITE']['Name']
            res['site_lat'] = config['SITE']['lat']
            res['site_lon'] = config['SITE']['lon']
            res['site_elev'] = float(config['SITE']['h'])

            res["auto_plot"] = config.getboolean('PLOT', 'auto_plot', fallback=False)
            res["plot_command"] = config.get('PLOT', 'plot_command', fallback='None')
            res["python"] = config.get('PLOT', 'python', fallback='None')

            res["obj_x"] = config.getfloat('OBJ_POS', 'OBJ_X', fallback=None)
            res["obj_y"] = config.getfloat('OBJ_POS', 'OBJ_Y', fallback=None)

            return res

        except Exception as E:
            print("Error in INI file\n", E)
            return False
    else:
        print("Error. Cant find config_sat.ini in " + conf_file)
        return False


def read_config_star_flux(conf_file):
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

            res['star_name'] = config['MAIN']['star_name']
            res['ra'] = config.get('MAIN',"RA", fallback=None)
            res['dec'] = config.get('MAIN',"DEC", fallback=None)

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


def calc_mag(flux, el, rg, zp, k, exp, min_mag=15, phase=None):
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
        return min_mag
    else:
        m_inst = -2.5 * math.log10(flux / exp)
        mz = 1 / (math.cos(math.pi / 2 - math.radians(el)))  # 1/cos(90-el)
        mr = -5 * math.log10(rg / 1000.0)
        mzr = k * mz

        m = zp + m_inst - mzr + mr

        # with open("tmp_mag.txt", "a+") as f:
        #     phase = math.degrees(phase)
        #     f.write(f"{m:5.3f}    {zp:5.3f}   {m_inst:5.3f}   {mzr:5.3f}  {mr:5.3f}  {el:5.3f}   {rg:5.3f}  {phase:5.3f}\n")
    return m


def calc_phase(sun, sat):
    # Calculate distance from observer to the sun, and to the satellite (length a and b of our triangle).
    sun_distance = (sun.earth_distance * ephem.meters_per_au) - ephem.earth_radius
    satellite_distance = sat.range

    # Calculate angle separating sun and satellite (angle C of our triangle)
    separation_angle = ephem.separation(sat, sun)

    # Calculate distance between sun and satellite (length c of our triangle)
    # c = sqrt(a*a + b*b - 2 * a * b * cos(C))
    sun_satellite_distance = math.sqrt((sun_distance * sun_distance) + (satellite_distance * satellite_distance) -
                                       (2 * sun_distance * satellite_distance * math.cos(separation_angle)))

    # Now we know the length of all sides of the triangle, calculate the phase angle (angle A of our triangle)
    # A = acos((b*b + c*c - a*a) / (2 * b * c))
    phase_angle = math.acos(((satellite_distance * satellite_distance) +
                             (sun_satellite_distance * sun_satellite_distance) - (sun_distance * sun_distance)) / (
                                        2 * satellite_distance * sun_satellite_distance))

    return phase_angle


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


def get_star_el(star_name, obs_lat, obs_lon, obs_elev, obs_date, star_ra_dec=None):
    """

    Parameters
    ----------
    star_name: name of star. Optional
    obs_lat, obs_lon, obs_elev: latitude and longitude of observatory (string). Example '52.5', '-1.91667', '223.3'
    obs_date: date of observation (datetime)
    star_ra_dec: ra/dec of star in a form of [ra, dec]. Example ['10:10:10', '+10:10:10']. Optional

    Returns
    -------
    az, el: star azimuth and elevation in degrees
    """
    if star_ra_dec is None:
        star_ra_dec = [0, 0]
        result = Vizier.query_object(str(star_name), catalog=['FK5'])[0]
        star_ra_dec[0] = result['RAJ2000'].value.data[0]  # str
        star_ra_dec[0] = star_ra_dec[0].replace(" ", ":")

        star_ra_dec[1] = result['DEJ2000'].value.data[0]
        star_ra_dec[1] = star_ra_dec[1].replace(" ", ":")

    star = ephem.FixedBody()
    star._ra = ephem.hours(star_ra_dec[0]) #('10:10:10')
    star._dec = ephem.degrees(star_ra_dec[1]) #('10:10:10')

    observer = ephem.Observer()
    observer.date = datetime.strptime(obs_date, "%Y-%m-%dT%H:%M:%S.%f0")

    # longitude = ephem.degrees('-1.91667')
    observer.lon = ephem.degrees(obs_lon)
    observer.lat = ephem.degrees(obs_lat)
    observer.elev = float(obs_elev)

    star.compute(observer)

    return star.az / ephem.degree, star.alt / ephem.degree  # float in degrees
