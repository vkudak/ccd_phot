import numpy as np
# import astropy.io.fits as fits
import os
# from astropy.stats import mad_std
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import ephem
from datetime import datetime
import math
from astropy import units as u
from astroquery.vizier import Vizier
import astropy.coordinates as coord
from tqdm import tqdm
# from sklearn.linear_model import LinearRegression
from astropy.modeling import Fittable2DModel, Parameter, models, fitting
from astropy.modeling.functional_models import Moffat2D, Gaussian2D
import warnings


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


def R2_calc(calc, observed):
    Zmm = calc
    Z3 = observed
    absError = Zmm - Z3
    SE = np.square(absError)  # squared errors
    MSE = np.mean(SE)  # mean squared errors
    RMSE = np.sqrt(MSE)  # Root Mean Squared Error, RMSE
    Rsquared = 1.0 - (np.var(absError) / np.var(Z3))
    return RMSE, Rsquared


def fit_lin_reg(x, y, yerr=None):
    # popt, pconv = curve_fit(linear, x, y, sigma=yerr)
    if yerr is not None:
        yerr = np.array(yerr)
        popt, pconv = np.polyfit(x, y, 1, cov=True, w=1.0 / yerr)
    else:
        popt, pconv = np.polyfit(x, y, 1, cov=True)
    a, b = popt
    ae, be = np.sqrt(pconv.diagonal())

    # residuals = y - linear(x, *popt)
    # ss_res = np.sum(residuals**2)
    # ss_tot = np.sum((y - np.mean(y))**2)
    # r2 = 1 - (ss_res / ss_tot)
    rmse, r2 = R2_calc(linear(x, *popt), y)

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
        var_table = Vizier.query_region(c, radius='0d0m3s', catalog=["GCVS"], column_filters=Filter)  # , catalog=["GCVS"]
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
    A = 1 / (2 * sigma**2)
    eq = amp * np.exp(-A * ((x - x0)**2 + (y - y0)**2))  # Gaussian
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

    params, pcov = curve_fit(gaussian, (xx.ravel(), yy.ravel()), image.ravel(), p0=guess, bounds=bounds, maxfev=mf)  # optimal fit.  Not sure what pcov is.

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
            plt.title("parameters= %s \nerrors=%s \ngate=%s  x,y=%2.3f, %2.3f" % (par, err, gate, x, y), pad=-50, fontsize=10)
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
            print("centring...")
        # ------------------------------------------------------------------ find brighter pixel and center it!!!!
        bx0, by0 = np.unravel_index(data_fit.argmax(), data_fit.shape)
        if not silent:
            print("bx, by=", bx0, by0)

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
    par, pcov = fit(data_fit)

    err = np.sqrt(np.diag(pcov))
    amp = par[-1]
    target = [x0 - gate + par[0], y0 - gate + par[1], err[0], err[1], amp]
    if debug:
        # plotting(data_fit, par, save=True, filename=fig_name, par=target, gate=gate)
        plotting(data_fit, par, save=True, filename=fig_name, par=par, err=err, tar=target, gate=gate)
    return target


def fit_moff(image, x0, y0, gate, debug, fig_name=None, centring=False, silent=False):
    data_fit = image[y0 - gate:y0 + gate, x0 - gate:x0 + gate]
    if centring:
        if not silent:
            print("centring...")
        # ------------------------------------------------------------------ find brighter pixel and center it!!!!
        bx0, by0 = np.unravel_index(data_fit.argmax(), data_fit.shape)
        if not silent:
            print("bx, by=", bx0, by0)

        sx = by0 - gate
        sy = bx0 - gate

        data_fit = image[y0 - gate + sy:y0 + gate + sy, x0 - gate + sx:x0 + gate + sx]
        # ----------------------------------------------------------------------------------------
    Z3 = data_fit
    X3 = np.arange(0, gate * 2, 1)
    Y3 = np.arange(0, gate * 2, 1)
    X3, Y3 = np.meshgrid(X3, Y3)

    m_init = Moffat2D(amplitude=np.max(Z3), x_0=gate, y_0=gate, gamma=1., alpha=1.0)
    fit_m = fitting.LevMarLSQFitter()
    with warnings.catch_warnings():
        # Ignore model linearity warning from the fitter
        warnings.simplefilter('ignore')
        m = fit_m(m_init, X3, Y3, Z3)
    # print("#####Moffat#########")
    # print(m.x_0.value)
    # print(m.y_0.value)

    Zmm = m(X3, Y3)
    RMSE, Rsquared = R2_calc(Zmm, Z3)
    # print("x=%2.3f  y=%2.3f   R^2=%2.4f" % (m.x_0.value + xs - gate, m.y_0.value + ys - gate, Rsquared))
    amp = np.max(Z3)
    target = [x0 - gate + m.x_0.value, y0 - gate + m.y_0.value, m.fwhm, Rsquared, amp]
    if debug:
        # plotting(data_fit, par, save=True, filename=fig_name, par=target, gate=gate)
        par = [m.x_0.value, m.y_0.value, m.fwhm]
        plotting(data_fit, par, save=True, filename=fig_name, par=par, err=[0., 0.], tar=target, gate=gate)
    return target


def fit_gaus(image, x0, y0, gate, debug, fig_name=None, centring=False, silent=False):
    data_fit = image[y0 - gate:y0 + gate, x0 - gate:x0 + gate]
    if centring:
        if not silent:
            print("centring...")
        # ------------------------------------------------------------------ find brighter pixel and center it!!!!
        bx0, by0 = np.unravel_index(data_fit.argmax(), data_fit.shape)
        if not silent:
            print("bx, by=", bx0, by0)

        sx = by0 - gate
        sy = bx0 - gate

        data_fit = image[y0 - gate + sy:y0 + gate + sy, x0 - gate + sx:x0 + gate + sx]
        # ----------------------------------------------------------------------------------------
    Z3 = data_fit
    X3 = np.arange(0, gate * 2, 1)
    Y3 = np.arange(0, gate * 2, 1)
    X3, Y3 = np.meshgrid(X3, Y3)

    sigma = np.std(Z3)
    g_init = Gaussian2D(amplitude=np.max(Z3), x_mean=gate, y_mean=gate)
    fit_g = fitting.LevMarLSQFitter()
    with warnings.catch_warnings():
        # Ignore model linearity warning from the fitter
        warnings.simplefilter('ignore')
        g = fit_g(g_init, X3, Y3, Z3)
    # print("#####Moffat#########")
    # print(m.x_0.value)
    # print(m.y_0.value)

    Zmm = g(X3, Y3)
    RMSE, Rsquared = R2_calc(Zmm, Z3)
    # print("x=%2.3f  y=%2.3f   R^2=%2.4f" % (m.x_0.value + xs - gate, m.y_0.value + ys - gate, Rsquared))
    amp = np.max(Z3)
    fwhm = (g.x_fwhm  + g.y_fwhm) /2.
    target = [x0 - gate + g.x_mean.value, y0 - gate + g.y_mean.value, fwhm, Rsquared, amp]
    if debug:
        # plotting(data_fit, par, save=True, filename=fig_name, par=target, gate=gate)
        par = [g.x_mean.value, g.y_mean.value, fwhm]
        plotting(data_fit, par, save=True, filename=fig_name, par=par, err=[0., 0.], tar=target, gate=gate)
    return target


def get_tle(tle_path):
    TLE_list = []
    ephf = open(tle_path, 'r')
    i = 1
    tle1 = []
    for l in ephf:
        if l != '':
            tle1.append(l.strip())
            i = i + 1
            if i > 3:
                TLE_list.append(tle1)
                tle1 = []
                i = 1
    ephf.close()
    return TLE_list


def calc_from_tle(TLE_list, date_time, COSPAR, NORAD, NAME):
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
                station.lat = '48.5635505'
                station.long = '22.453751'
                station.elevation = 231.1325
                sat = ephem.readtle(tle[0], tle[1], tle[2])
                try:
                    station.date = datetime.strptime(date_time[:-1], "%Y-%m-%dT%H:%M:%S.%f")
                except Exception:
                    print ("Error - Wrong date time format...")

                sat.compute(station)
                el = math.degrees(sat.alt)
                rg = sat.range / 1000.0  # km
                az = math.degrees(sat.az)
                return el, rg, az, n, no, c, tle
            else:
                print('Cant find TLE for such satellite!')
        except Exception as ex:
            print("Error occurred. Probably wrong TLE file. \n Error = " + ex.message)
