#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EOSSA FITS Converter v3.1.1
Конвертер фотометричних .phX файлів у FITS формат згідно стандарту EOSSA v3.1.1
"""

import os
import sys
import glob
import argparse
import numpy as np
from datetime import datetime
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, SkyCoord, get_sun
import astropy.units as u
import re

# === Constants === #
DISTANCE_NORM_KM = 1000.0
PLACEHOLDER_INT = -2147483648


# === Reading the header in phX === #
def read_header_to_dict(filepath):
    header_data = {}
    tle_lines = []
    in_tle_block = False
    key_value_pattern = re.compile(r'^#?\s*(\w+)\s*=\s*(.+)$')

    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            stripped = line.strip()
            if stripped == '':
                continue
            if stripped.startswith('# TLE:'):
                in_tle_block = True
                tle_lines = []
                continue
            if in_tle_block:
                tle_lines.append(stripped.lstrip('#').strip())
                if len(tle_lines) == 3:
                    header_data['TLE'] = tuple(tle_lines)
                    in_tle_block = False
                continue

            m = key_value_pattern.match(stripped)
            if m:
                key, val = m.group(1).strip(), m.group(2).strip()
                header_data[key.upper()] = val.replace('"', '')

    return header_data


# === Reading the table in phX === #
def read_data(fname):
    """
    Read table part of *.phX file, deleting char '#' from header
    Automatically determines the number of columns even with uneven spacing
    """
    with open(fname, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    # find first row that starts from 'Date' or 'DATE'
    header_line_idx = next((i for i, l in enumerate(lines) if re.match(r'#?\s*Date', l, re.IGNORECASE)), None)
    if header_line_idx is None:
        raise ValueError(f"Cant find header of the table in file {fname}")

    # Preparing a temporary file in memory for genfromtxt
    clean_lines = []
    for line in lines[header_line_idx:]:
        clean_lines.append(line.lstrip('#').strip())

    # Create structured array
    from io import StringIO
    buffer = StringIO('\n'.join(clean_lines))

    data = np.genfromtxt(
        buffer,
        names=True,
        dtype=None,
        encoding='utf-8',
        delimiter=None,  # any number of spaces
        autostrip=True
    )

    return data


# === Converting Az/El → RA/DEC === #
def azel_to_radec(az, el, obstime, location):
    altaz = AltAz(az=az*u.deg, alt=el*u.deg, obstime=obstime, location=location)
    sky = SkyCoord(altaz)
    return sky.icrs.ra.deg, sky.icrs.dec.deg


# === Calculating Sun phase angle === #
def compute_solar_phase_angle(obs_coord, sun_coord):
    return obs_coord.separation(sun_coord).deg


# ---- Helper function for searching columns by "fuzzy" name ----
def find_column_by_token(names, token):
    """
    names: iterable of column names (as returned by numpy recarray.dtype.names)
    token: normalized token, e.g. 'az', 'el', 'rg', 'date', 'ut', 'mag'
    returns: original name if found, else None
    """
    norm_map = {n: re.sub(r'[^0-9a-z]', '', n.lower()) for n in names}
    # 1) точний старт
    for orig, norm in norm_map.items():
        if norm.startswith(token):
            return orig
    # 2) містить
    for orig, norm in norm_map.items():
        if token in norm:
            # avoid matching 'mag_err' when token == 'mag' (we'll check separately)
            return orig
    return None


def convert_to_fits(input_file):
    print(f"→ Processing: {os.path.basename(input_file)}")

    header = read_header_to_dict(input_file)
    data = read_data(input_file)  # numpy recarray

    # Read filter from header of extension of the phX file
    filt = header.get('FILTER')
    if not filt:
        filt = os.path.splitext(input_file)[1][-1].upper()

    names = list(data.dtype.names)
    # find mag col: any col that starts from 'mag' but not ends on 'err'
    mag_col = None
    for n in names:
        nl = n.lower()
        if nl.startswith('mag') and not nl.endswith('err'):
            mag_col = n
            break
    if mag_col is None:
        # as a fallback, let's try to find columns containing 'mag' but not 'err'
        for n in names:
            nl = re.sub(r'[^0-9a-z]', '', n.lower())
            if 'mag' in nl and not nl.endswith('err'):
                mag_col = n
                break

    if mag_col is None:
        raise ValueError("Cannot find eny column of magnitude (mag*). Available cols: " + ", ".join(names))

    # find date / time
    date_col = find_column_by_token(names, 'date')
    time_col = find_column_by_token(names, 'ut') or find_column_by_token(names, 'time')

    if date_col is None or time_col is None:
        raise ValueError(f"Не знайдено колонки дати/часу. Доступні: {', '.join(names)}")

    # find az/el/rg
    az_col = find_column_by_token(names, 'az')
    el_col = find_column_by_token(names, 'el')
    rg_col = find_column_by_token(names, 'rg') or find_column_by_token(names, 'range')

    if az_col is None or el_col is None:
        raise ValueError(f"Не знайдено колонок az/el. Доступні: {', '.join(names)}")

    if rg_col is None:
        # Rg can be marked differently; we will put NaN, but it is better to warn
        print("[Warning] Cant find range loc (Rg). Tel_Obj_Range will be filled with NaN.")

    # Now we take the arrays from data by the found names
    mag_norm = data[mag_col]
    date_strs = data[date_col]
    time_strs = data[time_col]
    az = data[az_col].astype(float)
    el = data[el_col].astype(float)
    dist_m = np.array(
        [float(x) if rg_col and x != '' else np.nan for x in (data[rg_col] if rg_col else [np.nan] * len(data))])
    if rg_col:
        # If Rg is in Mm should convert it to meters
        # let's try to guess: if the value of the order 1000..2000 is Mm
        # If the values are large (>=1e5) — these are meters -> we don't multiply
        sample = dist_m[~np.isnan(dist_m)]
        if sample.size > 0:
            median = np.median(sample)
            if median < 1e5:  # most likely it is Mm
                dist_m = dist_m * 1e6
    # Time: reliably parse strings with fractional seconds
    utc_times = []
    for d, t in zip(date_strs, time_strs):
        s = f"{d} {t}"
        # some files may have microseconds truncated; try multiple formats
        tried = False
        for fmt in ('%Y-%m-%d %H:%M:%S.%f', '%Y-%m-%d %H:%M:%S'):
            try:
                dtobj = datetime.strptime(s, fmt)
                utc_times.append(dtobj)
                tried = True
                break
            except Exception:
                continue
        if not tried:
            # last resort: use astropy to parse
            try:
                utc_times.append(Time(s).to_datetime())
            except Exception:
                raise ValueError(f"Unable to parse date/time string: '{s}'")
    obstimes = Time(utc_times, scale='utc')

    # RA/DEC
    location = EarthLocation(lat=float(header.get('SITE_LAT', 0)) * u.deg,
                             lon=float(header.get('SITE_LON', 0)) * u.deg,
                             height=float(header.get('SITE_ELEV', 0)) * u.m)
    ra, dec = azel_to_radec(az, el, obstimes, location)

    # Mag_Exo_Atm (reduction to ephemeris distance)
    # if dist_m NaN -> Mag_Exo_Atm = NaN
    mag_exo = np.array([m + 5.0 * np.log10(d / 1e3 / DISTANCE_NORM_KM) if np.isfinite(d) and d > 0 else np.nan
                        for m, d in zip(mag_norm, dist_m)])

    # Sun, phase angle
    sun_coords = get_sun(obstimes)
    obs_coords = SkyCoord(ra * u.deg, dec * u.deg, frame='icrs')
    solar_phase = compute_solar_phase_angle(obs_coords, sun_coords)

    # Formation of the FITS table
    dt = float(header['DT'])
    cols = [
        fits.Column(name='UTC_Begin_exp', format='23A', array=[t.iso.replace(' ', 'T') for t in obstimes]),
        fits.Column(name='UTC_End_exp', format='23A', array=[(t + dt*u.s).iso.replace(' ', 'T') for t in obstimes]),
        fits.Column(name='JD_Mid_Exp', format='D', array=obstimes.jd + dt/(2*86400.0)),
        fits.Column(name='Exp_Duration', format='D', array=[dt]*len(obstimes)),
        fits.Column(name='Cur_Spec_Filt_Num', format='J', array=[1]*len(obstimes)),
        fits.Column(name='Cur_ND_Filt_Num', format='J', array=[PLACEHOLDER_INT]*len(obstimes)),
        fits.Column(name='Mag_Exo_Atm', format='D', array=mag_exo),
        fits.Column(name='Mag_Range_Norm', format='D', array=mag_norm),
        fits.Column(name='Tel_Obj_Range', format='D', array=dist_m),
        fits.Column(name='Eph_RA_DE', format='2D', array=np.column_stack((ra, dec))),
        fits.Column(name='Met_RA_DE', format='2D', array=np.column_stack((ra, dec))),
        fits.Column(name='Eph_Az_El', format='2D', array=np.column_stack((az, el))),
        fits.Column(name='Met_Az_El', format='2D', array=np.column_stack((az, el))),
        fits.Column(name='Solar_Phase_Ang', format='D', array=solar_phase)
    ]

    hdu = fits.BinTableHDU.from_columns(cols)

    # Primary Header
    phdr = fits.Header()
    phdr['SIMPLE'] = (True, 'FITS standard')
    phdr['BITPIX'] = (8, 'bits per data pixel')
    phdr['NAXIS'] = (0, 'no image data')
    phdr['EXTEND'] = (True, 'extensions are possible')
    phdr['CLASSIF'] = ('UNCLASS', 'Security classification level')
    phdr['EOSSAVER'] = ('3.1.1', 'EOSSA data standard version')
    phdr['OBSEPH'] = ('GROUND', 'Observation ephemeris source')
    phdr['FILTER'] = (filt, 'Observation filter')
    phdr['OBJECT'] = (header.get('NAME', 'Unknown'), 'Observed object')
    phdr['NORAD'] = (header.get('NORAD', 'Unknown'), 'NORAD ID')

    primary_hdu = fits.PrimaryHDU(header=phdr)

    # Output path
    wd = os.path.dirname(input_file)
    fname = os.path.splitext(os.path.basename(input_file))[0]
    eossa_path = os.path.join(wd, fname + '.eossa.fits')

    # Writing file
    hdul = fits.HDUList([primary_hdu, hdu])
    hdul.writeto(eossa_path, overwrite=True)

    print(f"✓ Created: {os.path.basename(eossa_path)}")


# === CLI === #
def main():
    parser = argparse.ArgumentParser(description='Конвертер EOSSA v3.1.1 для .phX файлів')
    parser.add_argument('input', nargs='+', help='Вхідні файли або шаблон (*.ph*)')
    args = parser.parse_args()

    all_files = []
    for pattern in args.input:
        if any(ch in pattern for ch in ['*', '?']):
            matched = glob.glob(pattern)
            all_files.extend(matched)
        elif os.path.isfile(pattern):
            all_files.append(pattern)

    if not all_files:
        print("[!] No files found to process.")
        sys.exit(1)

    print(f"→ Found {len(all_files)} files to process.\n")

    for file in all_files:
        try:
            convert_to_fits(file)
        except Exception as e:
            print(f"[Error] {os.path.basename(file)}: {e}")

    print("\nDone ✅")


if __name__ == '__main__':
    main()
