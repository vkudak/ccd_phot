#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EOSSA v3.1.1 compliant converter for .phX photometry files -> EOSSA FITS (binary table)

One-file, modular script. Functions:
 - read_phx_file(path) -> (header_dict, recarray)
 - build_eossa_table(header, table_rows) -> numpy structured array (EOSSA dtype)
 - write_eossa_fits(outpath, primary_meta, eossa_array)

Notes / assumptions:
 - Input .phX files contain header lines beginning with '#' as in example and a whitespace-separated table
   whose first non-comment line is the column names (e.g. "Date  UT  X  Y ... filename").
 - Az/El in degrees are present (columns named like 'Az(deg)' and 'El(deg)').
 - Rg column is in Mm (column name 'Rg(Mm)') as in example.
 - If RA/DEC are not present, they will be computed from Az/El using observer location from header.
 - Eph_* and Met_* will be identical (we only have ephemeris-based values), per user instruction.
 - Solar disk fraction is approximated geometrically (linear interpolation across solar radius).

Requires: astropy, numpy

"""

from __future__ import annotations

import os
import re
import math
from typing import Tuple, Dict, Any

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, SkyCoord, get_sun
import astropy.units as u


# --- Constants / placeholders ---
PLACEHOLDER_INT = -2147483648
PLACEHOLDER_FLOAT = -9999.0
EOSSA_VERSION = "3.1.1"
CLASSIF_VALUE = "UNCLASS"

# approximate solar angular radius (deg) — small variation with distance ignored here
SOLAR_ANGULAR_RADIUS_DEG = 0.2666

# Required ground-based K-2 columns (14 fields). We'll add a few related/error fields that are common.
EOSSA_COLUMNS = [
    ("UTC_Begin_exp", "U23"),      # 23A / string
    ("UTC_End_exp", "U23"),        # 23A
    ("JD_Mid_Exp", "f8"),          # D
    ("Exp_Duration", "f8"),        # D
    ("Cur_Spec_Filt_Num", "i4"),   # J
    ("Cur_ND_Filt_Num", "i4"),     # J
    ("Mag_Exo_Atm", "f8"),        # D
    ("Mag_Exo_Atm_Unc", "f8"),    # D (recommended)
    ("Mag_Range_Norm", "f8"),     # D
    ("Eph_RA_DE", ("f8", 2)),      # 2D
    ("Met_RA_DE", ("f8", 2)),      # 2D
    ("Eph_AZ_EL", ("f8", 2)),      # 2D
    ("Met_AZ_EL", ("f8", 2)),      # 2D
    ("Sun_AZ_EL", ("f8", 2)),      # 2D
    ("Tel_Obj_Range", "f8"),      # 1D meters
    ("Solar_Disk_Frac", "f8"),    # 1D (0..1)
]


# --------------------- I/O: read .phX ---------------------
def read_phx_file(path: str) -> Tuple[Dict[str, str], np.ndarray]:
    """Parse .phX file.

    Returns:
        header: dictionary of header key -> value (from lines starting with '#')
        data: numpy.recarray of table data (column names taken from first non-comment line after header)
    """
    header: Dict[str, str] = {}
    tle_block = []
    table_lines = []
    with open(path, 'r', encoding='utf-8') as fh:
        in_tle = False
        for raw in fh:
            line = raw.rstrip('\n')
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith('#'):
                body = stripped.lstrip('#').strip()
                # handle TLE block marker
                if body == 'TLE:':
                    in_tle = True
                    tle_block = []
                    continue
                if in_tle:
                    tle_block.append(body)
                    if len(tle_block) == 3:
                        header['TLE_LINE0'] = tle_block[0]
                        header['TLE_LINE1'] = tle_block[1]
                        header['TLE_LINE2'] = tle_block[2]
                        in_tle = False
                    continue
                # parse key = value if present
                if '=' in body:
                    parts = re.split(r'\s*=\s*', body, maxsplit=1)
                    if len(parts) == 2:
                        key, val = parts
                        header[key.strip()] = val.strip()
                        continue
                # other header-only lines (for example start/end times without '=')
                # attempt to detect ISO datetimes (e.g. 2025-08-18 19:12:57.981767)
                if re.match(r'^\d{4}-\d{2}-\d{2}\s+\d{2}:\d{2}:\d{2}', body):
                    # store sequentially as START_TIME and END_TIME if not present
                    if 'START_TIME' not in header:
                        header['START_TIME'] = body
                    elif 'END_TIME' not in header:
                        header['END_TIME'] = body
                    continue
                # otherwise skip
            else:
                # non-comment line -> part of table; include as-is
                table_lines.append(line)

    if not table_lines:
        raise ValueError(f"No data table found in {path}")

    # The first non-comment line is header (column names)
    # Normalize spaces and split -> column names
    header_line = table_lines[0].strip()
    colnames = re.split(r'\s+', header_line)

    # The remaining lines are data rows; use np.genfromtxt on a string buffer
    from io import StringIO
    data_text = '\n'.join(table_lines[1:])
    # Use genfromtxt with names=colnames
    try:
        data = np.genfromtxt(StringIO(data_text), dtype=None, names=colnames, encoding='utf-8')
    except Exception:
        # fallback: read as float columns and build recarray
        data = np.genfromtxt(StringIO(data_text), names=colnames, dtype=None, encoding='utf-8')

    return header, data


# --------------------- Utility: parse date/time & compute mid/JD ---------------------
def make_utc_strings_and_mid_jd(date_str: str, time_str: str, dt_seconds: float) -> Tuple[str, str, float]:
    """Given date 'YYYY-MM-DD' and time 'HH:MM:SS.sss', return UTC_begin_str, UTC_end_str (both ISO like yyyy-mm-ddThh:mm:ss.sss) and JD_mid (float).

    Note: we keep 6 fractional seconds digits (microseconds) where available.
    """
    # Try to parse time with optional fractional seconds
    iso_in = f"{date_str} {time_str}"
    # Some files may have 3 or 6 fractional digits; use astropy Time for reliable JD
    t_start = Time(iso_in, format='iso', scale='utc')
    t_end = Time(t_start.jd + (dt_seconds / 86400.0), format='jd', scale='utc')
    t_mid = Time((t_start.jd + t_end.jd) / 2.0, format='jd', scale='utc')

    # Format strings to 23 chars like 'yyyy-mm-ddThh:mm:ss.sss' (keep microseconds up to 6 digits)
    utc_begin = t_start.iso.replace(' ', 'T')[:23]
    utc_end = t_end.iso.replace(' ', 'T')[:23]
    jd_mid = t_mid.jd
    return utc_begin, utc_end, jd_mid


# --------------------- Compute RA/DEC from Az/El ---------------------
def azel_to_radec(az_deg: float, el_deg: float, obs_time: Time, site_lat: float, site_lon: float, site_elev_m: float) -> Tuple[float, float]:
    """Convert Az/El (deg) in local horizon frame to ICRS RA/Dec (deg) for given observer and time.

    Returns (ra_deg, dec_deg).
    """
    location = EarthLocation(lat=site_lat * u.deg, lon=site_lon * u.deg, height=site_elev_m * u.m)
    altaz = AltAz(obstime=obs_time, location=location)
    sc = SkyCoord(az=az_deg * u.deg, alt=el_deg * u.deg, frame=altaz)
    icrs = sc.transform_to('icrs')
    return float(icrs.ra.degree), float(icrs.dec.degree)


# --------------------- Build EOSSA structured array ---------------------
def build_eossa_array(header: Dict[str, str], data_rows: np.ndarray) -> np.ndarray:
    """Create a NumPy structured array conforming to the required EOSSA ground-based columns (K-2).

    header: parsed header dict (should contain SITE_LAT, SITE_LON, SITE_ELEV, dt, NORAD, NAME, etc.)
    data_rows: numpy.recarray with columns from file (expected: Date, UT, Az(deg), El(deg), Rg(Mm), magB or mag etc.)
    """
    n = len(data_rows)
    # Build dtype dynamically to include vector 2D where required
    dtype_list = []
    for name, fmt in EOSSA_COLUMNS:
        if isinstance(fmt, tuple):
            # ("f8", 2)
            dtype_list.append((name, np.float64, fmt[1]))
        else:
            dtype_list.append((name, np.float64 if fmt == 'f8' else (np.int32 if fmt == 'i4' else 'U23')))

    # But we want exact types: override mapping
    # Let's explicitly define dtype according to EOSSA_COLUMNS
    dtype = np.dtype([
        ('UTC_Begin_exp', 'U23'),
        ('UTC_End_exp', 'U23'),
        ('JD_Mid_Exp', np.float64),
        ('Exp_Duration', np.float64),
        ('Cur_Spec_Filt_Num', np.int32),
        ('Cur_ND_Filt_Num', np.int32),
        ('Mag_Exo_Atm', np.float64),
        ('Mag_Exo_Atm_Unc', np.float64),
        ('Mag_Range_Norm', np.float64),
        ('Eph_RA_DE', np.float64, (2,)),
        ('Met_RA_DE', np.float64, (2,)),
        ('Eph_AZ_EL', np.float64, (2,)),
        ('Met_AZ_EL', np.float64, (2,)),
        ('Sun_AZ_EL', np.float64, (2,)),
        ('Tel_Obj_Range', np.float64),
        ('Solar_Disk_Frac', np.float64),
    ])

    out = np.empty(n, dtype=dtype)

    # Observer location
    try:
        site_lat = float(header.get('SITE_LAT', header.get('TELLAT', '0')))
        site_lon = float(header.get('SITE_LON', header.get('TELLONG', '0')))
        site_elev = float(header.get('SITE_ELEV', header.get('TELALT', '0')))
    except Exception:
        site_lat = 0.0
        site_lon = 0.0
        site_elev = 0.0

    # exposure dt
    try:
        dt_val = float(header.get('dt', header.get('DT', header.get('Exp_Duration', 0.0))))
    except Exception:
        dt_val = 0.0

    # filter (from filename extension .phX)
    filt = None
    m = re.search(r"\.ph([A-Za-z0-9])$", os.path.basename(header.get('FILENAME', '')))
    if m:
        filt = m.group(1)
    # default spec filter index
    cur_spec_filt_num = 1

    # parse Rg column name possibilities
    rg_name = None
    for cand in ['Rg(Mm)', 'Rg', 'Rg_Mm', 'Rg(Mm)']:
        if cand in data_rows.dtype.names:
            rg_name = cand
            break
    # parse az/el names
    az_name = next((c for c in data_rows.dtype.names if 'Az' in c), None)
    el_name = next((c for c in data_rows.dtype.names if 'El' in c or 'el' in c), None)
    # parse mag name
    mag_name = next((c for c in data_rows.dtype.names if re.search(r"mag", c, flags=re.I)), None)
    mag_err_name = next((c for c in data_rows.dtype.names if re.search(r"mag_err|magerr|mag_err", c, flags=re.I)), None)

    # date/time names (common: Date and UT)
    date_col = next((c for c in data_rows.dtype.names if re.match(r'Date', c, flags=re.I)), None)
    time_col = next((c for c in data_rows.dtype.names if re.match(r'UT|Time', c, flags=re.I)), None)

    # If data_rows[time_col] contains microseconds in different format, astropy Time still parses

    for i in range(n):
        row = data_rows[i]
        date_str = str(row[date_col]) if date_col else ''
        time_str = str(row[time_col]) if time_col else ''
        # build begin/end/jd
        if date_str and time_str:
            try:
                utc_begin, utc_end, jd_mid = make_utc_strings_and_mid_jd(date_str, time_str, dt_val)
            except Exception:
                # fallback: try without fractional parts
                t_iso = f"{date_str} {time_str}"
                t0 = Time(t_iso, format='iso', scale='utc')
                t_end = Time(t0.jd + (dt_val / 86400.0), format='jd', scale='utc')
                t_mid = Time((t0.jd + t_end.jd) / 2.0, format='jd', scale='utc')
                utc_begin = t0.iso.replace(' ', 'T')[:23]
                utc_end = t_end.iso.replace(' ', 'T')[:23]
                jd_mid = t_mid.jd
        else:
            utc_begin = ''
            utc_end = ''
            jd_mid = np.nan

        out['UTC_Begin_exp'][i] = utc_begin
        out['UTC_End_exp'][i] = utc_end
        out['JD_Mid_Exp'][i] = jd_mid
        out['Exp_Duration'][i] = dt_val
        out['Cur_Spec_Filt_Num'][i] = cur_spec_filt_num
        out['Cur_ND_Filt_Num'][i] = PLACEHOLDER_INT

        # magnitude and range
        mag_val = float(row[mag_name]) if mag_name and row[mag_name] != '' else np.nan
        mag_unc = float(row[mag_err_name]) if mag_err_name and row[mag_err_name] != '' else np.nan
        rg_mm = float(row[rg_name]) if rg_name and row[rg_name] != '' else np.nan

        # Range-normalized magnitude (1000 km) per doc: Mag_Range_Norm = mag + 5*log10(d_km / 1000)
        # since d_km = rg_mm * 1000 => d_km/1000 = rg_mm
        if np.isfinite(mag_val) and np.isfinite(rg_mm) and rg_mm > 0:
            mag_range_norm = mag_val + 5.0 * math.log10(rg_mm)
        else:
            mag_range_norm = np.nan

        # Mag_Exo_Atm: without atmospheric correction available, set equal to Mag_Range_Norm
        mag_exo_atm = mag_range_norm

        out['Mag_Exo_Atm'][i] = mag_exo_atm if np.isfinite(mag_exo_atm) else PLACEHOLDER_FLOAT
        out['Mag_Exo_Atm_Unc'][i] = mag_unc if np.isfinite(mag_unc) else PLACEHOLDER_FLOAT
        out['Mag_Range_Norm'][i] = mag_range_norm if np.isfinite(mag_range_norm) else PLACEHOLDER_FLOAT

        # Az/El
        az = float(row[az_name]) if az_name and row[az_name] != '' else np.nan
        el = float(row[el_name]) if el_name and row[el_name] != '' else np.nan
        out['Eph_AZ_EL'][i, 0] = az if np.isfinite(az) else PLACEHOLDER_FLOAT
        out['Eph_AZ_EL'][i, 1] = el if np.isfinite(el) else PLACEHOLDER_FLOAT
        out['Met_AZ_EL'][i, 0] = out['Eph_AZ_EL'][i, 0]
        out['Met_AZ_EL'][i, 1] = out['Eph_AZ_EL'][i, 1]

        # Tel_Obj_Range (meters)
        out['Tel_Obj_Range'][i] = rg_mm * 1e6 if np.isfinite(rg_mm) else np.nan

        # Compute RA/DEC from Az/El (use mid-exposure time)
        try:
            obs_time = Time(out['JD_Mid_Exp'][i], format='jd', scale='utc')
            ra_deg, dec_deg = azel_to_radec(az, el, obs_time, site_lat, site_lon, site_elev)
            out['Eph_RA_DE'][i, 0] = ra_deg
            out['Eph_RA_DE'][i, 1] = dec_deg
            out['Met_RA_DE'][i, 0] = ra_deg
            out['Met_RA_DE'][i, 1] = dec_deg
        except Exception:
            out['Eph_RA_DE'][i, :] = PLACEHOLDER_FLOAT
            out['Met_RA_DE'][i, :] = PLACEHOLDER_FLOAT

        # Sun AZ/EL at mid time
        try:
            tmid = Time(out['JD_Mid_Exp'][i], format='jd', scale='utc')
            location = EarthLocation(lat=site_lat * u.deg, lon=site_lon * u.deg, height=site_elev * u.m)
            sun = get_sun(tmid)
            sun_altaz = sun.transform_to(AltAz(obstime=tmid, location=location))
            sun_az = float(sun_altaz.az.deg)
            sun_alt = float(sun_altaz.alt.deg)
            out['Sun_AZ_EL'][i, 0] = sun_az
            out['Sun_AZ_EL'][i, 1] = sun_alt
            # estimate solar disk fraction visible above horizon (simple linear approx):
            r = SOLAR_ANGULAR_RADIUS_DEG
            frac = (sun_alt + r) / (2.0 * r)
            frac = max(0.0, min(1.0, frac))
            out['Solar_Disk_Frac'][i] = frac
        except Exception:
            out['Sun_AZ_EL'][i, :] = PLACEHOLDER_FLOAT
            out['Solar_Disk_Frac'][i] = PLACEHOLDER_FLOAT

    return out


# --------------------- Write FITS ---------------------
def write_eossa_fits(out_path: str, header_meta: Dict[str, Any], eossa_array: np.ndarray):
    """Write primary header + EOSSA binary table to out_path (overwrite).

    header_meta: dictionary possibly containing TELESCOP, SITE_LAT, SITE_LON, SITE_ELEV, NORAD, NAME, etc.
    eossa_array: structured array matching EOSSA dtype
    """
    # Primary HDU
    prihdr = fits.Header()
    prihdr['SIMPLE'] = (True, 'standard FITS')
    prihdr['BITPIX'] = (8, 'bytes')
    prihdr['NAXIS'] = (0, 'no image')
    prihdr['EXTEND'] = (True, 'extensions follow')
    prihdu = fits.PrimaryHDU(header=prihdr)

    # Table HDU
    bintab = fits.BinTableHDU(data=eossa_array)

    # Fill required EOSSA header keywords for ground-based (Appendix K table K-1)
    # Use defaults if missing
    bintab.header['EXTNAME'] = os.path.basename(out_path)
    bintab.header['CLASSIF'] = header_meta.get('CLASSIF', CLASSIF_VALUE)
    bintab.header['VERS'] = EOSSA_VERSION
    bintab.header['OBSEPH'] = header_meta.get('OBSEPH', 'GROUND')
    bintab.header['TELESCOP'] = header_meta.get('SITE_NAME', header_meta.get('TELESCOP', 'UNKNOWN'))
    # lat/lon/elev
    try:
        bintab.header['TELLAT'] = float(header_meta.get('SITE_LAT', header_meta.get('TELLAT', 0.0)))
        bintab.header['TELLONG'] = float(header_meta.get('SITE_LON', header_meta.get('TELLONG', 0.0)))
        bintab.header['TELALT'] = float(header_meta.get('SITE_ELEV', header_meta.get('TELALT', 0.0)))
    except Exception:
        pass

    bintab.header['OBSNAME'] = header_meta.get('SITE_NAME', '')
    bintab.header['OBJEPH'] = header_meta.get('OBJEPH', 'TLE')
    bintab.header['OBJTYPE'] = header_meta.get('OBJTYPE', 'SCN')
    if 'NORAD' in header_meta:
        try:
            bintab.header['OBJNUM'] = int(header_meta['NORAD'])
        except Exception:
            bintab.header['OBJNUM'] = header_meta['NORAD']
    bintab.header['OBJECT'] = header_meta.get('NAME', header_meta.get('OBJECT', ''))

    # Add TLE lines if available
    if 'TLE_LINE0' in header_meta and 'TLE_LINE1' in header_meta:
        bintab.header['TLELN1'] = header_meta.get('TLE_LINE1', '')
        bintab.header['TLELN2'] = header_meta.get('TLE_LINE2', '')

    hdul = fits.HDUList([prihdu, bintab])
    hdul.writeto(out_path, overwrite=True)
    return out_path


# --------------------- High-level converter ---------------------
def convert_phx_to_eossa(phx_path: str, out_fits: str | None = None) -> str:
    if out_fits is None:
        base = os.path.splitext(os.path.basename(phx_path))[0]
        out_fits = base + '.eossa.fits'

    header, data = read_phx_file(phx_path)
    # add filename to header so build_eossa_array can see it if needed
    header['FILENAME'] = os.path.basename(phx_path)

    eossa_array = build_eossa_array(header, data)
    write_eossa_fits(out_fits, header, eossa_array)
    print(f"Wrote EOSSA file: {out_fits}")
    return out_fits


# --------------------- CLI ---------------------
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Convert .phX photometry file to EOSSA v3.1.1 FITS (binary table)')
    parser.add_argument('input', help='Input .phX file (or directory)')
    parser.add_argument('-o', '--out', help='Output FITS file (or directory for many)', default=None)
    args = parser.parse_args()

    inp = args.input
    out = args.out
    if os.path.isdir(inp):
        # process all .ph* files
        files = [os.path.join(inp, f) for f in os.listdir(inp) if re.match(r'.*\.ph.', f)]
        out_dir = out or inp
        for f in files:
            base = os.path.splitext(os.path.basename(f))[0]
            out_path = os.path.join(out_dir, base + '.eossa.fits')
            convert_phx_to_eossa(f, out_path)
    else:
        out_path = out
        convert_phx_to_eossa(inp, out_path)
