import numpy as np
import os
from astropy.io import fits
from astropy.time import Time
from astropy import units as u
from astropy import constants as const

# --- 1. КОНСТАНТИ EOSSA ТА ЗАГЛУШКИ ---

# Версія формату EOSSA
EOSSA_VERSION = '3.1.1'
# Значення-заглушки відповідно до вимог EOSSA [2]
D_NULL = -9999.0  # Double precision floating point (8 bytes) [2]
J_NULL = -2147483648  # 32-bit integer (4 bytes) [2]
A_NULL = 'NULLSTRING' # Character [2]

# Координати спостереження та дані цілі з джерела [6]
RAW_METADATA = {
    "TLE_LINE1_FULL": "1 23233U 94057A 25254.88090343 .00000005 00000-0 24291-4 0 9998",
    "TLE_LINE2_FULL": "2 23233 98.7142 306.8771 0010784 355.9423 4.1662 14.16019207602580",
    "NORAD": 23233,
    "OBJECT": "DMSP 5D-2 F12 (USA 106)",
    "SITE_NAME": "Derenivka",
    "SITE_LAT": 48.5635505,      # degrees N
    "SITE_LON": 22.453751,       # degrees E
    "SITE_ELEV": 231.1325,       # meters
    # Часові мітки для загального інтервалу спостереження [6]
    "OBS_UTC_START": "2025-09-12T19:23:12.404333",
    "OBS_UTC_END": "2025-09-12T19:28:47.102120",
    # Для цього прикладу використовуємо одиночний рядок даних [6]
    "DATA_ROW_TIME": "2025-09-12 19:23:12.404", # UTC час першого виміру
    "EXP_DURATION_SEC": 1.0, # Припускаємо тривалість експозиції 1.0 сек для цього кадру
    "magV": 6.587,
    "mag_err": 0.006,
    "Az_deg": 211.950,
    "El_deg": 39.089,
    "Range_Mm": 1231.248,
    "SPECTRAL_FILTER": 'V', # Припускаємо використання V-фільтра
}

# --- 2. ФУНКЦІЇ ОБРОБКИ ТА РОЗРАХУНКУ ---

def calculate_derived_fields(metadata):
    """Виконує необхідні астрономічні та фізичні розрахунки."""
    
    # Визначення часу початку та кінця експозиції для цього рядка даних
    t_start = Time(metadata["DATA_ROW_TIME"], format='iso', scale='utc')
    t_end = t_start + metadata["EXP_DURATION_SEC"] * u.second
    
    # 1. JD_Mid_Exp: Юліанська дата середнього часу експозиції [7]
    t_mid = t_start + (metadata["EXP_DURATION_SEC"] / 2.0) * u.second
    jd_mid_exp = t_mid.jd
    
    # 2. Tel_Obj_Range: Дальність у метрах [8]
    range_mm = metadata["Range_Mm"]
    range_m = range_mm * 1.0e6 

    # 3. Mag_Range_Norm: Нормована за дальністю магнітуда (еталон 1000 км) [9, 10]
    # d повинна бути в km, D = 1000 km.
    range_km = range_mm * 1000.0
    mag_exo_atm = metadata["magV"] # Приймаємо magV як Mag_Exo_Atm [11]
    # Mag_Exo_Atm = magV + 5 log_{10} (Tel_Obj_Range/1000)
    
    # M_r = m_r - 5 * log10(d / D) [10]
    if range_km > 0 and 1000.0 > 0:
        mag_range_norm = mag_exo_atm - 5.0 * np.log10(range_km / 1000.0)
    else:
        mag_range_norm = D_NULL

    # 4. Обчислення тривалості (для колонки)
    exp_duration = metadata["EXP_DURATION_SEC"]

    return {
        'UTC_Begin_Exp': t_start.isot,
        'UTC_End_Exp': t_end.isot,
        'JD_Mid_Exp': jd_mid_exp,
        'Exp_Duration': exp_duration,
        'Mag_Exo_Atm': mag_exo_atm,
        'Mag_Exo_Atm_Unc': metadata["mag_err"], # Пріоритет 28 [12]
        'Mag_Range_Norm': mag_range_norm,
        'Met_AZ_EL': np.array([metadata["Az_deg"], metadata["El_deg"]], dtype=np.float64), # 2D [8]
        'Tel_Obj_Range': range_m
    }

# --- 3. ВИЗНАЧЕННЯ СТРУКТУРИ ТАБЛИЦІ (FITS Columns) ---

def define_fits_columns(num_data_rows, derived_data):
    """
    Створює об'єкти fits.Column на основі обов'язкових полів EOSSA 
    для наземного сенсора [5, 8] та деяких важливих опціональних полів.
    """
    
    # Використовуємо astropy.io.fits.Column
    cols = [
        # REQUIRED FIELDS (ALL)
        fits.Column(name='UTC_Begin_Exp', format='24A', unit='yyyy-mm-ddThh:m:ss{.sssssss}', array=[derived_data['UTC_Begin_Exp']]), # 1 [5]
        fits.Column(name='UTC_End_Exp', format='24A', unit='yyyy-mm-ddThh:m:ss{.sssssss}', array=[derived_data['UTC_End_Exp']]),     # 2 [5]
        fits.Column(name='JD_Mid_Exp', format='1D', unit='days', array=[derived_data['JD_Mid_Exp']]),                                  # 3 [5]
        fits.Column(name='Exp_Duration', format='1D', unit='sec', array=[derived_data['Exp_Duration']]),                               # 4 [5]
        fits.Column(name='Cur_Spec_Filt_Num', format='1J', unit='num', array=[1]),                                                    # 5 (1=V-filter) [5]
        fits.Column(name='Mag_Exo_Atm', format='1D', unit='mag', array=[derived_data['Mag_Exo_Atm']]),                                 # 7 [5]
        fits.Column(name='Mag_Range_Norm', format='1D', unit='mag', array=[derived_data['Mag_Range_Norm']]),                            # 8 [5]
        fits.Column(name='Eph_RA_DE', format='2D', unit='deg', array=[np.array([D_NULL, D_NULL], dtype=np.float64)]),                   # 9 (Predicted RA/Dec) [5]
        fits.Column(name='Met_RA_DE', format='2D', unit='deg', array=[np.array([D_NULL, D_NULL], dtype=np.float64)]),                   # 10 (Measured RA/Dec - Невідомі, використовуємо заглушки) [5]
        fits.Column(name='Tel_Obj_Range', format='1D', unit='m', array=[derived_data['Tel_Obj_Range']]),                               # 14 [8]

        # REQUIRED FIELDS (GROUND-BASED - G)
        fits.Column(name='Eph_AZ_EL', format='2D', unit='deg', array=[np.array([D_NULL, D_NULL], dtype=np.float64)]),                   # 11 [8]
        fits.Column(name='Met_AZ_EL', format='2D', unit='deg', array=[derived_data['Met_AZ_EL']]),                                     # 12 (Measured Az/El) [8]
        fits.Column(name='Sun_AZ_EL', format='2D', unit='deg', array=[np.array([D_NULL, D_NULL], dtype=np.float64)]),                   # 13 [8]
        
        # OPTIONAL FIELDS (Prioritized 6 & 28)
        fits.Column(name='Cur_ND_Filt_Num', format='1J', unit='num', array=[J_NULL]),                                                 # 6 (ND Filter - заглушка) [13]
        fits.Column(name='UTC_Unc', format='1D', unit='sec', array=[1.0e-7]),                                                         # Uncertainty in time, Priority 6 [14]
        fits.Column(name='Mag_Exo_Atm_Unc', format='1D', unit='mag', array=[derived_data['Mag_Exo_Atm_Unc']]),                        # Uncertainty in magnitude, Priority 28 [12]
    ]
    
    # Перетворення списку fits.Column у ColDefs
    return fits.ColDefs(cols)

# --- 4. ФОРМУВАННЯ ЗАГОЛОВКІВ ---

def create_primary_hdu_header():
    """Створює Primary Header (HDU0)."""
    hdr = fits.Header()
    # Обов'язкові ключові слова Primary Header [15, 16]
    hdr['SIMPLE'] = (True, 'file does conform to FITS standard')
    hdr['BITPIX'] = (8, 'number of bits per data pixel')
    hdr['NAXIS'] = (0, 'number of data axes')
    hdr['EXTEND'] = (True, 'FITS dataset may contain extensions')
    hdr['CLASSIF'] = ('UNCLASS', 'Security classification level') [15]
    # END додається автоматично
    return fits.PrimaryHDU(header=hdr)

def create_binary_table_header(metadata, col_defs):
    """Створює Binary Table Extension Header (HDU1) з ключовими словами EOSSA."""
    
    # Визначення обов'язкових метаданих для наземного сенсора [4, 8]
    
    # 1. TLE: обрізаємо '1 ' та '2 ' на початку (повинно бути 67 символів) [17, 18]
    tleln1 = metadata['TLE_LINE1_FULL'][2:].strip()
    tleln2 = metadata['TLE_LINE2_FULL'][2:].strip()

    hdr = fits.Header()
    
    # Загальні ключові слова FITS
    hdr['XTENSION'] = ('BINTABLE', 'Binary table extension') [3]
    hdr['TFIELDS'] = (len(col_defs), 'number of fields in each row') [3]
    
    # Ключові слова EOSSA (Required for ALL / Ground-based - G)
    hdr['EXTNAME'] = (f"{metadata['NORAD']}_{metadata['SITE_NAME']}.eossa", 'Filename of the FITS binary extension table file.') [3]
    hdr['CLASSIF'] = ('UNCLASS', 'Security classification level of the data contained in the file.') [3]
    hdr['VERS'] = (EOSSA_VERSION, 'Version number of the EOSSA data format.') [3]
    hdr['OBSEPH'] = ('GROUND', 'Observer type, i.e., ground-based') [3]
    hdr['TELESCOP'] = (metadata['SITE_NAME'], 'Telescope site name or SSN sensor identifier.') [19]
    
    # Координати сенсора (Required G)
    hdr['TELLAT'] = (metadata['SITE_LAT'], 'Geographical latitude of the telescope (degrees N).') [19]
    hdr['TELLONG'] = (metadata['SITE_LON'], 'Geographical longitude of the telescope (degrees E).') [19]
    hdr['TELALT'] = (metadata['SITE_ELEV'], 'Distance above sea level of the telescope (m).') [19]
    hdr['OBSNAME'] = (metadata['SITE_NAME'], 'Telescope name.') [20]
    
    # Дані цілі (Required ALL)
    hdr['OBJEPH'] = ('TLE', 'Target object ephemeris source.') [17]
    hdr['OBJTYPE'] = ('SCN', 'If the target has a space catalog number, the string SCN is the value.') [18]
    hdr['OBJNUM'] = (metadata['NORAD'], 'Identification number of the object from the catalog in OBJTYPE.') [21]
    hdr['OBJECT'] = (metadata['OBJECT'], 'Common name of target object.') [17]
    hdr['TLELN1'] = (tleln1, 'Target Object Truncated TLE Line 1 (67 characters).') [17]
    hdr['TLELN2'] = (tleln2, 'Target Object Truncated TLE Line 2 (67 characters).') [18]
    
    # Фільтри (Required ALL)
    hdr['SPFNUM'] = (1, 'The value is the number of spectral filters used.') [18]
    hdr['SPFNAM1'] = (metadata['SPECTRAL_FILTER'], f"Spectral filter name (n=1).") [18]
    
    # Додаткові бажані поля (Priority 2) [22]
    hdr['PHOTTYP'] = ('Metric', 'Data product type.') [22]
    hdr['BASING'] = ('Ground', 'The type of sensor basing is identified here.') [22]
    
    # Створення HDU двійкової таблиці
    return fits.BinTableHDU.from_columns(col_defs, header=hdr)

# --- 5. ОСНОВНА ФУНКЦІЯ СКРИПТА ---

def convert_to_eossa_with_astropy(input_data_meta, output_filepath="output_eossa.fits"):
    
    # 1. Розрахунок похідних даних
    derived_data = calculate_derived_fields(input_data_meta)
    
    # 2. Визначення колонок FITS
    # Оскільки ми обробляємо один рядок даних [6], кількість рядків = 1
    col_defs = define_fits_columns(1, derived_data)
    
    # 3. Формування HDU
    primary_hdu = create_primary_hdu_header()
    binary_hdu = create_binary_table_header(input_data_meta, col_defs)

    # 4. Об'єднання та запис у файл
    hdul = fits.HDUList([primary_hdu, binary_hdu])
    
    # Перевірка на існування та перезапис
    if os.path.exists(output_filepath):
        os.remove(output_filepath)
        
    hdul.writeto(output_filepath, overwrite=True)
    hdul.close()

    print(f"--- Успішно створено файл EOSSA: {output_filepath} ---")
    print(f"Формат: FITS Binary Table (v{EOSSA_VERSION})")
    print("\nКлючові вимірювання:")
    print(f"JD_Mid_Exp: {derived_data['JD_Mid_Exp']:.8f} days [5]")
    print(f"Mag_Exo_Atm (V): {derived_data['Mag_Exo_Atm']:.3f} mag [5]")
    print(f"Mag_Range_Norm (1000km): {derived_data['Mag_Range_Norm']:.3f} mag [5]")
    print(f"Met_AZ_EL: Az={derived_data['Met_AZ_EL']:.3f} deg, El={derived_data['Met_AZ_EL'][1]:.3f} deg [8]")

# --- ВИКОНАННЯ ---
convert_to_eossa_with_astropy(RAW_METADATA)

# Висновок:
# Для перевірки створеного файлу FITS можна використати:
# fits.open('output_eossa.fits')[1].header
# fits.open('output_eossa.fits')[1].data