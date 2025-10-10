import sys
import os
import numpy as np
import glob

from datetime import datetime

from astropy.io import fits
from astropy.time import Time # Зовнішній інструмент для роботи з JD та UTC [4]

import re
from typing import Dict, Any, Union, List, Tuple


# Константа для стандартної дальності нормалізації (1000 км)
DISTANCE_NORM_KM = 1000.0
# Placeholder для 32-bit Integer (J) [7]
PLACEHOLDER_INT = -2147483648


# --- 1. Визначення Метаданих EOSSA ---

# Визначення обов'язкових полів для наземного сенсора (згідно з Таблицею K-2 [5])
# TTYPE: назва колонки, TFORM: формат FITS (наприклад, D для Double, A для ASCII, J для 32-bit Integer, 2D для пари Double)
REQUIRED_EOSSA_COLUMNS_METADATA = [
    # UTC & Time
    {'name': 'UTC_Begin_exp', 'format': '23A', 'unit': 'yyyy-mm-ddThh:mm:ss.sss'}, # [5]
    {'name': 'UTC_End_exp', 'format': '23A', 'unit': 'yyyy-mm-ddThh:mm:ss.sss'}, # [5]
    {'name': 'JD_Mid_Exp', 'format': 'D', 'unit': 'days'}, # [5]
    {'name': 'Exp_Duration', 'format': 'D', 'unit': 'sec'}, # [5]

    # Photometry & Filters
    {'name': 'Cur_Spec_Filt_Num', 'format': 'J', 'unit': 'num'}, # [5]
    {'name': 'Cur_ND_Filt_Num', 'format': 'J', 'unit': 'num'}, # [6]
    {'name': 'Mag_Exo_Atm', 'format': 'D', 'unit': 'mag'}, # Позаатмосферна зоряна величина (mr) [6]
    {'name': 'Mag_Range_Norm', 'format': 'D', 'unit': 'mag'}, # Нормована за дальністю (Mr) [6]

    # Astrometry & Geometry (приклад скорочено, але включено важливі поля)
    {'name': 'Tel_Obj_Range', 'format': 'D', 'unit': 'm'}, # Відстань до цілі [7]
    {'name': 'Eph_RA_DE', 'format': '2D', 'unit': 'deg'}, # Predicted RA/Dec (Ефемериди) [6]
    {'name': 'Met_RA_DE', 'format': '2D', 'unit': 'deg'}, # Measured RA/Dec [6]
    {'name': 'Eph_Az_El', 'format': '2D', 'unit': 'deg'}, # Predicted RA/Dec (Ефемериди) [6]
    {'name': 'Met_Az_El', 'format': '2D', 'unit': 'deg'}, # Predicted RA/Dec (Ефемериди) [6]
]


def transform_data_to_stream(
        date_time: List[str],
        flux: List[float],
        flux_err: List[float],
        mag: List[float],
        mag_err: List[float],
        Az: List[float],
        El: List[float],
        Rg: List[float],
        dt_value: float = 0.1  # Параметр dt має бути відомий з шапки файлу
    ) -> List[Dict[str, Any]]:
    """
    Перетворює паралельні списки даних на список словників.

    :param dt_value: Значення 'dt', зчитане з шапки файлу.
    :return: Список словників з структурованими даними.
    """

    # 1. Створюємо список кортежів, де кожен кортеж - це один "рядок" даних.
    #    Ми використовуємо 'zip' для поєднання елементів з усіх списків.
    #    Важливо: dt ігнорується, тому що це константне значення з шапки.
    zipped_data = zip(date_time, mag, Rg, Az, El)

    result_stream = []

    for row_data in zipped_data:
        # row_data - це кортеж: (date_time_str, magV, Rg, Az, El)
        # 3. Формуємо словник для поточного рядка
        data_row = {
            'Date UT': row_data[0].strftime("%Y-%m-%d"),
            'Time Start': row_data[0].strftime("%H:%M:%S.%f"),
            'dt': dt_value,  # Константне значення з шапки файлу
            'magV': row_data[1],
            'Rg(Mm)': row_data[2],
            'Az(deg)': row_data[3],
            'El(deg)': row_data[4],
            "RA(deg)": 0.00,
            "DEC(deg)": 0.00
        }
        result_stream.append(data_row)

    return result_stream


def read_header_to_dict(filepath: str) -> Dict[str, Any]:
    """
    Зчитує дані із шапки файлу (до початку табличних даних) і повертає їх у вигляді словника.

    Спеціально обробляє:
    1. Коментарі (починаються з '#').
    2. Блок TLE (три наступні рядки після '# TLE:').
    3. Два наступні рядки після TLE як час початку і кінця ('START_TIME', 'END_TIME').
    4. Пари 'ключ = значення' (ігноруючи пробіли).

    :param filepath: Шлях до файлу.
    :return: Словник з назвами та їх значеннями.
    """
    header_data: Dict[str, Any] = {}

    # Стан для TLE
    in_tle_block = False
    tle_lines = []

    # Стан для часу
    time_lines_read = 0

    # Регулярний вираз для пошуку пар 'ключ = значення'
    key_value_pattern = re.compile(r'^\s*(\w+)\s*=\s*(.+)$')

    # Список для тимчасового зберігання рядків часу, незалежно від TLE
    potential_time_lines = []

    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            for line in f:
                # Видаляємо зайві пробіли з кінців рядка
                stripped_line = line.strip()

                # Рядок, з якого видалені символи коментаря
                cleaned_line = stripped_line.lstrip('#').strip()

                # --- КРОК 1: Обробка блоку TLE ---
                if in_tle_block:
                    # Додаємо повний рядок (включно з '#') до TLE
                    tle_lines.append(stripped_line)
                    if len(tle_lines) == 3:
                        header_data['TLE'] = tuple(tle_lines)
                        in_tle_block = False
                    continue

                if stripped_line == '# TLE:':
                    in_tle_block = True
                    continue

                # --- КРОК 2: Збір потенційних рядків часу ---
                # Збираємо перші два НЕ порожні рядки, які не є TLE, і не є "ключ = значення"
                if not in_tle_block and time_lines_read < 2:
                    if cleaned_line and not key_value_pattern.match(cleaned_line):
                        # Перевіряємо, чи це не заголовок стовпців
                        if len(cleaned_line.split()) > 4 and (' ' in cleaned_line or '\t' in cleaned_line):
                            # Якщо це схоже на рядок заголовків таблиці, припиняємо зчитування шапки
                            break

                        # Зберігаємо чистий рядок (без '#')
                        potential_time_lines.append(cleaned_line)
                        time_lines_read += 1

                        # Якщо ми знайшли два рядки, додаємо їх у словник
                        if time_lines_read == 2:
                            header_data['START_TIME'] = potential_time_lines[0]
                            header_data['END_TIME'] = potential_time_lines[1]
                        continue

                # --- КРОК 3: Обробка 'ключ = значення' ---
                if cleaned_line:
                    match = key_value_pattern.match(cleaned_line)
                    if match:
                        key = match.group(1).strip()
                        value = match.group(2).strip()

                        # Спроба конвертувати числові значення
                        try:
                            if re.fullmatch(r'[-+]?\d+$', value):
                                header_data[key] = int(value)
                            elif re.fullmatch(r'[-+]?\d*\.\d+|\d+', value):
                                header_data[key] = float(value)
                            else:
                                header_data[key] = value.replace('"', '').replace("'", '')
                        except ValueError:
                            header_data[key] = value.replace('"', '').replace("'", '')

                        continue

                # --- КРОК 4: Припинення зчитування шапки ---
                # Якщо рядок не є порожнім, коментарем, TLE, часом, або 'ключ = значення',
                # це, ймовірно, початок табличних даних.
                if stripped_line and not stripped_line.startswith('#'):
                    break

    except FileNotFoundError:
        print(f"Помилка: Файл не знайдено за шляхом '{filepath}'")
        return {}

    return header_data


def read_data(fname):
    date_time = []
    date, time = np.genfromtxt(fname, unpack=True, skip_header=True, usecols=(0, 1), dtype=None, encoding="utf-8")
    x, y, xerr, yerr, flux, flux_err, mag, mag_err, Az, El, Rg = \
        np.genfromtxt(fname, skip_header=True,
                      usecols=(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,),
                      unpack=True)
    fit_file = np.genfromtxt(fname, unpack=True, skip_header=True, usecols=(-1, ), dtype=None, encoding="utf-8")

    for i in range(0, len(date)):
        date_time.append(datetime.strptime(date[i] + ' ' + time[i] + "000", "%Y-%m-%d %H:%M:%S.%f"))

    return date_time, x, y, xerr, yerr, flux, flux_err, mag, mag_err, Az, El, Rg, fit_file


# # --- ФУНКЦІЯ 1: Створення об'єктів Column ---
#
# def create_eossa_columns(columns_metadata: list) -> list:
#     """
#     Створює список об'єктів fits.Column на основі визначених метаданих EOSSA.
#     Ця функція відповідає за визначення структури таблиці FITS (TTYPEn, TFORMn, TUNITn).
#     """
#     column_objects = []
#     for col_def in columns_metadata:
#         # Створення об'єкта fits.Column згідно з вимогами FITS/EOSSA [1, 2]
#         column = fits.Column(
#             name=col_def['name'],  # TTYPEn
#             format=col_def['format'],  # TFORMn
#             unit=col_def.get('unit', '')  # TUNITn (рекомендовано) [1]
#         )
#         column_objects.append(column)
#
#     return column_objects
#
#
# # --- ФУНКЦІЯ 2: Заповнення даними та створення таблиці ---
#
# def populate_eossa_data_table(column_definitions: list, metadata: dict, observations: list) -> fits.BinTableHDU:
#     """
#     Заповнює створені об'єкти Column даними та створює об'єкт BinTableHDU.
#
#     Args:
#         column_definitions: Список fits.Column об'єктів (створених ФУНКЦІЄЮ 1).
#         metadata: Загальні метадані (наприклад, TLE, назви, час).
#         observations: Список фактичних рядків спостережень.
#
#     Returns:
#         Об'єкт astropy.io.fits.BinTableHDU.
#     """
#     num_rows = len(observations)
#     data_columns = {}
#
#     # 1. Ініціалізація порожніх масивів для кожного стовпця
#     for col_def in column_definitions:
#         # Визначення типу numpy на основі формату FITS TFORM [8]
#         fmt = col_def.format
#         if 'A' in fmt:
#             dtype = f'U{int(fmt.replace("A", "")) if fmt.replace("A", "").isdigit() else 23}'  # Приклад для A
#         elif 'D' in fmt:
#             dtype = np.float64
#         elif 'J' in fmt:
#             dtype = np.int32
#         else:
#             dtype = np.float64  # За замовчуванням
#
#         # Обробка векторних полів (наприклад, 2D)
#         count = int(fmt.replace('D', '')) if fmt.replace('D', '').isdigit() and 'D' in fmt else 1
#
#         # Якщо це вектор, створюємо масив з відповідними розмірами
#         if count > 1:
#             data_columns[col_def.name] = np.empty((num_rows, count), dtype=dtype)
#         else:
#             data_columns[col_def.name] = np.empty(num_rows, dtype=dtype)
#
#     # 2. Заповнення масивів даними
#     # Тут використовуємо дані з нашої історії розмови як приклад [9]
#
#     # ПРИКЛАД ВХІДНИХ ДАНИХ (для однієї точки спостереження)
#     obs_data = observations
#
#     # Перетворення часу
#     print(obs_data)
#     t_start = Time(obs_data['Date UT'] + obs_data['Time Start'], format='iso', scale='utc')
#     t_end = Time(obs_data['Date UT'] + obs_data['Time End'], format='iso', scale='utc')
#     t_mid = t_start + (t_end - t_start) / 2
#
#     # Заповнення колонками
#
#     # 1. UTC & Time
#     data_columns['UTC_Begin_exp'] = t_start.iso.replace(' ', 'T')[:23]
#     data_columns['UTC_End_exp'] = t_end.iso.replace(' ', 'T')[:23]
#     data_columns['JD_Mid_Exp'] = t_mid.jd
#     data_columns['Exp_Duration'] = obs_data['Exp_Duration']
#
#     # 2. Photometry & Filters (використовуємо приклади з попередньої розмови)
#     # Змінні Mag_Exo_Atm та Mag_Range_Norm:
#     mag_exo_atm = 6.602  # Припустиме значення після денормалізації/обчислень
#     mag_range_norm = obs_data['magV']  # 6.587 (якщо magV вважається вже нормованим) [9]
#
#     data_columns['Mag_Exo_Atm'] = mag_exo_atm
#     data_columns['Mag_Range_Norm'] = mag_range_norm
#
#     # Використовуємо заглушки або константи для відсутніх/фіксованих полів
#     data_columns['Cur_Spec_Filt_Num'] = 1  # Зазвичай 1, якщо не мультиспектральні [5]
#     data_columns['Cur_ND_Filt_Num'] = -2147483648  # Placeholder для 32-bit integer (J) [6, 8]
#
#     # 3. Geometry (використовуємо дані з прикладу)
#     # Tel_Obj_Range (1231.248 Mm = 1231248000 м)
#     range_m = obs_data['Rg(Mm)'] * 10 ** 6
#     data_columns['Tel_Obj_Range'] = range_m
#
#     # RA/DE, AZ/EL - для прикладу використовуємо нульові значення або дані з вхідного файлу
#     # Вхідні дані містять Az(deg) і El(deg), але не RA/Dec.
#     # Для демонстрації заповнення 2D-векторів:
#     # Припускаємо, що Met_AZ_EL — це (Az, El) [10].
#     data_columns['Met_RA_DE'] = [obs_data['RA'], obs_data['DE']]  # Якщо ці значення були обчислені
#     data_columns['Eph_RA_DE'] = [obs_data['RA_Eph'], obs_data['DE_Eph']]  # Якщо ці значення були обчислені
#
#     data_columns['Met_AZ_EL'] = [obs_data['Az(deg)'], obs_data['El(deg)']]  # [9]
#
#     # 3. Створення об'єкта таблиці FITS
#     c = fits.ColDefs(column_definitions)
#     hdu = fits.BinTableHDU.from_columns(c, nrows=num_rows)
#
#     return hdu


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


# Визначення dtype для NumPy/FITS_rec (використовується для створення структури)
# Це допомагає гарантувати, що кожен рядок має правильний формат перед об'єднанням.
def get_numpy_dtype_from_metadata(columns_metadata):
    """Створює структурований dtype для NumPy масиву на основі метаданих FITS."""
    np_dtypes = []
    for col_def in columns_metadata:
        name = col_def['name']
        fmt = col_def['format']

        # Визначення типу та розміру
        if 'A' in fmt:
            size = int(fmt.replace('A', '')) if fmt.replace('A', '').isdigit() else 23  # UTC size
            dtype = f'U{size}'
            count = 1
        elif 'D' in fmt:
            dtype = np.float64
            count = int(fmt.replace('D', '')) if fmt.replace('D', '').isdigit() else 1
        elif 'J' in fmt:
            dtype = np.int32
            count = 1
        else:
            dtype = np.float64
            count = 1

        if count > 1:
            np_dtypes.append((name, dtype, count))
        else:
            np_dtypes.append((name, dtype))
    return np.dtype(np_dtypes)


# Створення структури NumPy DTYPE
EOSSA_DTYPE = get_numpy_dtype_from_metadata(REQUIRED_EOSSA_COLUMNS_METADATA)


# --- ФУНКЦІЯ 1: Створення об'єктів Column ---
# Ця функція визначає лише заголовок/структуру, але не дані.
def create_eossa_columns(columns_metadata: list) -> fits.ColDefs:
    """Створює об'єкт fits.ColDefs, який описує структуру таблиці."""
    column_objects = []
    for col_def in columns_metadata:
        column = fits.Column(
            name=col_def['name'],
            format=col_def['format'],
            unit=col_def.get('unit', '')
        )
        column_objects.append(column)

    return fits.ColDefs(column_objects)


# --- ФУНКЦІЯ 2: Обробка ОДНОГО РЯДКА ---
# На вхід – один словник (рядок), на виході – один структурований елемент даних.

def process_single_observation(obs_data: dict) -> np.void:
    """
    Обробляє один рядок вхідних даних, виконує необхідні розрахунки
    (наприклад, денормалізацію) і повертає його у вигляді структурованого елемента.
    """
    # Створення пустого структурованого масиву на 1 рядок
    single_row = np.empty(1, dtype=EOSSA_DTYPE)

    # 1. Обчислення часу та тривалості
    # (Припускаємо, що в obs_data є 'Date UT' та інші поля часу)
    # print(obs_data['Date UT'] + obs_data['Time Start'])
    t_start = Time(obs_data['Date UT'] + " " + obs_data['Time Start'], format='iso', scale='utc', precision=6)
    t_end = Time(obs_data['Date UT'] + " " + obs_data['Time Start'], format='iso', scale='utc',
                 precision=6)  # Використовуємо початок, оскільки Exp_Duration = 0.1s
    t_end += obs_data['dt'] / 24.0 / 3600.0  # Припускаємо dt = 0.1s
    t_mid = t_start + (t_end - t_start) / 2

    single_row['UTC_Begin_exp'] = t_start.iso.replace(' ', 'T')[:23]
    single_row['UTC_End_exp'] = t_end.iso.replace(' ', 'T')[:23]
    single_row['JD_Mid_Exp'] = t_mid.jd
    single_row['Exp_Duration'] = obs_data['dt']  # 0.1 секунди, якщо dt = 0.1

    # 2. Фотометрія та нормалізація дальності

    # Mag_Range_Norm (Mr) - це вхідне значення magV (згідно з вашим описом)
    M_r = obs_data['magV']

    # Tel_Obj_Range (d) - відстань у Мегаметрах (Mm). Переводимо в метри.
    d_km = obs_data['Rg(Mm)'] * 1000.0

    # Денормалізація: m_r = M_r + 5 * log10(d / D)
    # Mag_Exo_Atm (mr) [9, 10]
    m_r = M_r + 5.0 * np.log10(d_km / DISTANCE_NORM_KM)

    single_row['Mag_Range_Norm'] = M_r
    single_row['Mag_Exo_Atm'] = m_r

    # 3. Геометрія
    single_row['Tel_Obj_Range'] = d_km * 1000.0  # м
    single_row['Cur_Spec_Filt_Num'] = 1
    single_row['Cur_ND_Filt_Num'] = PLACEHOLDER_INT  # Немає нейтрального фільтра [7]

    # Векторні поля (2D)
    # Eph_RA_DE, Met_RA_DE (Припускаємо, що RA/DE обчислюються)
    print(single_row)
    print(single_row.keys())
    single_row['Eph_RA_DE'] = np.array([0.0, 0.0])  # Заглушка, оскільки ці дані відсутні у вихідному файлі
    single_row['Met_RA_DE'] = np.array([0.0, 0.0])

    # Eph_AZ_EL, Met_AZ_EL (Вхідні дані містять Az(deg) і El(deg))
    single_row['Eph_AZ_EL'] = np.array([0.0, 0.0])  # Ефемеридні Азимут/Елевація
    single_row['Met_AZ_EL'] = np.array([obs_data['Az(deg)'], obs_data['El(deg)']])  # Виміряні Азимут/Елевація [8]

    # Solar_Phase_Ang (Припускаємо, що обчислюється)
    single_row['Solar_Phase_Ang'] = -9999.0  # Placeholder для D [7]

    return single_row  # Повертаємо єдиний структурований елемент


# --- ГОЛОВНА ЛОГІКА СКРИПТУ: ЗБІРКА ДАНИХ ТА ФІНАЛІЗАЦІЯ ---

def finalize_eossa_table(collected_rows: list, column_definitions: fits.ColDefs) -> fits.BinTableHDU:
    """
    Об'єднує всі зібрані рядки та створює фінальний об'єкт FITS BinTableHDU.
    """
    if not collected_rows:
        return None, 0

    # Об'єднання всіх окремих структурованих елементів у єдиний NumPy масив
    final_data_array = np.stack(collected_rows)

    # Створення об'єкта FITS HDU
    # astropy автоматично визначає NAXIS2 (кількість рядків)
    hdu = fits.BinTableHDU.from_columns(column_definitions, data=final_data_array)

    # Кількість рядків (data_len) буде автоматично записана в заголовок як NAXIS2
    data_len = len(final_data_array)

    # Оновлення заголовка (додавання обов'язкових ключових слів, які залежать від кількості полів)
    hdu.header['TFIELDS'] = len(column_definitions)  # TFIELDS (Number of fields) [11]
    hdu.header['NAXIS2'] = data_len

    return hdu, data_len


if __name__ == "__main__":
    filename = sys.argv[1]

    if os.path.isfile(filename):
        filelist = [filename]
    else:
        # we got path with mask
        filelist = glob.glob(filename)

    for filename in filelist:
        wd = os.path.dirname(filename)
        base = os.path.basename(filename)
        fname, ext = os.path.splitext(base)

        header = read_header_to_dict(filename)
        ymd = header['START_TIME'][11].strip('-')
        ut1 = (header['START_TIME'][11:])
        ut1 = ut1[:-4].strip(':')
        band = ext[2:]
        eossa_path = os.path.join(wd,
                                  "result_" + str(header['NORAD']) + "_" + ymd + "_UT" + ut1 + "_" + band + ".eossa")


        date_time, _, _, _, _, flux, flux_err, mag, mag_err, Az, El, Rg, _ = read_data(filename)

        # obs_data = list(zip(date_time, flux, flux_err, mag, mag_err, Az, El, Rg))
        obs_data = transform_data_to_stream(
            date_time,
            flux,
            flux_err,
            mag,
            mag_err,
            Az,
            El,
            Rg,
            dt_value=header['dt']  # Припустимо, що dt = 0.1
        )

        #TODO: calc (RA, DEC) to each row
        # print(obs_data[:2])

        # # Крок 1: Створення визначень колонок (Функція 1)
        # eossa_column_definitions = create_eossa_columns(REQUIRED_EOSSA_COLUMNS_METADATA)
        #
        # # Крок 2: Заповнення таблиці даними (Функція 2)
        # # Ми передаємо визначення колонок та фактичні дані для заповнення
        # binary_table_hdu = populate_eossa_data_table(eossa_column_definitions, {}, obs_data)
        #
        # primary_hdu = create_primary_hdu_header()
        # # 4. Об'єднання та запис у файл
        # hdul = fits.HDUList([primary_hdu, binary_table_hdu])
        #
        # # Перевірка на існування та перезапис
        # if os.path.exists(eossa_path):
        #     os.remove(eossa_path)
        #
        # try:
        #     hdul.writeto(eossa_path, overwrite=True)
        #     hdul.close()
        # except Exception as e:
        #     print(f"Помилка при записі FITS-файлу: {e}")




        # 1. Визначення структури FITS (Крок 1)
        eossa_column_definitions = create_eossa_columns(REQUIRED_EOSSA_COLUMNS_METADATA)

        # 2. Ініціалізація списку для збору даних
        collected_rows = []

        print("--- Обробка потоку даних по одному рядку ---")

        for data_row in obs_data:
            # Крок 2: Обробка одного рядка (виконується в процесі)
            processed_entry = process_single_observation(data_row)
            collected_rows.append(processed_entry)
            print(f"Оброблено запис JD: {processed_entry['JD_Mid_Exp']:.5f}")

        # 4. Фіналізація таблиці (виконується в кінці скрипту)
        final_hdu, data_len = finalize_eossa_table(collected_rows, eossa_column_definitions)

        if final_hdu:
            print("\n--- ФІНАЛІЗАЦІЯ ТАБЛИЦІ EOSSA ---")
            print(f"Загальна кількість рядків (NAXIS2): {data_len}")

            # Перевірка заповнення обов'язкового поля NAXIS2
            print(f"Значення NAXIS2 у заголовку: {final_hdu.header.get('NAXIS2')}")

            # Виведення перших 5 значень Mag_Exo_Atm для демонстрації
            print("\nПерші 5 значень Mag_Exo_Atm (позаатмосферна, денормована):")
            # Примітка: Mag_Exo_Atm має бути близькою до Mag_Range_Norm (magV), але трохи більшою
            # оскільки d (1231 км) > D (1000 км).
            print(final_hdu.data['Mag_Exo_Atm'][:5])

            # primary_hdu = create_primary_hdu_header()
            # 4. Об'єднання та запис у файл
            hdul = fits.HDUList(final_hdu)

            # Перевірка на існування та перезапис
            if os.path.exists(eossa_path):
                os.remove(eossa_path)

            try:
                hdul.writeto(eossa_path, overwrite=True)
                hdul.close()
            except Exception as e:
                print(f"Помилка при записі FITS-файлу: {e}")

