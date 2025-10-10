from astropy.io import fits
import numpy as np
from astropy.time import Time # Зовнішній інструмент для роботи з JD та UTC [4]

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
]

# --- ФУНКЦІЯ 1: Створення об'єктів Column ---

def create_eossa_columns(columns_metadata: list) -> list:
    """
    Створює список об'єктів fits.Column на основі визначених метаданих EOSSA.
    Ця функція відповідає за визначення структури таблиці FITS (TTYPEn, TFORMn, TUNITn).
    """
    column_objects = []
    for col_def in columns_metadata:
        # Створення об'єкта fits.Column згідно з вимогами FITS/EOSSA [1, 2]
        column = fits.Column(
            name=col_def['name'], # TTYPEn
            format=col_def['format'], # TFORMn
            unit=col_def.get('unit', '') # TUNITn (рекомендовано) [1]
        )
        column_objects.append(column)
        
    return column_objects

# --- ФУНКЦІЯ 2: Заповнення даними та створення таблиці ---

def populate_eossa_data_table(column_definitions: list, metadata: dict, observations: list) -> fits.BinTableHDU:
    """
    Заповнює створені об'єкти Column даними та створює об'єкт BinTableHDU.

    Args:
        column_definitions: Список fits.Column об'єктів (створених ФУНКЦІЄЮ 1).
        metadata: Загальні метадані (наприклад, TLE, назви, час).
        observations: Список фактичних рядків спостережень.

    Returns:
        Об'єкт astropy.io.fits.BinTableHDU.
    """
    num_rows = len(observations)
    data_columns = {}
    
    # 1. Ініціалізація порожніх масивів для кожного стовпця
    for col_def in column_definitions:
        # Визначення типу numpy на основі формату FITS TFORM [8]
        fmt = col_def.format
        if 'A' in fmt:
            dtype = f'U{int(fmt.replace("A", "")) if fmt.replace("A", "").isdigit() else 23}' # Приклад для A
        elif 'D' in fmt:
            dtype = np.float64
        elif 'J' in fmt:
            dtype = np.int32
        else:
            dtype = np.float64 # За замовчуванням
            
        # Обробка векторних полів (наприклад, 2D)
        count = int(fmt.replace('D', '')) if fmt.replace('D', '').isdigit() and 'D' in fmt else 1
        
        # Якщо це вектор, створюємо масив з відповідними розмірами
        if count > 1:
            data_columns[col_def.name] = np.empty((num_rows, count), dtype=dtype)
        else:
            data_columns[col_def.name] = np.empty(num_rows, dtype=dtype)


    # 2. Заповнення масивів даними
    # Тут використовуємо дані з нашої історії розмови як приклад [9]
    
    # ПРИКЛАД ВХІДНИХ ДАНИХ (для однієї точки спостереження)
    obs_data = observations
    
    # Перетворення часу
    t_start = Time(obs_data['Date UT'] + obs_data['Time Start'], format='iso', scale='utc')
    t_end = Time(obs_data['Date UT'] + obs_data['Time End'], format='iso', scale='utc')
    t_mid = t_start + (t_end - t_start) / 2
    
    # Заповнення колонками
    
    # 1. UTC & Time
    data_columns['UTC_Begin_exp'] = t_start.iso.replace(' ', 'T')[:23]
    data_columns['UTC_End_exp'] = t_end.iso.replace(' ', 'T')[:23]
    data_columns['JD_Mid_Exp'] = t_mid.jd
    data_columns['Exp_Duration'] = obs_data['Exp_Duration']

    # 2. Photometry & Filters (використовуємо приклади з попередньої розмови)
    # Змінні Mag_Exo_Atm та Mag_Range_Norm:
    mag_exo_atm = 6.602 # Припустиме значення після денормалізації/обчислень
    mag_range_norm = obs_data['magV'] # 6.587 (якщо magV вважається вже нормованим) [9]
    
    data_columns['Mag_Exo_Atm'] = mag_exo_atm
    data_columns['Mag_Range_Norm'] = mag_range_norm

    # Використовуємо заглушки або константи для відсутніх/фіксованих полів
    data_columns['Cur_Spec_Filt_Num'] = 1 # Зазвичай 1, якщо не мультиспектральні [5]
    data_columns['Cur_ND_Filt_Num'] = -2147483648 # Placeholder для 32-bit integer (J) [6, 8]

    # 3. Geometry (використовуємо дані з прикладу)
    # Tel_Obj_Range (1231.248 Mm = 1231248000 м)
    range_m = obs_data['Rg(Mm)'] * 10**6 
    data_columns['Tel_Obj_Range'] = range_m
    
    # RA/DE, AZ/EL - для прикладу використовуємо нульові значення або дані з вхідного файлу
    # Вхідні дані містять Az(deg) і El(deg), але не RA/Dec.
    # Для демонстрації заповнення 2D-векторів:
    # Припускаємо, що Met_AZ_EL — це (Az, El) [10].
    data_columns['Met_RA_DE'] = [obs_data['RA'], obs_data['DE']] # Якщо ці значення були обчислені
    data_columns['Eph_RA_DE'] = [obs_data['RA_Eph'], obs_data['DE_Eph']] # Якщо ці значення були обчислені
    
    data_columns['Met_AZ_EL'] = [obs_data['Az(deg)'], obs_data['El(deg)']] # [9]

    
    # 3. Створення об'єкта таблиці FITS
    c = fits.ColDefs(column_definitions)
    hdu = fits.BinTableHDU.from_columns(c, nrows=num_rows)
    
    return hdu

# --- Виконання прикладу ---

# 1. Дані (імітація одного рядка з джерела [9])
EXAMPLE_OBSERVATION_DATA = [{
    'Date UT': '2025-09-12 ',
    'Time Start': '19:23:12.404',
    'Time End': '19:28:47.102',
    'Exp_Duration': 334.698, # Різниця часу початку/кінця
    'magV': 6.587,
    'Rg(Mm)': 1231.248, # [9]
    'Az(deg)': 211.950, 
    'El(deg)': 39.089,
    'RA': 75.31, # Уявні виміряні RA
    'DE': -5.201, # Уявні виміряні Dec
    'RA_Eph': 75.33, # Уявні ефемеридні RA
    'DE_Eph': -5.001, # Уявні ефемеридні Dec
}]


# Крок 1: Створення визначень колонок (Функція 1)
eossa_column_definitions = create_eossa_columns(REQUIRED_EOSSA_COLUMNS_METADATA)

# Крок 2: Заповнення таблиці даними (Функція 2)
# Ми передаємо визначення колонок та фактичні дані для заповнення
binary_table_hdu = populate_eossa_data_table(eossa_column_definitions, {}, EXAMPLE_OBSERVATION_DATA)

# Виведення результату
print("--- Створена FITS Таблиця (HDU) ---")
print(binary_table_hdu.header)
print("\n--- Перевірка даних Mag_Exo_Atm ---")
print(binary_table_hdu.data['Mag_Exo_Atm'])