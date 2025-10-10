from astropy.io import fits
import numpy as np
from astropy.time import Time
import logging

# Константа для стандартної дальності нормалізації (1000 км)
DISTANCE_NORM_KM = 1000.0
# Placeholder для 32-bit Integer (J) [7]
PLACEHOLDER_INT = -2147483648

# --- 1. ВИЗНАЧЕННЯ МЕТАДАНИХ (СТРУКТУРА КОЛОНОК) ---

REQUIRED_EOSSA_COLUMNS_METADATA = [
    # TTYPE | Format (TFORM) | Unit (TUNIT) | Notes
    {'name': 'UTC_Begin_exp', 'format': '23A', 'unit': 'yyyy-mm-ddThh:mm:ss.sss'},
    {'name': 'UTC_End_exp', 'format': '23A', 'unit': 'yyyy-mm-ddThh:mm:ss.sss'},
    {'name': 'JD_Mid_Exp', 'format': 'D', 'unit': 'days'},
    {'name': 'Exp_Duration', 'format': 'D', 'unit': 'sec'},
    {'name': 'Cur_Spec_Filt_Num', 'format': 'J', 'unit': 'num'},
    {'name': 'Cur_ND_Filt_Num', 'format': 'J', 'unit': 'num'},
    {'name': 'Mag_Exo_Atm', 'format': 'D', 'unit': 'mag'},
    {'name': 'Mag_Range_Norm', 'format': 'D', 'unit': 'mag'},
    {'name': 'Eph_RA_DE', 'format': '2D', 'unit': 'deg'},
    {'name': 'Met_RA_DE', 'format': '2D', 'unit': 'deg'},
    {'name': 'Eph_AZ_EL', 'format': '2D', 'unit': 'deg'}, # Required for Ground-based [8]
    {'name': 'Met_AZ_EL', 'format': '2D', 'unit': 'deg'}, # Required for Ground-based [8]
    {'name': 'Tel_Obj_Range', 'format': 'D', 'unit': 'm'},
    {'name': 'Solar_Phase_Ang', 'format': 'D', 'unit': 'deg'},
]

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
            size = int(fmt.replace('A', '')) if fmt.replace('A', '').isdigit() else 23 # UTC size
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
    t_start = Time(obs_data['Date UT'] + obs_data['Time Start'], format='iso', scale='utc', precision=6)
    t_end = Time(obs_data['Date UT'] + obs_data['Time Start'], format='iso', scale='utc', precision=6) # Використовуємо початок, оскільки Exp_Duration = 0.1s
    t_end += obs_data['dt'] / 24.0 / 3600.0 # Припускаємо dt = 0.1s
    t_mid = t_start + (t_end - t_start) / 2
    
    single_row['UTC_Begin_exp'] = t_start.iso.replace(' ', 'T')[:23]
    single_row['UTC_End_exp'] = t_end.iso.replace(' ', 'T')[:23]
    single_row['JD_Mid_Exp'] = t_mid.jd
    single_row['Exp_Duration'] = obs_data['dt'] # 0.1 секунди, якщо dt = 0.1

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
    single_row['Tel_Obj_Range'] = d_km * 1000.0 # м
    single_row['Cur_Spec_Filt_Num'] = 1 
    single_row['Cur_ND_Filt_Num'] = PLACEHOLDER_INT # Немає нейтрального фільтра [7]
    
    # Векторні поля (2D)
    # Eph_RA_DE, Met_RA_DE (Припускаємо, що RA/DE обчислюються)
    single_row['Eph_RA_DE'] = np.array([0.0, 0.0]) # Заглушка, оскільки ці дані відсутні у вихідному файлі
    single_row['Met_RA_DE'] = np.array([0.0, 0.0])
    
    # Eph_AZ_EL, Met_AZ_EL (Вхідні дані містять Az(deg) і El(deg))
    single_row['Eph_AZ_EL'] = np.array([0.0, 0.0]) # Ефемеридні Азимут/Елевація
    single_row['Met_AZ_EL'] = np.array([obs_data['Az(deg)'], obs_data['El(deg)']]) # Виміряні Азимут/Елевація [8]

    # Solar_Phase_Ang (Припускаємо, що обчислюється)
    single_row['Solar_Phase_Ang'] = -9999.0 # Placeholder для D [7]

    return single_row # Повертаємо єдиний структурований елемент

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
    hdu.header['TFIELDS'] = len(column_definitions) # TFIELDS (Number of fields) [11]
    hdu.header['NAXIS2'] = data_len
    
    return hdu, data_len

# --- ПРИКЛАД ОБРОБКИ ПОТОКУ ДАНИХ ---

# 1. Визначення структури FITS (Крок 1)
eossa_column_definitions = create_eossa_columns(REQUIRED_EOSSA_COLUMNS_METADATA)

# 2. Ініціалізація списку для збору даних
collected_rows = []

# 3. Імітація обробки кадрів по одному (потік даних)
# Використовуємо дані з вашого файлу "Спостереження супутника DMSP 5D-2 F12" [12-14]
example_data_stream = [
    {'Date UT': '2025-09-12', 'Time Start': '19:23:12.404', 'dt': 0.1, 'magV': 6.587, 'Rg(Mm)': 1231.248, 'Az(deg)': 211.950, 'El(deg)': 39.089},
    {'Date UT': '2025-09-12', 'Time Start': '19:23:12.504', 'dt': 0.1, 'magV': 6.540, 'Rg(Mm)': 1230.869, 'Az(deg)': 211.981, 'El(deg)': 39.107},
    {'Date UT': '2025-09-12', 'Time Start': '19:23:12.604', 'dt': 0.1, 'magV': 6.692, 'Rg(Mm)': 1230.491, 'Az(deg)': 212.011, 'El(deg)': 39.126},
    # ... тут можуть бути ще сотні рядків ...
    {'Date UT': '2025-09-12', 'Time Start': '19:23:14.707', 'dt': 0.1, 'magV': 6.276, 'Rg(Mm)': 1222.620, 'Az(deg)': 212.665, 'El(deg)': 39.521},
    {'Date UT': '2025-09-12', 'Time Start': '19:23:14.807', 'dt': 0.1, 'magV': 6.345, 'Rg(Mm)': 1222.249, 'Az(deg)': 212.696, 'El(deg)': 39.540},
]

print("--- Обробка потоку даних по одному рядку ---")

for data_row in example_data_stream:
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


# Пояснення реалізації
# 1. Функція process_single_observation:
#     ◦ Вона приймає лише один словник (obs_data) [Запит].
#     ◦ Вона використовує попередньо визначений EOSSA_DTYPE (структура NumPy), щоб створити контейнер для одного рядка даних (np.empty(1, dtype=EOSSA_DTYPE)).
#     ◦ Всі необхідні розрахунки (перетворення часу, денормалізація $M_{r} \rightarrow m_{r}$) виконуються всередині цієї функції для цього одного рядка.
#     ◦ Вона повертає готовий структурований елемент, який відповідає одному рядку в майбутній таблиці FITS.
# 2. Збір даних (collected_rows):
#     ◦ Записи, повернені process_single_observation, накопичуються у звичайному Python-списку (collected_rows).
# 3. Функція finalize_eossa_table (Фіналізація):
#     ◦ Крок 1: Об'єднання. Список collected_rows перетворюється на єдиний, фінальний, двовимірний NumPy масив (final_data_array) за допомогою np.stack().
#     ◦ Крок 2: Створення HDU та $NAXIS2$. При виклику fits.BinTableHDU.from_columns(..., data=final_data_array) бібліотека astropy автоматично обчислює загальну кількість рядків (len(final_data_array)) і записує це значення в заголовок як $NAXIS2$, а також записує загальну кількість стовпців як $TFIELDS$. Таким чином, це задовольняє вашій вимозі "заповнити поля типу data_len, попередньо порахувавши їх" в кінці скрипту.
