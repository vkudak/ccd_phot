import os
import sys
import argparse
from datetime import datetime
from astropy.io import fits
from tqdm import tqdm


# 1. СПИСОК ПОЛІВ ДЛЯ ЗАМІНИ
# Додайте сюди всі ключі заголовока (keywords), які потрібно перевірити
KEYWORDS_TO_UPDATE = ['DATE-OBS', 'MJD-OBS', 'DATE', 'GPS_ST', 'GPS_ET', 'GPS_NT', 'DATE-END', 'DATE-OB2', 'DATE-AVG']

# Нова дата, яку потрібно встановити (формат: РРРР-ММ-ДД)
# new_date_str = "2026-05-28" 


def modify_date_keep_time(current_value, new_date_iso):
    """
    Аналізує поточне значення ключа, виділяє час і підставляє нову дату.
    """
    if not isinstance(current_value, str):
        return None
    
    # Спроба 1: Стандартний ISO формат (YYYY-MM-DDTHH:MM:SS.xxx)
    if 'T' in current_value:
        try:
            date_part, time_part = current_value.split('T', 1)
            # Перевіряємо, чи перша частина схожа на дату
            datetime.strptime(date_part, "%Y-%m-%d") 
            return f"{new_date_iso}T{time_part}"
        except ValueError:
            pass

    # Спроба 2: Старий формат FITS (DD/MM/YY) - якщо раптом час записаний через пробіл або його немає
    if '/' in current_value and ' ' in current_value:
        try:
            parts = current_value.split(' ', 1)
            time_part = parts[1]
            # Конвертуємо нову дату у старий формат для збереження стилю, або залишаємо ISO
            # За стандартом краще переводити в ISO, але змініть за потреби
            new_date_obj = datetime.strptime(new_date_iso, "%Y-%m-%d")
            new_date_old_fmt = new_date_obj.strftime("%d/%m/%y")
            return f"{new_date_old_fmt} {time_part}"
        except ValueError:
            pass

    return None


def process_fits_directory(directory_path, new_date_str):
    if not os.path.isdir(directory_path):
        print(f"Помилка: Директорія '{directory_path}' не існує.")
        sys.exit(1)

    print(f"Start processing of FITS files in: {directory_path}")
    print(f"New date: {new_date_str}")
    print(f"Field wre we change date: {KEYWORDS_TO_UPDATE}\n" + "-"*40)


    # Збираємо всі потрібні файли в список, щоб tqdm знав загальну кількість для progress bar
    fits_files = [
        f for f in os.listdir(directory_path) 
        if f.lower().endswith(('.fits', '.fit'))
    ]

    if not fits_files:
        print("No FITS files in selected DIR.")
        return

    print(f"Find {len(fits_files)} files. Date update...")

    # Ініціалізуємо tqdm прогрес-бар
    for filename in tqdm(fits_files, desc="Обробка FITS", unit="file"):

    # Шукаємо всі .fits та .fit файли
    # for filename in os.listdir(directory_path):
        if filename.lower().endswith(('.fits', '.fit')):
            file_path = os.path.join(directory_path, filename)
            file_changed = False

            try:
                # Відкриваем файл у режимі оновлення ('update')
                with fits.open(file_path, mode='update') as hdul:
                    # Проходимо по всіх HDU (Header Data Units)
                    for hdu in hdul:
                        header = hdu.header
                        
                        for key in KEYWORDS_TO_UPDATE:
                            if key in header:
                                old_value = header[key]
                                new_value = modify_date_keep_time(old_value, new_date_str)
                                
                                if new_value and old_value != new_value:
                                    header[key] = new_value
                                    # print(f"[{filename}] Оновлено {key}: '{old_value}' -> '{new_value}'")
                                    file_changed = True
                    
                    # Якщо зміни були, вони автоматично збережуться при закритті контекстного менеджера
                    if file_changed:
                        hdul.flush()  # Примусовий запис змін
                        
            except Exception as e:
                print(f"Error in processing file {filename}: {e}")

    print("-"*40 + "\nDone.")

if __name__ == "__main__":
    # Налаштування аргументів командного рядка
    parser = argparse.ArgumentParser(description="Change date in FITS files with time saving.")
    parser.add_argument("-nd", "--new_date", type=str, help="New date in format YYYY-MM-DD")
    parser.add_argument("-d", "--dir", type=str, help="Path to dir with FITS files")
    
    args = parser.parse_args()

    process_fits_directory(args.dir, args.new_date)