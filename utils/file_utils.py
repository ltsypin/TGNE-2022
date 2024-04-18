import os
from csv import DictWriter

# file_utils

def write_to_csv(csv_file_path, data_item, header):
    # Check if the CSV file exists and write header if it doesn't
    if not os.path.isfile(csv_file_path):
        with open(csv_file_path, 'w', newline='') as file:
            writer = DictWriter(file, fieldnames=header)
            writer.writeheader()

    with open(csv_file_path, 'a', newline='') as file:
        writer = DictWriter(file, fieldnames=header)
        writer.writerow(data_item)

def remove_file(file_path):
    if os.path.exists(file_path):
        os.remove(file_path)

def create_directories(dirPathString):
    if not os.path.exists(dirPathString):
        os.makedirs(dirPathString)

def generate_uniquely_numbered_export_path(export_dir: str, file_name, file_ext, tags=[], start_num=1, num_step=1):
    num = start_num

    create_directories(export_dir)

    files = os.listdir(export_dir)

    while any(file.startswith(f'{file_name}{num}') for file in files):
        num += num_step

    return f'{export_dir}{file_name}{num}_{"_".join(tags)}{file_ext}'