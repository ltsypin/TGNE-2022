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
        