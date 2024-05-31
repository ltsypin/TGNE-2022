import os
from csv import DictWriter
import zipfile
import requests
import tarfile
import re

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

def unzip_file(zip_file: str, extract_to=None):
    if extract_to is None:
        extract_to = zip_file[: len(zip_file) - len('.zip')]
    with zipfile.ZipFile(zip_file, 'r') as zip_ref:
        zip_ref.extractall(extract_to)

def download_file_chunks(url: str, output_file_path: str, chunk_size=128):
    r = requests.get(url, stream=True)
    with open(output_file_path, 'wb') as f:
        for chunk in r.iter_content(chunk_size=chunk_size):
                f.write(chunk)
    print(f'FILE SAVED: {os.path.abspath(output_file_path)}')

def extract_tar(tar_file_path: str):
    tar_folder = tarfile.open(tar_file_path)
    output_dir = tar_file_path[:len(tar_file_path) - len('.tar')]
    tar_folder.extractall(path=output_dir)
    tar_folder.close()

def download_geo_data_file(url: str, accession: str):
    r = requests.get(url, stream=True)
    d = r.headers['content-disposition']
    fname = re.findall("filename=(.+)", d)[0]
    accession = re.search(r'GSE.*(tar|gz)', fname).group()

    download_file_chunks(url, f'./new_raw_data/{accession}')
    
    if '.tar' in accession:
        extract_tar(f'./new_raw_data/{accession}')

def move_file(file_path: str, target_dir: str):
    name = os.path.basename(file_path)
    os.rename(file_path, os.path.join(target_dir, name))

def move_files(file_paths: list, target_dir: str):
    for f in file_paths:
        move_file(f, target_dir)
