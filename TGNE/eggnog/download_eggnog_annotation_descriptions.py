import sys
import os
import subprocess

file_dir = os.path.dirname(os.path.abspath(__file__))

download_scripts = [
    'download_annotation_descriptions_go.py',
    'download_annotation_descriptions_pfam.py',
    'download_annotation_descriptions_kegg_ec.py',
]

if not all([
        os.path.exists(os.path.join(file_dir, '../../active_files/ec_annotations.csv')), 
        os.path.exists(os.path.join(file_dir, '../../active_files/kegg_annotations.csv')), 
        os.path.exists(os.path.join(file_dir, '../../active_files/go_annotations.csv')), 
        os.path.exists(os.path.join(file_dir, '../../active_files/pfam_annotations.csv')), 
        ]):
    
    for script in download_scripts:
        print(f'EXECUTING: {script}')
        subprocess.run([sys.executable, os.path.join(file_dir, script)])