import pickle
from multiprocessing import Pool
from tqdm import tqdm
import os

# HELPER FUNCTION FOR BUILDING A TETRA_ID TO TTHERM_IDS LIST DICT
# FOR EACH TETRA_ID, LOOP THROUGH ALL THE CORRESPONDING PROBE_IDS IN THE TETRA_ID TO PROBE_IDS LIST DICT
    # FOR EACH PROBE_ID, LOOP THROUGH ALL OF THE PROBE_IDS LISTS IN THE TTHERM_ID TO PROBE_IDS LIST DICT
        # IF THE PROBE_ID IS FOUND IN THE PROBE_IDS LIST, APPEND THE CORRESPONDING TTHERM_ID TO THE TTHERM_IDS LIST FOR THAT TETRA_ID (IF IT IS NOT ALREADY PRESENT (PREVENTS DUPLICATES))
def process_tetra_id_wrapper(tetra_id):
    file_dir = os.path.dirname(os.path.abspath(__file__))

    tetra_ttherm_dict = {}

    with open(f'{file_dir}/seq_probe_list_dict.pkl', 'rb') as file:
        seq_probe_list_dict = pickle.load(file)

    with open(f'{file_dir}/align_list_dict.pkl', 'rb') as file:
        align_list_dict = pickle.load(file)

    align_list_dict_keys = list(align_list_dict.keys())

    lenC = len(align_list_dict_keys)

    lenB = len(seq_probe_list_dict[tetra_id])
    idxB = 0
    while idxB < lenB:

        probe_id = seq_probe_list_dict[tetra_id][idxB]
        
        idxC = 0
        while idxC < lenC:
            
            ttherm_id = align_list_dict_keys[idxC]

            if probe_id in align_list_dict[ttherm_id]:
                if tetra_id not in tetra_ttherm_dict:
                    tetra_ttherm_dict[tetra_id] = []
                # filter out duplicates
                if ttherm_id not in tetra_ttherm_dict[tetra_id]:
                    tetra_ttherm_dict[tetra_id].append(ttherm_id)

            idxC += 1
        
        idxB += 1
        
    return tetra_ttherm_dict

if __name__ == '__main__':

    file_dir = os.path.dirname(os.path.abspath(__file__))

    with open(f'{file_dir}/seq_probe_list_dict.pkl', 'rb') as file:
        seq_probe_list_dict = pickle.load(file)

    seq_probe_list_dict_keys = list(seq_probe_list_dict.keys())

    # CREATE A LIST OF ALL OF THE TETRA_IDS
    tetra_ids = seq_probe_list_dict_keys

    # PROCESS EACH OF THE TETRA_IDS IN A DIFFERENT PYTHON PROCESS AT THE SAME TIME
    with Pool() as pool:
        tetra_ttherm_dict_list = list(tqdm(pool.imap(process_tetra_id_wrapper, tetra_ids), total=len(tetra_ids)))

    # COMBINE EACH OF THE TETRA_ID TO TTHERM_IDS LIST DICTS INTO A SINGLE DICTIONARY
    tetra_ttherm_dict = {}
    for result_dict in tetra_ttherm_dict_list:
        tetra_ttherm_dict.update(result_dict)

    # SERIALIZE THE COMPLETE TETRA_ID TO TTHERM_IDS LIST DICT
    with open(f'{file_dir}/tetra_ttherm_dict.pkl', 'wb') as file:
        pickle.dump(tetra_ttherm_dict, file)
