import pandas as pd
import numpy as np
import scipy.stats as st

# rna_seq_utils

# Z-SCORE
def normalizer(array):
    log_a = [np.log10(tpm + 1) for tpm in array]
    normalized = st.zscore(log_a)
    return normalized

def normalize_expression_per_gene(expression_df):
    if 'TTHERM_ID' in expression_df.columns:
        ttids = expression_df['TTHERM_ID'].values
        data = expression_df[list(expression_df.columns)[1:]]
        
        norm_expression_df = pd.DataFrame(data.apply(normalizer, axis=1).tolist(), columns=data.columns)

        norm_expression_df['TTHERM_ID'] = ttids

        columns = norm_expression_df.columns.tolist()

        rearrangement = columns[-1:] + columns[:-1]

        norm_expression_df = norm_expression_df[rearrangement]
    else:
        norm_expression_df = pd.DataFrame(expression_df.apply(normalizer, axis=1).tolist(), columns=expression_df.columns)

    return norm_expression_df

def geo_mean_df_of_duplicates(df: pd.DataFrame):
    return_df = df.loc[:, df.columns[0:1]]
    num_df = df.loc[:, df.columns[1:]]

    idxa = 0
    idxb = 1

    num_df_cols = list(num_df.columns)

    while idxb < len(num_df_cols):
        col_a_split = num_df_cols[idxa].split('_')
        col_name = '_'.join(col_a_split[:len(col_a_split)-1])
        return_df[col_name] = np.sqrt(num_df[num_df_cols[idxa]] * num_df[num_df_cols[idxb]])

        idxa += 1
        idxb += 1
    
    return return_df

def ari_mean_df_of_duplicates(df: pd.DataFrame):
    return_df = df.loc[:, df.columns[0:1]]
    num_df = df.loc[:, df.columns[1:]]

    idxa = 0
    idxb = 1

    num_df_cols = list(num_df.columns)

    while idxb < len(num_df_cols):
        col_a_split = num_df_cols[idxa].split('_')
        col_name = '_'.join(col_a_split[:len(col_a_split)-1])
        return_df[col_name] = num_df[num_df_cols[idxa]] + num_df[num_df_cols[idxb]] / 2

        idxa += 1
        idxb += 1
    
    return return_df

# # FOLD CHANGE FROM MEAN EXPRESSION
# def normalizer(array):
#     array_mean = np.mean(array)
#     normalized = [np.log2(tpm) for tpm in ((array + 1) / (array_mean + 1))]
#     return normalized

# def normalize_expression_per_gene(expression_df):
#     if 'TTHERM_ID' in expression_df.columns:
#         ttids = expression_df['TTHERM_ID'].values
#         data = expression_df[list(expression_df.columns)[1:]]
        
#         norm_expression_df = pd.DataFrame(data.apply(normalizer, axis=1).tolist(), columns=data.columns)

#         norm_expression_df['TTHERM_ID'] = ttids

#         columns = norm_expression_df.columns.tolist()

#         rearrangement = columns[-1:] + columns[:-1]

#         norm_expression_df = norm_expression_df[rearrangement]
#     else:
#         norm_expression_df = pd.DataFrame(expression_df.apply(normalizer, axis=1).tolist(), columns=expression_df.columns)

#     return norm_expression_df

# # MIN-MAX
# def normalizer(array):
#     """
#     Normalizes the values of an array to range from zero to one
#     """
    
#     a = np.array(array)
    
#     normalized = (array - np.min(array)) / (np.max(array) - np.min(array))
    
#     return normalized

# def normalize_expression_per_gene(expression_df):
#     """
#     Function to normalize all gene expression to range from zero to one.
#     """
#     if 'TTHERM_ID' in expression_df.columns:
#         ttids = expression_df['TTHERM_ID'].values
#         data = expression_df[list(expression_df.columns)[1:]]
        
#         norm_expression_df = data.apply(lambda row: normalizer(row), axis=1)
#         norm_expression_df['TTHERM_ID'] = ttids
        
#         columns = norm_expression_df.columns.tolist()
        
#         rearrangment = columns[-1:] + columns[:-1]
        
#         norm_expression_df = norm_expression_df[rearrangment]
        
#     else:
#         norm_expression_df = expression_df.apply(lambda row: normalizer(row), axis=1)
    
#     return norm_expression_df

