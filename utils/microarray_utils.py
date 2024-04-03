import pandas as pd
import numpy as np
import scipy.stats as st

# microarray_utils

def get_geom_mean_expression(expression_df):
    """
    Function to take an expression dataframe from the microarrays and collapse it into the means of
    all replicate chips.
    """
    # C2 and S12 got removed during quality control
    x = [
        'Ll', 
        'Lm', 
        'Lh', 
        'S0', 
        'S3', 
        'S6', 
        'S9', 
        # 'S12', 
        'S15', 
        'S24', 
        'C0', 
        # 'C2', 
        'C4', 
        'C6', 
        'C8', 
        'C10', 
        'C12', 
        'C14', 
        'C16', 
        'C18']
    
    # cols = expression_df.columns[1:]
    # x = [c for c in x if c in cols]
    
    condition_expr_dict = {c.split("_")[0]: [] for c in expression_df.columns[1:]}
    
    for c in list(expression_df.columns)[1:]:
        
        cond = c.split('_')[0]
        if cond in condition_expr_dict.keys():
            expr_list = condition_expr_dict.get(cond, [])

            # Need to avoid true zeros
            expr_list.append(expression_df[c].values)
            condition_expr_dict[cond] = expr_list
        
    condition_mean_dict = {c: (st.mstats.gmean(np.array(condition_expr_dict[c]) + 1, 0) - 1) for c in condition_expr_dict.keys() if c in x}
    
    mean_expr_df = pd.DataFrame(condition_mean_dict)
    mean_expr_df['TTHERM_ID'] = expression_df['TTHERM_ID'].values
    cols = list(mean_expr_df.columns)
    reorder = cols[-1:] + cols[:-1]
    mean_expr_df = mean_expr_df[reorder]
    
    return mean_expr_df

def normalizer(array):
    """
    Normalizes the values of an array to range from zero to one
    """
    
    a = np.array(array)
    
    normalized = (array - np.min(array)) / (np.max(array) - np.min(array))
    
    return normalized

def normalize_expression_per_gene(expression_df):
    """
    Function to normalize all gene expression to range from zero to one.
    """
    if 'TTHERM_ID' in expression_df.columns:
        ttids = expression_df['TTHERM_ID'].values
        data = expression_df[list(expression_df.columns)[1:]]
        
        norm_expression_df = data.apply(lambda row: normalizer(row), axis=1)
        norm_expression_df['TTHERM_ID'] = ttids
        
        columns = norm_expression_df.columns.tolist()
        
        rearrangment = columns[-1:] + columns[:-1]
        
        norm_expression_df = norm_expression_df[rearrangment]
        
    else:
        norm_expression_df = expression_df.apply(lambda row: normalizer(row), axis=1)
    
    return norm_expression_df
