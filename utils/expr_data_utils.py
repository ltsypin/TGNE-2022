import pandas as pd
import numpy as np
import scipy.stats as st
from functools import partial

def mean_df_of_duplicates(df: pd.DataFrame, mean_type: str):
    """_summary_

    Parameters
    ----------
    df : pd.DataFrame
        _description_
    mean_type : str
        'geometric' or 'arithmetic'

    Returns
    -------
    _type_
        _description_
    """
    return_df = df.loc[:, list(df.columns)[0:1]]
    num_df = df.loc[:, list(df.columns)[1:]]

    num_df_cols = list(num_df.columns)

    col_splits = [col.split('_') for col in num_df_cols]

    all_col_prefixes = ['_'.join(col_split[:len(col_split)-1]) for col_split in col_splits]

    indexes = np.unique(all_col_prefixes, return_index=True)[1]

    col_prefixes = [all_col_prefixes[index] for index in sorted(indexes)]

    for col_prefix in col_prefixes:
        cols = [num_df[col].values for col in num_df_cols if col.startswith(col_prefix)]

        if mean_type == 'geometric':
            return_df[col_prefix] = st.mstats.gmean(np.array(cols) + 1, axis=0) - 1
        elif mean_type == 'arithmetic':
            return_df[col_prefix] = np.mean(np.array(cols), axis=0)
        else:
            raise ValueError(f'Invalid mean_type: {mean_type}')
    
    return return_df


# Z-SCORE
def z_score_normalizer(array, add_scalar=0):
    log_a = np.log10(array + add_scalar)
    normalized = st.zscore(log_a)
    return normalized

# MIN-MAX
def min_max_normalizer(array):    
    normalized = np.divide((array - np.min(array)), (np.max(array) - np.min(array)))
    return normalized


def normalize_expression_per_gene(expression_df: pd.DataFrame, norm_type: str, add_scalar=0):
    """_summary_

    Parameters
    ----------
    expression_df : pd.DataFrame
        _description_
    norm_type : str
        'z_score' or 'min_max'
    add_scalar : int, optional
        Only used if 'z_score' is chosen as the norm_type, by default 0

    Returns
    -------
    _type_
        _description_
    """

    z_score_normalizer_with_scalar = partial(z_score_normalizer, add_scalar=add_scalar)

    if expression_df.dtypes.to_list()[0] == object:
        id_col = list(expression_df.columns)[0]
        val_cols = list(expression_df.columns)[1:]
        ids = expression_df[id_col].values
        data = expression_df[val_cols]

        if norm_type == 'z_score':
            norm_expression_df = data.apply(z_score_normalizer_with_scalar, axis=1)
        elif norm_type == 'min_max':
            norm_expression_df = data.apply(min_max_normalizer, axis=1)
        else:
            raise ValueError(f'Invalid norm_type: {norm_type}')
        
        norm_expression_df[id_col] = ids
        
        columns = norm_expression_df.columns.tolist()
        
        rearrangment = columns[-1:] + columns[:-1]
        
        norm_expression_df = norm_expression_df[rearrangment]
        
    else:
        if norm_type == 'z_score':
            norm_expression_df = expression_df.apply(z_score_normalizer_with_scalar, axis=1)
        elif norm_type == 'min_max':
            norm_expression_df = expression_df.apply(min_max_normalizer, axis=1)
        else:
            raise ValueError(f'Invalid norm_type: {norm_type}')
        
    return norm_expression_df

