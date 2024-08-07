import sys
import os
import pandas as pd
import numpy as np
import scipy.stats as st
# import tqdm
import sys
import pickle
from multiprocessing import Pool

import annotation_info

# ENRICHMENT_ESSENTIAL_FINAL

main_dir = '../'

complete_annotation = pd.read_csv(os.path.join(main_dir, 'active_files/complete_annotation.csv'))

background_annotation = complete_annotation
term_columns = ['COG_category', 'GOs', 'KEGG_ko', 'EC', 'PFAMs', 'InterPro']

if 'mucocysts' in list(complete_annotation.columns):
    term_columns.append('mucocysts')

def term_count_dict_from_annotation_df(annot_df, term_column):
    
    column = annot_df[term_column].values
    
    funct_terms = []
    for entry in column:
        if entry != '-':
            if term_column == 'COG_category':
                # terms = [f'{e}: {COG_dict[e]}' for e in entry]
                terms = list(entry)
            else:
                terms = entry.split(',')
            for t in terms:
                funct_terms.append(t)

#     len(terms)
    
    term_count_dict = {}
    
    for t in funct_terms:
        count = term_count_dict.get(t, 0)
        count += 1
        term_count_dict[t] = count
        
    return term_count_dict

def enrichment_analysis(module, label_df, background_annotation, term_column, report_counts=False):
    
    # module_ttids = np.array([id for id in list(label_df.loc[label_df[f'label'] == module]['TTHERM_ID'].values) if "TTHERM" in id])
    module_ttids = label_df.loc[label_df[f'label'] == module]['TTHERM_ID'].values

    module_annotation = background_annotation.loc[background_annotation['TTHERM_ID'].isin(module_ttids)]
    
    background_term_dict = term_count_dict_from_annotation_df(background_annotation, term_column)
    module_term_dict = term_count_dict_from_annotation_df(module_annotation, term_column)
    
    bs = []
    ps = []
    folds = []
    terms = []

    if report_counts:
        mcs = []
        mss = []
    
    for t, module_count in module_term_dict.items():
        
        background_count = background_term_dict[t]
        module_size = len(module_annotation)
        background_size = len(background_annotation)
        
        standard_contingency_table = [
                                [module_count, background_count - module_count], 
                                [module_size - module_count, background_size - module_size - (background_count - module_count)]
                            ]
        
        # The -1 and +1 make this more conservative (see explanation from the DAVID database: 
        # https://david.ncifcrf.gov/helps/functional_annotation.html#geneenrich)
        conservative_contingency_table = [
                                [module_count - 1, background_count - module_count + 1], 
                                [module_size - module_count, background_size - module_size - (background_count - module_count)]
                            ]
        
        
        odds, p_standard = st.fisher_exact(standard_contingency_table, 'greater')
        odds, p_conservative = st.fisher_exact(conservative_contingency_table, 'greater')
        
        p_reasonable = np.mean([p_standard, p_conservative])
        
        bonferroni  = p_reasonable * len(module_term_dict)

        fold_enrichment = (module_count/module_size) / (background_count/background_size)

        if bonferroni <= 0.05:
            
            ps.append(p_reasonable)
            bs.append(bonferroni)
            folds.append(fold_enrichment)
            terms.append(t)
            if report_counts:
                mcs.append(module_count)
                mss.append(module_size)
            
#         else:
#             ps.append('')
#             bs.append('')
#             folds.append('')
#             terms.append('')
    if report_counts:
        return ps, bs, folds, terms, mcs, mss
            
    return ps, bs, folds, terms
        
# ENRICHMENT_ESSENTIAL_FINAL END
def init_pool(data):
    global lldf
    lldf = data
    global REPORT_COUNTS
    REPORT_COUNTS = True # FIXME

def process_module(m):
    term_dfs = []
    
    for tc in term_columns:

        if REPORT_COUNTS:
            ps, bs, folds, terms, mcs, mss = enrichment_analysis(m, lldf, background_annotation, tc, report_counts=REPORT_COUNTS)
        else:
            ps, bs, folds, terms = enrichment_analysis(m, lldf, background_annotation, tc, report_counts=REPORT_COUNTS)

        info = []

        if tc == 'GOs':
            for t in terms:
                name, definition, obsolete = annotation_info.get_GO_info(t)
                info.append(f'{name.capitalize()}: {definition}')
                    
        elif tc == 'COG_category':
            for t in terms:
                info.append(annotation_info.COG_dict[t])
                            
        elif tc == 'KEGG_ko':
            for t in terms:
                info.append(annotation_info.get_KEGG_info(t))
                
        elif tc == 'EC':
            for t in terms:
                info.append(annotation_info.get_EC_info(t))

        elif tc == 'PFAMs':
            for t in terms:
                info.append(annotation_info.get_PFAM_info(t))

        elif tc == 'InterPro':
            for t in terms:
                info.append(annotation_info.get_InterPro_info(t))

        elif tc == 'mucocysts':
            for t in terms:
                if t == 'DE':
                    info.append('differentially expressed in mucocyst regranulation dataset')
                elif t == 'EV':
                    info.append('experimentally validated mucocyst-associated gene')
                
        term_df = pd.DataFrame({'module': [m]*len(terms),
                                'term': terms,
                                'info': info,
                                'fold_change': folds,
                                'bonferroni': bs,
                                'term_count': mcs,
                                'module_size': mss
                                })
        
        term_dfs.append(term_df)
    
    module_df = pd.concat(term_dfs)
    return module_df

if __name__ == '__main__':
    
    data_bytes = sys.stdin.buffer.read()
    lldf = pickle.loads(data_bytes)
    
    # Process data in parallel
    with Pool(initializer=init_pool, initargs=(lldf, )) as pool:
        # module_dfs = list(tqdm.tqdm(pool.imap(process_module, sorted(lldf['label'].unique())), total=len(lldf['label'].unique())))
        module_dfs = list(pool.imap(process_module, sorted(lldf['label'].unique())))

    # Concatenate results
    all_enrichment_df = pd.concat(module_dfs)

    output_bytes = pickle.dumps(all_enrichment_df)
    sys.stdout.buffer.write(output_bytes)