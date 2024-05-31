import pandas as pd
import numpy as np
import sqlite3
import scipy.stats as st

# dataframe_utils

def get_hypercube_sample(dimensions, samples):
    sampler = st.qmc.LatinHypercube(d=dimensions)
    # sampler = st.qmc.Sobol(d=dimensions)

    hypercube_sample = sampler.random(n=samples)

    # print(st.qmc.discrepancy(hypercube_sample, workers=-1))

    return pd.DataFrame(hypercube_sample)

def shuffle_row(row):
    shuffled_row = row.values.copy()
    np.random.shuffle(shuffled_row)
    return pd.Series(shuffled_row, index=row.index)

def shuffle_rows(df):
    columns_to_shuffle = df.columns[1:]
    df[columns_to_shuffle] = df[columns_to_shuffle].apply(shuffle_row, axis=1)
    return df

def sql_query_df(dataframes: dict, query: str):
    # Create an in-memory SQLite database connection
    connection = sqlite3.connect(':memory:')

    try:
        # Iterate over each dataframe in the list
        for table_name, dataframe in dataframes.items():
            # Convert the dataframe to an SQLite table
            dataframe.to_sql(table_name, con=connection, if_exists="replace", index=False)

        # Execute the provided SQL query on the transformed data
        cursor = connection.cursor()
        cursor.execute(query)

        # Fetch all the results from the query
        results = cursor.fetchall()

        # Get the column names from the query results
        column_names = [desc[0] for desc in cursor.description]

        # Create a new DataFrame from the query results
        transformed_data = pd.DataFrame(results, columns=column_names)

        return transformed_data

    finally:
        # Close the cursor and connection
        cursor.close()
        connection.close()

def csv_files_to_df(files: list):

    combined_df = None

    for f in files:

        if combined_df is None:
            combined_df = pd.read_csv(f)
            continue

        curr_df = pd.read_csv(f)

        combined_df = pd.concat([combined_df, curr_df], ignore_index=True)

    return combined_df