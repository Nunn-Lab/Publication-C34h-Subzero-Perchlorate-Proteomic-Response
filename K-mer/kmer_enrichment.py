from scipy.stats import chisquare, chi2_contingency
import statsmodels.stats.multicomp
import numpy as np
from matplotlib import pyplot
import matplotlib.pyplot as plt
import random
import pandas as pd
import math
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency

def calculate_chi2_enrichment(df_control, df_experiment, features):
    """
    change this so that it takes an array of features (kmers) and does a chi2 analysis on each one. 
    """
    # Merge the two DataFrames on the shared feature column
    merged_df = df_control.merge(df_experiment, on='feature', suffixes=('_control', '_experiment'))

    # Filter the merged DataFrame to include only the specified features
    merged_df = merged_df[merged_df['feature'].isin(features)]

    # Calculate chi-squared enrichment for each feature
    chi2_results = []
    for index, row in merged_df.iterrows():
        contingency_table = np.array([[row['count_control'], row['count_experiment']],
                                      [row['total_count_control'] - row['count_control'],
                                       row['total_count_experiment'] - row['count_experiment']]])
        chi2, p, _, _ = chi2_contingency(contingency_table, correction=False)
        chi2_results.append({'feature': row['feature'], 'chi2': chi2, 'p-value': p})

    return pd.DataFrame(chi2_results)

