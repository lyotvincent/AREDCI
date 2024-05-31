
import os
import pandas as pd
import numpy as np
from scipy.stats import norm

def norm_sf(infile, outfile, group):

    # * read data
    df = pd.read_csv(infile, index_col=0, header=0)
    print(f"[INFO] read data: {df.shape}")
    print(df[:4])
    pets = df.to_numpy()

    # * compute D
    # mean of first 4 col, mean of last col
    d_list = pets[:, :4].mean(axis=1) - pets[:, 4:].mean(axis=1)
    print(f"[INFO] compute D: {d_list.shape}")

    # * normalize m
    d_list = ( d_list - np.mean(d_list) ) / np.std(d_list)

    # * compute pvalue
    pvalues = norm.sf(d_list, loc=0, scale=1)

    pvalues[pvalues > 0.5] = 1 - pvalues[pvalues > 0.5]
    pvalues *= 2

    # * save pvalues
    pvalues = pd.DataFrame({
        "logFC": [0]*len(pvalues),
        "logCPM": [0]*len(pvalues),
        "LR": [0]*len(pvalues),
        "PValue": pvalues,
        "FDR": [0]*len(pvalues)
    }, index=df.index)
    pvalues.to_csv(outfile)
    print(f"[INFO] save pvalues to {outfile}")
    


if __name__ == '__main__':
    ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    K562_MCF7_EXP_FILE = os.path.join(ROOT_DIR, 'data', 'DNase_seq', 'dnaseseq_k562_mcf7_signalValue_matrix.csv')
    K562_MCF7_EXP_NORM_SF_FILE = os.path.join(ROOT_DIR, 'data', 'DNase_seq', 'K562_MCF7_pvalues.csv')
    group = ["K562_REP1", "K562_REP2", "K562_REP3", "K562_REP4", "MCF7_REP1", "MCF7_REP2", "MCF7_REP3", "MCF7_REP4"]
    norm_sf(K562_MCF7_EXP_FILE, K562_MCF7_EXP_NORM_SF_FILE, group)
    
