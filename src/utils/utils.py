import numpy as np
import os
import pandas as pd


def get_PETs_from_loopinfo(data):
    """
    Get PETs from loopinfo data
    """
    PETs = np.empty(shape=(len(data[0]), len(data)))
    for i in range(len(data)):
        for j in range(len(data[i])):
            loop = data[i][j]
            PETs[j][i] = loop.pet
    return PETs

def count_anchor_frequency_in_a_sample(loopinfos):
    """
    Count the frequency of each anchor in a sample
    """
    anchor_frequency = {}
    for loopinfo in loopinfos:
        for anchor in (loopinfo.loop.anchor1.anchor_id, loopinfo.loop.anchor2.anchor_id):
            if anchor == 302:
                print(loopinfo)
            if anchor in anchor_frequency:
                anchor_frequency[anchor] += 1
            else:
                anchor_frequency[anchor] = 1
    return anchor_frequency

def print_counts(loopinfos_matrix, path):
    counts_matrix = np.zeros((len(loopinfos_matrix), len(loopinfos_matrix[0])))
    for i in range(len(loopinfos_matrix)):
        for j in range(len(loopinfos_matrix[i])):
            counts_matrix[i][j] = loopinfos_matrix[i][j].pet
    counts_matrix = counts_matrix.T
    np.savetxt(path, counts_matrix, delimiter=",", fmt="%d")


def print_my_dci_reuslt(data, pvalue, fdr, output_path):
    with open(f"{output_path}/my_dci_result.csv", "w") as f:
        for i in range(len(pvalue)):
            anchor1_start, anchor1_end = data[0][i].loop.anchor1.start, data[0][i].loop.anchor1.end
            anchor2_start, anchor2_end = data[0][i].loop.anchor2.start, data[0][i].loop.anchor2.end
            f.write(f"{data[0][i].loop.chr}\t{anchor1_start}\t{anchor1_end}\t{anchor2_start}\t{anchor2_end}\t{pvalue[i]}\t{fdr[i]}\n")

def print_my_dci_data_in_narrowPeak_format(data, samples, output_path):
    for i in range(len(data)):
        with open(f"{output_path}/{samples[i]}_dci_data.narrowPeak", "w") as f:
            for j in range(len(data[i])):
                anchor1_start, anchor1_end = data[i][j].loop.anchor1.start, data[i][j].loop.anchor1.end
                anchor2_start, anchor2_end = data[i][j].loop.anchor2.start, data[i][j].loop.anchor2.end
                f.write(f"{data[i][j].loop.chr}\t{anchor1_start}\t{anchor2_end}\t{data[i][j].get_loop_id()}\t1000\t.\t{data[i][j].pet}\t{data[i][j].pvalue}\t{data[i][j].fdr}\t-1\n")
    
    # with open(f"{output_path}/peak_list_for_idr.bed", "w") as f:
    #     for j in range(len(data[i])):
    #         anchor1_start, anchor1_end = data[i][j].loop.anchor1.start, data[i][j].loop.anchor1.end
    #         anchor2_start, anchor2_end = data[i][j].loop.anchor2.start, data[i][j].loop.anchor2.end
    #         f.write(f"{data[i][j].loop.chr}\t{anchor1_start}\t{anchor2_end}\t{data[i][j].get_loop_id()}\n")
    

def get_sim_counts(suffix):
    ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    SIM_COUNTS_DIR = os.path.join(ROOT_DIR, "data", "simulated_data", f"sim_counts.{suffix}.csv")
    df = pd.read_csv(SIM_COUNTS_DIR, index_col=0)
    df = df.to_numpy(dtype=int)
    # print(df.shape)
    return df

def print_my_sim_result(pvalue, fdr, output_path):
    with open(output_path, "w") as f:
        for i in range(len(pvalue)):
            f.write(f"{pvalue[i]}\t{fdr[i]}\n")

