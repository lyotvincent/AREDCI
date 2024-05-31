from preprocessing.quality_controller import quality_control
from preprocessing.ma_normalizer import normalize as ma_normalize
from preprocessing.md_normalizer import normalize as md_normalize
from preprocessing.tmm_normalizer import tmm_normalize
from preprocessing.ssrd_normalizer import normalize as ssrd_normalize
from reproducibility.reproducibility import compute_reproducibility, compute_reproducibility_chr
from reproducibility.process_idr import gen_idr_for_interactions
from dci.kde_dci import run_dci_by_kde_for_two_groups
from dci.neighbor_dci import run_dci_by_neighbors
from dci.md_dci import run_dci_by_md_for_two_groups
from dci.ssrd_dci import run_dci_by_ssrd_for_two_groups
from dci.dci_utils import compute_difference_between_group_mean
from utils.utils import get_PETs_from_loopinfo, count_anchor_frequency_in_a_sample, print_counts, print_my_dci_reuslt, print_my_dci_data_in_narrowPeak_format, get_sim_counts, print_my_sim_result

import numpy as np
import os

def run_dci(files,
            samples,
            groups,
            blacklist_path,
            gap,
            remove_loops_in_blacklist,
            remove_self_ligated_loops,
            fdr_threshold,
            loop_len_threshold,
            intra_only,
            chr_filter,
            pet_threshold,
            norm_method,
            fig_demand,
            reproducibility_checkbox,
            idr_checkbox,
            dci_method,
            output_path):

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    # * Quality control
    if not os.path.exists(os.path.join(output_path, "loops_filtered_in_QC")):
        os.mkdir(os.path.join(output_path, "loops_filtered_in_QC"))
    data = quality_control(files=files,
                           samples=samples,
                           groups=groups,
                           blacklist_path=blacklist_path,
                           gap=gap,
                           remove_loops_in_blacklist=remove_loops_in_blacklist,
                           remove_self_ligated_loops=remove_self_ligated_loops,
                           fdr_threshold=fdr_threshold,
                           loop_len_threshold=loop_len_threshold,
                           intra_only=intra_only,
                           chr_filter=chr_filter,
                           pet_threshold=pet_threshold,
                           output_path=output_path)

    # * Normalization
    if not os.path.exists(f"{output_path}/figures"):
        os.mkdir(f"{output_path}/figures")

    if norm_method == "ma":
        data = ma_normalize(data=data,
                            samples=samples,
                            output_path=output_path,
                            fig_demand=fig_demand)
    elif norm_method == "md":
        data = md_normalize(data=data,
                            samples=samples,
                            output_path=output_path,
                            fig_demand=fig_demand)
    elif norm_method == "tmm":
        data = tmm_normalize(data=data,
                            samples=samples,
                            output_path=output_path,
                            fig_demand=fig_demand)
    elif norm_method == "ssrd":
        data = ssrd_normalize(data=data,
                            samples=samples,
                            output_path=output_path,
                            fig_demand=fig_demand)


    # * Reproducibility
    ## * reproducibility for interactions
    if idr_checkbox:
        if not os.path.exists(os.path.join(output_path, "reproducibility")):
            os.mkdir(os.path.join(output_path, "reproducibility"))
        print_my_dci_data_in_narrowPeak_format(data, samples, output_path)
        gen_idr_for_interactions(samples, output_path)
    ## * sample-level reproducibility
    if reproducibility_checkbox:
        BLOCK_SIZE = 3E6 # selfish default: 3,000,000, 3E6。对于人类，3亿bp大约是染色体的总长度，所以40kb的分辨率大约是7500个bin。所以block=7500/100=75bin=75*40kb=3,000,000bp.
        # BLOCK_SIZE = 5E6
        # the weight of reproducibility score, selfish default is 5 for 40kb resolution Hi-C data. 
        # 但是我这个是ChIA-PET，数据丰度不一样，所以这个权值没办法设成一样的，并且需要设的比selfish小很多。?为什么我写要小很多？
        # C = 1E-5 
        C = 1E-3
        # compute_reproducibility(data=data,
        #                         samples=samples,
        #                         block_size=BLOCK_SIZE,
        #                         stride=BLOCK_SIZE//2,
        #                         c=C,
        #                         output_path=output_path)
        reproducibility_result = compute_reproducibility_chr(data=data,
                                    samples=samples,
                                    block_size=BLOCK_SIZE,
                                    stride=BLOCK_SIZE//2,
                                    c=C,
                                    output_path=output_path)

    # * DCI
    # pets_matrix = get_sim_counts("outlier_0.5")
    # * DCI by KDE
    if dci_method == "kde":
        pets_matrix = get_PETs_from_loopinfo(data)
        pvalues_2, fdrs_2 = run_dci_by_kde_for_two_groups(pets_matrix=pets_matrix, groups=groups)
        print_my_dci_reuslt(data, pvalues_2, fdrs_2, output_path)
    # * DCI by neighbors
    # 在k562和mcf7两个group中, anchor的数量是11974，anchor的频率是1-69
    elif dci_method == "neighbor":
        R = 5
        pvalues_3, fdrs_3, has_neighbor = run_dci_by_neighbors(loopinfos=data, groups=groups, r=R)
        print_my_dci_reuslt(data, pvalues_3, fdrs_3, output_path)
    # * DCI by MD plot
    elif dci_method == "md":
        pets_matrix = get_PETs_from_loopinfo(data)
        pvalues_md, fdrs_md = run_dci_by_md_for_two_groups(pets_matrix=pets_matrix,
                                                        samples=samples,
                                                        groups=groups,
                                                        a_threshold=5,
                                                        should_fit=True)
        print_my_dci_reuslt(data, pvalues_md, fdrs_md, output_path)
    # * DCI by SSRD
    elif dci_method == "ssrd":
        pets_matrix = get_PETs_from_loopinfo(data)
        pvalues_ssrd, fdrs_ssrd = run_dci_by_ssrd_for_two_groups(pets_matrix=pets_matrix, groups=groups, should_fit=True)
        print_my_dci_reuslt(data, pvalues_ssrd, fdrs_ssrd, output_path)

    # * DCI combined
    # print(f"<FUNC> DCI combined")
    # pvalues_merged = list()
    # fdrs_merged = list()
    # from scipy.stats import combine_pvalues
    # from statsmodels.stats.multitest import fdrcorrection
    # for i in range(len(has_neighbor)):
    #     p_values = [pvalues_2[i], pvalues_3[i], pvalues_md[i]]
    #     _, merged_p_value = combine_pvalues(p_values, method='mudholkar_george')
    #     if np.nan in p_values or np.isnan(merged_p_value):
    #         print("[WARN]", p_values, merged_p_value)
    #     pvalues_merged.append(merged_p_value)
    # pvalues_merged = np.array(pvalues_merged)
    # # fdrs = np.array(fdrs)
    # print(pvalues_merged.shape)
    # print(pvalues_merged[:5])
    # _, fdrs_merged = fdrcorrection(pvalues_merged, alpha=0.05, method='indep', is_sorted=False)

    # * metrics
    # print_my_dci_reuslt(data, pvalues_2, fdrs_2, output_path)
    # print_my_dci_reuslt(data, pvalues_3, fdrs_3, output_path)
    # print_my_dci_reuslt(data, pvalues_md, fdrs_md, output_path)
    # print_my_dci_reuslt(data, pvalues_ssrd, fdrs_ssrd, output_path)
    # print_my_dci_reuslt(data, pvalues_merged, fdrs_merged, output_path)

    # print_my_sim_result(pvalues_2, fdrs_2, os.path.join(ROOT_PATH, "test", "sim_test", "my_sim_result.csv"))

    # print("=============================")
    # print(f"has_neighbor num: {sum(has_neighbor)}/{len(has_neighbor)}")
    # print("=============================")
    # count_fdr_and_diff(fdrs_2, pets_matrix, pvalues_2)
    # print("=============================")
    # count_fdr_and_diff(fdrs_3, pets_matrix, pvalues_3)
    # print("=============================")
    # count_fdr_and_diff(fdrs_md, pets_matrix, pvalues_md)
    # print("=============================")
    # count_fdr_and_diff(fdrs_ssrd, pets_matrix, pvalues_ssrd)
    # print("=============================")
    # count_fdr_and_diff(fdrs_merged, pets_matrix, pvalues_merged)

    if reproducibility_checkbox:
        return reproducibility_result
    else:
        return None

def count_fdr_and_diff(fdrs, pets_matrix, pvalues):
    print(f"fdr.shape: {fdrs.shape}")
    print(fdrs[:5])
    print("fdrs < 0.05", sum(fdrs < 0.05))
    # get index when fdr < 0.05
    index_fdr = np.where(fdrs < 0.05)
    diffs = compute_difference_between_group_mean(pets_matrix, groups)
    # get topmax 10 indexes of diffs
    index_dif = np.argsort(diffs)
    # print(f"topmax indexes of diffs: {sorted(index)}")
    # print(f"topmax diffs: {diffs[index]}")
    # print(f"topmax tmm_pets: {tmm_pets[index]}")
    # print(f"topmax fdr: {fdr[index]}")
    # print(f"topmax cumulative_prob: {cumulative_prob[index]}")
    # print(f"topmax fitted_values[index]: {fitted_values[index]}")
    # print(f"topmax fitted_values0[index]: {fitted_values0[index]}")
    index_dif = index_dif[-sum(fdrs < 0.05):]
    print(len(index_dif))
    print(len(index_fdr[0]))
    print("len(set(index_dif).intersection(set(index_fdr[0])))")
    print(len(set(index_dif).intersection(set(index_fdr[0]))))
    print("X30331")
    print(fdrs[30331])
    print(pvalues[30331])
    print(pets_matrix[30331])
    print("X30332")
    print(fdrs[30332])
    print(pvalues[30332])
    print(pets_matrix[30332])
    print("X30333")
    print(fdrs[30333])
    print(pvalues[30333])
    print(pets_matrix[30333])

if __name__ == "__main__":
    np.random.seed(9523)
    ROOT_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    GM12878_REP1 = os.path.join(ROOT_PATH, "data", "ChIA-PET2_result", "GM12878_ChIA-PET2_result", "replicate1", "GM12878.interactions.MICC")
    GM12878_REP2 = os.path.join(ROOT_PATH, "data", "ChIA-PET2_result", "GM12878_ChIA-PET2_result", "replicate2", "GM12878.interactions.MICC")
    K562_REP1 = os.path.join(ROOT_PATH, "data", "ChIA-PET2_result", "K562_ChIA-PET2_result", "replicate1", "K562.interactions.MICC")
    K562_REP2 = os.path.join(ROOT_PATH, "data", "ChIA-PET2_result", "K562_ChIA-PET2_result", "replicate2", "K562.interactions.MICC")
    MCF7_REP1 = os.path.join(ROOT_PATH, "data", "ChIA-PET2_result", "MCF-7_ChIA-PET2_result", "replicate1", "MCF-7.interactions.MICC")
    MCF7_REP2 = os.path.join(ROOT_PATH, "data", "ChIA-PET2_result", "MCF-7_ChIA-PET2_result", "replicate2", "MCF-7.interactions.MICC")
    files = [K562_REP1, K562_REP2, MCF7_REP1, MCF7_REP2]
    samples = ["k562_rep1", "k562_rep2", "mcf7_rep1", "mcf7_rep2"]
    groups = ["k562", "k562", "mcf7", "mcf7"]

    ESC_REP1 = os.path.join(ROOT_PATH, "test", "mouse_neuronal_differentiation", "ESC_rep1.csv")
    ESC_REP2 = os.path.join(ROOT_PATH, "test", "mouse_neuronal_differentiation", "ESC_rep2.csv")
    NSC_REP1 = os.path.join(ROOT_PATH, "test", "mouse_neuronal_differentiation", "NSC_rep1.csv")
    NSC_REP2 = os.path.join(ROOT_PATH, "test", "mouse_neuronal_differentiation", "NSC_rep2.csv")
    NPC_REP1 = os.path.join(ROOT_PATH, "test", "mouse_neuronal_differentiation", "NPC_rep1.csv")
    # files = [ESC_REP1, ESC_REP2, NSC_REP1, NSC_REP2]
    # samples = ["ESC_rep1", "ESC_rep2", "NSC_rep1", "NSC_rep2"]
    # groups = ["ESC", "ESC", "NSC", "NSC"]
    # files = [ESC_REP1, ESC_REP2, NPC_REP1]
    # samples = ["ESC_rep1", "ESC_rep2", "NPC_rep1"]
    # groups = ["ESC", "ESC", "NPC"]
    # files = [NSC_REP1, NSC_REP2, NPC_REP1]
    # samples = ["NSC_rep1", "NSC_rep2", "NPC_rep1"]
    # groups = ["NSC", "NSC", "NPC"]

    # * default parameters
    gap = 500
    blacklist_path = ROOT_PATH+"/data/blacklist/hg38-blacklist.v2.bed"
    remove_loops_in_blacklist=True
    remove_self_ligated_loops=True
    fdr_threshold=0.01
    loop_len_threshold=5000
    intra_only=True
    chr_filter=["chrM", "chrX", "chrUn"]
    pet_threshold=5
    norm_method='tmm'
    fig_demand=True
    reproducibility_checkbox=True
    idr_checkbox=True
    dci_method='kde'
    output_path = ROOT_PATH+"/data/my_test_data/"
    # output_path = os.path.join(ROOT_PATH, "test", "mouse_neuronal_differentiation", "ESC_vs_NSC_result")
    # output_path = os.path.join(ROOT_PATH, "test", "mouse_neuronal_differentiation", "ESC_vs_NPC_result")
    # output_path = os.path.join(ROOT_PATH, "test", "mouse_neuronal_differentiation", "NSC_vs_NPC_result")

    # run_dci(files=files,
    #         samples=samples,
    #         groups=groups,
    #         blacklist_path=blacklist_path,
    #         gap=gap,
    #         remove_loops_in_blacklist=remove_loops_in_blacklist,
    #         remove_self_ligated_loops=remove_self_ligated_loops,
    #         fdr_threshold=fdr_threshold,
    #         loop_len_threshold=loop_len_threshold,
    #         intra_only=intra_only,
    #         chr_filter=chr_filter,
    #         pet_threshold=pet_threshold,
    #         output_path=output_path)

    run_dci(files=files,
            samples=samples,
            groups=groups,
            blacklist_path=blacklist_path,
            gap=gap,
            remove_loops_in_blacklist=remove_loops_in_blacklist,
            remove_self_ligated_loops=remove_self_ligated_loops,
            fdr_threshold=fdr_threshold,
            loop_len_threshold=loop_len_threshold,
            intra_only=intra_only,
            chr_filter=chr_filter,
            pet_threshold=pet_threshold,
            norm_method=norm_method,
            fig_demand=fig_demand,
            reproducibility_checkbox=reproducibility_checkbox,
            idr_checkbox=idr_checkbox,
            dci_method=dci_method,
            output_path=output_path)

