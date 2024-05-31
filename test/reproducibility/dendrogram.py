
import os, sys
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt

ROOT_PATH = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(ROOT_PATH, "src"))
from preprocessing.quality_controller import quality_control
from reproducibility.reproducibility import compute_reproducibility_chr

def run_reproducibility(files,
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
                        output_path):

    # * Quality control
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

    # * Reproducibility
    BLOCK_SIZE = 3E6 # selfish default: 3,000,000, 3E6。对于人类，3亿bp大约是染色体的总长度，所以40kb的分辨率大约是7500个bin。所以block=7500/100=75bin=75*40kb=3,000,000bp.
    # the weight of reproducibility score, selfish default is 5 for 40kb resolution Hi-C data. 
    # 但是我这个是ChIA-PET，数据丰度不一样，所以这个权值没办法设成一样的，并且需要设的比selfish小很多。?为什么我写要小很多？
    # C = 1E-5 
    C = 1E-3
    reproducibility_list = compute_reproducibility_chr(data=data,
                                samples=samples,
                                block_size=BLOCK_SIZE,
                                stride=BLOCK_SIZE//2,
                                c=C,
                                output_path=output_path)
    return reproducibility_list[0][1]

if __name__ == "__main__":
    np.random.seed(9523)

    GM12878_REP1 = os.path.join(ROOT_PATH, "data", "ChIA-PET2_result", "GM12878_ChIA-PET2_result", "replicate1", "GM12878.interactions.MICC")
    GM12878_REP2 = os.path.join(ROOT_PATH, "data", "ChIA-PET2_result", "GM12878_ChIA-PET2_result", "replicate2", "GM12878.interactions.MICC")
    K562_REP1 = os.path.join(ROOT_PATH, "data", "ChIA-PET2_result", "K562_ChIA-PET2_result", "replicate1", "K562.interactions.MICC")
    K562_REP2 = os.path.join(ROOT_PATH, "data", "ChIA-PET2_result", "K562_ChIA-PET2_result", "replicate2", "K562.interactions.MICC")
    MCF7_REP1 = os.path.join(ROOT_PATH, "data", "ChIA-PET2_result", "MCF-7_ChIA-PET2_result", "replicate1", "MCF-7.interactions.MICC")
    MCF7_REP2 = os.path.join(ROOT_PATH, "data", "ChIA-PET2_result", "MCF-7_ChIA-PET2_result", "replicate2", "MCF-7.interactions.MICC")
    HCT116_REP1 = os.path.join(ROOT_PATH, "data", "ChIA-PET2_result", "HCT116_ChIA-PET2_result", "replicate1", "HCT116.interactions.MICC")
    HCT116_REP2 = os.path.join(ROOT_PATH, "data", "ChIA-PET2_result", "HCT116_ChIA-PET2_result", "replicate2", "HCT116.interactions.MICC")
    ESC_REP1 = os.path.join(ROOT_PATH, "data", "ChIA-PET2_result", "ESC_ChIA-PET2_result", "replicate1", "ESC.interactions.MICC")
    ESC_REP2 = os.path.join(ROOT_PATH, "data", "ChIA-PET2_result", "ESC_ChIA-PET2_result", "replicate2", "ESC.interactions.MICC")
    SPLEEN_REP1 = os.path.join(ROOT_PATH, "data", "ChIA-PET2_result", "spleen_ChIA-PET2_result", "replicate1", "spleen.interactions.MICC")

    files_list = [GM12878_REP1, GM12878_REP2, K562_REP1, K562_REP2, MCF7_REP1, MCF7_REP2, HCT116_REP1, HCT116_REP2, ESC_REP1, ESC_REP2, SPLEEN_REP1]
    samples = ["gm12878_rep1", "gm12878_rep2", "k562_rep1", "k562_rep2", "mcf7_rep1", "mcf7_rep2", "hct116_rep1", "hct116_rep2", "esc_rep1", "esc_rep2", "spleen_rep1"]
    groups = ["gm12878", "gm12878", "k562", "k562", "mcf7", "mcf7", "hct116", "hct116", "esc", "esc", "spleen"]

    files_list = [GM12878_REP1, GM12878_REP2, K562_REP1, K562_REP2, MCF7_REP1, MCF7_REP2, HCT116_REP1, HCT116_REP2, ESC_REP1, ESC_REP2]
    samples = ["gm12878_rep1", "gm12878_rep2", "k562_rep1", "k562_rep2", "mcf7_rep1", "mcf7_rep2", "hct116_rep1", "hct116_rep2", "esc_rep1", "esc_rep2"]
    groups = ["gm12878", "gm12878", "k562", "k562", "mcf7", "mcf7", "hct116", "hct116", "esc", "esc"]


    # * default parameters
    gap = 500
    blacklist_path = ROOT_PATH+"/data/blacklist/hg38-blacklist.v2.bed"
    remove_loops_in_blacklist=False
    remove_self_ligated_loops=True
    fdr_threshold=0.01
    loop_len_threshold=5000
    intra_only=True
    chr_filter=["chrM", "chrX", "chrUn"]
    pet_threshold=5
    output_path = os.path.join(ROOT_PATH, "test", "reproducibility", "temp_output")

    reproducibility_matrix = np.zeros((len(files_list), len(files_list)))
    for i in range(0, len(files_list)):
        for j in range(i+1, len(files_list)):
            reproducibility = run_reproducibility(files=[files_list[i], files_list[j]],
                                                samples=[samples[i], samples[j]],
                                                groups=[groups[i], groups[j]],
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
            reproducibility_matrix[i][j] = reproducibility
            reproducibility_matrix[j][i] = reproducibility

    reproducibility_matrix = 1 - reproducibility_matrix
    np.fill_diagonal(reproducibility_matrix, 0)

    np.save(os.path.join(ROOT_PATH, "test", "reproducibility", "temp_output", "reproducibility_matrix.npy"), reproducibility_matrix)
    reproducibility_matrix = np.load(os.path.join(ROOT_PATH, "test", "reproducibility", "temp_output", "reproducibility_matrix.npy"))

    # 数值=距离，所以越大越远，越小越近
    # reproducibility_matrix = np.array([[0.,         0.09569471, 0.18749197, 0.16470426],
    #                                 [0.09569471, 0.,         0.19149057, 0.1651398 ],
    #                                 [0.18749197, 0.19149057, 0.,         0.11607364],
    #                                 [0.16470426, 0.1651398,  0.11607364, 0.        ]])
    print(reproducibility_matrix)

    # 使用scipy的linkage函数来计算层次聚类
    linked = linkage(reproducibility_matrix, 'single')
    print(linked)

    # 使用scipy的dendrogram函数来绘制树状图
    plt.figure(dpi=300, figsize=(8, 6))
    FONTSIZE = 10
    plt.rc("font", family="Times New Roman")
    params = {"axes.titlesize": FONTSIZE,
              "legend.fontsize": FONTSIZE,
              "axes.labelsize": FONTSIZE,
              "xtick.labelsize": FONTSIZE,
              "ytick.labelsize": FONTSIZE,
              "figure.titlesize": FONTSIZE,
              "font.size": FONTSIZE}
    plt.rcParams.update(params)
    dendrogram(linked,
                orientation='top',  # 树状图的方向，'top'表示从上到下
                labels=samples,
                distance_sort=True,  # 是否根据距离排序
                show_leaf_counts=False  # 是否显示每个叶节点的计数
                )

    plt.xticks(rotation=45)
    # 显示树状图
    plt.title('Dendrogram using Single Linkage')
    # plt.xlabel('Data Points')
    # plt.ylabel('Distance')
    plt.tight_layout()
    plt.savefig(f"{output_path}/tree.png", format="png", bbox_inches="tight")
    # plt.show()
    plt.close()

