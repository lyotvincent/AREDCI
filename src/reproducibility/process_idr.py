
import os, math
import numpy as np
from scipy.stats import rankdata
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import subprocess

def draw_rep_pet_and_p(file_path, output_path):
    """
    Draw pet and pvalue for rep
    """
    pet_list = list()
    pvalue_list = list()
    with open(file_path, "r") as f:
        line = f.readline()
        while line:
            fields = line.strip().split()
            pet = float(fields[6])
            pvalue = float(fields[7])
            if pvalue == math.inf:
                line = f.readline()
                continue
            pet_list.append(pet)
            pvalue_list.append(pvalue)
            line = f.readline()

    pet_sort_index = np.argsort(pet_list)[::-1]
    pet_rank = np.zeros_like(pet_sort_index)
    pet_rank[pet_sort_index] = np.arange(1, len(pet_sort_index)+1)

    pvalue_sort_index = np.argsort(pvalue_list)[::-1]
    pvalue_rank = np.zeros_like(pvalue_sort_index)
    pvalue_rank[pvalue_sort_index] = np.arange(1, len(pvalue_sort_index)+1)

    import matplotlib.pyplot as plt
    plt.figure(dpi=300, figsize=(8, 6))
    FONTSIZE = 12
    plt.rc("font", family="Times New Roman")
    params = {"axes.titlesize": FONTSIZE,
              "legend.fontsize": FONTSIZE,
              "axes.labelsize": FONTSIZE,
              "xtick.labelsize": FONTSIZE,
              "ytick.labelsize": FONTSIZE,
              "figure.titlesize": FONTSIZE,
              "font.size": FONTSIZE}
    plt.rcParams.update(params)

    plt.scatter(pet_rank, pvalue_rank, s=2, alpha=0.75)
    plt.xlabel("PET in replicate")
    plt.ylabel("P value in replicate")

    plt.tight_layout()
    plt.savefig(output_path, format="png", bbox_inches="tight")
    # plt.show()
    plt.close()


def load_dci_data(file_path, rank_indicator="pvalue"):
    """
    Load dci data from file, pet, pvalue, fdr
    """
    if rank_indicator not in ["pet", "pvalue", "fdr"]:
        raise ValueError("rank_indicator should be one of ['pet', 'pvalue', 'fdr']")
    if rank_indicator == "pet":
        rank_indicator = 6
    elif rank_indicator == "pvalue":
        rank_indicator = 7
    elif rank_indicator == "fdr":
        rank_indicator = 8
    # we just name the indicator is pvalue
    pvalues = list()
    pos_list = list()
    with open(file_path, "r") as f:
        line = f.readline()
        while line:
            fields = line.strip().split()
            pvalues.append(float(fields[rank_indicator]))
            pos_list.append("_".join(fields[:3]))
            line = f.readline()
    return np.array(pvalues), pos_list

def load_idrs(idr_path, pos_list):
    with open(idr_path, "r") as f:
        line = f.readline()
        idr_pos_dict = dict()
        while line:
            fields = line.strip().split()
            idr_pos_dict["_".join(fields[:3])] = np.power(10, -float(fields[11]))
            line = f.readline() 
    idrs = list()
    for pos in pos_list:
        idrs.append(idr_pos_dict[pos])
    return np.array(idrs)

def get_rank_for_two_reps(rep1_pvalues, rep2_pvalues, idrs, rank_indicator="pvalue"):
    """
    Get rank for two reps
    """
    # filter pvalue==-1 in k562_rep1_pvalues or k562_rep2_pvalues
    if rank_indicator not in ["pet", "pvalue", "fdr"]:
        raise ValueError("rank_indicator should be one of ['pet', 'pvalue', 'fdr']")
    if rank_indicator == "pet":
        indexes1 = np.where(rep1_pvalues != 0)
        indexes2 = np.where(rep2_pvalues != 0)
    elif rank_indicator == "pvalue":
        indexes1 = np.where(rep1_pvalues != math.inf)
        indexes2 = np.where(rep2_pvalues != math.inf)
    indexes = np.intersect1d(indexes1, indexes2)
    rep1_pvalues = rep1_pvalues[indexes]
    rep2_pvalues = rep2_pvalues[indexes]
    idrs = idrs[indexes]

    # ! old code, the same value will not be the same rank
    # 获得k562_rep1的每个rank_indicator在k562_rep1中的排序，pvalue的话，升序（从小到大）排序
    # rep1_pvalues_sort_index = np.argsort(rep1_pvalues)
    # if rank_indicator == "pet": # pet的话，降序（从大到小）排序
    #     rep1_pvalues_sort_index = rep1_pvalues_sort_index[::-1]
    # rep1_pvalues_rank = np.zeros_like(rep1_pvalues_sort_index)
    # rep1_pvalues_rank[rep1_pvalues_sort_index] = np.arange(1, len(rep1_pvalues_sort_index)+1)
    # ! new code, the same value will be the same rank
    if rank_indicator == "pvalue": # pvalue的话，升序（从小到大）排序
        rep1_pvalues_rank = rankdata(rep1_pvalues, method='min')
    elif rank_indicator == "pet": # pet的话，降序（从大到小）排序
        rep1_pvalues_rank = rankdata(-rep1_pvalues, method='min')

    # 获得k562_rep2的每个rank_indicator在k562_rep2中的排序，pvalue的话，升序（从小到大）排序
    # rep2_pvalues_sort_index = np.argsort(rep2_pvalues)
    # if rank_indicator == "pet": # pet的话，降序（从大到小）排序
    #     rep2_pvalues_sort_index = rep2_pvalues_sort_index[::-1]
    # rep2_pvalues_rank = np.zeros_like(rep2_pvalues_sort_index)
    # rep2_pvalues_rank[rep2_pvalues_sort_index] = np.arange(1, len(rep2_pvalues_sort_index)+1)
    if rank_indicator == "pvalue": # pvalue的话，升序（从小到大）排序
        rep2_pvalues_rank = rankdata(rep2_pvalues, method='min')
    elif rank_indicator == "pet": # pet的话，降序（从大到小）排序
        rep2_pvalues_rank = rankdata(-rep2_pvalues, method='min')

    return rep1_pvalues_rank, rep2_pvalues_rank, idrs

def draw_rank_scatter(rep1_pvalues_rank, rep2_pvalues_rank, idrs, output_path):
    """
    Draw rank scatter
    """
    # normalize idrs to [0, 1]
    # idrs = np.log(idrs + 1.001)
    # idrs = (idrs - np.min(idrs)) / (np.max(idrs) - np.min(idrs))

    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker

    plt.figure(dpi=300, figsize=(8, 8))
    FONTSIZE = 12
    plt.rc("font", family="Times New Roman")
    params = {"axes.titlesize": FONTSIZE,
              "legend.fontsize": FONTSIZE,
              "axes.labelsize": 40,
              "xtick.labelsize": 32,
              "ytick.labelsize": 32,
              "figure.titlesize": FONTSIZE,
              "font.size": FONTSIZE}
    plt.rcParams.update(params)
    ax = plt.gca()

    # color from idrs
    cmap = plt.get_cmap("YlOrRd")
    colors = cmap(np.arange(cmap.N))
    # 从黄色开始，即从颜色列表的1/4开始
    start = int(cmap.N / 4)
    new_cmap = ListedColormap(colors[start:,:-1])

    # # 定义颜色列表
    # colors = ["red", "orange", "yellow", "green", "blue", "indigo", "violet"]
    # # 创建自定义颜色映射
    # new_cmap = LinearSegmentedColormap.from_list("mycmap", colors)

    # 使用颜色映射
    # plt.pcolormesh(np.random.rand(10,10), cmap=new_cmap)
    # plt.colorbar(shrink=0.2)
    # clist = [abs(m) for m in idrs]

    plt.scatter(rep1_pvalues_rank, rep2_pvalues_rank, s=3, alpha=0.75, c=idrs, cmap=new_cmap)

    # 创建一个自定义的刻度格式器
    def custom_formatter(x, pos):
        return f"{int(x/1000):.0f}" if x % 1000 == 0 else ''

    # 应用自定义刻度格式器到y轴
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(custom_formatter))
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(custom_formatter))


    plt.xlabel("rank in replicate1 (kilo)")
    plt.ylabel("rank in replicate2 (kilo)")
    # plt.xticks([])
    # plt.yticks([])
    plt.yticks(rotation=90)

    plt.tight_layout()
    plt.savefig(output_path+'.png', format="png", bbox_inches="tight")
    # plt.show()
    plt.close()

def process_idr_for_two_reps(rep1_path, rep2_path, idr_path, fig_path, rank_indicator="pvalue"):
    """
    Process idr for two reps
    """
    rep1_pvalues, pos_list = load_dci_data(rep1_path, rank_indicator=rank_indicator)
    rep2_pvalues, _ = load_dci_data(rep2_path, rank_indicator=rank_indicator)
    idrs = load_idrs(idr_path, pos_list)
    print(max(idrs), min(idrs))
    rep1_pvalues_rank, rep2_pvalues_rank, idrs = get_rank_for_two_reps(rep1_pvalues, rep2_pvalues, idrs, rank_indicator=rank_indicator)
    print(max(idrs), min(idrs))
    draw_rank_scatter(rep1_pvalues_rank, rep2_pvalues_rank, idrs, fig_path)

def main():
    ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    # print(ROOT_DIR)
    K562_REP1_DCI_PATH = os.path.join(ROOT_DIR, "data", "my_test_data", "K562_rep1_dci_data.narrowPeak")
    K562_REP2_DCI_PATH = os.path.join(ROOT_DIR, "data", "my_test_data", "K562_rep2_dci_data.narrowPeak")
    MCF7_REP1_DCI_PATH = os.path.join(ROOT_DIR, "data", "my_test_data", "MCF7_rep1_dci_data.narrowPeak")
    MCF7_REP2_DCI_PATH = os.path.join(ROOT_DIR, "data", "my_test_data", "MCF7_rep2_dci_data.narrowPeak")

    # * draw consistence within a replicate
    K562_REP1_FIG_PATH = os.path.join(ROOT_DIR, "data", "my_test_data", "reproducibility", "k562_rep1.consistence.png")
    K562_REP2_FIG_PATH = os.path.join(ROOT_DIR, "data", "my_test_data", "reproducibility", "k562_rep2.consistence.png")
    MCF7_REP1_FIG_PATH = os.path.join(ROOT_DIR, "data", "my_test_data", "reproducibility", "mcf7_rep1.consistence.png")
    MCF7_REP2_FIG_PATH = os.path.join(ROOT_DIR, "data", "my_test_data", "reproducibility", "mcf7_rep2.consistence.png")
    draw_rep_pet_and_p(K562_REP1_DCI_PATH, K562_REP1_FIG_PATH)
    draw_rep_pet_and_p(K562_REP2_DCI_PATH, K562_REP2_FIG_PATH)
    draw_rep_pet_and_p(MCF7_REP1_DCI_PATH, MCF7_REP1_FIG_PATH)
    draw_rep_pet_and_p(MCF7_REP2_DCI_PATH, MCF7_REP2_FIG_PATH)


    K562_REP1_K562_REP2_IDR_PATH = os.path.join(ROOT_DIR, "data", "my_test_data", "reproducibility", "k562_rep1_K562_rep2_idr")
    K562_REP1_MCF7_REP1_IDR_PATH = os.path.join(ROOT_DIR, "data", "my_test_data", "reproducibility", "k562_rep1_mcf7_rep1_idr")
    K562_REP1_MCF7_REP2_IDR_PATH = os.path.join(ROOT_DIR, "data", "my_test_data", "reproducibility", "k562_rep1_mcf7_rep2_idr")
    K562_REP2_MCF7_REP1_IDR_PATH = os.path.join(ROOT_DIR, "data", "my_test_data", "reproducibility", "k562_rep2_mcf7_rep1_idr")
    K562_REP2_MCF7_REP2_IDR_PATH = os.path.join(ROOT_DIR, "data", "my_test_data", "reproducibility", "k562_rep2_mcf7_rep2_idr")
    MCF7_REP1_MCF7_REP2_IDR_PATH = os.path.join(ROOT_DIR, "data", "my_test_data", "reproducibility", "mcf7_rep1_mcf7_rep2_idr")

    # * draw pvalue axis
    K562_REP1_K562_REP2_FIG_PATH = os.path.join(ROOT_DIR, "data", "my_test_data", "reproducibility", "k562_rep1_K562_rep2_idr.rank")
    K562_REP1_MCF7_REP1_FIG_PATH = os.path.join(ROOT_DIR, "data", "my_test_data", "reproducibility", "k562_rep1_mcf7_rep1_idr.rank")
    K562_REP1_MCF7_REP2_FIG_PATH = os.path.join(ROOT_DIR, "data", "my_test_data", "reproducibility", "k562_rep1_mcf7_rep2_idr.rank")
    K562_REP2_MCF7_REP1_FIG_PATH = os.path.join(ROOT_DIR, "data", "my_test_data", "reproducibility", "k562_rep2_mcf7_rep1_idr.rank")
    K562_REP2_MCF7_REP2_FIG_PATH = os.path.join(ROOT_DIR, "data", "my_test_data", "reproducibility", "k562_rep2_mcf7_rep2_idr.rank")
    MCF7_REP1_MCF7_REP2_FIG_PATH = os.path.join(ROOT_DIR, "data", "my_test_data", "reproducibility", "mcf7_rep1_mcf7_rep2_idr.rank")

    process_idr_for_two_reps(K562_REP1_DCI_PATH, K562_REP2_DCI_PATH, K562_REP1_K562_REP2_IDR_PATH, K562_REP1_K562_REP2_FIG_PATH, rank_indicator="pvalue")
    process_idr_for_two_reps(K562_REP1_DCI_PATH, MCF7_REP1_DCI_PATH, K562_REP1_MCF7_REP1_IDR_PATH, K562_REP1_MCF7_REP1_FIG_PATH, rank_indicator="pvalue")
    process_idr_for_two_reps(K562_REP1_DCI_PATH, MCF7_REP2_DCI_PATH, K562_REP1_MCF7_REP2_IDR_PATH, K562_REP1_MCF7_REP2_FIG_PATH, rank_indicator="pvalue")
    process_idr_for_two_reps(K562_REP2_DCI_PATH, MCF7_REP1_DCI_PATH, K562_REP2_MCF7_REP1_IDR_PATH, K562_REP2_MCF7_REP1_FIG_PATH, rank_indicator="pvalue")
    process_idr_for_two_reps(K562_REP2_DCI_PATH, MCF7_REP2_DCI_PATH, K562_REP2_MCF7_REP2_IDR_PATH, K562_REP2_MCF7_REP2_FIG_PATH, rank_indicator="pvalue")
    process_idr_for_two_reps(MCF7_REP1_DCI_PATH, MCF7_REP2_DCI_PATH, MCF7_REP1_MCF7_REP2_IDR_PATH, MCF7_REP1_MCF7_REP2_FIG_PATH, rank_indicator="pvalue")

    # * draw pet axis
    K562_REP1_K562_REP2_FIG_PATH = os.path.join(ROOT_DIR, "data", "my_test_data", "reproducibility", "k562_rep1_K562_rep2_idr.value")
    K562_REP1_MCF7_REP1_FIG_PATH = os.path.join(ROOT_DIR, "data", "my_test_data", "reproducibility", "k562_rep1_mcf7_rep1_idr.value")
    K562_REP1_MCF7_REP2_FIG_PATH = os.path.join(ROOT_DIR, "data", "my_test_data", "reproducibility", "k562_rep1_mcf7_rep2_idr.value")
    K562_REP2_MCF7_REP1_FIG_PATH = os.path.join(ROOT_DIR, "data", "my_test_data", "reproducibility", "k562_rep2_mcf7_rep1_idr.value")
    K562_REP2_MCF7_REP2_FIG_PATH = os.path.join(ROOT_DIR, "data", "my_test_data", "reproducibility", "k562_rep2_mcf7_rep2_idr.value")
    MCF7_REP1_MCF7_REP2_FIG_PATH = os.path.join(ROOT_DIR, "data", "my_test_data", "reproducibility", "mcf7_rep1_mcf7_rep2_idr.value")

    process_idr_for_two_reps(K562_REP1_DCI_PATH, K562_REP2_DCI_PATH, K562_REP1_K562_REP2_IDR_PATH, K562_REP1_K562_REP2_FIG_PATH, rank_indicator="pet")
    process_idr_for_two_reps(K562_REP1_DCI_PATH, MCF7_REP1_DCI_PATH, K562_REP1_MCF7_REP1_IDR_PATH, K562_REP1_MCF7_REP1_FIG_PATH, rank_indicator="pet")
    process_idr_for_two_reps(K562_REP1_DCI_PATH, MCF7_REP2_DCI_PATH, K562_REP1_MCF7_REP2_IDR_PATH, K562_REP1_MCF7_REP2_FIG_PATH, rank_indicator="pet")
    process_idr_for_two_reps(K562_REP2_DCI_PATH, MCF7_REP1_DCI_PATH, K562_REP2_MCF7_REP1_IDR_PATH, K562_REP2_MCF7_REP1_FIG_PATH, rank_indicator="pet")
    process_idr_for_two_reps(K562_REP2_DCI_PATH, MCF7_REP2_DCI_PATH, K562_REP2_MCF7_REP2_IDR_PATH, K562_REP2_MCF7_REP2_FIG_PATH, rank_indicator="pet")
    process_idr_for_two_reps(MCF7_REP1_DCI_PATH, MCF7_REP2_DCI_PATH, MCF7_REP1_MCF7_REP2_IDR_PATH, MCF7_REP1_MCF7_REP2_FIG_PATH, rank_indicator="pet")

def gen_idr_for_interactions(samples, output_path):
    # idr --samples ./data/my_test_data/k562_rep1_dci_data.narrowPeak ./data/my_test_data/k562_rep2_dci_data.narrowPeak --peak-list ./data/my_test_data/k562_rep1_dci_data.narrowPeak --input-file-type narrowPeak --rank signal.value --output-file ./data/my_test_data/reproducibility/k562_rep1_k562_rep2_idr --plot --log-output-file ./data/my_test_data/reproducibility/k562_rep1_k562_rep2_idr.log
    reproducibility_path = os.path.join(output_path, "reproducibility")
    for spl1 in samples:
        for spl2 in samples:
            if spl1 == spl2:
                continue
            spl1_path = os.path.join(output_path, f"{spl1}_dci_data.narrowPeak")
            spl2_path = os.path.join(output_path, f"{spl2}_dci_data.narrowPeak")
            idr_path = os.path.join(reproducibility_path, f"{spl1}_{spl2}_idr")
            print("idr --samples", spl1_path, spl2_path, "--peak-list", spl1_path, "--input-file-type narrowPeak --rank signal.value --output-file", idr_path, "--plot --log-output-file", idr_path+".log")
            subprocess.run(f"idr --samples {spl1_path} {spl2_path} --peak-list {spl1_path} --input-file-type narrowPeak --rank signal.value --output-file {idr_path} --plot --log-output-file {idr_path}.log", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8")


if __name__ == "__main__":
    main()

    # a = np.array([33, 22, 55, 22, 11, 66])
    # k562_rep1_pvalues_sort_index = np.argsort(a)
    # k562_rep1_pvalues_rank = np.zeros_like(k562_rep1_pvalues_sort_index)
    # k562_rep1_pvalues_rank[k562_rep1_pvalues_sort_index] = np.arange(1, len(k562_rep1_pvalues_sort_index)+1)

    # k562_rep1_pvalues_rank = rankdata(a, method='min')
    # print(k562_rep1_pvalues_rank)
    # k562_rep1_pvalues_rank = rankdata(-a, method='min')
    # print(k562_rep1_pvalues_rank)
