import math, time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import statsmodels.api as sm


def get_common_loops_for_two_samples(loop_infos_1, loop_infos_2, func_depth=0):
    # get the common loops
    time1 = time.time()
    print('\t'*func_depth+"<FUNC> get_common_loops_for_two_samples")
    i = 0
    while i < len(loop_infos_1) and i < len(loop_infos_2):
        loop_id_1 = loop_infos_1[i].get_loop_id()
        loop_id_2 = loop_infos_2[i].get_loop_id()
        if loop_id_1 == loop_id_2:
            i += 1
        elif loop_id_1 < loop_id_2:
            loop_infos_1.pop(i)
        else:
            loop_infos_2.pop(i)
    if len(loop_infos_1) > len(loop_infos_2):
        loop_infos_1 = loop_infos_1[:i]
    else:
        loop_infos_2 = loop_infos_2[:i]
    print('\t'*func_depth+f"<TIME> get_common_loops_for_two_samples cost {time.time()-time1} seconds")
    return loop_infos_1, loop_infos_2

def filter_zero_pet_loops(pets_1, pets_2, func_depth=0):
    time1 = time.time()
    print('\t'*func_depth+"<FUNC> filter_zero_pet_loops")
    pets_1 = np.array(pets_1)
    pets_2 = np.array(pets_2)
    # get index of zero pet loops
    no_zero_pets_index_1 = np.where(pets_1 != 0)[0]
    no_zero_pets_index_2 = np.where(pets_2 != 0)[0]
    no_zero_pets_index = np.intersect1d(no_zero_pets_index_1, no_zero_pets_index_2)
    # remove zero pet loops
    # pets_1 = np.delete(pets_1, zero_pets_index)
    # pets_2 = np.delete(pets_2, zero_pets_index)
    print('\t'*func_depth+f"<TIME> filter_zero_pet_loops cost {time.time()-time1} seconds")
    return no_zero_pets_index

def compute_D_for_two_samples(loop_infos_1, func_depth=0):
    """
    run after get_common_loops_for_two_samples
    Compute the D for two samples, the two samples should have the same set of loops(but may different in PETs)
    D=distance between two interacting regions
    """
    time1 = time.time()
    print('\t'*func_depth+"<FUNC> compute_D_for_two_samples")
    d_list = list()
    for loop_info_1 in loop_infos_1:
        d = loop_info_1.loop.get_loop_length()
        d_list.append(d)
    print('\t'*func_depth+f"<TIME> compute_D_for_two_samples cost {time.time()-time1} seconds")
    return d_list

def compute_M_for_two_samples(loop_infos_1, loop_infos_2, func_depth=0):
    """
    run after get_common_loops_for_two_samples
    Compute the MD for two samples, the two samples should have the same set of loops(but may different in PETs)
    M=log2(IF2/IF1)
    """
    time1 = time.time()
    print('\t'*func_depth+"<FUNC> compute_M_for_two_samples")
    m_list = list()
    for loop_info_1, loop_info_2 in zip(loop_infos_1, loop_infos_2):
        m = math.log2(loop_info_1.pet/loop_info_2.pet)
        m_list.append(m)
    print('\t'*func_depth+f"<TIME> compute_M_for_two_samples cost {time.time()-time1} seconds")
    return m_list

def get_pet_matrix_from_loop_infos(data, func_depth=0):
    time1 = time.time()
    print('\t'*func_depth+"<FUNC> get_pet_matrix_from_loop_infos")
    pet_matrix = np.empty((len(data), len(data[0])), np.float32)
    for i in range(len(data)):
        for j in range(len(data[i])):
            pet_matrix[i][j] = data[i][j].pet
    print('\t'*func_depth+f"<TIME> get_pet_matrix_from_loop_infos cost {time.time()-time1} seconds")
    return pet_matrix

def gen_lowess_curve(x, y, it=3, delta=0, func_depth=0):
    time1 = time.time()
    print('\t'*func_depth+f"<FUNC> gen_lowess_curve, data num:{len(x)}")
    lowess = sm.nonparametric.lowess
    smoothed = lowess(y, x, frac=0.1, it=it, delta=delta, return_sorted=False)
    print('\t'*func_depth+f"<TIME> gen_lowess_curve cost {time.time()-time1} seconds")
    return smoothed

def draw_MD_scatter_subplot(samples_num, i, j, m_list, d_list, smoothed, title, func_depth=0):
    # x-axis: D, y-axis: M
    print('\t'*func_depth+"<FUNC> draw_MD_scatter_subplot")
    plt.subplot(samples_num, samples_num, samples_num*i+j+1)
    # the points color range from yellow to red according to the value of M
    cmap = plt.get_cmap("YlOrRd")
    colors = cmap(np.arange(cmap.N))

    # 从黄色开始，即从颜色列表的一半开始
    start = int(cmap.N / 4)
    new_cmap = ListedColormap(colors[start:,:-1])

    clist = [abs(m) for m in m_list]

    sc = plt.scatter(d_list, m_list, s=0.5, c=clist, cmap=new_cmap)
    # plt.colorbar(sc)
    # draw lowess curve, merge x and smoothed and sort by x
    lowess_xy = np.dstack((d_list, smoothed)) # its shape = [1, len(d_list), 2]
    lowess_xy.sort(axis=1)
    plt.plot(lowess_xy[0][:, 0], lowess_xy[0][:, 1], c='#00A2E8', lw=1, label='LOESS Fit')
    # draw x axis
    plt.axhline(0, c="gray", ls="--", lw=0.5)
    plt.xlim(0, np.percentile(d_list, 98))
    # compute mean and std of M, and add text at right-top in plt
    m_list = np.array(m_list)
    m_mean = np.mean(m_list)
    m_std = np.std(m_list)
    plt.text(0.65, 0.85, f"mean={m_mean:.2f}\nstd={m_std:.2f}", transform=plt.gca().transAxes)
    plt.title(title)
    return m_mean, m_std

def normalize_two_samples(loop_infos_1, loop_infos_2, smoothed, func_depth=0):
    time1 = time.time()
    print('\t'*func_depth+"<FUNC> normalize_two_samples")
    for i in range(len(loop_infos_1)):
        loop_infos_1[i].pet = math.exp2(math.log2(loop_infos_1[i].pet)-smoothed[i]/2)
        loop_infos_2[i].pet = math.exp2(math.log2(loop_infos_2[i].pet)+smoothed[i]/2)
    print('\t'*func_depth+f"<TIME> normalize_two_samples cost {time.time()-time1} seconds")
    return loop_infos_1, loop_infos_2

def normalize_two_samples_from_pets(pets_1, pets_2, smoothed, func_depth=0):
    time1 = time.time()
    print('\t'*func_depth+"<FUNC> normalize_two_samples_from_pets")
    pets_1 = np.array(pets_1)
    pets_2 = np.array(pets_2)
    pets_1 = np.exp2(np.log2(pets_1)-smoothed/2)
    pets_2 = np.exp2(np.log2(pets_2)+smoothed/2)
    pets_1 = np.maximum(pets_1, 1E-6)
    pets_2 = np.maximum(pets_2, 1E-6)
    print('\t'*func_depth+f"<TIME> normalize_two_samples_from_pets cost {time.time()-time1} seconds")
    return pets_1, pets_2

def draw_MD_scatter_plot(data, samples, output_path, func_depth=0):
    time1 = time.time()
    next_func_depth = func_depth + 1
    print('\t'*func_depth+"<FUNC> draw_MD_scatter_plot")
    FONTSIZE = 10
    plt.figure(dpi=300, figsize=(12, 12))
    plt.rc("font", family="Times New Roman")
    params = {"axes.titlesize": FONTSIZE,
              "legend.fontsize": FONTSIZE,
              "axes.labelsize": FONTSIZE,
              "xtick.labelsize": FONTSIZE,
              "ytick.labelsize": FONTSIZE,
              "figure.titlesize": FONTSIZE,
              "font.size": FONTSIZE}
    plt.rcParams.update(params)
    samples_num = len(data)
    m_mean_list, m_std_list = list(), list()
    for i in range(samples_num):
        for j in range(samples_num):
            if i >= j:
                # plt.subplot(samples_num, samples_num, samples_num*i+j+1)
                continue
            loop_infos_1 = data[i].copy()
            loop_infos_2 = data[j].copy()
            loop_infos_1, loop_infos_2 = get_common_loops_for_two_samples(loop_infos_1, loop_infos_2, next_func_depth)
            # * compute D
            d_list = compute_D_for_two_samples(loop_infos_1, next_func_depth)
            # * compute M
            m_list = compute_M_for_two_samples(loop_infos_1, loop_infos_2, next_func_depth)
            # * draw lowess curve
            smoothed = gen_lowess_curve(d_list, m_list, func_depth=next_func_depth)
            # * draw scatter plot
            m_mean, m_std = draw_MD_scatter_subplot(samples_num, i, j, m_list, d_list, smoothed, title=f"{samples[i]} vs {samples[j]}", func_depth=next_func_depth)
            m_mean_list.append(m_mean)
            m_std_list.append(m_std)
    plt.savefig(f"{output_path}.png", format="png", bbox_inches="tight")
    plt.close()
    print('\t'*func_depth+f"<TIME> draw_MD_scatter_plot cost {time.time()-time1} seconds")
    return m_mean_list, m_std_list

def draw_MD_scatter_plot_v2(pet_matrix, no_zero_pets_index_dict, d_list_all, samples, output_path, func_depth=0):
    time1 = time.time()
    next_func_depth = func_depth + 1
    print('\t'*func_depth+"<FUNC> draw_MD_scatter_plot")
    FONTSIZE = 10
    plt.figure(dpi=300, figsize=(12, 12))
    plt.rc("font", family="Times New Roman")
    params = {"axes.titlesize": FONTSIZE,
              "legend.fontsize": FONTSIZE,
              "axes.labelsize": FONTSIZE,
              "xtick.labelsize": FONTSIZE,
              "ytick.labelsize": FONTSIZE,
              "figure.titlesize": FONTSIZE,
              "font.size": FONTSIZE}
    plt.rcParams.update(params)
    samples_num = len(pet_matrix)
    m_mean_list, m_std_list = list(), list()
    for i in range(samples_num):
        for j in range(samples_num):
            if i >= j:
                continue
            # * filter zero pet loops
            pets_1 = pet_matrix[i][no_zero_pets_index_dict[(i, j)]]
            pets_2 = pet_matrix[j][no_zero_pets_index_dict[(i, j)]]
            # * compute D
            # d_list = compute_D_for_two_samples(loop_infos_1, next_func_depth)
            d_list = d_list_all[no_zero_pets_index_dict[(i, j)]]
            # * compute M
            # m_list = compute_M_for_two_samples(loop_infos_1, loop_infos_2, next_func_depth)
            m_list = np.log2(pets_1/pets_2)
            # * draw lowess curve
            smoothed = gen_lowess_curve(d_list, m_list, it=1, delta=0.1, func_depth=next_func_depth)
            # * draw scatter plot
            m_mean, m_std = draw_MD_scatter_subplot(samples_num, i, j, m_list, d_list, smoothed, title=f"{samples[i]} vs {samples[j]}", func_depth=next_func_depth)
            m_mean_list.append(m_mean)
            m_std_list.append(m_std)
    plt.savefig(f"{output_path}.png", format="png", bbox_inches="tight")
    plt.close()
    print('\t'*func_depth+f"<TIME> draw_MD_scatter_plot cost {time.time()-time1} seconds")
    return m_mean_list, m_std_list

class EarlyStop:
    def __init__(self, patience=5, mean_threshold=0.01, func_depth=0):
        self.patience = patience
        self.best = math.inf
        self.current_patience = 0
        self.mean_threshold = mean_threshold
        self.func_depth = func_depth

    def update(self, m_mean_list):
        # if m_mean change less than mean_threshold "patience" times, return True
        if np.mean(m_mean_list) - self.best < self.mean_threshold:
            self.current_patience += 1
            print('\t'*self.func_depth+f"[EARLYSTOP] patience: {self.current_patience}")
        else:
            self.current_patience = 0
        if self.current_patience >= self.patience:
            return True
        return False

def verify_coverage(data, func_depth=0):
    time1 = time.time()
    next_func_depth = func_depth + 1
    print('\t'*func_depth+"<FUNC> verify_coverage")
    samples_num = len(data)
    m_mean_list, m_std_list = list(), list()
    for i in range(samples_num):
        for j in range(samples_num):
            if i >= j:
                # plt.subplot(samples_num, samples_num, samples_num*i+j+1)
                continue
            # // loop_infos_1 = data[i].copy()
            # // loop_infos_2 = data[j].copy()
            # // loop_infos_1, loop_infos_2 = get_common_loops_for_two_samples(loop_infos_1, loop_infos_2, next_func_depth)
            loop_infos_1 = data[i]
            loop_infos_2 = data[j]
            # * compute M
            m_list = compute_M_for_two_samples(loop_infos_1, loop_infos_2, next_func_depth)
            # * get mean and std of M
            m_list = np.array(m_list)
            m_mean = np.mean(m_list)
            m_std = np.std(m_list)
            m_mean_list.append(m_mean)
            m_std_list.append(m_std)
    print('\t'*func_depth+f"<TIME> verify_coverage cost {time.time()-time1} seconds")
    return m_mean_list, m_std_list

def verify_coverage_v2(pet_matrix, no_zero_pets_index_dict, func_depth=0):
    time1 = time.time()
    print('\t'*func_depth+"<FUNC> verify_coverage_v2")
    samples_num = len(pet_matrix)
    m_mean_list, m_std_list = list(), list()
    for i in range(samples_num):
        for j in range(samples_num):
            if i >= j:
                continue
            # * filter zero pet loops
            pets_1 = pet_matrix[i][no_zero_pets_index_dict[(i, j)]]
            pets_2 = pet_matrix[j][no_zero_pets_index_dict[(i, j)]]
            # * compute M
            m_list = np.log2(pets_1/pets_2)
            m_mean = np.mean(m_list)
            m_std = np.std(m_list)
            m_mean_list.append(m_mean)
            m_std_list.append(m_std)
    print('\t'*func_depth+f"<TIME> verify_coverage_v2 cost {time.time()-time1} seconds")
    return m_mean_list, m_std_list

def normalize(data, samples, output_path, fig_demand=True, func_depth=0, complement_zero_pet_loops=False):
    """
    用 MD plot 来 normalize PETs values. 只修改了PETs的值，没有修改loop的信息。
    """
    time1 = time.time()
    print('\t'*func_depth+"<FUNC> normalize")
    next_func_depth = func_depth + 1
    early_stop_obj = EarlyStop(patience=5, mean_threshold=0.01, func_depth=next_func_depth)

    # * filter zero pet loops
    samples_num = len(data)
    no_zero_pets_index_dict = dict()
    for i in range(samples_num):
        for j in range(samples_num):
            if i >= j:
                continue
            pets_1 = [k.pet for k in data[i]]
            pets_2 = [k.pet for k in data[j]]
            no_zero_pets_index = filter_zero_pet_loops(pets_1, pets_2, next_func_depth)
            no_zero_pets_index_dict[(i, j)] = no_zero_pets_index

    # * compute D
    d_list_all = compute_D_for_two_samples(data[0], next_func_depth)
    d_list_all = np.array(d_list_all)

    pet_matrix = get_pet_matrix_from_loop_infos(data, next_func_depth)

    if fig_demand:
        # * draw origin MD scatter plot
        draw_MD_scatter_plot_v2(pet_matrix, no_zero_pets_index_dict, d_list_all, samples, f"{output_path}/figures/MD_scatter_plot_0", next_func_depth)
        # * draw origin MA scatter plot
        draw_MA_scatter_plot_v2(pet_matrix, no_zero_pets_index_dict, samples, f"{output_path}/figures/MA_scatter_plot_0", next_func_depth)

    # * normalize
    epoch_num = 0
    is_stop = False
    while is_stop == False:
        epoch_num += 1
        for i in range(samples_num):
            for j in range(samples_num):
                if i >= j:
                    continue
                # * get common loops
                # // loop_infos_1 = data[i].copy()
                # // loop_infos_2 = data[j].copy()
                # // loop_infos_1, loop_infos_2 = get_common_loops_for_two_samples(loop_infos_1, loop_infos_2, next_func_depth)
                # * filter zero pet loops
                pets_1 = pet_matrix[i][no_zero_pets_index_dict[(i, j)]]
                pets_2 = pet_matrix[j][no_zero_pets_index_dict[(i, j)]]
                # # * compute D
                # d_list = compute_D_for_two_samples(loop_infos_1, next_func_depth)
                d_list = d_list_all[no_zero_pets_index_dict[(i, j)]]
                # * compute M
                # m_list = compute_M_for_two_samples(loop_infos_1, loop_infos_2, next_func_depth)
                m_list = np.log2(pets_1/pets_2)
                # * draw lowess curve
                smoothed = gen_lowess_curve(d_list, m_list, func_depth=next_func_depth)
                # * normalize two samples
                # loop_infos_1, loop_infos_2 = normalize_two_samples(loop_infos_1, loop_infos_2, smoothed, next_func_depth)
                pets_1, pets_2 = normalize_two_samples_from_pets(pets_1, pets_2, smoothed, next_func_depth)
                # * update data[i]
                # // index_i, index_j = 0, 0
                # // while index_i < len(loop_infos_1):
                # //     loop_id_1 = loop_infos_1[index_i].get_loop_id()
                # //     loop_id_data = data[i][index_j].get_loop_id()
                # //     if loop_id_1 == loop_id_data:
                # //         data[i][index_j].pet = loop_infos_1[index_i].pet
                # //         index_i += 1
                # //         index_j += 1
                # //     elif loop_id_1 < loop_id_data:
                # //         index_i += 1
                # //     else:
                # //         index_j += 1
                for k, data_index in enumerate(no_zero_pets_index_dict[(i, j)]):
                    # data[i][data_index].pet = pets_1[k]
                    pet_matrix[i][data_index] = pets_1[k]
                # * update data[j]
                # // index_i, index_j = 0, 0
                # // while index_i < len(loop_infos_2):
                # //     loop_id_2 = loop_infos_2[index_i].get_loop_id()
                # //     loop_id_data = data[j][index_j].get_loop_id()
                # //     if loop_id_2 == loop_id_data:
                # //         data[j][index_j].pet = loop_infos_2[index_i].pet
                # //         index_i += 1
                # //         index_j += 1
                # //     elif loop_id_2 < loop_id_data:
                # //         index_i += 1
                # //     else:
                # //         index_j += 1
                for k, data_index in enumerate(no_zero_pets_index_dict[(i, j)]):
                    # data[j][data_index].pet = pets_2[k]
                    pet_matrix[j][data_index] = pets_2[k]

        if fig_demand:
            # * draw normalized MD scatter plot
            m_mean_list, _ = draw_MD_scatter_plot_v2(pet_matrix, no_zero_pets_index_dict, d_list_all, samples, f"{output_path}/figures/MD_scatter_plot_(md_norm)_epoch_{epoch_num}", next_func_depth)
            # * draw normalized MA scatter plot
            draw_MA_scatter_plot_v2(pet_matrix, no_zero_pets_index_dict, samples, f"{output_path}/figures/MA_scatter_plot_(md_norm)_epoch_{epoch_num}", next_func_depth)
        else:
            # m_mean_list, _ = verify_coverage(data, next_func_depth)
            m_mean_list, _ = verify_coverage_v2(pet_matrix, no_zero_pets_index_dict, next_func_depth)
        is_stop = early_stop_obj.update(m_mean_list)

    # * restore pets to data
    for i in range(samples_num):
        for j in range(len(data[i])):
            data[i][j].pet = pet_matrix[i][j]
    if complement_zero_pet_loops:
        from preprocessing.preprocessing_utils import complement_zero_pet_loops_for_samples
        data = complement_zero_pet_loops_for_samples(data, func_depth=next_func_depth)
    print('\t'*func_depth+f"<TIME> normalize cost {time.time()-time1} seconds")
    return data

def compute_A_for_two_samples(loop_infos_1, loop_infos_2, func_depth=0):
    """
    A <- (log2(obs/libsize.obs) + log2(ref/libsize.ref))/2
    """
    time1 = time.time()
    print('\t'*func_depth+"<FUNC> compute_A_for_two_samples")
    a_list = list()
    for loop_info_1, loop_info_2 in zip(loop_infos_1, loop_infos_2):
        a = (math.log2(loop_info_1.pet)+math.log2(loop_info_2.pet))/2
        a_list.append(a)
    print('\t'*func_depth+f"<TIME> compute_A_for_two_samples cost {time.time()-time1} seconds")
    return a_list

def draw_MA_scatter_subplot(samples_num, i, j, m_list, d_list, smoothed, title, func_depth=0):
    # x-axis: D, y-axis: M
    print('\t'*func_depth+"<FUNC> draw_MA_scatter_subplot")
    plt.subplot(samples_num, samples_num, samples_num*i+j+1)
    # the points color range from yellow to red according to the value of M
    cmap = plt.get_cmap("YlOrRd")
    colors = cmap(np.arange(cmap.N))

    # 从黄色开始，即从颜色列表的一半开始
    start = int(cmap.N / 4)
    new_cmap = ListedColormap(colors[start:,:-1])

    clist = [abs(m) for m in m_list]

    sc = plt.scatter(d_list, m_list, s=0.5, c=clist, cmap=new_cmap)
    # plt.colorbar(sc)
    # draw lowess curve, merge x and smoothed and sort by x
    lowess_xy = np.dstack((d_list, smoothed)) # its shape = [1, len(d_list), 2]
    lowess_xy.sort(axis=1)
    plt.plot(lowess_xy[0][:, 0], lowess_xy[0][:, 1], c='#00A2E8', lw=1, label='LOESS Fit')
    # draw x axis
    plt.axhline(0, c="gray", ls="--", lw=0.5)
    # plt.xlim(0, np.percentile(d_list, 98))
    # compute mean and std of M, and add text at right-top in plt
    m_list = np.array(m_list)
    m_mean = np.mean(m_list)
    m_std = np.std(m_list)
    plt.text(0.65, 0.85, f"mean={m_mean:.2f}\nstd={m_std:.2f}", transform=plt.gca().transAxes)
    plt.title(title)

def draw_MA_scatter_plot(data, samples, output_path, func_depth=0):
    """
    MA normalization uses a similar concept of representing
    measures from two datasets on a single plot, except
    it uses the Average chromatin interaction frequency as
    the X-axis instead of the Distance
    """
    time1 = time.time()
    next_func_depth = func_depth + 1
    print('\t'*func_depth+"<FUNC> draw_MA_scatter_plot")
    FONTSIZE = 10
    plt.figure(dpi=300, figsize=(12, 12))
    plt.rc("font", family="Times New Roman")
    params = {"axes.titlesize": FONTSIZE,
              "legend.fontsize": FONTSIZE,
              "axes.labelsize": FONTSIZE,
              "xtick.labelsize": FONTSIZE,
              "ytick.labelsize": FONTSIZE,
              "figure.titlesize": FONTSIZE,
              "font.size": FONTSIZE}
    plt.rcParams.update(params)
    samples_num = len(data)
    for i in range(samples_num):
        for j in range(samples_num):
            if i >= j:
                # plt.subplot(samples_num, samples_num, samples_num*i+j+1)
                continue
            loop_infos_1 = data[i].copy()
            loop_infos_2 = data[j].copy()
            loop_infos_1, loop_infos_2 = get_common_loops_for_two_samples(loop_infos_1, loop_infos_2, next_func_depth)
            # * compute A
            a_list = compute_A_for_two_samples(loop_infos_1, loop_infos_2, next_func_depth)
            # * compute M
            m_list = compute_M_for_two_samples(loop_infos_1, loop_infos_2, next_func_depth)
            # * draw lowess curve
            smoothed = gen_lowess_curve(a_list, m_list, func_depth=next_func_depth)
            # * draw scatter plot
            draw_MA_scatter_subplot(samples_num, i, j, m_list, a_list, smoothed, title=f"{samples[i]} vs {samples[j]}", func_depth=next_func_depth)
    plt.savefig(f"{output_path}.png", format="png", bbox_inches="tight")
    plt.close()
    print('\t'*func_depth+f"<TIME> draw_MA_scatter_plot cost {time.time()-time1} seconds")


def draw_MA_scatter_plot_v2(pet_matrix, no_zero_pets_index_dict, samples, output_path, func_depth=0):
    """
    MA normalization uses a similar concept of representing
    measures from two datasets on a single plot, except
    it uses the Average chromatin interaction frequency as
    the X-axis instead of the Distance
    """
    time1 = time.time()
    next_func_depth = func_depth + 1
    print('\t'*func_depth+"<FUNC> draw_MA_scatter_plot")
    FONTSIZE = 10
    plt.figure(dpi=300, figsize=(12, 12))
    plt.rc("font", family="Times New Roman")
    params = {"axes.titlesize": FONTSIZE,
              "legend.fontsize": FONTSIZE,
              "axes.labelsize": FONTSIZE,
              "xtick.labelsize": FONTSIZE,
              "ytick.labelsize": FONTSIZE,
              "figure.titlesize": FONTSIZE,
              "font.size": FONTSIZE}
    plt.rcParams.update(params)
    samples_num = len(pet_matrix)
    for i in range(samples_num):
        for j in range(samples_num):
            if i >= j:
                # plt.subplot(samples_num, samples_num, samples_num*i+j+1)
                continue
            # // loop_infos_1 = data[i].copy()
            # // loop_infos_2 = data[j].copy()
            # // loop_infos_1, loop_infos_2 = get_common_loops_for_two_samples(loop_infos_1, loop_infos_2, next_func_depth)
            pets_1 = pet_matrix[i][no_zero_pets_index_dict[(i, j)]]
            pets_2 = pet_matrix[j][no_zero_pets_index_dict[(i, j)]]
            # * compute A
            # a_list = compute_A_for_two_samples(loop_infos_1, loop_infos_2, next_func_depth)
            a_list = (np.log2(pets_1)+np.log2(pets_2))/2
            # * compute M
            # m_list = compute_M_for_two_samples(loop_infos_1, loop_infos_2, next_func_depth)
            m_list = np.log2(pets_1/pets_2)
            # * draw lowess curve
            smoothed = gen_lowess_curve(a_list, m_list, it=1, delta=0.1, func_depth=next_func_depth)
            # * draw scatter plot
            draw_MA_scatter_subplot(samples_num, i, j, m_list, a_list, smoothed, title=f"{samples[i]} vs {samples[j]}", func_depth=next_func_depth)
    plt.savefig(f"{output_path}.png", format="png", bbox_inches="tight")
    plt.close()
    print('\t'*func_depth+f"<TIME> draw_MA_scatter_plot cost {time.time()-time1} seconds")



