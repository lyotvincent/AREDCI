
import time, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from preprocessing.md_normalizer import EarlyStop, get_pet_matrix_from_loop_infos, compute_D_for_two_samples, filter_zero_pet_loops, gen_lowess_curve, draw_MD_scatter_plot_v2, draw_MA_scatter_plot_v2

def ssrd_normalize_two_samples(pets_1, pets_2, smoothed, func_depth=0):
    """
    用 SSRD 来 normalize PETs values. 只修改了PETs的值，没有修改loop的信息。
    """
    # time1 = time.time()
    print('\t'*func_depth+"<FUNC> ssrd_normalize_two_samples")
    # * update samples
    eps = np.finfo(pets_1.dtype).eps
    denominator = np.where(pets_1 != pets_2, pets_1 - pets_2, eps)
    alpha = (pets_1 + pets_2) / (2 * denominator)
    new_pets_1 = pets_1 - np.sign(pets_1-pets_2) * alpha * smoothed
    new_pets_1 = np.maximum(new_pets_1, 0)
    new_pets_2 = pets_2 + np.sign(pets_1-pets_2) * alpha * smoothed
    new_pets_2 = np.maximum(new_pets_2, 0)
    # print('\t'*func_depth+f"<TIME> ssrd_normalize_two_samples cost {time.time()-time1} seconds")
    return new_pets_1, new_pets_2

def ssrd_normalize_two_samples_v2(pets_1, pets_2, smoothed, func_depth=0):
    """
    用 SSRD 来 normalize PETs values. 只修改了PETs的值，没有修改loop的信息。
    保持PETs_1和PETs_2的值是整数
    """
    # time1 = time.time()
    print('\t'*func_depth+"<FUNC> ssrd_normalize_two_samples_v2")
    # * update samples
    denominator = np.where(pets_1 != pets_2, pets_1 - pets_2, 1.)
    alpha = (pets_1 + pets_2) / (2 * denominator)
    new_pets_1 = pets_1 - np.sign(pets_1-pets_2) * alpha * smoothed
    new_pets_1 = np.round(np.maximum(new_pets_1, 0.))
    new_pets_2 = pets_2 + np.sign(pets_1-pets_2) * alpha * smoothed
    new_pets_2 = np.round(np.maximum(new_pets_2, 0.))
    # print('\t'*func_depth+f"<TIME> ssrd_normalize_two_samples_v2 cost {time.time()-time1} seconds")
    return new_pets_1, new_pets_2

def verify_coverage(pet_matrix, func_depth=0):
    time1 = time.time()
    print('\t'*func_depth+"<FUNC> verify_coverage_v2")
    samples_num = len(pet_matrix)
    ssrd_mean_list, ssrd_std_list = list(), list()
    for i in range(samples_num):
        for j in range(samples_num):
            if i >= j:
                continue
            # * filter zero pet loops
            pets_1 = pet_matrix[i]
            pets_2 = pet_matrix[j]
            # * compute SSRD
            eps = np.finfo(np.float32).eps
            ssrd_list = np.sign(pets_1-pets_2)*np.square(pets_1-pets_2)/(pets_1+pets_2+eps)
            ssrd_mean = np.mean(ssrd_list)
            ssrd_std = np.std(ssrd_list)
            ssrd_mean_list.append(ssrd_mean)
            ssrd_std_list.append(ssrd_std)
    print('\t'*func_depth+f"<TIME> verify_coverage_v2 cost {time.time()-time1} seconds")
    return ssrd_mean_list, ssrd_std_list

def draw_SSRD_scatter_subplot(samples_num, i, j, ssrd_list, d_list, smoothed, title, func_depth=0):
    # x-axis: D, y-axis: M
    print('\t'*func_depth+"<FUNC> draw_MD_scatter_subplot")
    plt.subplot(samples_num, samples_num, samples_num*i+j+1)
    # the points color range from yellow to red according to the value of M
    cmap = plt.get_cmap("YlOrRd")
    colors = cmap(np.arange(cmap.N))

    # 从黄色开始，即从颜色列表的一半开始
    start = int(cmap.N / 4)
    new_cmap = ListedColormap(colors[start:,:-1])

    clist = [abs(m) for m in ssrd_list]

    sc = plt.scatter(d_list, ssrd_list, s=0.5, c=clist, cmap=new_cmap)
    # plt.colorbar(sc)
    # draw lowess curve, merge x and smoothed and sort by x
    lowess_xy = np.dstack((d_list, smoothed)) # its shape = [1, len(d_list), 2]
    lowess_xy.sort(axis=1)
    plt.plot(lowess_xy[0][:, 0], lowess_xy[0][:, 1], c='#00A2E8', lw=1, label='LOESS Fit')
    # draw x axis
    plt.axhline(0, c="gray", ls="--", lw=0.5)
    plt.xlim(0, np.percentile(d_list, 98))
    # compute mean and std of M, and add text at right-top in plt
    ssrd_list = np.array(ssrd_list)
    ssrd_mean = np.mean(ssrd_list)
    ssrd_std = np.std(ssrd_list)
    plt.text(0.65, 0.85, f"mean={ssrd_mean:.2f}\nstd={ssrd_std:.2f}", transform=plt.gca().transAxes)
    plt.title(title)


def draw_SSRD_scatter_plot(pet_matrix, d_list, samples, output_path, func_depth=0):
    time1 = time.time()
    next_func_depth = func_depth + 1
    print('\t'*func_depth+"<FUNC> draw_SSRD_scatter_plot")
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
                continue
            pets_1 = pet_matrix[i]
            pets_2 = pet_matrix[j]
            # * compute SSRD
            eps = np.finfo(np.float32).eps
            ssrd_list = np.sign(pets_1-pets_2)*np.square(pets_1-pets_2)/(pets_1+pets_2+eps)
            # * draw lowess curve
            smoothed = gen_lowess_curve(d_list, ssrd_list, it=1, delta=0.1, func_depth=next_func_depth)
            # * draw scatter plot
            draw_SSRD_scatter_subplot(samples_num, i, j, ssrd_list, d_list, smoothed, title=f"{samples[i]} vs {samples[j]}", func_depth=next_func_depth)
    plt.savefig(f"{output_path}.png", format="png", bbox_inches="tight")
    plt.close()
    print('\t'*func_depth+f"<TIME> draw_SSRD_scatter_plot cost {time.time()-time1} seconds")


def normalize(data, samples, output_path, fig_demand=True, func_depth=0):
    """
    用 SSRD (Signed Squared Relative Difference) plot 来 normalize PETs values. 只修改了PETs的值，没有修改loop的信息。
    """
    time1 = time.time()
    print('\t'*func_depth+"<FUNC> normalize")
    next_func_depth = func_depth + 1
    early_stop_obj = EarlyStop(patience=5, mean_threshold=0.01, func_depth=next_func_depth)
    samples_num = len(data)

    # * compute D
    d_list = compute_D_for_two_samples(data[0], next_func_depth)
    d_list = np.array(d_list)

    pet_matrix = get_pet_matrix_from_loop_infos(data, next_func_depth)

    if fig_demand:
        no_zero_pets_index_dict = dict()
        for i in range(samples_num):
            for j in range(samples_num):
                if i >= j:
                    continue
                pets_1 = [k.pet for k in data[i]]
                pets_2 = [k.pet for k in data[j]]
                no_zero_pets_index = filter_zero_pet_loops(pets_1, pets_2, next_func_depth)
                no_zero_pets_index_dict[(i, j)] = no_zero_pets_index
        if os.path.exists(f"{output_path}/figures/MD_scatter_plot_0") == False:
            # * draw origin MD scatter plot
            draw_MD_scatter_plot_v2(pet_matrix, no_zero_pets_index_dict, d_list, samples, f"{output_path}/figures/MD_scatter_plot_0", next_func_depth)
            # * draw origin MA scatter plot
            draw_MA_scatter_plot_v2(pet_matrix, no_zero_pets_index_dict, samples, f"{output_path}/figures/MA_scatter_plot_0", next_func_depth)
            # * draw origin SSRD scatter plot
            draw_SSRD_scatter_plot(pet_matrix, d_list, samples, f"{output_path}/figures/SSRD_scatter_plot_0", next_func_depth)

    # * normalize
    epoch_num = 0
    is_stop = False
    while is_stop == False:
        epoch_num += 1
        for i in range(samples_num):
            for j in range(samples_num):
                if i >= j:
                    continue
                pets_1 = pet_matrix[i]
                pets_2 = pet_matrix[j]
                # * compute SSRD
                eps = np.finfo(np.float32).eps
                ssrd_list = np.sign(pets_1-pets_2)*np.square(pets_1-pets_2)/(pets_1+pets_2+eps)
                # * draw lowess curve
                smoothed = gen_lowess_curve(d_list, ssrd_list, func_depth=next_func_depth)
                # * normalize two samples
                pets_1, pets_2 = ssrd_normalize_two_samples_v2(pets_1, pets_2, smoothed, next_func_depth)
                # * update data[i] & data[j]
                for k in range(len(pets_1)):
                    pet_matrix[i][k] = pets_1[k]
                    pet_matrix[j][k] = pets_2[k]

        if fig_demand:
            no_zero_pets_index_dict = dict()
            for i in range(samples_num):
                for j in range(samples_num):
                    if i >= j:
                        continue
                    pets_1 = pet_matrix[i]
                    pets_2 = pet_matrix[j]
                    no_zero_pets_index = filter_zero_pet_loops(pets_1, pets_2, next_func_depth)
                    no_zero_pets_index_dict[(i, j)] = no_zero_pets_index
            # * draw normalized MD scatter plot
            draw_MD_scatter_plot_v2(pet_matrix, no_zero_pets_index_dict, d_list, samples, f"{output_path}/figures/MD_scatter_plot_(ssrd_norm)_epoch_{epoch_num}", next_func_depth)
            # * draw normalized MA scatter plot
            draw_MA_scatter_plot_v2(pet_matrix, no_zero_pets_index_dict, samples, f"{output_path}/figures/MA_scatter_plot_(ssrd_norm)_epoch_{epoch_num}", next_func_depth)
            # * draw normalized SSRD scatter plot
            draw_SSRD_scatter_plot(pet_matrix, d_list, samples, f"{output_path}/figures/SSRD_scatter_plot_(ssrd_norm)_epoch_{epoch_num}", next_func_depth)

        ssrd_mean_list, _ = verify_coverage(pet_matrix, next_func_depth)
        is_stop = early_stop_obj.update(ssrd_mean_list)
    
    # * restore pets to data
    for i in range(samples_num):
        for j in range(len(data[i])):
            data[i][j].pet = pet_matrix[i][j]

    print('\t'*func_depth+f"<TIME> normalize cost {time.time()-time1} seconds")
    return data




