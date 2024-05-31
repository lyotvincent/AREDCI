
import time
import numpy as np


from preprocessing.md_normalizer import EarlyStop, get_pet_matrix_from_loop_infos, filter_zero_pet_loops, gen_lowess_curve, draw_MD_scatter_plot_v2, draw_MA_scatter_plot_v2, verify_coverage_v2


def normalize_two_samples_from_pets(a_list, m_correct, func_depth=0):
    time1 = time.time()
    print('\t'*func_depth+"<FUNC> normalize_two_samples_from_pets")
    a_list = np.array(a_list)
    pets_1 = np.power(2, (a_list+0.5*m_correct))
    pets_2 = np.power(2, (a_list-0.5*m_correct))
    pets_1 = np.maximum(pets_1, 1E-6)
    pets_2 = np.maximum(pets_2, 1E-6)
    print('\t'*func_depth+f"<TIME> normalize_two_samples_from_pets cost {time.time()-time1} seconds")
    return pets_1, pets_2

def normalize(data, samples, output_path, fig_demand=True, func_depth=0, complement_zero_pet_loops=False):
    """
    用 MA plot 来 normalize PETs values. 只修改了PETs的值，没有修改loop的信息。
    """
    time1 = time.time()
    print('\t'*func_depth+"<FUNC> ma normalize")
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

    pet_matrix = get_pet_matrix_from_loop_infos(data, next_func_depth)

    if fig_demand:
        from preprocessing.md_normalizer import compute_D_for_two_samples
        from preprocessing.ssrd_normalizer import draw_SSRD_scatter_plot
        # * compute D
        d_list_all = compute_D_for_two_samples(data[0], next_func_depth)
        d_list_all = np.array(d_list_all)
        # * draw origin MD scatter plot
        draw_MD_scatter_plot_v2(pet_matrix, no_zero_pets_index_dict, d_list_all, samples, f"{output_path}/figures/MD_scatter_plot_0", next_func_depth)
        # * draw origin MA scatter plot
        draw_MA_scatter_plot_v2(pet_matrix, no_zero_pets_index_dict, samples, f"{output_path}/figures/MA_scatter_plot_0", next_func_depth)
        # * draw normalized SSRD scatter plot
        draw_SSRD_scatter_plot(pet_matrix, d_list_all, samples, f"{output_path}/figures/SSRD_scatter_plot_0", next_func_depth)


    # * normalize
    epoch_num = 0
    is_stop = False
    while is_stop == False:
        epoch_num += 1
        for i in range(samples_num):
            for j in range(samples_num):
                if i >= j:
                    continue
                # * filter zero pet loops
                pets_1 = pet_matrix[i][no_zero_pets_index_dict[(i, j)]]
                pets_2 = pet_matrix[j][no_zero_pets_index_dict[(i, j)]]
                # * compute A
                a_list = np.log2(np.multiply(pets_1,pets_2))/2
                # * compute M
                m_list = np.log2(pets_1/pets_2)
                # * draw lowess curve
                smoothed = gen_lowess_curve(a_list, m_list, func_depth=next_func_depth)
                # * normalize two samples
                pets_1, pets_2 = normalize_two_samples_from_pets(a_list, m_list-smoothed, next_func_depth)
                # * update data[i]
                for k, data_index in enumerate(no_zero_pets_index_dict[(i, j)]):
                    pet_matrix[i][data_index] = pets_1[k]
                # * update data[j]
                for k, data_index in enumerate(no_zero_pets_index_dict[(i, j)]):
                    pet_matrix[j][data_index] = pets_2[k]

        if fig_demand:
            # * compute D
            d_list_all = compute_D_for_two_samples(data[0], next_func_depth)
            d_list_all = np.array(d_list_all)
            # * draw normalized MD scatter plot
            m_mean_list, _ = draw_MD_scatter_plot_v2(pet_matrix, no_zero_pets_index_dict, d_list_all, samples, f"{output_path}/figures/MD_scatter_plot_(ma_norm)_epoch_{epoch_num}", next_func_depth)
            # * draw normalized MA scatter plot
            draw_MA_scatter_plot_v2(pet_matrix, no_zero_pets_index_dict, samples, f"{output_path}/figures/MA_scatter_plot_(ma_norm)_epoch_{epoch_num}", next_func_depth)
            # * draw normalized SSRD scatter plot
            draw_SSRD_scatter_plot(pet_matrix, d_list_all, samples, f"{output_path}/figures/SSRD_scatter_plot_(ma_norm)_epoch_{epoch_num}", next_func_depth)

        m_mean_list, _ = verify_coverage_v2(pet_matrix, no_zero_pets_index_dict, next_func_depth)
        is_stop = early_stop_obj.update(m_mean_list)
    
    # * restore pets to data
    for i in range(samples_num):
        for j in range(len(data[i])):
            data[i][j].pet = pet_matrix[i][j]
    print('\t'*func_depth+f"<TIME> ma normalize cost {time.time()-time1} seconds")
    return data
