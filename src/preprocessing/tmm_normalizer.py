import time, math
import numpy as np
from preprocessing.preprocessing_utils import complement_zero_pet_loops_for_samples
import warnings


def find_reference(pet_matrix, func_depth=0):
    # return the index of reference sample
    print('\t'*func_depth+"<FUNC> find_reference")
    samples_num = len(pet_matrix)
    pet_75_list = [np.percentile(pet_matrix[i], 75) for i in range(samples_num)]
    mean_pet_75 = np.mean(pet_75_list)
    # find the sample whose pet_75 is closest to mean_pet_75
    min_diff = math.inf
    reference_sample_index = 0
    for i in range(samples_num):
        cur = abs(pet_75_list[i] - mean_pet_75)
        if cur < min_diff:
            min_diff = cur
            reference_sample_index = i
    return reference_sample_index

def tmm_normalize(data, samples, output_path, fig_demand=True, func_depth=0):
    # warnings.filterwarnings("ignore", category=RuntimeWarning)
    time1 = time.time()
    print('\t'*func_depth+"<FUNC> tmm_normalize")
    LOGRATIOTRIM = 0.3
    SUMTRIM = 0.05
    ACUTOFF = -1E10
    next_func_depth = func_depth + 1
    samples_num = len(data)

    # * Step 0: 补全pet为0的样本中的0 PET loop
    # data_full = data.copy()
    '''
    data.copy() 是浅拷贝，它只会拷贝 data 的最外层元素，如果 data 是一个嵌套的列表（即列表的元素还是列表），那么 data.copy() 只会拷贝最外层的列表，而不会拷贝内部的子列表。所以，如果你修改了 data 的子列表，那么 data.copy() 得到的结果也会被修改。
    [s.copy() for s in data] 是深拷贝，它会拷贝 data 的每一个元素，包括所有的子列表。所以，如果你修改了 data 的子列表，那么 [s.copy() for s in data] 得到的结果不会被修改。
    '''
    data_full = [s.copy() for s in data]
    # data_full = complement_zero_pet_loops_for_samples(data_full, func_depth=next_func_depth)

    # * Step 1: 选择参考样本
    """
    TMM normalization首先需要选择一个参考样本，以它为基准进行校正。
    默认下，参考样本的选择是通过比较每个样本的CPM (counts per million)的上四分位数与所有样本CPM的平均上四分位数之间的差值，找出差值最小的样本作为参考样本。
    """
    # TODO 在edgeR::calcNormFactors()中，也可以通过refColumn=参数指定参考样本
    # total has samples_num columns, every column is a sample containing pet list of loops
    pet_matrix = [list() for _ in range(samples_num)]
    for i in range(samples_num):
        for j in range(len(data_full[i])):
            pet_matrix[i].append(data_full[i][j].pet)
    pet_matrix = np.array(pet_matrix)
    libsize = np.sum(pet_matrix, axis=1)
    cpm_pet_matrix = pet_matrix/libsize[:, np.newaxis]*1E6 # CPM
    reference_index = find_reference(cpm_pet_matrix, func_depth=next_func_depth)
    print('\t'*func_depth+f"<INFO> reference sample: {samples[reference_index]}")

    # * Step 2: calculation of sample-reference pairwise M and A
    # 接着，在参考样本和非参考样本间两两计算校正因子（normalization factors）。
    # 我们首先需要计算参考样本和非参考样本间的Fold change (M)和平均表达量 (A)
    factors = np.zeros(samples_num)
    for i in range(samples_num):
        obs = pet_matrix[i] # obs: 非参考样本的原始count
        ref = pet_matrix[reference_index] # ref: 参考样本的原始count
        libsize_obs = np.sum(obs) # libsize.obs: 非参考样本的原始文库大小
        libsize_ref = np.sum(ref) # libsize.ref: 参考样本的原始文库大小
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            # M value: log ratio of expression, accounting for library size
            m = np.log2((obs/libsize_obs)/(ref/libsize_ref))
            # A value:absolute expression
            a = (np.log2(obs/libsize_obs) + np.log2(ref/libsize_ref))/2
            # estimated asymptotic variance
            # v <- (libsize.obs - obs)/libsize.obs/obs + (libsize.ref - ref)/libsize.ref/ref
            v = (libsize_obs - obs)/libsize_obs/obs + (libsize_ref - ref)/libsize_ref/ref
        # remove infinite values, cutoff based on A
        # fin <- is.finite(M) & is.finite(A) & (A > Acutoff)
        # a: 在obs和ref中pet都为0的会过滤
        # m: 在obs或ref中pet都为0的会过滤
        fin = np.isfinite(m) & np.isfinite(a) & (a > ACUTOFF)

        # filter out infinite values by fin
        m = m[fin]
        a = a[fin]
        v = v[fin]

        # * Step 3: trimmed mean of M values
        # Double trim the upper and lower percentages of the data
        # trim M values by 30% and A values by 5%
        n = len(m)
        loM = np.floor(n * LOGRATIOTRIM) + 1
        hiM = n - loM + 1
        loA = np.floor(n * SUMTRIM) + 1
        hiA = n - loA + 1

        # keep <- (rank(M)>=loM & rank(M)<=hiM) & (rank(A)>=loA & rank(A)<=hiA)
        keep = (np.argsort(np.argsort(m)) >= loM) & (np.argsort(np.argsort(a)) >= loA) & (np.argsort(np.argsort(m)) <= hiM) & (np.argsort(np.argsort(a)) <= hiA)

        # Weighted mean of M after trimming
        factors[i] = np.sum(m[keep]/v[keep])/np.sum(1/v[keep])
        factors[i] = 2**factors[i]

    # Factors should multiple to one, 让所有样本的校正因子乘积为1
    # f <- f/exp(mean(log(f)))
    factors = factors/np.exp(np.mean(np.log(factors)))
    
    # * Step 4: TMM normalized counts
    # 原方法是用counts/(libsize*factors)，但是这里的counts是pet，不太一样，此处暂不用libsize
    for i in range(samples_num):
        sized_libsize = libsize[i]*factors[i]
        for j in range(len(data[i])):
            data[i][j].pet = data[i][j].pet/sized_libsize * 1E6
            data[i][j].pet = round(data[i][j].pet, 4)

    if fig_demand:
        from preprocessing.md_normalizer import filter_zero_pet_loops, compute_D_for_two_samples, get_pet_matrix_from_loop_infos, draw_MD_scatter_plot_v2, draw_MA_scatter_plot_v2
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
        d_list_all = compute_D_for_two_samples(data[i], next_func_depth)
        d_list_all = np.array(d_list_all)

        pet_matrix = get_pet_matrix_from_loop_infos(data, next_func_depth)
        # * draw TMM MD scatter plot
        draw_MD_scatter_plot_v2(pet_matrix, no_zero_pets_index_dict, d_list_all, samples, f"{output_path}/figures/MD_scatter_plot_tmm", next_func_depth)
        draw_MA_scatter_plot_v2(pet_matrix, no_zero_pets_index_dict, samples, f"{output_path}/figures/MA_scatter_plot_tmm", next_func_depth)

    print('\t'*func_depth+f"<TIME> tmm_normalize cost {time.time()-time1} seconds")
    # np.save(f"{output_path}/tmm_pets.npy", pet_matrix)
    # np.savetxt(f'{output_path}/tmm_pets.csv', pet_matrix.T, fmt='%d', delimiter=',')
    return data

# if __name__ == "__main__":
#     a= [[1,2,3,4,5,6,7,8,9,10],[11,12,13,14,15,16,17,18,19,20],[21,22,23,24,25,26,27,28,29,30],[31,32,33,34,35,36,37,38,39,40]]
#     a = np.array(a)
#     libsize = np.sum(a, axis=1)
#     print(libsize)
#     print(libsize[:, np.newaxis])
#     # every row divide by libsize
#     b = a/libsize[:, np.newaxis]
#     print(b)
#     print(np.array([1,2,3,4,5,6,7,8,9,10])/55)


