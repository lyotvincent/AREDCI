
import time
import numpy as np
from scipy.stats import norm
from statsmodels.stats.multitest import fdrcorrection

from dci.kde_dci import fit_predict



def run_dci_by_md_for_two_groups(pets_matrix, samples, groups, a_threshold=1, should_fit=True, func_depth=0):
    time1 = time.time()
    print('\t'*func_depth+"<FUNC> run_dci_by_md_for_two_groups")
    next_func_depth = func_depth + 1
    groups_set = list(set(groups))
    if should_fit:
        # * fit 1st group
        print('\t'*func_depth+f"[INFO] fit 1st group {groups_set[0]} alternative model")
        x = list()
        y = list()
        y_index = [i for i,g in enumerate(groups) if g == groups_set[0]]
        for i in range(len(pets_matrix)):
            # x.append([i+1, 1, 0])
            # x.append([i+1, 0, 1])
            # x.append([i+1])
            # x.append([i+1])
            for _ in range(len(y_index)):
                x.append([i+1, 1, 0])
            # x.append([i+1, 1, 0])
            # y.append(pets_matrix[i][0])
            # y.append(pets_matrix[i][1])
        # 把groups里属于groups_set[0]的样本的列提取出来放入y
        y = pets_matrix[:, y_index].flatten()
        x = np.array(x)
        # y = np.array(y)
        exog = np.insert(x, 0, 1, axis=1) # 在第0列，添加常数项1, 用于计算截距, same as sm.add_constant(x)
        fitted_values1 = fit_predict(exog, y)
        fitted_values1 = fitted_values1.reshape((pets_matrix.shape[0], len(y_index)))[:, 0]

        # * fit 2nd group
        print('\t'*func_depth+f"[INFO] fit 2nd group {groups_set[1]} alternative model")
        x = list()
        y = list()
        y_index = [i for i,g in enumerate(groups) if g == groups_set[1]]
        for i in range(len(pets_matrix)):
            # x.append([i+1, 1, 0])
            # x.append([i+1, 0, 1])
            # x.append([i+1])
            # x.append([i+1])
            for _ in range(len(y_index)):
                x.append([i+1, 0, 1])
            # x.append([i+1, 0, 1])
            # y.append(pets_matrix[i][2])
            # y.append(pets_matrix[i][3])
        y = pets_matrix[:, y_index].flatten()
        x = np.array(x)
        # y = np.array(y)
        exog = np.insert(x, 0, 1, axis=1) # 在第0列，添加常数项1, 用于计算截距, same as sm.add_constant(x)
        fitted_values2 = fit_predict(exog, y)
        fitted_values2 = fitted_values2.reshape((pets_matrix.shape[0], len(y_index)))[:, 0]

        # * validate/count 0 in fitted value
        zero_count1 = sum(fitted_values1 == 0)
        zero_count2 = sum(fitted_values2 == 0)
        print('\t'*func_depth+f"<INFO> zero count in fitted_values1: {zero_count1}")
        print('\t'*func_depth+f"<INFO> zero count in fitted_values2: {zero_count2}")
    else:
        y_index = [i for i,g in enumerate(groups) if g == groups_set[0]]
        fitted_values1 = np.maximum(np.mean(pets_matrix[:, y_index], axis=1), 1E-6)
        y_index = [i for i,g in enumerate(groups) if g == groups_set[1]]
        fitted_values2 = np.maximum(np.mean(pets_matrix[:, y_index], axis=1), 1E-6)

    # * compute M, A
    m_list = np.log2(fitted_values1/fitted_values2)
    # a_list = (np.log2(fitted_values1)+np.log2(fitted_values2))/2
    a_list = (fitted_values1+fitted_values2)/2

    # * set m = 0 where a < a_threshold
    m_list[a_list < a_threshold] = 0

    # * normalize m
    m_list = ( m_list - np.mean(m_list) ) / np.std(m_list)

    # * compute pvalue
    pvalues = norm.sf(m_list, loc=0, scale=1)

    pvalues[pvalues > 0.5] = 1 - pvalues[pvalues > 0.5]
    pvalues *= 2

    _, fdr = fdrcorrection(pvalues, alpha=0.05, method='indep', is_sorted=False)
    print('\t'*func_depth+f"<TIME> run_dci_by_md_for_two_groups cost {time.time()-time1} seconds")
    return pvalues, fdr


