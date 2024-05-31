
"""
Signed Squared Relative Difference (SSRD) DCI
eps = np.finfo(np.float32).eps
ssrd = np.sign(x-y)*np.square(x-y)/(x+y+eps)
"""

import time
import numpy as np
from scipy.stats import norm
from statsmodels.stats.multitest import fdrcorrection

from dci.kde_dci import fit_predict

def run_dci_by_ssrd_for_two_groups(pets_matrix, groups, should_fit=True, func_depth=0):
    time1 = time.time()
    print('\t'*func_depth+"<FUNC> run_dci_by_ssrd_for_two_groups")
    next_func_depth = func_depth + 1
    groups_set = list(set(groups))
    if should_fit:
        # * fit 1st group
        print('\t'*func_depth+"[INFO] fit 1st group alternative model")
        x = list()
        y = list()
        y_index_1 = [i for i,g in enumerate(groups) if g == groups_set[0]]
        for i in range(len(pets_matrix)):
            for _ in range(len(y_index_1)):
                x.append([i+1, 1, 0])
        y = pets_matrix[:, y_index_1].flatten()
        x = np.array(x)
        exog = np.insert(x, 0, 1, axis=1) # 在第0列，添加常数项1, 用于计算截距, same as sm.add_constant(x)
        fitted_values1 = fit_predict(exog, y)
        fitted_values1 = fitted_values1.reshape((pets_matrix.shape[0], len(y_index_1)))[:, 0]

        # * fit 2nd group
        print('\t'*func_depth+"[INFO] fit 2nd group alternative model")
        x = list()
        y = list()
        y_index_2 = [i for i,g in enumerate(groups) if g == groups_set[1]]
        for i in range(len(pets_matrix)):
            for _ in range(len(y_index_2)):
                x.append([i+1, 0, 1])
        y = pets_matrix[:, y_index_2].flatten()
        x = np.array(x)
        exog = np.insert(x, 0, 1, axis=1) # 在第0列，添加常数项1, 用于计算截距, same as sm.add_constant(x)
        fitted_values2 = fit_predict(exog, y)
        fitted_values2 = fitted_values2.reshape((pets_matrix.shape[0], len(y_index_2)))[:, 0]
    else:
        y_index = [i for i,g in enumerate(groups) if g == groups_set[0]]
        fitted_values1 = np.mean(pets_matrix[:, y_index], axis=1)
        y_index = [i for i,g in enumerate(groups) if g == groups_set[1]]
        fitted_values2 = np.mean(pets_matrix[:, y_index], axis=1)

    # * compute SSRD
    eps = np.finfo(np.float32).eps
    ssrd_list = np.sign(fitted_values1-fitted_values2)*np.square(fitted_values1-fitted_values2)/(fitted_values1+fitted_values2+eps)

    # * normalize SSRD
    ssrd_list = ( ssrd_list - np.mean(ssrd_list) ) / np.std(ssrd_list)

    # * compute pvalue
    pvalues = norm.sf(ssrd_list, loc=0, scale=1)

    pvalues[pvalues > 0.5] = 1 - pvalues[pvalues > 0.5]
    pvalues *= 2

    _, fdr = fdrcorrection(pvalues, alpha=0.05, method='indep', is_sorted=False)
    print('\t'*func_depth+f"<TIME> run_dci_by_ssrd_for_two_groups cost {time.time()-time1} seconds")
    return pvalues, fdr

