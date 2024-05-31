
import numpy as np
# import statsmodels.api as sm
from statsmodels.stats.multitest import fdrcorrection
# from sklearn.ensemble import RandomForestRegressor
# from sklearn.metrics import mean_squared_error, mean_absolute_error, accuracy_score, roc_auc_score, precision_score, recall_score, f1_score, fbeta_score, roc_curve, auc, average_precision_score, precision_recall_curve, PrecisionRecallDisplay
from sklearn.metrics import mean_squared_error
from xgboost import XGBRegressor
from scipy.stats import chi2, fisher_exact, chisquare
import time

# from kde_skl import compute_P_value_by_kde_skl
from dci.kde_scipy import compute_P_value_by_kde_scipy

import warnings, os

def fit_predict(
        exog, # x, exogenous
        endog # y, endogenous 
    ):
    # 创建随机森林回归模型
    # model = RandomForestRegressor()
    # or sm.nonparametric.lowess?
    # model = XGBRegressor(max_depth=3, n_estimators=128, tree_method='gpu_hist', gpu_id=0, seed=9523)
    model = XGBRegressor(max_depth=8, n_estimators=1024, n_jobs=8, seed=9523)
    # 拟合模型
    model.fit(exog, endog)
    # 预测拟合值
    fit_results_array = model.predict(exog)
    return fit_results_array

def compute_rmse(pets, fitted_values):
    return [np.sqrt(mean_squared_error(pets[i], fitted_values[i])) for i in range(pets.shape[0])]

def run_dci_by_fisher_for_two_groups(pets_matrix, group):
    # OUTPUT_PATH = "./data/my_test_data/"
    # tmm_pets = np.loadtxt(OUTPUT_PATH+"/tmm_pets.csv", delimiter=",", dtype=np.float32)
    # group = ("k562", "k562", "mcf7", "mcf7")

    # * alternative model
    # * fit 1st group
    x = list()
    y = list()
    for i in range(len(pets_matrix)):
        x.append([i+1, 1, 0])
        x.append([i+1, 1, 0])
        y.append(pets_matrix[i][0])
        y.append(pets_matrix[i][1])
    x = np.array(x)
    y = np.array(y)
    exog = np.insert(x, 0, 1, axis=1) # 在第0列，添加常数项1, 用于计算截距, same as sm.add_constant(x)
    fitted_values1 = fit_predict(exog, y)
    fitted_values1 = fitted_values1.reshape((pets_matrix.shape[0], 2))

    # * fit 2nd group
    x = list()
    y = list()
    for i in range(len(pets_matrix)):
        x.append([i+1, 0, 1])
        x.append([i+1, 0, 1])
        y.append(pets_matrix[i][2])
        y.append(pets_matrix[i][3])
    x = np.array(x)
    y = np.array(y)
    exog = np.insert(x, 0, 1, axis=1) # 在第0列，添加常数项1, 用于计算截距, same as sm.add_constant(x)
    fitted_values2 = fit_predict(exog, y)
    fitted_values2 = fitted_values2.reshape((pets_matrix.shape[0], 2))

    # compute mean
    mean_1_2 = np.mean(fitted_values1, axis=1)
    mean_3_4 = np.mean(fitted_values2, axis=1)
    # stack two means
    fitted_values = np.stack((mean_1_2, mean_3_4), axis=1)
    print(fitted_values.shape)

    p_values = []
    loop_num = fitted_values.shape[0] - 1
    s1 = sum(fitted_values[:, 0])
    s2 = sum(fitted_values[:, 1])
    for gene_data in fitted_values:
        # 执行确切检验
        contingency_table = [gene_data,
                            [s1-gene_data[0], s2-gene_data[1]]]
        # ? 和减gene[0]再除以loop_num，p值全是1.0，可能是因为gene_data[0]相对于loop_num太小了。
        # contingency_table = [gene_data,
        #                     [(s1-gene_data[0])/loop_num, (s2-gene_data[1])/loop_num]]

        # res = chi2_contingency(contingency_table)
        _, p_value = fisher_exact(contingency_table, alternative='two-sided')
        # bers = barnard_exact(contingency_table, alternative='two-sided')
        # bers = boschloo_exact(contingency_table, alternative='two-sided')
        # 存储p值
        # p_values.append(res[1])
        p_values.append(p_value)
    p_values = np.array(p_values)
    # fdr = calculate_fdr(p_values)
    _, fdr = fdrcorrection(p_values, alpha=0.05, method='indep', is_sorted=False)

def run_dci_by_kde_for_two_groups(pets_matrix, groups, func_depth=0):
    time1 = time.time()
    print('\t'*func_depth+"<FUNC> run_dci_by_kde_for_two_groups")
    next_func_depth = func_depth + 1
    groups_set = list(set(groups))
    # * alternative model
    # * fit 1st group
    print('\t'*func_depth+"[INFO] fit 1st group alternative model")
    x = list()
    y = list()
    y_index_1 = [i for i,g in enumerate(groups) if g == groups_set[0]]
    for i in range(len(pets_matrix)):
        # x.append([i+1, 1, 0])
        # x.append([i+1, 0, 1])
        # x.append([i+1])
        # x.append([i+1])
        # x.append([i+1, 1, 0])
        # x.append([i+1, 1, 0])
        # y.append(pets_matrix[i][0])
        # y.append(pets_matrix[i][1])
        for _ in range(len(y_index_1)):
            x.append([i+1, 1, 0])
    y = pets_matrix[:, y_index_1].flatten()
    x = np.array(x)
    exog = np.insert(x, 0, 1, axis=1) # 在第0列，添加常数项1, 用于计算截距, same as sm.add_constant(x)
    fitted_values1 = fit_predict(exog, y)
    fitted_values1 = fitted_values1.reshape((pets_matrix.shape[0], len(y_index_1)))

    # * fit 2nd group
    print('\t'*func_depth+"[INFO] fit 2nd group alternative model")
    x = list()
    y = list()
    y_index_2 = [i for i,g in enumerate(groups) if g == groups_set[1]]
    for i in range(len(pets_matrix)):
        # x.append([i+1, 1, 0])
        # x.append([i+1, 0, 1])
        # x.append([i+1])
        # x.append([i+1])
        # x.append([i+1, 0, 1])
        # x.append([i+1, 0, 1])
        # y.append(pets_matrix[i][2])
        # y.append(pets_matrix[i][3])
        for _ in range(len(y_index_2)):
            x.append([i+1, 0, 1])
    y = pets_matrix[:, y_index_2].flatten()
    x = np.array(x)
    exog = np.insert(x, 0, 1, axis=1) # 在第0列，添加常数项1, 用于计算截距, same as sm.add_constant(x)
    fitted_values2 = fit_predict(exog, y)
    fitted_values2 = fitted_values2.reshape((pets_matrix.shape[0], len(y_index_2)))

    # cat two fitted_values
    fitted_values = np.concatenate((fitted_values1, fitted_values2), axis=1)

    print('\t'*func_depth+"[INFO] compute alternative model deviance")
    deviance = compute_rmse(pets_matrix, fitted_values)
    deviance = np.array(deviance)
    # print(deviance[:5])

    # * null model
    print('\t'*func_depth+"[INFO] fit null model")
    x0 = list()
    y0 = list()
    for i in range(len(pets_matrix)):
        # x0.append([i+1, 1, 0, 0, 0])
        # x0.append([i+1, 0, 1, 0, 0])
        # x0.append([i+1, 0, 0, 1, 0])
        # x0.append([i+1, 0, 0, 0, 1])
        # x0.append([i+1])
        # x0.append([i+1])
        # x0.append([i+1])
        # x0.append([i+1])
        # x0.append([i+1, 1, 1])
        # x0.append([i+1, 1, 1])
        # x0.append([i+1, 1, 1])
        # x0.append([i+1, 1, 1])
        # y0.append(pets_matrix[i][0])
        # y0.append(pets_matrix[i][1])
        # y0.append(pets_matrix[i][2])
        # y0.append(pets_matrix[i][3])
        for _ in range(len(groups)):
            x0.append([i+1, 1, 1])
    y_group1 = pets_matrix[:, y_index_1]
    y_group2 = pets_matrix[:, y_index_2]
    y0 = np.concatenate((y_group1, y_group2), axis=1).flatten()

    x0 = np.array(x0)

    exog0 = np.insert(x0, 0, 1, axis=1) # 在第0列，添加常数项1, 用于计算截距, same as sm.add_constant(x)

    fitted_values0 = fit_predict(exog0, y0)
    fitted_values0 = fitted_values0.reshape(pets_matrix.shape)
    # print(f"null fit value = alternative fit value? ", sum(fitted_values == fitted_values0))

    print('\t'*func_depth+"[INFO] compute null model deviance")
    deviance0 = compute_rmse(pets_matrix, fitted_values0)
    deviance0 = np.array(deviance0)
    # print(deviance0[:5])

    likelihood_ratio = deviance0/deviance
    # print(f"number of likelihood_ratio > 0: {sum(likelihood_ratio > 1)}/{len(likelihood_ratio)}")

    # likelihood_ratio = deviance0 - deviance
    # print(f"number of likelihood_ratio > 0: {sum(likelihood_ratio > 0)}/{len(likelihood_ratio)}")
    # likelihood_ratio = np.abs(deviance0 - deviance)

    print('\t'*func_depth+"[INFO] compute P values and FDRs")
    # cumulative_prob = chi2.sf(x=likelihood_ratio, df=degree_freedom_test) # Survival function (also defined as 1 - cdf, but sf is sometimes more accurate).
    cumulative_prob = compute_P_value_by_kde_scipy(likelihood_ratio, next_func_depth)
    cumulative_prob = np.where(cumulative_prob > 1, 1., cumulative_prob)

    # p值小的差异显著
    # print(f"cumulative_prob.shape: {cumulative_prob.shape}")
    # print(cumulative_prob[:5])
    # print(sum(cumulative_prob))
    # print(sum(cumulative_prob < 0.05))

    _, fdr = fdrcorrection(cumulative_prob, alpha=0.05, method='indep', is_sorted=False)
    print('\t'*func_depth+f"<TIME> run_dci_by_kde_for_two_groups cost {time.time()-time1} seconds")
    return cumulative_prob, fdr


if __name__ == "__main__":
    ROOT_PATH = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    OUTPUT_PATH = ROOT_PATH+"/data/my_test_data/"
    SIM_PATH = ROOT_PATH+"/data/simulated_data/"
    # run_dci_by_fisher()
    # run_dci_temp()


