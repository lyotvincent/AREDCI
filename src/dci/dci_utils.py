import numpy as np


@DeprecationWarning
def plot_pet_and_fitted_value(obs, fit):
    import matplotlib.pyplot as plt
    bins_range = (min(obs), max(obs))
    print(bins_range)
    plt.hist(obs, bins=np.arange(bins_range[0], bins_range[1]), alpha=0.5, label='Observed')

    # # Plot the distribution of fitted values
    plt.hist(fit, bins=np.arange(bins_range[0], bins_range[1], 0.2), alpha=0.5, label='Fitted')

    # # Add a legend
    plt.legend()
    # # x轴标签
    # # plt.xticks(np.arange(1, 11, 1))
    # # plt.xlim(0.5, 10.5)
    # # Show the plot
    plt.show()



@DeprecationWarning
def compute_unit_deviance(observed_value, fitted_value, dispersion):
    ONE_TENTHOUSANDTH = 1E-4
    MILDLY_LOW_VALUE = 1E-8
    ONE_MILLION = 1E6
    fitted_value = max(fitted_value, MILDLY_LOW_VALUE)

    observed_value += MILDLY_LOW_VALUE
    fitted_value += MILDLY_LOW_VALUE
    # unit_deviance = 2 * (observed_value * np.log(observed_value / fitted_value) - (observed_value - fitted_value))
    # # unit_deviance = unit_deviance / dispersion /2
    # return unit_deviance
    if dispersion < ONE_TENTHOUSANDTH:
        resid = observed_value - fitted_value
        return 2 * (observed_value * np.log(observed_value / fitted_value) - resid - 0.5 * resid * resid * dispersion * (1 + dispersion * (2 / 3 * resid - observed_value)))
    else:
        product = fitted_value * dispersion
        if product > ONE_MILLION:
            return 2 * ((observed_value - fitted_value) / fitted_value - np.log(observed_value / fitted_value)) * fitted_value / (1 + product)
        else:
            invphi = 1 / dispersion
            return 2 * (observed_value * np.log(observed_value / fitted_value) + (observed_value + invphi) * np.log((fitted_value + invphi) / (observed_value + invphi)))

@DeprecationWarning
def nbinomDeviance(observed_value, fitted_value, dispersion=0.2, weights=1):
    """
    Calculate the deviance value for a given observed value, fitted value, and dispersion.

    Parameters:
    observed_value (numpy.ndarray): The observed value.
    fitted_value (numpy.ndarray): The fitted value.
    dispersion (float): The dispersion parameter. Default is 0.2.

    Returns:
    float: The deviance value.
    """
    # * https://bioconductor.org/packages/release/bioc/manuals/edgeR/man/edgeR.pdf
    # * 里说了是unit deviance的加权和，而weights默认全是1。
    # * 而对于每个unit deviance，dispersion都是一样的，比如本例子里的0.2。每个unit deviance都是0.2。
    # deviance1 = compute_unit_deviance(observed_value[0], fitted_value[0], dispersion)
    # deviance2 = compute_unit_deviance(observed_value[1], fitted_value[1], dispersion)
    # deviance_sum = np.dot([deviance1, deviance2], weights)
    deviance_sum = 0.
    for i in range(len(observed_value)):
        deviance = compute_unit_deviance(observed_value[i], fitted_value[i], dispersion)
        # print(deviance)
        deviance_sum += deviance * weights[i]
    return deviance_sum

@DeprecationWarning
def compute_deviance(pets, fitted_values, dispersion=0.2, weights=np.array([1,1,1,1])):
    return [nbinomDeviance(pets[i], fitted_values[i], dispersion=dispersion, weights=weights) for i in range(pets.shape[0])]

@DeprecationWarning
def compute_residual(sample_num, group_num):
    return sample_num - group_num

def compute_dispersion(values):
    mean = np.mean(values)
    variance = np.var(values)
    dispersion = (variance - mean) / (mean**2)
    return dispersion

def compute_difference_between_group_mean(pets, group):
    """
    compute mean of each group and then compute difference between group mean
    """
    group_set = list(set(group))
    diffs = list()
    for interac in pets:
        group_dict = {g:[] for g in group_set}
        for i in range(len(interac)):
            group_dict[group[i]].append(interac[i])
        group_dict = {g:np.mean(group_dict[g]) for g in group_set}
        diffs.append(abs(group_dict[group_set[0]] - group_dict[group_set[1]]))
    return np.array(diffs)

def plot_difference_between_fit_and_obs(pets, fitted_pets):
    # reshape to (-1, 1)
    pets = pets.reshape((-1))
    fitted_pets = fitted_pets.reshape((-1))
    diff = pets - fitted_pets
    # normalize
    diff = (diff - np.mean(diff)) / np.std(diff)
    print(f"number of diff > 0: {sum(diff > 0)}")
    print(f"number of diff < 0: {sum(diff < 0)}")
    import matplotlib.pyplot as plt
    bins_range = (min(diff), max(diff))
    print(bins_range)
    plt.hist(diff, bins=np.arange(bins_range[0], bins_range[1], 0.2), alpha=0.5, label='Difference between fit and obs')
    plt.show()

def plot_mse(mse):
    import matplotlib.pyplot as plt
    bins_range = (min(mse), max(mse))
    print(bins_range)

    # # Plot the distribution of fitted values
    plt.hist(mse, bins=np.arange(bins_range[0], bins_range[1], 0.2), alpha=0.5, label='mse')

    # # Add a legend
    plt.legend()
    # # x轴标签
    # # plt.xticks(np.arange(1, 11, 1))
    # # plt.xlim(0.5, 10.5)
    # # Show the plot
    plt.show()

@DeprecationWarning
def goodness_of_fit_test(observed_data):
    # 假设你有一组连续的观测数据
    observed_data = np.random.chisquare(5, size=1000)  # 生成服从自由度为5的卡方分布的随机样本

    # 确定区间（bins）的数量和边界
    bin_edges = np.linspace(0, 50, num=50 + 1)  # 根据实际情况调整区间数量

    # 计算每个区间的频数
    observed_counts, _ = np.histogram(observed_data, bins=bin_edges)

    # 计算每个区间在给定自由度下的理论期望频数
    df = 1  # 自由度，这里假设你知道或猜测的自由度
    bin_widths = np.diff(bin_edges)
    expected_counts = [chi2.pdf((bin_edges[i] + bin_edges[i+1]) / 2, df) * len(observed_data) * bin_widths[i] for i in range(len(bin_edges[:-1]))]

    # 执行卡方拟合优度检验
    chi2_statistic, p_value = chisquare(observed_counts, expected_counts)

    print("卡方统计量:", chi2_statistic)
    print("p值:", p_value)

    alpha = 0.05
    if p_value < alpha:
        print(f"在显著性水平{alpha}下，拒绝数据服从自由度为{df}的卡方分布的原假设")
    else:
        print(f"在显著性水平{alpha}下，无法拒绝数据可能服从自由度为{df}的卡方分布的原假设")

# Calculate the dispersion for each sample
# 在负二项分布中，方差等于均值加上均值的平方乘以离散度，所以反过来可以求离散度alpha，离散度和r成反比
# https://zhuanlan.zhihu.com/p/111632687
# edgeR的公式也是这样的
def compute_dispersion(values):
    mean = np.mean(values)
    variance = np.var(values)
    dispersion = (variance - mean) / (mean**2)
    return dispersion

@DeprecationWarning
def run_dci_temp():
    # tmm_pets = np.load("D:/workspace/dci/data/my_test_data/tmm_pets.npy")
    tmm_pets = np.loadtxt(OUTPUT_PATH+"/tmm_pets.csv", delimiter=",", dtype=np.float32)
    # print(f"tmm_pets.shape: {tmm_pets.shape}") # (40373, 4)
    # print(tmm_pets[:5])
    group = ("k562", "k562", "mcf7", "mcf7")

    # dispersion = compute_dispersion(np.concatenate(tmm_pets))
    # print(f"dispersion: {dispersion}")

    # * alternative model
    # * fit 1st group
    x = list()
    y = list()
    for i in range(len(tmm_pets)):
        # x.append([i+1, 1, 0])
        # x.append([i+1, 0, 1])
        # x.append([i+1])
        # x.append([i+1])
        x.append([i+1, 1, 0])
        x.append([i+1, 1, 0])
        y.append(tmm_pets[i][0])
        y.append(tmm_pets[i][1])
    x = np.array(x)
    y = np.array(y)
    exog = np.insert(x, 0, 1, axis=1) # 在第0列，添加常数项1, 用于计算截距, same as sm.add_constant(x)
    # print(f"exog.shape: {exog.shape}")
    fitted_values1 = fit_predict(exog, y)
    # print(f"fitted_values.shape: {fitted_values1.shape}")
    fitted_values1 = fitted_values1.reshape((tmm_pets.shape[0], 2))

    # * fit 2nd group
    x = list()
    y = list()
    for i in range(len(tmm_pets)):
        # x.append([i+1, 1, 0])
        # x.append([i+1, 0, 1])
        # x.append([i+1])
        # x.append([i+1])
        x.append([i+1, 0, 1])
        x.append([i+1, 0, 1])
        y.append(tmm_pets[i][2])
        y.append(tmm_pets[i][3])
    x = np.array(x)
    y = np.array(y)
    exog = np.insert(x, 0, 1, axis=1) # 在第0列，添加常数项1, 用于计算截距, same as sm.add_constant(x)
    # print(f"exog.shape: {exog.shape}")
    fitted_values2 = fit_predict(exog, y)
    # print(f"fitted_values.shape: {fitted_values2.shape}")
    fitted_values2 = fitted_values2.reshape((tmm_pets.shape[0], 2))

    # cat two fitted_values
    fitted_values = np.concatenate((fitted_values1, fitted_values2), axis=1)
    # print(f"fitted_values.shape after reshape: {fitted_values.shape}")
    # print(fitted_values[:5])

    # plot_pet_and_fitted_value(y, fitted_values)
    # plot_difference_between_fit_and_obs(tmm_pets, fitted_values)


    # deviance = compute_deviance(tmm_pets, fitted_values, dispersion)
    deviance = compute_rmse(tmm_pets, fitted_values)
    deviance = np.array(deviance)
    # print(f"deviance.shape: {deviance.shape}")
    print(deviance[:5])

    # plot_mse(deviance)
    # goodness_of_fit_test(deviance)
    # with open(OUTPUT_PATH+"/deviance.csv", "w") as f:
    #     for d in deviance:
    #         f.write(f"{d}\n")


    residual = compute_residual(tmm_pets.shape[1], 2)
    residual = np.array([residual]*tmm_pets.shape[0])

    # * null model
    x0 = list()
    y0 = list()
    for i in range(len(tmm_pets)):
        # x0.append([i+1, 1, 0, 0, 0])
        # x0.append([i+1, 0, 1, 0, 0])
        # x0.append([i+1, 0, 0, 1, 0])
        # x0.append([i+1, 0, 0, 0, 1])
        # x0.append([i+1])
        # x0.append([i+1])
        # x0.append([i+1])
        # x0.append([i+1])
        x0.append([i+1, 1, 1])
        x0.append([i+1, 1, 1])
        x0.append([i+1, 1, 1])
        x0.append([i+1, 1, 1])
        y0.append(tmm_pets[i][0])
        y0.append(tmm_pets[i][1])
        y0.append(tmm_pets[i][2])
        y0.append(tmm_pets[i][3])

    x0 = np.array(x0)
    y0 = np.array(y0)

    exog0 = np.insert(x0, 0, 1, axis=1) # 在第0列，添加常数项1, 用于计算截距, same as sm.add_constant(x)
    # print(f"exog0.shape: {exog0.shape}")

    fitted_values0 = fit_predict(exog0, y0)
    # print(f"fitted_values0.shape: {fitted_values0.shape}")
    fitted_values0 = fitted_values0.reshape(tmm_pets.shape)
    # print(f"fitted_values0.shape after reshape: {fitted_values0.shape}")
    # print(fitted_values0[:5])
    print(f"null fit value = alternative fit value? ", sum(fitted_values == fitted_values0))

    # plot_difference_between_fit_and_obs(tmm_pets, fitted_values0)

    # deviance0 = compute_deviance(tmm_pets, fitted_values0, dispersion)
    deviance0 = compute_rmse(tmm_pets, fitted_values0)
    deviance0 = np.array(deviance0)
    # print(f"deviance0.shape: {deviance0.shape}")
    print(deviance0[:5])

    # plot_mse(deviance0)

    residual0 = compute_residual(tmm_pets.shape[1], 1)
    residual0 = np.array([residual0]*tmm_pets.shape[0])

    # likelihood_ratio = deviance0 - deviance # "Likelihood Ratio"（似然比）
    # likelihood_ratio = np.abs(deviance0 - deviance)
    likelihood_ratio = deviance0/deviance
    # likelihood_ratio = -2 * np.log(deviance0/deviance)
    # likelihood_ratio = compute_mse(fitted_values, fitted_values0)
    # plot_difference_between_fit_and_obs(fitted_values, fitted_values0)
    # plot_mse(likelihood_ratio)
    np.savetxt(f'{OUTPUT_PATH}/likelihood_ratio.csv', likelihood_ratio.T, fmt='%f', delimiter=',')

    # with open(OUTPUT_PATH+"/likelihood_ratio.csv", "w") as f:
    #     for d in likelihood_ratio:
    #         f.write(f"{d}\n")
    print("number of likelihood_ratio > 0: ", sum(likelihood_ratio > 0))

    degree_freedom_test = residual0 - residual # 似然比检验中用到的自由度通常是指两个模型之间参数个数的差异

    # cumulative_prob = chi2.sf(x=likelihood_ratio, df=degree_freedom_test) # Survival function (also defined as 1 - cdf, but sf is sometimes more accurate).
    cumulative_prob = compute_P_value_by_kde_scipy(likelihood_ratio)

    # p值小的差异显著
    print(f"cumulative_prob.shape: {cumulative_prob.shape}")
    print(cumulative_prob[:5])
    print(sum(cumulative_prob))
    print(sum(cumulative_prob < 0.05))
    np.savetxt(f'{OUTPUT_PATH}/cumulative_prob.csv', cumulative_prob.T, fmt='%f', delimiter=',')

    # cumulative_prob = np.loadtxt(OUTPUT_PATH+"/cumulative_prob.csv", delimiter=",", dtype=np.float32)
    # fdr小的差异显著
    # fdr = calculate_fdr(cumulative_prob)
    _, fdr = fdrcorrection(cumulative_prob, alpha=0.05, method='indep', is_sorted=False)
    np.savetxt(f'{OUTPUT_PATH}/fdr.csv', fdr.T, fmt='%f', delimiter=',')
    print(f"fdr.shape: {fdr.shape}")
    print(fdr[:5])
    print(sum(fdr < 0.05))
    # get index when fdr < 0.05
    index_fdr = np.where(fdr < 0.05)
    # print(f"fdr < 0.05 index: {index_fdr}")
    # print(tmm_pets[index])
    # print(f"fdr tmm_pets: {tmm_pets[index]}")
    # print(f"fdr fdr: {fdr[index]}")
    # print(f"fdr cumulative_prob: {cumulative_prob[index]}")
    # print(f"fdr fitted_values[index]: {fitted_values[index]}")
    # print(f"fdr fitted_values0[index]: {fitted_values0[index]}")

    # difference between group mean
    diffs = compute_difference_between_group_mean(tmm_pets, group)
    # get topmax 10 indexes of diffs
    index_dif = np.argsort(diffs)
    # print(f"topmax indexes of diffs: {sorted(index)}")
    # print(f"topmax diffs: {diffs[index]}")
    # print(f"topmax tmm_pets: {tmm_pets[index]}")
    # print(f"topmax fdr: {fdr[index]}")
    # print(f"topmax cumulative_prob: {cumulative_prob[index]}")
    # print(f"topmax fitted_values[index]: {fitted_values[index]}")
    # print(f"topmax fitted_values0[index]: {fitted_values0[index]}")
    index_dif = index_dif[-sum(fdr < 0.05):]



    print(len(index_dif))
    print(len(index_fdr[0]))
    print("len(set(index_dif).intersection(set(index_fdr[0])))")
    print(len(set(index_dif).intersection(set(index_fdr[0]))))

    print("X11013")
    print(fdr[11013])
    print(cumulative_prob[11013])
    print(tmm_pets[11013])
    print(fitted_values[11013])
    print(fitted_values0[11013])
    print(likelihood_ratio[11013])
    print(deviance[11013])
    print(deviance0[11013])
    print("X11014")
    print(fdr[11014])
    print(cumulative_prob[11014])
    print(tmm_pets[11014])
    print(fitted_values[11014])
    print(fitted_values0[11014])
    print(likelihood_ratio[11014])
    print(deviance[11014])
    print(deviance0[11014])
    print("X11015")
    print(fdr[11015])
    print(cumulative_prob[11015])
    print(tmm_pets[11015])
    print(fitted_values[11015])
    print(fitted_values0[11015])
    print(likelihood_ratio[11015])
    print(deviance[11015])
    print(deviance0[11015])


def load_sim_data(x_path, y_path):
    f = open(x_path, "r")
    lines = f.readlines()[1:]
    f.close()
    sim_data = list()
    for line in lines:
        line = [int(i) for i in line.strip().split(",")[1:]]
        sim_data.append(line)

    f = open(y_path, "r")
    lines = f.readlines()[1:]
    f.close()
    y = list()
    for line in lines:
        line = [1. if i == "TRUE" else 0. for i in line.strip().split(",")[1:]]
        assert sum(line) == 4 or sum(line) == 0
        line = 1 if sum(line) == 4 else 0
        y.append(line)
    return np.array(sim_data), np.array(y)

def compute_metrics(pvalues, fdrs, labels):
    print(f"sum(fdrs < 0.05): {sum(fdrs < 0.05)}")
    # make prediction=1 when fdrs < 0.05
    predictions = np.zeros(fdrs.shape)
    predictions[fdrs < 0.05] = 1
    # compute metrics
    acc = accuracy_score(labels, predictions)
    precision = precision_score(labels, predictions)
    recall = recall_score(labels, predictions)
    f1 = f1_score(labels, predictions)
    fbeta = fbeta_score(y_true=labels, y_pred=predictions, beta=2)

    auroc = roc_auc_score(labels, 1-fdrs)
    fpr, tpr, thresholds = roc_curve(labels, 1-fdrs)
    auroc_by_auc = auc(fpr, tpr)

    auprc = average_precision_score(labels, 1-fdrs)
    precision_point, recall_point, thresholds = precision_recall_curve(labels, 1-fdrs, pos_label=1)
    # precision_point[(recall_point==0)] = 1.0
    auprc_by_auc = auc(recall_point, precision_point)
    print(f"accuracy: {acc}, precision: {precision}, recall: {recall}, f1: {f1}, fbeta: {fbeta}, auroc: {auroc}, auroc_by_auc: {auroc_by_auc}, auprc: {auprc}, auprc_by_auc: {auprc_by_auc}")
    return fpr, tpr, auroc_by_auc, precision_point, recall_point, auprc_by_auc

def load_edgeR_sim_results():
    f = open("D:/workspace/dci/data/simulated_data/edgeR_sim_results.csv", "r")
    lines = f.readlines()
    f.close()
    pvalues = list()
    fdrs = list()
    for line in lines[1:]:
        line = line.strip().split(",")
        pvalues.append(float(line[4]))
        fdrs.append(float(line[5]))
    return np.array(pvalues), np.array(fdrs)

def plot_sim_results():
    import matplotlib.pyplot as plt
    import seaborn as sns
    # ax = plt.gca()
    my_sim_results = np.loadtxt(OUTPUT_PATH+"/my_sim_results.csv", delimiter=",", dtype=np.float32)
    pvalues, fdrs1 = my_sim_results[:, 0], my_sim_results[:, 1]
    fpr1, tpr1, auroc_by_auc1, precision_point1, recall_point1, auprc_by_auc1 = compute_metrics(pvalues, fdrs1, sim_labels)
    # plt.plot(precision_point, recall_point, label=r"%s (AUC=%.2f)" % ("mine", auprc_by_auc), color="red")
    # pr_display = PrecisionRecallDisplay(precision=precision_point, recall=recall_point)
    # pr_display.plot(ax=ax, name="mine")

    pvalues, fdrs2 = load_edgeR_sim_results()
    fpr2, tpr2, auroc_by_auc2, precision_point2, recall_point2, auprc_by_auc2 = compute_metrics(pvalues, fdrs2, sim_labels)
    # plt.plot(precision_point, recall_point, label=r"%s (AUC=%.2f)" % ("edgeR", auprc_by_auc), color="green")
    # plt.plot(mean_fpr, mean_tpr, color=color, label=r'%s (AUC=%0.4f$\pm$%0.4f)' % ("dci", mean_auc, std_auc), lw=2, alpha=.8)

    plt.subplot(3, 2, 1)
    ax = plt.gca()
    pr_display1 = PrecisionRecallDisplay(precision=precision_point1, recall=recall_point1)
    pr_display1.plot(ax=ax, name="mine", color="red")
    pr_display2 = PrecisionRecallDisplay(precision=precision_point2, recall=recall_point2)
    pr_display2.plot(ax=ax, name="edgeR", color="green")
    # pr_display1 = PrecisionRecallDisplay.from_predictions(y_true=sim_labels, y_pred=1-fdrs1, name="mine", ax=ax, color="red")
    # pr_display2 = PrecisionRecallDisplay.from_predictions(y_true=sim_labels, y_pred=1-fdrs2, name="edgeR", ax=ax, color="green")
    plt.title('PR Curve (PrecisionRecallDisplay)')
    plt.legend()
    # 创建第二个子图
    plt.subplot(3, 2, 2)  # 1行2列的子图中的第2个
    plt.plot(precision_point1, recall_point1, label=r"%s (AUC=%.2f)" % ("mine", auprc_by_auc1), color="red")
    plt.plot(precision_point2, recall_point2, label=r"%s (AUC=%.2f)" % ("edgeR", auprc_by_auc2), color="green")
    plt.title('PR Curve (plt.plot)')
    plt.legend()

    plt.subplot(3, 2, 3)
    fnr1 = 1 - tpr1
    plt.plot(fpr1, fnr1, label="mine", color="red")
    fnr2 = 1 - tpr2
    plt.plot(fpr2, fnr2, label="edgeR", color="green")
    plt.title('Detection Error Tradeoff curve')
    plt.legend(loc="best")

    plt.subplot(3, 2, 5)
    plt.hist(fdrs1, bins=30, alpha=0.5, label='Group 1')
    plt.hist(fdrs2, bins=30, alpha=0.5, label='Group 2')
    plt.legend(loc='upper right')
    plt.title('Histogram of Predictions')

    # 创建密度图
    plt.subplot(3, 2, 6)
    sns.kdeplot(fdrs1, label='Group 1')
    sns.kdeplot(fdrs2, label='Group 2')
    plt.title('Density Plot of Predictions')

    plt.tight_layout()
    plt.show()

def calculate_fdr(p_values):
    # https://zhuanlan.zhihu.com/p/328926817
    # Step 1: Sort the p-values
    print(p_values[0])
    print(p_values[7])
    print(p_values[3301])
    sorted_p_indices = np.argsort(-p_values)
    sorted_p_values = p_values[sorted_p_indices]
    print(f"the index0 p: {np.where(sorted_p_values == 0.579017)}")
    print(f"the index7 p: {np.where(sorted_p_values == 0.018143)}")
    print(f"the index3301 p: {np.where(sorted_p_values == 0.000923)}")

    # Step 2: Calculate initial FDR values
    p_count = len(p_values)
    print(f"p_count: {p_count}")
    ranks = np.arange(1, p_count + 1)[::-1]
    print(f"ranks: {ranks[:5]}")
    initial_fdr_values = sorted_p_values * p_count / ranks
    print(f"14416 index: {initial_fdr_values[14416]}")
    print(f"39647 index: {initial_fdr_values[39647]}")
    print(f"40335 index: {initial_fdr_values[40335]}")

    # Step 3: Ensure FDR is non-decreasing
    fdr_values = np.minimum.accumulate(initial_fdr_values)

    # Put the FDR values back in the original order
    fdr_values_original_order = np.empty(p_count)
    fdr_values_original_order[sorted_p_indices] = fdr_values

    return fdr_values_original_order

def compute_mse(pets, fitted_values):
    return [mean_squared_error(pets[i], fitted_values[i]) for i in range(pets.shape[0])]

def compute_mae(pets, fitted_values):
    return [mean_absolute_error(pets[i], fitted_values[i]) for i in range(pets.shape[0])]

if __name__ == "__main__":
    
    sim_pets, sim_labels = load_sim_data(SIM_PATH+"/sim_counts.csv", SIM_PATH+"/sim_DE_labels.csv")
    group = ("k562", "k562", "mcf7", "mcf7")
    # pvalues, fdrs = run_dci_for_two_groups(sim_pets, group)
    plot_sim_results()


