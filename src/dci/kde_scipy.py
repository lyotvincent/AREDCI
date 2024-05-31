from sklearn.neighbors import KernelDensity
from scipy.stats import gaussian_kde
from scipy.integrate import quad
import numpy as np
import time

# 计算SF
def _evaluate_sf(
        pdf: gaussian_kde, # probability density function
        data_point_a,
        data_point_b
    ):
    # upper_bound = np.inf  # KDE的上界
    # 调整积分精度：quad函数的epsabs和epsrel参数可以用来控制积分的精度。减小这些值可以提高积分速度，但可能会降低结果的精度。
    sf_integral, _ = quad(func=lambda x: pdf(x),
                          a=data_point_a,
                          b=data_point_b,
                          epsabs=1e-6, # 设置绝对误差容忍度
                          epsrel=1e-4, # 设置相对误差容忍度，把epsabs和epsrel都设大一点，可以加快积分速度，比如这个设置节约了1/3的时间
                          limit=64) # Time: 33.829s -> Time: 20.684s
    return sf_integral

def _fit_kde_by_scipy(data, low_threshold=0.2, high_threshold=0.8, func_depth=0):
    print(func_depth*"\t"+f"<FUNC> _fit_kde, fitting KDE to {data.shape[0]} observations")
    # 按数值截掉数据最大的5%和最小的5%，以避免极端值对KDE的影响
    data = np.sort(data)
    data = data[int(len(data) * low_threshold):int(len(data) * high_threshold)]
    # 随机不重复抽样10%的data
    data = np.random.choice(data, size=int(len(data) * 0.1), replace=False)
    kde_scipy = gaussian_kde(data,
                             bw_method='silverman') # "scott", "silverman", or a scalar
    return kde_scipy

def _plot_kde(kde: gaussian_kde, min_x, max_x):
    # 生成网格以评估密度
    grid = np.linspace(min_x, max_x, 500)
    density = kde.evaluate(grid)

    # 绘制 KDE 曲线
    import matplotlib.pyplot as plt
    plt.plot(grid, density)
    plt.show()

def _compute_P_value(pdf, data, func_depth=0):
    print(func_depth*"\t"+f"<FUNC> _compute_P_value")
    # 降序排列数据
    sorted_data_indices = np.argsort(-data)
    sorted_data = data[sorted_data_indices]

    # 创建一个列表，用于存储每个数据点的累积概率密度sf，避免重复计算。
    # 第一列值是降序排序后的数据点，第二列值是sf值survival function
    sf_sorted_record = list()
    sf_sorted_record.append([np.inf, 0])
    i = 0
    time1 = time.time()
    for point in sorted_data:
        i += 1
        temp_sf = _evaluate_sf(pdf, point, sf_sorted_record[-1][0])
        if i % 1000 == 0:
            print(func_depth*"\t"+f"compute_P_value: {i}, time: {time.time() - time1:.3f}s")
            time1 = time.time()
        sf_sorted_record.append([point, temp_sf + sf_sorted_record[-1][1]])
    sf_sorted_record = np.array(sf_sorted_record)

    # 恢复正常顺序的数据，并将对应的P值存储在sf_list的对应位置中
    sf_list = np.empty(len(sf_sorted_record)-1)
    sf_list[sorted_data_indices] = sf_sorted_record[:, 1][1:]

    # threshold = np.percentile(sf_list, 5)  # 设置一个例如第95百分位的阈值
    # print(f"Threshold: {threshold:.3f}")
    # threshold = 0.05
    # significant_indices = np.where(sf_list < threshold)[0]

    # # 输出显著差异的索引及其对应的密度值
    # for idx in significant_indices:
    #     print(f"Index {idx}, {data[idx]} has a density value of {sf_list[idx]:.3f}, suggesting a potentially significant difference.")

    # 注意：实际应用中，显著性判断通常与假设检验相关联，而不仅仅是
    return sf_list

def compute_P_value_by_kde_scipy(data, func_depth=0):
    print(func_depth*"\t"+"<FUNC> compute_P_value_by_kde_scipy")
    # print(data.shape)
    time1 = time.time()
    next_func_depth = func_depth + 1

    pdf = _fit_kde_by_scipy(data, func_depth=next_func_depth)
    # _plot_kde(pdf, min(data), max(data))
    p_values = _compute_P_value(pdf, data, next_func_depth)
    print(func_depth*"\t"+f"<TIME> compute_P_value_by_kde_scipy cost {time.time()-time1} seconds")
    return p_values


if __name__ == "__main__":
    import time
    time1 = time.time()
    # 假设我们有这样一个数组a，表示差异
    np.random.seed(42)  # 设置随机种子以重现结果
    n = 1000  # 示例数据长度
    a = np.abs(np.random.normal(size=n, scale=1, loc=0))  # 示例差异数组

    compute_P_value_by_kde_scipy(a)

    print(f"Time: {time.time() - time1:.3f}s")

