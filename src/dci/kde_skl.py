from sklearn.neighbors import KernelDensity
from scipy.stats import gaussian_kde
from scipy.integrate import quad
import numpy as np
from functools import partial

# 计算整个样本空间的概率密度函数
def _evaluate_density(
        pdf: KernelDensity, # probability density function
        data_point
    ):
    # global pdf
    if np.isscalar(data_point):
        data_point = np.array([data_point])
    log_pd = pdf.score_samples(data_point[:, np.newaxis])
    pd = np.exp(log_pd)
    return pd

# 计算SF
def _evaluate_sf(
        pdf, # probability density function
        data_point
    ):
    upper_bound = np.inf  # KDE的上界
    partial_func = partial(_evaluate_density, pdf)
    # 调整积分精度：quad函数的epsabs和epsrel参数可以用来控制积分的精度。减小这些值可以提高积分速度，但可能会降低结果的精度。
    sf_integral, _ = quad(func=partial_func,
                          a=data_point,
                          b=upper_bound,
                          epsabs=1e-6, # 设置绝对误差容忍度
                          epsrel=1e-4, # 设置相对误差容忍度，把epsabs和epsrel都设大一点，可以加快积分速度，比如这个设置节约了1/3的时间
                          limit=1000) # Time: 33.829s -> Time: 20.684s
    return sf_integral

# 计算CDF
@DeprecationWarning
def _evaluate_cdf(data_point):
    lower_bound = -np.inf  # KDE的下界
    cdf_integral, _ = quad(_evaluate_density, lower_bound, data_point)
    return cdf_integral

# 使用sklearn的KernelDensity进行拟合，并计算P值。
def _fit_kde_by_skl(data):
    print(f"<FUNC> _fit_kde, fitting KDE to {data.shape[0]} observations")
    # kde_skl = KernelDensity(bandwidth=0.2, kernel='gaussian')
    kde_skl = KernelDensity(bandwidth='scott', kernel='gaussian')
    kde_skl.fit(data)  # KDE需要一维数据转为二维（样本x特征）
    print(f"[INFO] bandwidth: {kde_skl.bandwidth}")
    return kde_skl

def _compute_P_value(pdf, a):
    print(f"<FUNC> _compute_P_value")
    # 对于数组a中的每个点，计算其概率密度
    # a = a[:, np.newaxis]  # KDE需要一维数据转为二维（样本x特征）
    # densities = [evaluate_density(point)[0] for point in a]
    densities = []
    i = 0
    for point in a:
        i += 1
        # print(f"compute_P_value: {i}")
        # if i % 100 == 0:
        #     print(f"compute_P_value: {i}")
        densities.append(_evaluate_sf(pdf, point))
    densities = np.array(densities)

    # 虽然不能直接计算P值，但可以将每个点的密度与所有点的密度分布相比，
    # 寻找低于某个阈值的点，这些点可能代表显著的极端值
    threshold = np.percentile(densities, 5)  # 设置一个例如第95百分位的阈值
    print(f"Threshold: {threshold:.3f}")
    threshold = 0.05
    significant_indices = np.where(densities < threshold)[0]

    # 输出显著差异的索引及其对应的密度值
    for idx in significant_indices:
        print(f"Index {idx}, {a[idx]} has a density value of {densities[idx]:.3f}, suggesting a potentially significant difference.")

    # 注意：实际应用中，显著性判断通常与假设检验相关联，而不仅仅是比较密度大小。
    # 在这里仅作为一个直观展示KDE如何用于发现潜在异常点的例子。
    # 上述代码首先对数组 a 进行了KDE拟合，然后计算了每个数据点的概率密度。通过设定一个密度阈值（这里是第95个百分位数），我们可以找出那些在密度上相对较低、因此在统计意义上可能是显著差异的数据点。不过，请注意这种方法是简化的，并且没有严格遵循统计学意义上的假设检验流程来得出精确的P值。如果要进行更严格的显著性测试，您可能需要结合具体的统计假设检验方法来设计算法。

def compute_P_value_by_kde_skl(data):
    # 已知a是数组，如果a是一维的，把a变成二维的
    if data.ndim == 1:
        data = data[:, np.newaxis]
    print(data.shape)

    pdf = _fit_kde_by_skl(data)
    # pdf = _fit_kde_by_scipy(data)
    _compute_P_value(pdf, data)



if __name__ == "__main__":
    import time
    time1 = time.time()
    # 假设我们有这样一个数组a，表示差异
    np.random.seed(42)  # 设置随机种子以重现结果
    n = 1000  # 示例数据长度
    a = np.abs(np.random.normal(size=n, scale=1, loc=0))  # 示例差异数组

    compute_P_value_by_kde_skl(a)
    print(f"Time: {time.time() - time1:.3f}s")

