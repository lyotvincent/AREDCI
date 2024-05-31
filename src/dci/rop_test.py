

def find_neighbor(loopinfos, row_index, neighbor_num) -> list:
    """
    find the neighbor of a loop
    loopinfos: LoopInfos of a sample, the loop information of every sample are same
    row_index: the index of the loop we want to find the neighbor
    neighbor_num: the number of neighbor of each anchor
    """
    #

def rop_test(loopinfos, pvalues):
    """
    rth ordered P-value (rOP) statistic
    loopinfos: matrix of LoopInfo, each row is a loop, each column is a sample
    pvalues: list of pvalues, each pvalue is a float corresponding to a row in loopinfos
    1. 从一个环的每个锚点两边各找一个环，加上本体就是5个环，和对比的sample形成5v5。会得到5个P值。
    2. 我们设定一个阈值p ∈ (0, 1]，它代表了被认为是DCI所需的、在两种条件下显著变化的KNN的比例。
    3. 然后计算r，它是小于或等于p × W2的最大整数值。
    4. 首先对所有2W2个P值进行从小到大的排序，得到一组顺序统计量。
    5. 然后选取第r个顺序统计量。
    6. 接着用这个选定的顺序统计量去与具有α = r 和 β = 2W2 − r + 1自由度的beta分布进行比较或者拟合，以评估在零假设下观察到该顺序统计量的可能性。如果观察到的顺序统计量在β分布的尾部区域，根据显著性水平，则可能拒绝零假设（假设无差异），意味着存在足够的证据支持所设定的生物学条件间的差异。
    """
    # * 1. 从一个环的每个锚点两边各找一个环，加上本体就是5个环，和对比的sample形成5v5。会得到5个P值。



