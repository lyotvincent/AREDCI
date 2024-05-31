import networkx as nx
from collections import deque, defaultdict
from scipy.stats import norm, combine_pvalues
import numpy as np
from statsmodels.stats.multitest import fdrcorrection
import time

def build_graph(loopinfos, func_depth=0):
    print('\t'*func_depth+"<FUNC> build_graph")
    time1 = time.time()
    graph = nx.Graph()
    for loopinfo in loopinfos:
        graph.add_edge(loopinfo.loop.anchor1.anchor_id, loopinfo.loop.anchor2.anchor_id, weight=loopinfo.pet)
    print('\t'*func_depth+f"<TIME> build_graph cost {time.time()-time1} seconds")
    return graph

def bfs(graph, source1, source2, r, func_depth=0):
    """
    计算每个边的距离为1的邻居的weight，距离为2的邻居的weight，直到距离为r的邻居的weight。每个边的这r个距离数组构成一个字典。
    """
    # print('\t'*func_depth+"<FUNC> bfs")
    # time1 = time.time()
    neighbor_edge_dict = dict() # key radius, value list of edges
    for i in range(1, r + 1):
        neighbor_edge_dict[i] = list()
    visited = {source1, source2}
    queue = deque([(source1, 0), (source2, 0)])
    while queue:
        node, distance = queue.popleft()
        if distance < r:
            for neighbor in graph.neighbors(node):
                if neighbor not in visited:
                    if node < neighbor:
                        neighbor_edge_dict[distance+1].append((node, neighbor))
                    else:
                        neighbor_edge_dict[distance+1].append((neighbor, node))
                    visited.add(neighbor)
                    queue.append((neighbor, distance + 1))
    # print('\t'*func_depth+f"<TIME> bfs cost {time.time()-time1} seconds")
    return neighbor_edge_dict

def get_gaussian_kernel(radius):
    truncate = 4.0
    # sigma = 2.0  # 你想要的标准差

    # # 计算所需窗口大小，确保包含大部分高斯分布能量
    # window_size = int(truncate * sigma + 0.5)
    # print(window_size)
    # r = 4
    gaussian_filter_dict = dict()
    for r in range(1, radius+2):
        sigma = r/truncate
        # 创建一维均匀网格，并计算对应的高斯值
        x = np.linspace(-r, r, 2*r)
        # get half of x
        x = x[r:]
        # print(x)
        gaussian_kernel_1d = np.exp(-(x ** 2 / (2 * sigma ** 2)))
        # print("gaussian_kernel_1d", gaussian_kernel_1d, np.sum(gaussian_kernel_1d))
        # print("norm: ", gaussian_kernel_1d/np.sum(gaussian_kernel_1d))
        gaussian_filter_dict[r-1] = gaussian_kernel_1d/np.sum(gaussian_kernel_1d)
    return gaussian_filter_dict

def transfer_edge_values(loopinfos, graph: nx.Graph, r: int, func_depth=0):
    print('\t'*func_depth+"<FUNC> transfer_edge_values")
    time1 = time.time()
    # next_func_depth = func_depth + 1
    edge_derivatives = []
    # 直接用graph.edges()，序号会乱，和原数据的interaction不对应
    # for u, v, d in graph.edges(data=True):
    for loopinfo in loopinfos:
        u, v = loopinfo.loop.anchor1.anchor_id, loopinfo.loop.anchor2.anchor_id
        # * 找到两个anchors的邻居，然后合并作为边的邻居
        neighbor_edge_dict = bfs(graph, u, v, r) # 对于一个edge的邻居字典
        # * 计算一个边的不同半径r临域的平均值
        averages = [graph[u][v]['weight']]
        for i in range(1, r + 1):
            if neighbor_edge_dict[i]:
                weights = [graph[node_a][node_b]['weight'] for node_a, node_b in neighbor_edge_dict[i]]
                averages.append(sum(weights) / len(weights))
            else:
                averages.append(0)
        # * 求不同r之间的导数，并存储
        derivatives = [averages[i] - averages[i+1] for i in range(0, len(averages)-1)]
        edge_derivatives.append(derivatives)
    print('\t'*func_depth+f"<TIME> transfer_edge_values cost {time.time()-time1} seconds")
    # * 返回所有边的不同半径r临域的平均值们， 形状为[边数, r]
    return edge_derivatives

def transfer_edge_values_with_gaussian(loopinfos, graph: nx.Graph, r: int, func_depth=0):
    print('\t'*func_depth+"<FUNC> transfer_edge_values")
    time1 = time.time()
    gaussian_kernel_dict = get_gaussian_kernel(r)
    # next_func_depth = func_depth + 1
    edge_derivatives = list()
    has_neighbor = list()
    # 直接用graph.edges()，序号会乱，和原数据的interaction不对应
    # for u, v, d in graph.edges(data=True):
    for loopinfo in loopinfos:
        u, v = loopinfo.loop.anchor1.anchor_id, loopinfo.loop.anchor2.anchor_id
        # * 找到两个anchors的邻居，然后合并作为边的邻居
        neighbor_edge_dict = bfs(graph, u, v, r) # 对于一个edge的邻居字典
        # * 计算一个边的不同半径r临域的平均值
        averages = [graph[u][v]['weight']]
        for i in range(1, r + 1):
            if neighbor_edge_dict[i]:
                weights = [graph[node_a][node_b]['weight'] for node_a, node_b in neighbor_edge_dict[i]]
                averages.append(sum(weights) / len(weights))
            else:
                averages.append(0)
        averages = np.array(averages)
        # * 判断是否有邻居, 有则True, 无则False, 用于后续判断是否使用其他方法的P值和FDR
        if np.sum(averages[1:]) == 0:
            has_neighbor.append(False)
        else:
            has_neighbor.append(True)
        # * 求不同r之间的导数，并存储
        gaussian_pets = np.zeros(r+1)
        for i in range(0, r+1):
            # 点乘
            gaussian_pets[i] = np.sum(np.multiply(averages[:i+1], gaussian_kernel_dict[i]))

        derivatives = [gaussian_pets[i] - gaussian_pets[i+1] for i in range(0, len(gaussian_pets)-1)]
        edge_derivatives.append(derivatives)
    print('\t'*func_depth+f"<TIME> transfer_edge_values cost {time.time()-time1} seconds")
    # * 返回所有边的不同半径r临域的平均值们， 形状为[边数, r]
    return edge_derivatives, has_neighbor

def run_dci_by_neighbors(loopinfos, groups, r, func_depth=0):
    """
    loopinfos: list of list of LoopInfo [samples_num, loops_num]
    """
    print('\t'*func_depth+"<FUNC> run_dci_by_neighbors")
    time1 = time.time()
    next_func_depth = func_depth + 1
    # * 计算每个sample 每个环 临域的 不同r的 导数
    derivatives_for_samples = list()
    # TODO delete the "has_neighbor" variables
    has_neighbor_matrix = list()
    for loopinfos_one_sample in loopinfos:
        # loopinfos_one_sample: loopinfos in one sample
        graph = build_graph(loopinfos_one_sample, next_func_depth)
        # edge_derivatives = transfer_edge_values(loopinfos_one_sample, graph, r, next_func_depth)
        edge_derivatives, has_neighbor = transfer_edge_values_with_gaussian(loopinfos_one_sample, graph, r, next_func_depth)
        derivatives_for_samples.append(edge_derivatives)
        has_neighbor_matrix.append(has_neighbor)
    # derivatives_for_samples shape: [samples_num, loops_num, r]
    derivatives_for_samples = np.array(derivatives_for_samples)
    # has_neighbor_matrix shape: [samples_num, loops_num]
    has_neighbor_matrix = np.array(has_neighbor_matrix)
    # 如果至少有一个sample有邻居，就是有邻居
    has_neighbor = np.sum(has_neighbor_matrix, axis=0) > 0
    # * count the index of sample in groups
    group_indices = defaultdict(list)
    for i, group in enumerate(groups):
        group_indices[group].append(i)
    group_set = list(set(groups))
    # * 计算各个半径r的导数的均值和标准差
    # // means = list()
    # // stds = list()
    # // for i in range(r):
    # //     means.append(np.mean(derivatives_for_samples[:, :, i]))
    # //     stds.append(np.std(derivatives_for_samples[:, :, i]))
    # // MINSTD = 1e-8
    # // stds = np.maximum(stds, MINSTD)
    means = list()
    stds = list()
    for rr in range(r):
        derivative_diffs = list()
        for i in group_indices[group_set[0]]:
            for j in group_indices[group_set[1]]:
                derivative_diffs.append(derivatives_for_samples[i, :, rr]-derivatives_for_samples[j, :, rr])
        derivative_diffs = np.array(derivative_diffs)
        means.append(np.mean(derivative_diffs))
        stds.append(np.std(derivative_diffs))
    MINSTD = 1e-8
    stds = np.maximum(stds, MINSTD)
    # * 计算不同group的sample间的对应的半径r之间的正态分布sf
    # 比如k562和mcf7各有n,m个rep，则每个loop要计算n*m*r个sf
    interaction_pvalues = list() # shape: [loops_num] every interaction has a final pvalue
    # i,j is index of sample; k is index of loop/interaction; rr is index of radius r
    for k in range(derivatives_for_samples.shape[1]):
        pvalues_for_an_interaction = list()
        for i in group_indices[group_set[0]]:
            for j in group_indices[group_set[1]]:
                for rr in range(r):
                    sf = norm.sf(derivatives_for_samples[i, k, rr]-derivatives_for_samples[j, k, rr], loc=means[rr], scale=stds[rr])
                    pvalues_for_an_interaction.append(sf)
        pvalues_for_an_interaction = np.array(pvalues_for_an_interaction)
        np.nan_to_num(pvalues_for_an_interaction, copy=False, posinf=1, neginf=1, nan=1)
        # 由于算的是差异，所以正态分布的左边那半边要反转过来，然后总体数值*2.
        pvalues_for_an_interaction[pvalues_for_an_interaction > 0.5] = 1 - pvalues_for_an_interaction[pvalues_for_an_interaction > 0.5]
        pvalues_for_an_interaction *= 2
        # interaction_pvalues.append(min(pvalues_for_an_interaction))
        combine_method = 'fisher'
        _, merged_p_value = combine_pvalues(pvalues_for_an_interaction, method=combine_method)
        interaction_pvalues.append(merged_p_value)

    _, fdrs = fdrcorrection(interaction_pvalues, alpha=0.05, method='indep', is_sorted=False)
    print('\t'*func_depth+f"<TIME> run_dci_by_neighbors cost {time.time()-time1} seconds")
    return interaction_pvalues, fdrs, has_neighbor

# @DeprecationWarning
# def test():
#     import sys
#     sys.path.insert(0, 'D:\\workspace\\dci/src/utils')
#     # print(sys.path)
#     from loop_object import LoopInfo, Loop, Anchor

#     anchor1 = Anchor(1, 2)
#     anchor2 = Anchor(3, 4)
#     anchor3 = Anchor(5, 6)
#     anchor4 = Anchor(7, 8)

#     loop0 = Loop("chr1", anchor1, anchor2)
#     loop1 = Loop("chr1", anchor2, anchor3)
#     loop2 = Loop("chr1", anchor3, anchor4)

#     loopinfo0 = LoopInfo(loop0, 0, 1)
#     loopinfo1 = LoopInfo(loop1, 1, 1)
#     loopinfo2 = LoopInfo(loop2, 2, 1)
#     loopinfos1 = [loopinfo1, loopinfo0, loopinfo2]
#     # graph = build_graph(loopinfos)
#     # print(graph.edges(data=True))
#     # edge_derivatives = transfer_edge_values(graph, r=3)
#     # print(edge_derivatives)

#     loopinfo0 = LoopInfo(loop0, 2, 1)
#     loopinfo1 = LoopInfo(loop1, 3, 1)
#     loopinfo2 = LoopInfo(loop2, 4, 1)
#     loopinfos2 = [loopinfo1, loopinfo0, loopinfo2]

#     # pvalues, fdrs = run_dci_by_neighbors([loopinfos1, loopinfos2], ["k562", "mcf7"], r=3)
#     # print(pvalues, fdrs)

#     # print(gaussian_filter1d([1.0, 1.0, 1.0, 1.0, 1.0], 2))
#     # print(gaussian_filter1d([1.0, 1.0, 1.0, 1.0, 1.0, 1.0], 2))
#     # print(gaussian_filter1d([1.0, 2.0, 3.0, 4.0, 5.0], 1))


