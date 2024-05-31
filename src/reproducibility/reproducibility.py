# from utils.loop_object import LoopInfo, Loop, Anchor
# import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import time, random, math, os

def sort_loop_infos_by_chr_start_end(loop_infos, func_depth=0):
    time1 = time.time()
    print('\t'*func_depth+"<FUNC> sort_loop_infos_by_chr_start_end")
    def custom_sort(x):
        # sort by int < X < Y
        x = x[3:]
        try:
            return int(x), 0
        except ValueError:
            letter_order = {'X': 1, 'Y': 2}  # Add more if needed
            return float('inf'), letter_order.get(x, 0)
    print('\t'*func_depth+f"<TIME> sort_loop_infos_by_chr_start_end cost {time.time()-time1} seconds")
    return sorted(loop_infos, key=lambda x: (custom_sort(x.loop.chr), x.loop.anchor1.start, x.loop.anchor1.end, x.loop.anchor2.start, x.loop.anchor2.end))

def get_nodes_and_edges_of_a_sample(loop_infos, func_depth=0):
    time1 = time.time()
    next_func_depth = func_depth + 1
    print('\t'*func_depth+"<FUNC> get_nodes_and_edges_of_a_sample")
    # 1. LoopInfo is node;
    # 2. LoopInfo1 and LoopInfo2 are connected if they have overlapped regions;
    # 3. edge weight is the absolute pet difference of the LoopInfo1 and LoopInfo2;
    nodes = sort_loop_infos_by_chr_start_end(loop_infos, func_depth=next_func_depth)
    # print(nodes[:10])
    # print(nodes[-10:])
    edges = list()
    nodes_num = len(nodes)
    for i in range(nodes_num):
        for j in range(i+1, nodes_num):
            if nodes[i].loop.chr != nodes[j].loop.chr:
                break
            if nodes[i].overlap(nodes[j]):
                edges.append((nodes[i], nodes[j], abs(nodes[i].pet - nodes[j].pet)))
    print('\t'*func_depth+f"<TIME> get_nodes_and_edges_of_a_sample cost {time.time()-time1} seconds")
    return nodes, edges

@DeprecationWarning
def draw_circlr_plot(nodes, edges, output_path, visible_proportion=0.1, func_depth=0):
    time1 = time.time()
    print('\t'*func_depth+"<FUNC> draw_circlr_plot")
    nodes_num = len(nodes)
    # random remove visible_proportion edges
    edges = [edge for edge in edges if random.random()<visible_proportion]
    FONTSIZE = 10
    plt.figure(dpi=300, figsize=(12, 12))
    plt.rc("font", family="Times New Roman")
    params = {"axes.titlesize": FONTSIZE,
              "legend.fontsize": FONTSIZE,
              "axes.labelsize": FONTSIZE,
              "xtick.labelsize": FONTSIZE,
              "ytick.labelsize": FONTSIZE,
              "figure.titlesize": FONTSIZE,
              "font.size": FONTSIZE}
    plt.rcParams.update(params)
    color_list = ['#000000', '#7F7F7F', '#880015', '#ED1C24', '#F08784', '#FF7F27', '#FFF200', '#22B14C', '#377D22', '#00A2E8', '#3F48CC', '#EE8AF8', '#EA3680', '#A349A4',
                  '#C3C3C3', '#B97A57', '#FFAEC9', '#FFC90E', '#EFE4B0', '#75F94D', '#B5E61D', '#99D9EA', '#73FBFD', '#7092BE', '#7E84F7', '#C8BFE7']
    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_weighted_edges_from(edges)
    UG = G.to_undirected()
    # 定义圆形布局
    pos = {nodes[i]: (-np.cos(np.pi/2 + 2*np.pi*i/nodes_num), np.sin(np.pi/2 + 2*np.pi*i/nodes_num)) for i in range(nodes_num)}
    options = {
        "node_color": "black",
        "node_size": 1,
        "edge_color": "gray",
        "linewidths": 0,
        "width": 0.1,
        "alpha": 0.1,
        "arrows": True,
        "connectionstyle": 'arc3,rad=1.0'
    }
    nx.draw(UG, pos, **options)
    # plt.tight_layout()
    plt.savefig(f"{output_path}/figures/reproducibility.png", format="png", bbox_inches="tight")
    # plt.show()
    print('\t'*func_depth+f"<TIME> draw_circlr_plot cost {time.time()-time1} seconds")

@DeprecationWarning
def compute_reproducibility_by_graph(data, output_path, func_depth=0):
    time1 = time.time()
    next_func_depth = func_depth + 1
    print('\t'*func_depth+"<FUNC> compute_reproducibility_by_graph")
    # * convert each data[i] to a graph
    # 1. LoopInfo is node;
    # 2. LoopInfo1 and LoopInfo2 are connected if they have overlapped regions;
    # 3. edge weight is the absolute pet difference of the LoopInfo1 and LoopInfo2;
    # for loopinfos in data:
    # nodes and edges are already sorted by (chr, start1, end1, start2, end2)
    # nodes, edges = get_nodes_and_edges_of_a_sample(data[1], func_depth=next_func_depth)
    # print("nodes num:", len(nodes))
    # print("edges num:", len(edges))
    # connected components
    # connected_components = nx.connected_components(UG)
    # print("connected_components:", len(list(connected_components)))
    # draw_circlr_plot(nodes, edges, output_path, visible_proportion=0.1, func_depth=next_func_depth)

    print('\t'*func_depth+f"<TIME> compute_reproducibility_by_graph cost {time.time()-time1} seconds")

def compute_reproducibility(data, samples, block_size, stride, c, output_path, func_depth=0):
    time1 = time.time()
    next_func_depth = func_depth + 1
    print('\t'*func_depth+"<FUNC> compute_reproducibility")
    samples_num = len(data)
    # * sort data by chr, start1, end1, start2, end2
    for i in range(samples_num):
        data[i] = sort_loop_infos_by_chr_start_end(data[i], func_depth=next_func_depth)
    # * compute each chromosome range among all samples
    chr_ranges = dict()
    start = math.inf
    end = -math.inf
    for i in range(samples_num):
        current_chr = data[i][0].loop.chr
        for loop_info in data[i]:
            if loop_info.loop.chr != current_chr:
                if current_chr not in chr_ranges:
                    chr_ranges[current_chr] = [start, end]
                else:
                    chr_ranges[current_chr][0] = min(chr_ranges[current_chr][0], start)
                    chr_ranges[current_chr][1] = max(chr_ranges[current_chr][1], end)
                current_chr = loop_info.loop.chr
                start = math.inf
                end = -math.inf
            start = min(start, loop_info.loop.anchor1.start)
            end = max(end, loop_info.loop.anchor2.end)
        if current_chr not in chr_ranges:
            chr_ranges[current_chr] = [start, end]
        else:
            chr_ranges[current_chr][0] = min(chr_ranges[current_chr][0], start)
            chr_ranges[current_chr][1] = max(chr_ranges[current_chr][1], end)
    # print(chr_ranges)
    # * get every block positions around all chrs
    block_positions = dict()
    for chr in chr_ranges:
        block_positions[chr] = list()
        window_start = chr_ranges[chr][0]
        window_end = window_start + block_size
        while window_start < chr_ranges[chr][1]:
            block_positions[chr].append((window_start, window_end))
            window_start += stride
            window_end = window_start + block_size

    # * scan blocks of each sample, & compute mean pet of each block
    samples_mean_pets = list() # samples_mean_pets[i][j] is the mean pet of the jth block of the ith sample, shape: (samples_num, blocks_num)
    for i in range(samples_num):
        loops_num = len(data[i]) # loops number of a sample
        mean_pets = list() # mean PET of each block for a sample
        pets_in_window = list() # PETs in a window/block, used for computing mean pets
        loop_index = 0
        current_chr = data[i][0].loop.chr
        current_block_positions = block_positions[current_chr] # a list containing current chr's blocks start & end
        window_start, window_end = current_block_positions[0] # get the first chr's start
        next_loop_index = None
        current_block_index = 0
        while loop_index < loops_num:
            # * jump to next chr, we need complement 上一个 chr's pets, and get new window/block for next chr.
            if data[i][loop_index].loop.chr != current_chr:
                if len(pets_in_window) > 0:
                    mean_pets.append(np.mean(pets_in_window))
                else:
                    mean_pets.append(0.)
                for _ in range(current_block_index+1, len(current_block_positions)):
                    mean_pets.append(0.)
                pets_in_window = list()
                current_chr = data[i][loop_index].loop.chr
                current_block_positions = block_positions[current_chr]
                window_start, window_end = current_block_positions[0] # get the current chr's start
                next_loop_index = None
                current_block_index = 0
            # * if current window/block still has remained loops, add next loop's pet.
            if data[i][loop_index].loop.anchor1.start <= window_end:
                pets_in_window.append(data[i][loop_index].pet)
                # 找到下一个block开始的loop的索引，找到loop的右边界与下个block的start重叠 的loop
                if current_block_index < len(current_block_positions)-1 and next_loop_index is None and data[i][loop_index].loop.anchor2.end >= current_block_positions[current_block_index+1][0]:
                    next_loop_index = loop_index
                loop_index += 1
            # * if current window/block's loops are already searched, jump into next window/block.
            else:
                if len(pets_in_window) > 0:
                    mean_pets.append(np.mean(pets_in_window))
                else:
                    mean_pets.append(0.)
                pets_in_window = list()
                if next_loop_index != None:
                    loop_index = next_loop_index
                window_start, window_end = current_block_positions[current_block_index+1]
                next_loop_index = None
                current_block_index += 1
        # * all chrs are searched, we need complement pets for 最后的 chr.
        if len(pets_in_window) > 0:
            mean_pets.append(np.mean(pets_in_window))
        else:
            mean_pets.append(0.)
        for _ in range(current_block_index+1, len(current_block_positions)):
            mean_pets.append(0.)
        samples_mean_pets.append(mean_pets)

    # for i in range(samples_num):
    #     print(len(samples_mean_pets[i]), samples_mean_pets[i][:10], samples_mean_pets[i][-10:])

    # * compute self-similarity indicator scores
    indicator_scores = [list() for _ in range(samples_num)]
    for i in range(samples_num):
        for j in range(len(samples_mean_pets[i])):
            for k in range(j+1, len(samples_mean_pets[i])):
                if samples_mean_pets[i][j] >= samples_mean_pets[i][k]:
                    indicator_scores[i].append(1)
                else:
                    indicator_scores[i].append(0)

    # for i in range(samples_num):
    #     print(len(indicator_scores[i]), indicator_scores[i][:10], indicator_scores[i][-10:])

    # * compute reproducibility
    def compute_reproducibility_of_two_samples(indicator_scores1, indicator_scores2):
        # s = 0
        # for i in range(len(indicator_scores1)):
        #     s += (indicator_scores1[i] - indicator_scores2[i])
        # s = s**2
        # s = math.sqrt(s)
        indicator_scores1 = np.array(indicator_scores1)
        indicator_scores2 = np.array(indicator_scores2)
        s = np.power(np.sum(np.power(indicator_scores1-indicator_scores2, 2)), 0.5)
        return math.exp(-c*s), s

    reproducibility_result = list()
    for i in range(samples_num):
        for j in range(i+1, samples_num):
            reproducibility, s = compute_reproducibility_of_two_samples(indicator_scores[i], indicator_scores[j])
            print(f"Reproducibility of {samples[i]} and {samples[j]} is {reproducibility}, {s}")
            reproducibility_result.append((samples[i], samples[j], reproducibility, s))

    print('\t'*func_depth+f"<TIME> compute_reproducibility cost {time.time()-time1} seconds")
    return reproducibility_result


def compute_reproducibility_chr(data, samples, block_size, stride, c, output_path, func_depth=0):
    time1 = time.time()
    next_func_depth = func_depth + 1
    print('\t'*func_depth+"<FUNC> compute_reproducibility")
    samples_num = len(data)
    # * sort data by chr, start1, end1, start2, end2
    for i in range(samples_num):
        data[i] = sort_loop_infos_by_chr_start_end(data[i], func_depth=next_func_depth)
    # * compute each chromosome range among all samples
    chr_ranges = dict()
    start = math.inf
    end = -math.inf
    for i in range(samples_num):
        current_chr = data[i][0].loop.chr
        for loop_info in data[i]:
            if loop_info.loop.chr != current_chr:
                if current_chr not in chr_ranges:
                    chr_ranges[current_chr] = [start, end]
                else:
                    chr_ranges[current_chr][0] = min(chr_ranges[current_chr][0], start)
                    chr_ranges[current_chr][1] = max(chr_ranges[current_chr][1], end)
                current_chr = loop_info.loop.chr
                start = math.inf
                end = -math.inf
            start = min(start, loop_info.loop.anchor1.start)
            end = max(end, loop_info.loop.anchor2.end)
        if current_chr not in chr_ranges:
            chr_ranges[current_chr] = [start, end]
        else:
            chr_ranges[current_chr][0] = min(chr_ranges[current_chr][0], start)
            chr_ranges[current_chr][1] = max(chr_ranges[current_chr][1], end)
    # print(chr_ranges)
    # * get every block positions around all chrs
    block_positions = dict()
    for chr in chr_ranges:
        block_positions[chr] = list()
        window_start = chr_ranges[chr][0]
        window_end = window_start + block_size
        while window_start < chr_ranges[chr][1]:
            block_positions[chr].append((window_start, window_end))
            window_start += stride
            window_end = window_start + block_size

    # * scan blocks of each sample, & compute mean pet of each block
    # samples_mean_pets[i]["chr"][j] is the mean pet of the jth block of "chr" in the ith sample, shape: (samples_num, chr_num, blocks_num)
    samples_mean_pets = list()
    for i in range(samples_num):
        loops_num = len(data[i]) # loops number of a sample
        mean_pets = {chr: list() for chr in chr_ranges} # mean PET of each block for a sample
        pets_in_window = list() # PETs in a window/block, used for computing mean pets
        loop_index = 0
        current_chr = data[i][0].loop.chr
        current_block_positions = block_positions[current_chr] # a list containing current chr's blocks start & end
        window_start, window_end = current_block_positions[0] # get the first chr's start
        next_loop_index = None
        current_block_index = 0
        while loop_index < loops_num:
            # * jump to next chr, we need complement 上一个 chr's pets, and get new window/block for next chr.
            if data[i][loop_index].loop.chr != current_chr:
                if len(pets_in_window) > 0:
                    mean_pets[current_chr].append(np.mean(pets_in_window))
                else:
                    mean_pets[current_chr].append(0.)
                for _ in range(current_block_index+1, len(current_block_positions)):
                    mean_pets[current_chr].append(0.)
                pets_in_window = list()
                current_chr = data[i][loop_index].loop.chr
                current_block_positions = block_positions[current_chr]
                window_start, window_end = current_block_positions[0] # get the current chr's start
                next_loop_index = None
                current_block_index = 0
            # * if current window/block still has remained loops, add next loop's pet.
            if data[i][loop_index].loop.anchor1.start <= window_end:
                pets_in_window.append(data[i][loop_index].pet)
                # 找到下一个block开始的loop的索引，找到loop的右边界与下个block的start重叠 的loop
                if current_block_index < len(current_block_positions)-1 and next_loop_index is None and data[i][loop_index].loop.anchor2.end >= current_block_positions[current_block_index+1][0]:
                    next_loop_index = loop_index
                loop_index += 1
            # * if current window/block's loops are already searched, jump into next window/block.
            else:
                if len(pets_in_window) > 0:
                    mean_pets[current_chr].append(np.mean(pets_in_window))
                else:
                    mean_pets[current_chr].append(0.)
                pets_in_window = list()
                if next_loop_index != None:
                    loop_index = next_loop_index
                window_start, window_end = current_block_positions[current_block_index+1]
                next_loop_index = None
                current_block_index += 1
        # * all chrs are searched, we need complement pets for 最后的 chr.
        if len(pets_in_window) > 0:
            mean_pets[current_chr].append(np.mean(pets_in_window))
        else:
            mean_pets[current_chr].append(0.)
        for _ in range(current_block_index+1, len(current_block_positions)):
            mean_pets[current_chr].append(0.)
        samples_mean_pets.append(mean_pets)

    # * compute self-similarity indicator scores
    indicator_scores = [list() for _ in range(samples_num)]
    for i in range(samples_num):
        for chr in chr_ranges:
            for j in range(len(samples_mean_pets[i][chr])):
                for k in range(j+1, len(samples_mean_pets[i][chr])):
                    if samples_mean_pets[i][chr][j] >= samples_mean_pets[i][chr][k]:
                        indicator_scores[i].append(1)
                    else:
                        indicator_scores[i].append(0)

    # * compute reproducibility
    def compute_reproducibility_of_two_samples(indicator_scores1, indicator_scores2):
        # s = 0
        # for i in range(len(indicator_scores1)):
        #     s += (indicator_scores1[i] - indicator_scores2[i])
        # s = s**2
        # s = math.sqrt(s)
        indicator_scores1 = np.array(indicator_scores1)
        indicator_scores2 = np.array(indicator_scores2)
        s = np.power(np.sum(np.power(indicator_scores1-indicator_scores2, 2)), 0.5)
        return math.exp(-c*s), s

    reproducibility_path = os.path.join(output_path, "reproducibility", "reproducibility_result.txt")
    reproducibility_result = list()
    with open(reproducibility_path, "w") as f:
        for i in range(samples_num):
            for j in range(i+1, samples_num):
                reproducibility, s = compute_reproducibility_of_two_samples(indicator_scores[i], indicator_scores[j])
                print(f"Reproducibility of {samples[i]} and {samples[j]} is {reproducibility}, {s}")
                f.write(f"Reproducibility of {samples[i]} and {samples[j]} is {reproducibility}\n")
                reproducibility_result.append((samples[i], samples[j], reproducibility, s))


    print('\t'*func_depth+f"<TIME> compute_reproducibility cost {time.time()-time1} seconds")
    return reproducibility_result


# if __name__ == "__main__":
#     anchor1 = Anchor(1, 2)
#     anchor2 = Anchor(3, 4)
#     anchor3 = Anchor(5, 6)
#     anchor4 = Anchor(7, 8)
#     anchor5 = Anchor(9, 10)
#     anchor6 = Anchor(11, 12)
#     anchor7 = Anchor(13, 14)
#     anchor8 = Anchor(15, 16)

#     anchor9 = Anchor(117, 118)
#     anchor10 = Anchor(119, 120)
#     anchor11 = Anchor(121, 122)
#     anchor12 = Anchor(123, 124)
#     anchor13 = Anchor(125, 126)
#     anchor14 = Anchor(127, 128)
#     anchor15 = Anchor(129, 130)
#     anchor16 = Anchor(131, 132)

#     anchor17 = Anchor(233, 234)
#     anchor18 = Anchor(235, 236)

#     loop1 = Loop("chr1", anchor1, anchor2)
#     loop2 = Loop("chr1", anchor3, anchor4)
#     loop3 = Loop("chr1", anchor5, anchor6)
#     loop4 = Loop("chr1", anchor7, anchor8)

#     loop5 = Loop("chr1", anchor9, anchor10)
#     loop6 = Loop("chr1", anchor11, anchor12)
#     loop7 = Loop("chr1", anchor13, anchor14)
#     loop8 = Loop("chr1", anchor15, anchor16)

#     loop9 = Loop("chr1", anchor17, anchor18)

#     loop_info1 = LoopInfo(loop1, 10, 0.01)
#     loop_info2 = LoopInfo(loop2, 8, 0.21)
#     loop_info3 = LoopInfo(loop3, 14, 0.31)
#     loop_info4 = LoopInfo(loop4, 21, 0.31)

#     loop_info5 = LoopInfo(loop5, 4, 0.1)
#     loop_info6 = LoopInfo(loop6, 2, 0.1)
#     loop_info7 = LoopInfo(loop7, 23, 0.3)
#     loop_info8 = LoopInfo(loop8, 12, 0.11)

#     loop_info9 = LoopInfo(loop9, 233, 0.131)

#     nodes = [loop_info1, loop_info2, loop_info3, loop_info4, loop_info5, loop_info6, loop_info7, loop_info8, loop_info9]
#     # edges = [(loop1, loop2, 1), (loop2, loop3, 2), (loop3, loop4, 3),
#     #          (loop5, loop6, 4), (loop6, loop7, 5), (loop7, loop8, 6)]
#     edges = [(loop_info1, loop_info2, 1), (loop_info2, loop_info3, 2), (loop_info3, loop_info4, 3),
#             (loop_info5, loop_info6, 4), (loop_info6, loop_info7, 5), (loop_info7, loop_info8, 6)]

#     G = nx.Graph()
#     G.add_nodes_from(nodes)
#     G.add_weighted_edges_from(edges)

#     # connected components
#     connected_components = nx.connected_components(G)
#     for component in connected_components:
#         print("Connected Component:", component)
#         # 根据节点集合创建一个子图
#         subgraph = G.subgraph(component)
        
#         # 获取子图中的所有边
#         edges_in_subgraph = subgraph.edges()
#         print("Edges in Connected Component:", edges_in_subgraph)
        
#         # 获取子图中的所有边及其权重, edge[0]和edge[1]是两个节点，G[edge[0]][edge[1]]['weight']是两个节点之间edge的权重
#         edges_with_weights = [(edge[0], edge[1], G[edge[0]][edge[1]]['weight']) for edge in subgraph.edges()]
#         print("Edges with weights in Connected Component:", edges_with_weights)

#     print([i for i in range(10,1)])


#     FONTSIZE = 10
#     plt.figure(dpi=300, figsize=(12, 12))
#     plt.rc("font", family="Times New Roman")
#     params = {"axes.titlesize": FONTSIZE,
#               "legend.fontsize": FONTSIZE,
#               "axes.labelsize": FONTSIZE,
#               "xtick.labelsize": FONTSIZE,
#               "ytick.labelsize": FONTSIZE,
#               "figure.titlesize": FONTSIZE,
#               "font.size": FONTSIZE}
#     plt.rcParams.update(params)
#     connected_components = nx.connected_components(G)
#     nx.draw(G, node_size=1, width=0.1, alpha=0.5)
#     plt.tight_layout()
#     plt.show()