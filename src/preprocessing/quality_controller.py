import os, time
import numpy as np
from collections import defaultdict

from utils.loop_object import Anchor, Loop, LoopInfo
from preprocessing.preprocessing_utils import complement_zero_pet_loops_for_samples, complement_zero_pet_loops_from_LoopClass

# combine anchors which distance is less than parameter 'gap' (default: 500)
# for example, in followed file, column 1 and 2 are anchor1 start and end, column 4 and 5 are anchor2 start and end
# the column 6 and 7 are names of anchor1 and anchor2, seperately. the column 6,7 should be changed after combined
# the function loads all input files, and combine anchors less than 'gap' in all files, and change column 6,7 and save all files in new files
def combine_anchors(files, samples, gap, func_depth=0):
    print('\t'*func_depth+"<FUNC> combine_anchors")
    assert len(files) > 0
    assert len(files) == len(samples)
    assert gap >= 0
    # samples2files = {samples[i]: files[i] for i in range(len(files))}
    files2samples = {files[i]: samples[i] for i in range(len(files))}
    # anchors_list: [sampleName_anchorName, chr, start, end]
    anchors_list = list()
    # * load anchors
    for file in files:
        data = open(file, 'r')
        lines = data.readlines()
        data.close()
        if "start" in lines[0]:
            lines = lines[1:]
        # save already recorded anchors
        already_saved_anchors = set()
        for line in lines:
            line = line.strip().split()
            # peak_1 in MICC file
            sample_anchor_name = f"{files2samples[file]}_{line[6]}"
            if sample_anchor_name not in already_saved_anchors:
                chr = line[0] if "chr" not in line[0] else line[0][3:]
                anchor = [sample_anchor_name, chr, int(line[1]), int(line[2])]
                anchors_list.append(anchor)
                already_saved_anchors.add(sample_anchor_name)
            # peak_2 in MICC file
            sample_anchor_name = f"{files2samples[file]}_{line[7]}"
            if sample_anchor_name not in already_saved_anchors:
                chr = line[3] if "chr" not in line[3] else line[3][3:]
                anchor = [sample_anchor_name, chr, int(line[4]), int(line[5])]
                anchors_list.append(anchor)
                already_saved_anchors.add(sample_anchor_name)
    # print(f"[TIME] time1 - load anchors: {time.time() - start_time}")
    # * sort anchors by chromosome, start and end
    def custom_sort(x):
        # sort by int < X < Y
        letter_order = {'X': 1, 'Y': 2}  # Add more if needed
        try:
            return int(x), 0
        except ValueError:
            return float('inf'), letter_order.get(x, 0)
    anchors_list = sorted(anchors_list, key=lambda x: (custom_sort(x[1]), x[2], x[3]))
    # print(len(anchors_list))
    # print(anchors_list[36327:])
    # print(f"[TIME] time2 - sort anchors: {time.time() - start_time}")
    # chr_list = list()
    # for anchor in anchors_list:
    #     if anchor[1] not in chr_list:
    #         chr_list.append(anchor[1])
    # print(chr_list)
    # * combine anchors dict
    def judge_overlap(start1, end1, start2, end2):
        if start1 <= start2:
            return start2 - end1 <= gap
        else:
            return start1 - end2 <= gap
    def add_to_list(current_list, input):
        if isinstance(current_list, str):
            current_list = [current_list]
        if isinstance(input, list):
            current_list.extend(input)
        elif isinstance(input, str):
            current_list.append(input)
        return current_list

    combined_anchors = []
    current_anchor = anchors_list[0]
    for anchor in anchors_list[1:]:
        # anchor[1] is chr
        if anchor[1] == current_anchor[1] and judge_overlap(current_anchor[2], current_anchor[3], anchor[2], anchor[3]):
            current_anchor[2] = min(current_anchor[2], anchor[2])
            current_anchor[3] = max(current_anchor[3], anchor[3])
            # combine first parameter peakname into list
            current_anchor[0] = add_to_list(current_anchor[0], anchor[0])
        else:
            combined_anchors.append(current_anchor)
            current_anchor = anchor
    # the combined_anchors is already ordered, because combine dont change the order
    combined_anchors.append(current_anchor)
    # print(len(combined_anchors))
    # * look how many anchors are combined
    # combined_statistics = collections.defaultdict(int)
    # for anchor in combined_anchors:
    #     if isinstance(anchor[0], list):
    #         combined_statistics[len(anchor[0])] += 1
    #     else:
    #         combined_statistics[1] += 1
    # print(combined_statistics)
    # * verify the order combined anchors
    # for i in range(len(combined_anchors)-1):
    #     if combined_anchors[i][1] == combined_anchors[i+1][1]:
    #         if combined_anchors[i][3] > combined_anchors[i+1][2]:
    #             print(combined_anchors[i], combined_anchors[i+1])
    # print(f"[TIME] time3 - combine anchors dict: {time.time() - start_time}")

    # * build new map from combined peaks to new anchors
    combined_peaks2new_anchors = dict()
    for i, anchor in enumerate(combined_anchors):
        if isinstance(anchor[0], list):
            for peak in anchor[0]:
                # map from combined peaks to Anchor(new_anchor_name, new_start, new_end)
                combined_peaks2new_anchors[peak] = Anchor(anchor[2], anchor[3])
        else:
            combined_peaks2new_anchors[anchor[0]] = Anchor(anchor[2], anchor[3])
    # print(len(combined_peaks2new_anchors))
    print('\t'*func_depth+"[INFO] old peaks num:", len(set(combined_peaks2new_anchors.keys())), ", combined anchors num:", len(set(combined_peaks2new_anchors.values())))
    # print(f"[TIME] time4 - build new map: {time.time() - start_time}")
    # * 返回的anchor_id已经根据start和end排好序，并且是跨samples的全局的排序。
    # * return: keys:   sample_anchor_name, corresponding to column 6,7 in MICC file,
    # *         values: [new_anchor_name, new_start, new_end]
    return combined_peaks2new_anchors

# [QC]1. remove inter-chromosome loops [default: True]
# [QC]2. remove self-ligated loops [default: True]
# [QC]3. remove loops with FDR > 0.01 [default: 0.01]
# [QC]4. remove loops in chr filter [default: ["chrM", "chrX"]]
# [QC]5. remove loops in blacklist [default: True]
# 1. replace old anchors with new anchors in MICC file
# 2. merge loops with same anchors.
def qc_and_replace_anchors_in_MICC_file(file, sample, combined_peaks2new_anchors, blacklist,
                                 remove_loops_in_blacklist, remove_self_ligated_loops,
                                 loop_len_threshold, intra_only, chr_filter, output_path, func_depth=0):

    print('\t'*func_depth+f"<FUNC> replace anchors in {file}")
    data = open(file, 'r')
    loops = data.readlines()
    data.close()
    # remove header
    if "start" in loops[0]:
        loops = loops[1:]
    new_loops = dict()
    # * QC metrics record
    inter_num, inter_loop_list = 0, list()
    chr_filter_num, loops_in_chr_filter = 0, list()
    self_ligated_num, self_ligated_loops = 0, list()
    # fdr_gt_threshold_num, loops_gt_fdr_threshold = 0, list()
    loop_len_lt_threshold_num, loops_lt_loop_len_threshold = 0, list()
    blacklist_loop_num, loops_in_black_region = 0, list()
    for line in loops:
        fields = line.strip().split()
        # * intra only
        if intra_only and fields[0] != fields[3]:
            # print(f"[WARNING] inter-chromosome loop: {line}")
            inter_num += 1
            inter_loop_list.append(line)
            continue
        # * chr filter
        if fields[0][:5] in chr_filter or fields[3][:5] in chr_filter:
            # print(f"[WARNING] chr filter: {line}")
            chr_filter_num += 1
            loops_in_chr_filter.append(line)
            continue
        # * diffloop去除了mango输出参数FDR>0.01的环
        # if float(fields[12]) > fdr_threshold:
        #     # print(f"[WARNING] FDR > {fdr_threshold}: {line}")
        #     fdr_gt_threshold_num += 1
        #     loops_gt_fdr_threshold.append(line)
        #     continue
        # * remove loop length < 5000
        if abs(int(fields[5]) - int(fields[1])) < loop_len_threshold:
            loop_len_lt_threshold_num += 1
            loops_lt_loop_len_threshold.append(line)
            continue
        # * combine anchors may cause self-ligated loops, remove self-ligated loops
        # peak_1 in MICC file
        peak_name_1 = f"{sample}_{fields[6]}"
        # peak_2 in MICC file
        peak_name_2 = f"{sample}_{fields[7]}"
        if remove_self_ligated_loops and combined_peaks2new_anchors[peak_name_1] == combined_peaks2new_anchors[peak_name_2]:
            # print(f"[WARNING] self-ligated loop: {line}")
            self_ligated_num += 1
            self_ligated_loops.append(line)
            continue
        # * remove CNV regions (blacklist region)
        if remove_loops_in_blacklist:
            start1 = int(fields[1])
            end1 = int(fields[2])
            start2 = int(fields[4])
            end2 = int(fields[5])
            loop_left = min(start1, end1, start2, end2)
            loop_right = max(start1, end1, start2, end2)
            overlapping = False
            for entry_b in blacklist[fields[0]]:
                start_b, end_b = entry_b
                # 如果有重叠，则标记并退出循环
                if not (loop_right < start_b or loop_left > end_b):
                    overlapping = True
                    break
            if overlapping:
                blacklist_loop_num += 1
                loops_in_black_region.append(line)
                continue
        # * 1. replace old "peak" with new "anchor"; 2. merge loops with same anchors.
        try:
            # old exists, merge loops with same anchors.
            key = (combined_peaks2new_anchors[peak_name_1], combined_peaks2new_anchors[peak_name_2])
            loop_info = new_loops[key]
            loop_info.pet += int(fields[10]) # idr的合并规则：Default: 'sum' for signal/score/column indexes, 'min' for p/q-value.
            pvalue = 1-np.power(10, -float(fields[11]))
            loop_info.pvalue = min(loop_info.pvalue, pvalue)
            loop_info.fdr = min(loop_info.fdr, float(fields[12])) # ? 为什么要取最小值
        except KeyError:
            # create new, replace old peak with new anchor
            pvalue = 1-np.power(10, -float(fields[11]))
            new_loops[key] = LoopInfo(loop=Loop(fields[0],
                                                combined_peaks2new_anchors[peak_name_1],
                                                combined_peaks2new_anchors[peak_name_2]),
                                      pet=int(fields[10]),
                                      pvalue=pvalue,
                                      fdr=float(fields[12]))
    # * QC metrics record
    print('\t'*func_depth+f"[INFO] inter-chromosome loop num: {inter_num}")
    f = open(f"{output_path}/loops_filtered_in_QC/{sample}_inter_chromosome_loops.MICC", 'w')
    for line in inter_loop_list:
        f.write(line)
    f.close()
    print('\t'*func_depth+f"[INFO] chr filter loop num: {chr_filter_num}")
    f = open(f"{output_path}/loops_filtered_in_QC/{sample}_chr_filter_loops.MICC", 'w')
    for line in loops_in_chr_filter:
        f.write(line)
    f.close()
    print('\t'*func_depth+f"[INFO] self-ligated loop num: {self_ligated_num}")
    f = open(f"{output_path}/loops_filtered_in_QC/{sample}_self_ligated_loops.MICC", 'w')
    for line in self_ligated_loops:
        f.write(line)
    f.close()
    # print('\t'*func_depth+f"[INFO] FDR > {fdr_threshold} loop num: {fdr_gt_threshold_num}")
    # f = open(f"{output_path}/loops_filtered_in_QC/{sample}_FDR_gt_{fdr_threshold}_loops.MICC", 'w')
    # for line in loops_gt_fdr_threshold:
    #     f.write(line)
    # f.close()
    print('\t'*func_depth+f"[INFO] loop length > {loop_len_threshold} loop num: {loop_len_lt_threshold_num}")
    f = open(f"{output_path}/loops_filtered_in_QC/{sample}_loop_len_lt_{loop_len_threshold}_loops.MICC", 'w')
    for line in loops_lt_loop_len_threshold:
        f.write(line)
    f.close()
    print('\t'*func_depth+f"[INFO] loop in blacklist num: {blacklist_loop_num}")
    f = open(f"{output_path}/loops_filtered_in_QC/{sample}_loops_in_blacklist.MICC", 'w')
    for line in loops_in_black_region:
        f.write(line)
    f.close()
    # * sort for filter PET
    new_loops = list(new_loops.values())
    new_loops = sorted(new_loops, key=lambda x: x.loop.loop_id)
    # verify that each loop is unique
    loop_num = len(set([loop_info.loop.loop_id for loop_info in new_loops]))
    assert loop_num == len(new_loops)
    # * print new loops to file
    data = open(f"{output_path}/{sample}_replace_anchors_test.MICC", 'w')
    data.write("<[loop_id,chr,(anchor_id,start,end),(anchor_id,start,end)],PET,FDR>\n")
    for line in new_loops:
        data.write(str(line)+"\n")
    data.close()
    # print(f"[TIME] time5 - replace anchors in MICC file: {time.time() - start_time}")
    return new_loops


# combine anchors may cause self-ligated loops, remove them
# def remove_self_ligated_loops():
#     pass


# def remove_cnv_regions():
#     pass

def get_blacklist(blacklist_path, func_depth=0):
    print('\t'*func_depth+"<FUNC> get_blacklist")
    # the file is bed format, chr start end annotation
    data = open(blacklist_path, 'r')
    lines = data.readlines()
    data.close()
    # key: chr; value: [start, end]
    blacklist = defaultdict(list)
    for line in lines:
        line = line.strip().split()
        # chr = line[0]
        # start = int(line[1])
        # end = int(line[2])
        blacklist[line[0]].append([int(line[1]), int(line[2])])
    for k in blacklist:
        blacklist[k] = sorted(blacklist[k], key=lambda x: (x[0], x[1]))
    return blacklist


# diffloop去除了mango输出参数FDR>0.01的环
# def remove_by_FDR_greater_than(fdr_threshold=0.01):
#     pass


# 去除了一组replicates之间，在一个replicate出现PETs>=5在另一个replicate出现PETs==0的loops。
def remove_by_inconsistent_PETs_between_replicates(data, samples, groups, pet_threshold, output_path, func_depth=0):
    print('\t'*func_depth+"<FUNC> remove_by_inconsistent_PETs_between_replicates")

    def judge_top_loops(top_loops):
        j = 0
        while j < len(top_loops) and top_loops[j] is None:
            j += 1
        min_order_loops = [top_loops[j]]
        min_order_loop_indexes = [j]
        for i in range(j+1, len(top_loops)):
            if top_loops[i] is None:
                continue
            if top_loops[i].get_loop_id() == min_order_loops[0].get_loop_id():
                min_order_loops.append(top_loops[i])
                min_order_loop_indexes.append(i)
            elif top_loops[i].get_loop_id() < min_order_loops[0].get_loop_id():
                min_order_loops = [top_loops[i]]
                min_order_loop_indexes = [i]
        return min_order_loop_indexes

    new_data = [None for _ in groups]
    for group_name in list(set(groups)):
        print('\t'*func_depth+"[INFO] remove inconsistent PETs between replicates in group:", group_name)
        # * get current group indexes
        current_group_index = list()
        for i in range(len(groups)):
            if group_name == groups[i]:
                current_group_index.append(i)
        # print(current_group_index, group_name)
        # * get current group replicates
        replicates = [data[i] for i in current_group_index]
        print('\t'*func_depth+f"[INFO] {len(replicates)} replicates in group: {group_name}")
        # 每个replicate都是已经按照loop_name排序过的，
        # 下面要类似归并排序的方式，从每个replicate中取出队列前面的loop，找到同位置的loop，然后比较PETs数量
        # 然后看是否出现在一个replicate出现PETs>=5在另一个replicate出现PETs==0，
        # 出现的话，在本group中所有replicates中去掉这个loop。
        # every replicates are already sorted by loop_id, get the first loop in each replicate
        rep_num = len(replicates)
        new_replicates = [list() for i in range(rep_num)]
        is_remain = True # if all replicates are empty, then stop
        remove_inconsistent_PETs_num = 0
        while is_remain:
            top_loops = list()
            for i in range(rep_num):
                if len(replicates[i]) > 0:
                    top_loops.append(replicates[i][0])
                else:
                    top_loops.append(None)
            min_order_loop_indexes = judge_top_loops(top_loops)
            # 在多个replicates中，如果至少有一个replicates没有包含min_order_loop，就要检查是否有replicates中的该loop的PETs>=5，有的话去掉所有replicates中的这个loop
            # 如果所有replicates都包含min_order_loop，就直接通过。
            if len(min_order_loop_indexes) == rep_num:
                for i in min_order_loop_indexes:
                    new_replicates[i].append(replicates[i][0])
                    replicates[i] = replicates[i][1:]
            else:
                is_remove = False # if at least one replicate has PETs >= 5, then remove all replicates
                for i in min_order_loop_indexes:
                    if top_loops[i].pet >= pet_threshold:
                        is_remove = True
                        remove_inconsistent_PETs_num += 1
                        break
                if is_remove:
                    for i in min_order_loop_indexes:
                        replicates[i] = replicates[i][1:]
                else:
                    for i in min_order_loop_indexes:
                        new_replicates[i].append(replicates[i][0])
                        replicates[i] = replicates[i][1:]
            is_remain = False
            for i in range(rep_num):
                if len(replicates[i]) > 0:
                    is_remain = True
                    break
        print('\t'*func_depth+f"[INFO] remove inconsistent PETs num: {remove_inconsistent_PETs_num} in group: {group_name}")
        # * print new loops to file
        for i, rep in enumerate(new_replicates):
            f = open(f"{output_path}/{samples[current_group_index[i]]}_after_remove_inconsistent_PETs.MICC", 'w')
            f.write("<[loop_id,chr,(anchor_id,start,end),(anchor_id,start,end)],PET,FDR>\n")
            for line in rep:
                f.write(str(line)+"\n")
            f.close()
        for i in range(len(new_replicates)):
            new_data[current_group_index[i]] = new_replicates[i]
        # print(f"[TIME] time5 - replace anchors in MICC file: {time.time() - start_time}")
    return new_data

def remove_by_inconsistent_PETs_and_fdr_thres(data, samples, groups, pet_threshold, fdr_threshold, output_path, func_depth=0):
    print('\t'*func_depth+"<FUNC> remove_by_inconsistent_PETs_between_replicates")
    new_data = [list() for _ in groups]
    rep_num = range(len(groups))
    group_dict = defaultdict(list)
    for i, group_name in enumerate(groups):
        group_dict[group_name].append(i)
    number_of_inconsistent_interactions = 0
    fdr_gt_threshold_num = 0
    for i in range(len(data[0])): # i is interactions number
        # * filter inconsistent PETs
        # 对于不同rep的同一个loop，如果有一个rep的PETs==0，一个rep的PETs>=5，就去掉这个loop
        is_exceed_threshold = False
        for group_name in group_dict:
            pets_of_replicates = [data[j][i].pet for j in group_dict[group_name]]
            # if at least one replicate has 0 PET and least one has PETs >= 5, then remove all replicates
            if 0 in pets_of_replicates and max(pets_of_replicates) >= pet_threshold:
                is_exceed_threshold = True
                break
        if is_exceed_threshold:
            number_of_inconsistent_interactions += 1
            continue
        # * filter FDR > threshold
        fdrs_of_replicates = [data[j][i].fdr > fdr_threshold for j in rep_num]
        if all(fdrs_of_replicates):
            fdr_gt_threshold_num += 1
            continue
        for j in rep_num:
            new_data[j].append(data[j][i])
    print('\t'*func_depth+f"[INFO] remove inconsistent PETs num: {number_of_inconsistent_interactions}")
    print('\t'*func_depth+f"[INFO] remove FDR > {fdr_threshold} num: {fdr_gt_threshold_num}")
        # print(f"[TIME] time5 - replace anchors in MICC file: {time.time() - start_time}")
    return new_data

def quality_control(files, samples, groups, blacklist_path, output_path,
                    gap=500, remove_loops_in_blacklist=True,
                    remove_self_ligated_loops=True, fdr_threshold=0.01, loop_len_threshold=5000,
                    intra_only=True, chr_filter=["chrM", "chrX"], pet_threshold=5, func_depth=0):

    print('\t'*func_depth+"<FUNC> quality_control")
    next_func_depth = func_depth + 1
    # * combine overlapped anchors
    combined_peaks2new_anchors = combine_anchors(files, samples, gap, next_func_depth)

    # replace_anchors_in_MICC_file(GM12878_REP1, "GM12878_rep1", combined_peaks2new_anchors)

    blacklist = get_blacklist(blacklist_path, next_func_depth)

    # * [QC]1. remove inter-chromosome loops [default: True]
    # * [QC]2. remove self-ligated loops [default: True]
    # * [QC]3. remove loops in chr filter [default: ["chrM", "chrX"]]
    # * [QC]4. remove loops in blacklist [default: True]
    # * [Merge]1. replace old anchors with new anchors in MICC file
    # * [Merge]2. merge loops with same anchors.
    data = list()
    if not os.path.exists(f"{output_path}/loops_filtered_in_QC"):
        os.mkdir(f"{output_path}/loops_filtered_in_QC")
    for i in range(len(files)):
        data.append(qc_and_replace_anchors_in_MICC_file(files[i],
                                                        samples[i],
                                                        combined_peaks2new_anchors,
                                                        blacklist,
                                                        remove_loops_in_blacklist,
                                                        remove_self_ligated_loops,
                                                        loop_len_threshold, intra_only, chr_filter, output_path, next_func_depth))

    data = complement_zero_pet_loops_from_LoopClass(data, next_func_depth)

    # * [QC]5. remove loops with inconsistent PETs between replicates, which one PET==0 and another PET>=threshold
    # * [QC]6. remove loops with FDR > 0.01 [default: 0.01]
    # data = remove_by_inconsistent_PETs_between_replicates(data, samples, groups, pet_threshold, output_path, next_func_depth)
    data = remove_by_inconsistent_PETs_and_fdr_thres(data, samples, groups, pet_threshold, fdr_threshold, output_path, next_func_depth)

    # return: list of list of LoopInfo
    # the index of data is samples index, data[0] means the first input sample, which is consistent with "samples" and "groups" parameters
    # the index of data[0] is LoopInfo index, data[0][0] means the first LoopInfo in first input sample
    return data



