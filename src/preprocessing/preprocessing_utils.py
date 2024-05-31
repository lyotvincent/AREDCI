import time, math
from utils.loop_object import LoopInfo, Loop

def complement_zero_pet_loops_for_samples(data_full, func_depth=0):
    time1 = time.time()
    print('\t'*func_depth+"<FUNC> complement_zero_pet_loops_for_samples")
    samples_num = len(data_full)
    i = 0
    while i < max([len(data_full[j]) for j in range(samples_num)]):
        is_zero_appear = False
        # * 初始化一个最小的loop_id和对应的loop
        for j in range(samples_num):
            if len(data_full[j]) > i:
                min_loopid = data_full[j][i].get_loop_id()
                min_loopid_loop = data_full[j][i]
                break
        # * 找到最小的loop_id和对应的loop
        for j in range(samples_num):
            if len(data_full[j]) <= i:
                is_zero_appear = True
                continue
            current_loop_id = data_full[j][i].get_loop_id()
            if current_loop_id < min_loopid:
                min_loopid = current_loop_id
                min_loopid_loop = data_full[j][i]
                is_zero_appear = True
            elif current_loop_id > min_loopid:
                is_zero_appear = True
        # if there is no "zero pet loop", then move to next loop. Otherwise, complement zero pet loop for all samples.
        # 只是为了加快速度，如果没有0 pet的loop，就不用补全了
        if is_zero_appear == False:
            i += 1
            continue
        # print(min_loopid_loop)
        # print(data[0][i], data[1][i], data[2][i], data[3][i])
        # * 补全0 pet的loop
        for j in range(samples_num):
            if len(data_full[j]) <= i or data_full[j][i].get_loop_id() > min_loopid:
                data_full[j].insert(i, LoopInfo(loop=min_loopid_loop.loop,
                                           pet=0,
                                           fdr=math.inf))
        assert data_full[0][i] == data_full[1][i] == data_full[2][i] == data_full[3][i]
        i += 1
    print('\t'*func_depth+f"<TIME> complement_zero_pet_loops_for_samples cost {time.time()-time1} seconds")
    return data_full


def complement_zero_pet_loops_from_LoopClass(data, func_depth=0):
    time1 = time.time()
    print('\t'*func_depth+f"<FUNC> complement_zero_pet_loops_from_LoopClass, Loop num: {len(Loop._instances)}")
    samples_num = len(data)
    new_data = [list() for i in range(samples_num)]
    for i in range(samples_num):
        j = 0
        cur_sample_len = len(data[i])
        for loop in Loop._instances.values():
            if j >= cur_sample_len or loop.loop_id < data[i][j].loop.loop_id:
                new_data[i].append(LoopInfo(loop=loop,
                                            pet=0,
                                            pvalue=math.inf,
                                            fdr=math.inf))
            elif loop.loop_id == data[i][j].loop.loop_id:
                new_data[i].append(data[i][j])
                j += 1
    print('\t'*func_depth+f"<TIME> complement_zero_pet_loops_from_LoopClass cost {time.time()-time1} seconds")
    return new_data

