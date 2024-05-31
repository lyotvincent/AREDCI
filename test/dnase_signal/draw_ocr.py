
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import LinearSegmentedColormap


ROOT_PATH = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))



# k562_rep1 = pd.read_csv("./K562/ENCFF221SKA.bed", sep="\t", header=None, names=["chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"])
# k562_rep2 = pd.read_csv("./K562/ENCFF582SNT.bed", sep="\t", header=None, names=["chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"])
# k562_rep3 = pd.read_csv("./K562/ENCFF660GHM.bed", sep="\t", header=None, names=["chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"])
# k562_rep4 = pd.read_csv("./K562/ENCFF736NYC.bed", sep="\t", header=None, names=["chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"])
# k562_rep5 = pd.read_csv("./K562/ENCFF769AUF.bed", sep="\t", header=None, names=["chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"])
# mcf7_rep1 = pd.read_csv("./MCF7/ENCFF138LHE.bed", sep="\t", header=None, names=["chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"])
# mcf7_rep2 = pd.read_csv("./MCF7/ENCFF163JHE.bed", sep="\t", header=None, names=["chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"])
# mcf7_rep3 = pd.read_csv("./MCF7/ENCFF278FNP.bed", sep="\t", header=None, names=["chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"])
# mcf7_rep4 = pd.read_csv("./MCF7/ENCFF542DEP.bed", sep="\t", header=None, names=["chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"])
# mcf7_rep5 = pd.read_csv("./MCF7/ENCFF596MGE.bed", sep="\t", header=None, names=["chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"])

k562_rep1_path = os.path.join(ROOT_PATH, 'data', 'DNase_seq', 'K562', 'ENCFF185XRG.bed')
k562_rep2_path = os.path.join(ROOT_PATH, 'data', 'DNase_seq', 'K562', 'ENCFF274YGF.bed')
k562_rep3_path = os.path.join(ROOT_PATH, 'data', 'DNase_seq', 'K562', 'ENCFF628UDH.bed')
k562_rep4_path = os.path.join(ROOT_PATH, 'data', 'DNase_seq', 'K562', 'ENCFF926HBH.bed')
mcf7_rep1_path = os.path.join(ROOT_PATH, 'data', 'DNase_seq', 'MCF7', 'ENCFF438LQM.bed')
mcf7_rep2_path = os.path.join(ROOT_PATH, 'data', 'DNase_seq', 'MCF7', 'ENCFF522NDW.bed')
mcf7_rep3_path = os.path.join(ROOT_PATH, 'data', 'DNase_seq', 'MCF7', 'ENCFF774JMS.bed')
mcf7_rep4_path = os.path.join(ROOT_PATH, 'data', 'DNase_seq', 'MCF7', 'ENCFF835KCG.bed')

# chr1
x_min, x_max = 113.7*1E6, 116.7*1E6

# all int between x_min and x_max
x = np.arange(int(x_min), int(x_max), 1)

y_k562 = [list() for _ in range(len(x))]
y_mcf7 = [list() for _ in range(len(x))]

def judge_overlap(start1, end1):
    start2, end2 = int(x_min), int(x_max)
    if start1 > end2 or start2 > end1:
        return False
    return True

def add_ocr_to_y_k562(file_path):
    global y_k562
    with open(file_path, 'r') as f:
        line = f.readline()
        while line:
            fields = line.strip().split()
            if fields[0] != 'chr1':
                line = f.readline()
                continue
            if not judge_overlap(int(fields[1]), int(fields[2])):
                line = f.readline()
                continue
            start = max(int(fields[1]), int(x_min))
            end = min(int(fields[2]), int(x_max))
            signal_value = float(fields[6])
            for i in range(start, end):
                y_k562[i-int(x_min)].append(signal_value)
            line = f.readline()

def add_ocr_to_y_mcf7(file_path):
    global y_mcf7
    with open(file_path, 'r') as f:
        line = f.readline()
        while line:
            fields = line.strip().split()
            if fields[0] != 'chr1':
                line = f.readline()
                continue
            if not judge_overlap(int(fields[1]), int(fields[2])):
                line = f.readline()
                continue
            start = max(int(fields[1]), int(x_min))
            end = min(int(fields[2]), int(x_max))
            signal_value = float(fields[6])
            for i in range(start, end):
                y_mcf7[i-int(x_min)].append(signal_value)
            line = f.readline()

def draw_bar(x,y,name,format):
    OUT_PATH = os.path.join(ROOT_PATH, 'test', 'dnase_signal', f'{name}.{format}')
    plt.figure(dpi=300, figsize=(20, 4))
    FONTSIZE = 20
    plt.rc("font", family="Times New Roman")
    params = {"axes.titlesize": FONTSIZE,
                "legend.fontsize": FONTSIZE,
                "axes.labelsize": FONTSIZE,
                "xtick.labelsize": FONTSIZE,
                "ytick.labelsize": FONTSIZE,
                "figure.titlesize": FONTSIZE,
                "font.size": FONTSIZE}
    plt.rcParams.update(params)
    ax = plt.gca()
    # 设置x轴和y轴的刻度
    plt.xlim(x_min, x_max)
    plt.ylim(0, 30)

    # y0_indexes = np.where(y != 0)[0]
    # x_new = x[y0_indexes]
    # y_new = y[y0_indexes]
    # c_list_new = c_list[y0_indexes]
    print("x len:", len(x))
    # print(x[np.where(y>400)[0]])
    # plt.bar(x_new, y_new, width=1, color=c_list_new)
    # plt.bar(x, y, width=1, color=c_list)
    for i in range(len(x)-1):
         plt.plot(x[i:i+2], y[i:i+2], color=c_list[i])
    # plt.plot(x, y)

    def custom_formatter(x, pos):
        return f"{float(x/1E6):.1f}"
    # 应用自定义刻度格式器到y轴
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(custom_formatter))

    plt.xticks([])
    plt.yticks(rotation=90)
    plt.tight_layout()
    # 显示图形
    plt.savefig(OUT_PATH, format=format, bbox_inches='tight')
    # plt.show()
    plt.close()

add_ocr_to_y_k562(k562_rep1_path)
add_ocr_to_y_k562(k562_rep2_path)
add_ocr_to_y_k562(k562_rep3_path)
add_ocr_to_y_k562(k562_rep4_path)

y_k562_new = list()
for i in y_k562:
    if i:
        y_k562_new.append(np.mean(i))
    else:
        y_k562_new.append(0.)
y_k562 = np.array(y_k562_new)

add_ocr_to_y_mcf7(mcf7_rep1_path)
add_ocr_to_y_mcf7(mcf7_rep2_path)
add_ocr_to_y_mcf7(mcf7_rep3_path)
add_ocr_to_y_mcf7(mcf7_rep4_path)

y_mcf7_new = list()
for i in y_mcf7:
    if i:
        y_mcf7_new.append(np.mean(i))
    else:
        y_mcf7_new.append(0.)
y_mcf7 = np.array(y_mcf7_new)

y_dif = np.abs(y_k562 - y_mcf7)

c_list = np.abs(y_k562 - y_mcf7)
# c_list = (c_list - np.min(c_list) ) / (np.max(c_list) - np.min(c_list))
c_list[c_list > 25] = 25
c_list = c_list / 25
c_list = [[i, 0, 0, 1] for i in c_list]
c_list = np.array(c_list, dtype=np.float32)

# # * increase speed of plt.plot drawing, filter 0 ocr
y_k562_no0 = np.where(y_k562 != 0)
y_mcf7_no0 = np.where(y_mcf7 != 0)
# y_k562_no0 和 y_mcf7_no0 的并集
any_no0 = np.union1d(y_k562_no0, y_mcf7_no0)
any_no0.sort()
print(len(any_no0))
# y_k562_no0 和 y_mcf7_no0 的交
# both_no0 = np.intersect1d(y_k562_no0, y_mcf7_no0)
# 去掉0的y_k562和y_mcf7和x
for i in range(len(any_no0)-1, 0, -1):
    # print(x.shape, x[:any_no0[i-1]+1].shape, x[any_no0[i]-1:].shape)
    if any_no0[i]-any_no0[i-1] < 6:
        continue
    x = np.concatenate((x[:any_no0[i-1]+2], x[any_no0[i]-2:]))
    y_k562_no0 = np.concatenate((y_k562[:any_no0[i-1]+2], y_k562[any_no0[i]-2:]))
    y_mcf7_no0 = np.concatenate((y_mcf7[:any_no0[i-1]+2], y_mcf7[any_no0[i]-2:]))
    c_list = np.concatenate((c_list[:any_no0[i-1]+2], c_list[any_no0[i]-2:]))
    y_dif = np.concatenate((y_dif[:any_no0[i-1]+2], y_dif[any_no0[i]-2:]))

print("len(x):", len(x))
print("len(y_k562_no0):", len(y_k562_no0))
print("len(y_mcf7_no0):", len(y_mcf7_no0))
print("len(c_list):", len(c_list))
print("len(y_dif):", len(y_dif))
print(max(y_dif))

# draw_bar(x, y_k562_no0, 'k562_chip', 'png')
# draw_bar(x, y_mcf7_no0, 'mcf7_chip', 'png')


draw_bar(x, y_dif, 'k562_mcf7_ocr_diff', 'pdf')

