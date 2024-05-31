import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as ticker
import numpy as np

# 输入矩阵, 基因组示例图，是画的my_interactions_nonorm_kde.csv的225行到256。
# anchors = [
#     [1, 2],
#     [3, 4],
#     [2, 4],
#     # ... 添加更多的锚点对
# ]
ROOT_PATH = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# print(ROOT_PATH)
mine = os.path.join(ROOT_PATH, 'test', 'chip_signal', 'my_dci_result', 'my_interactions_nonorm_kde.csv')
anchors = list()
fdrs = list()
with open(mine, 'r') as f:
    mine_line = f.readline()
    i = 0
    while mine_line:
        if 225 <= i <= 256:
            fields = mine_line.strip().split()
            anchors.append([(int(fields[1])+int(fields[2]))/2, (int(fields[3])+int(fields[4]))/2])
            fdrs.append(float(fields[6]))
        mine_line = f.readline()
        i += 1
fdrs = np.clip(fdrs, 0, 1)

# 设置图形大小
plt.figure(dpi=300, figsize=(10, 2))
FONTSIZE = 12
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
# x轴的范围
# x_min, x_max = min([a[0] for a in anchors]), max([a[1] for a in anchors])
# x_min -= 2.5E5
# x_max += 2.5E5
x_min, x_max = 113.7*1E6, 116.7*1E6

# y轴的范围
y_min, y_max = 0, 1

# 设置x轴和y轴的刻度
plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)

widths = [np.abs(a[0] - a[1]) for a in anchors]
width_min = min(widths)
width_max = max(widths)

# 绘制锚点
for a in anchors:
    if a[0] < 115.12*1E6:
        plt.plot(a[0], 0, 'o', color="#00A2E8")
    else:
        plt.plot(a[0], 0, 'o', color="#22B14C")
    if a[1] < 115.12*1E6:
        plt.plot(a[1], 0, 'o', color="#00A2E8")
    else:
        plt.plot(a[1], 0, 'o', color="#22B14C")

# # 定义颜色列表
# colors = ["red", "black"]
# # 创建自定义颜色映射
# new_cmap = LinearSegmentedColormap.from_list("mycmap", colors)

# # 使用颜色映射
# plt.pcolormesh(np.random.rand(10,10), cmap=new_cmap)
# plt.colorbar(shrink=1)

# 绘制关联曲线
for i in range(len(anchors)):
    a1, a2 = anchors[i]
    print(a1, a2)
    # 计算半圆的半径
    radius = np.abs(a1 - a2) / 2
    mid_point = (a1+a2)/2
    width = np.abs(a1 - a2)
    height = (width-width_min)/(width_max-width_min) + 1 - 1E-2
    print(width, height)
    # 创建半圆的路径
    circle = patches.Ellipse((mid_point, 0), width, height, angle=180, color=[1-fdrs[i],0,0,1], fill=False)
    # 添加到图形中
    plt.gca().add_artist(circle)

# 创建一个自定义的刻度格式器
def custom_formatter(x, pos):
    return f"{float(x/1E6):.1f}"

# 应用自定义刻度格式器到y轴
ax.xaxis.set_major_formatter(ticker.FuncFormatter(custom_formatter))
plt.yticks([])
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
# 显示图形
plt.savefig(os.path.join(ROOT_PATH, 'test', 'chip_signal', 'genome_example.pdf'), format='pdf', bbox_inches='tight')
# plt.show()
plt.close()


