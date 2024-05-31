
import os
import numpy as np

ROOT_PATH = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
f_path = os.path.join(ROOT_PATH, "test", "reproducibility", "temp_output", "reproducibility_vs_dif_block_size.txt")

with open(f_path, 'r') as f:
    lines = f.readlines()

block_sizes = np.array([1E6, 2E6, 3E6, 4E6, 5E6], dtype=int)
reproducibility_dict = dict()
reproducibility_dict["k562_rep1_vs_k562_rep2"] = list()
reproducibility_dict["k562_rep1_vs_mcf7_rep1"] = list()
reproducibility_dict["k562_rep1_vs_mcf7_rep2"] = list()
reproducibility_dict["k562_rep2_vs_mcf7_rep1"] = list()
reproducibility_dict["k562_rep2_vs_mcf7_rep2"] = list()
reproducibility_dict["mcf7_rep1_vs_mcf7_rep2"] = list()

for line in lines:
    if not line.startswith("#"):
        fields = line.strip().split()
        reproducibility_dict[fields[0]].append(float(fields[1]))

# print(reproducibility_dict)

import matplotlib.pyplot as plt

plt.figure(dpi=300, figsize=(8, 6))
FONTSIZE = 10
plt.rc("font", family="Times New Roman")
params = {"axes.titlesize": FONTSIZE,
            "legend.fontsize": 18,
            "axes.labelsize": 26,
            "xtick.labelsize": 26,
            "ytick.labelsize": 26,
            "figure.titlesize": FONTSIZE,
            "font.size": FONTSIZE}
plt.rcParams.update(params)

block_sizes = block_sizes / 1E6
plt.plot(block_sizes, reproducibility_dict["k562_rep1_vs_k562_rep2"], linewidth=3, marker="*", markersize=15, c="#00A0E9", label="k562_rep1_vs_k562_rep2")
plt.plot(block_sizes, reproducibility_dict["k562_rep1_vs_mcf7_rep1"], marker="+", c="#22B14C", label="k562_rep1_vs_mcf7_rep1")
plt.plot(block_sizes, reproducibility_dict["k562_rep1_vs_mcf7_rep2"], marker="+", c="#A349A4", label="k562_rep1_vs_mcf7_rep2")
plt.plot(block_sizes, reproducibility_dict["k562_rep2_vs_mcf7_rep1"], marker="+", c="#ED1C24", label="k562_rep2_vs_mcf7_rep1")
plt.plot(block_sizes, reproducibility_dict["k562_rep2_vs_mcf7_rep2"], marker="+", c="#FFAEC9", label="k562_rep2_vs_mcf7_rep2")
plt.plot(block_sizes, reproducibility_dict["mcf7_rep1_vs_mcf7_rep2"], linewidth=3, marker="D", markersize=10, c="#FFE200", label="mcf7_rep1_vs_mcf7_rep2")

plt.xticks(block_sizes)
# 显示树状图
# plt.title('Dendrogram using Single Linkage')
plt.xlabel('Block size (1E6 bp)')
plt.ylabel('Reproducibility')
plt.legend(loc='best')
plt.tight_layout()
output_path = os.path.join(ROOT_PATH, "test", "reproducibility", "temp_output")
plt.savefig(f"{output_path}/reproducibility_vs_dif_block_size.png", format="png", bbox_inches="tight")
# plt.show()
plt.close()



