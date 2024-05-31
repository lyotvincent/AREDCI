import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from brokenaxes import brokenaxes
import matplotlib.gridspec as gridspec

CUR_PATH = os.path.dirname(os.path.abspath(__file__))

EXCEL_PATH = os.path.join(CUR_PATH, "result_records.xlsx")

df = pd.read_excel(EXCEL_PATH, sheet_name="sim_test_v2", header=None, index_col=None)

# print(df)
fig = plt.figure(dpi=300, figsize=(8, 10))
FONTSIZE = 10
plt.rc("font", family="Times New Roman")
params = {"axes.titlesize": FONTSIZE,
            "legend.fontsize": FONTSIZE,
            "axes.labelsize": FONTSIZE,
            "xtick.labelsize": FONTSIZE,
            "ytick.labelsize": FONTSIZE,
            "figure.titlesize": FONTSIZE,
            "font.size": FONTSIZE}
plt.rcParams.update(params)

color_list = ["#00A2E8", "#ED1C24", "#22B14C", "#FFC90E"]

# accuracy, start from row 11
x = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
edgeR_accuracy = df.iloc[11:17, 1].values
kde_accuracy = df.iloc[11:17, 2].values
md_accuracy = df.iloc[11:17, 4].values
ssrd_accuracy = df.iloc[11:17, 6].values
# plt.subplot(3, 2, 1)
# plt.plot(x, edgeR_accuracy, label="edgeR", color=color_list[0], marker="x", linestyle="--")
# plt.plot(x, kde_accuracy, label="kde", color=color_list[1], marker="*")
# plt.plot(x, md_accuracy, label="md", color=color_list[2], marker="^")
# plt.plot(x, ssrd_accuracy, label="ssrd", color=color_list[3], marker="o")
# plt.ylabel("Accuracy")
# plt.xlabel("pOutlier")
# plt.xticks(x)
# plt.yticks([0.7, 0.75, 0.8, 0.85, 0.9])
# plt.legend(loc="best")

gs = gridspec.GridSpec(3, 2)  # Change this to match your grid size
subplot_spec = gs[0, 0]  # Change this to match your desired subplot location
bax = brokenaxes(fig=fig, subplot_spec=subplot_spec, ylims=((0.73, 0.76), (0.88, 0.91)), hspace=.2, despine=False, diag_color="r", d=0.)
bax.plot(x, edgeR_accuracy, label="edgeR", color=color_list[0], marker="x", linestyle="--")
bax.plot(x, kde_accuracy, label="kde", color=color_list[1], marker="*")
bax.plot(x, md_accuracy, label="md", color=color_list[2], marker="^")
bax.plot(x, ssrd_accuracy, label="ssrd", color=color_list[3], marker="o")

bax.set_ylabel("Accuracy")
bax.set_xlabel("pOutlier")
bax.legend(loc=(0.7, 0.3))


# precision, start from row 20
edgeR_precision = df.iloc[20:26, 1].values
kde_precision = df.iloc[20:26, 2].values
md_precision = df.iloc[20:26, 4].values
ssrd_precision = df.iloc[20:26, 6].values
plt.subplot(3, 2, 2)
plt.plot(x, edgeR_precision, label="edgeR", color=color_list[0], marker="x", linestyle="--")
plt.plot(x, kde_precision, label="kde", color=color_list[1], marker="*")
plt.plot(x, md_precision, label="md", color=color_list[2], marker="^")
plt.plot(x, ssrd_precision, label="ssrd", color=color_list[3], marker="o")

plt.ylabel("Precision")
plt.xlabel("pOutlier")
plt.xticks(x)
plt.legend(loc=(0.05, 0.1))

# recall, start from row 29
edgeR_recall = df.iloc[29:35, 1].values
kde_recall = df.iloc[29:35, 2].values
md_recall = df.iloc[29:35, 4].values
ssrd_recall = df.iloc[29:35, 6].values
plt.subplot(3, 2, 3)
plt.plot(x, edgeR_recall, label="edgeR", color=color_list[0], marker="x", linestyle="--")
plt.plot(x, kde_recall, label="kde", color=color_list[1], marker="*")
plt.plot(x, md_recall, label="md", color=color_list[2], marker="^")
plt.plot(x, ssrd_recall, label="ssrd", color=color_list[3], marker="o")

plt.ylabel("Recall")
plt.xlabel("pOutlier")
plt.xticks(x)
plt.legend(loc="best")

# f1, start from row 38
edgeR_f1 = df.iloc[38:44, 1].values 
kde_f1 = df.iloc[38:44, 2].values
md_f1 = df.iloc[38:44, 4].values
ssrd_f1 = df.iloc[38:44, 6].values
plt.subplot(3, 2, 4)
plt.plot(x, edgeR_f1, label="edgeR", color=color_list[0], marker="x", linestyle="--")
plt.plot(x, kde_f1, label="kde", color=color_list[1], marker="*")
plt.plot(x, md_f1, label="md", color=color_list[2], marker="^")
plt.plot(x, ssrd_f1, label="ssrd", color=color_list[3], marker="o")

plt.ylabel("F1-Score")
plt.xlabel("pOutlier")
plt.xticks(x)
plt.yticks(np.arange(0, 0.21, 0.03))
plt.legend(loc="best")

# auroc, start from row 47
edgeR_auroc = df.iloc[47:53, 1].values
kde_auroc = df.iloc[47:53, 2].values
md_auroc = df.iloc[47:53, 4].values
ssrd_auroc = df.iloc[47:53, 6].values
plt.subplot(3, 2, 5)
plt.plot(x, edgeR_auroc, label="edgeR", color=color_list[0], marker="x", linestyle="--")
plt.plot(x, kde_auroc, label="kde", color=color_list[1], marker="*")
plt.plot(x, md_auroc, label="md", color=color_list[2], marker="^")
plt.plot(x, ssrd_auroc, label="ssrd", color=color_list[3], marker="o")

plt.ylabel("AUROC")
plt.xlabel("pOutlier")
plt.xticks(x)
plt.legend(loc=(0.7, 0.45))

# auprc, start from row 56
edgeR_auroc = df.iloc[56:62, 1].values
kde_auroc = df.iloc[56:62, 2].values
md_auroc = df.iloc[56:62, 4].values
ssrd_auroc = df.iloc[56:62, 6].values
plt.subplot(3, 2, 6)
plt.plot(x, edgeR_auroc, label="edgeR", color=color_list[0], marker="x", linestyle="--")
plt.plot(x, kde_auroc, label="kde", color=color_list[1], marker="*")
plt.plot(x, md_auroc, label="md", color=color_list[2], marker="^")
plt.plot(x, ssrd_auroc, label="ssrd", color=color_list[3], marker="o")

plt.ylabel("AUPRC")
plt.xlabel("pOutlier")
plt.xticks(x)
plt.yticks([0.1, 0.11, 0.12, 0.13])
plt.legend(loc=(0.7, 0.4))

# plt.xticks(rotation=45)
# plt.title('Dendrogram using Single Linkage')
# plt.xlabel('Data Points')
# plt.ylabel('Distance')
plt.tight_layout()
plt.savefig(os.path.join(CUR_PATH, "sim_metrics.png"), format="png", bbox_inches="tight")
plt.savefig(os.path.join(CUR_PATH, "sim_metrics.pdf"), format="pdf", bbox_inches="tight")
# plt.show()
plt.close()
