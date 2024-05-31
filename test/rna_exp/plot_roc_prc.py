
import os
import numpy as np
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import combine_pvalues
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, fbeta_score, roc_auc_score, roc_curve, auc, average_precision_score, precision_recall_curve, RocCurveDisplay, PrecisionRecallDisplay


def metrics_compute(filename, fdr_threshold=0.05):
    ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    TEST_DIR = os.path.join(ROOT_DIR, 'test', 'rna_exp')
    print(filename)
    with open(os.path.join(TEST_DIR, filename), 'r') as f:
        lines = f.readlines()
    pred_score = list()
    labels = list()
    merged_p_values = list()
    for i in lines:
        i = i.strip().split('\t')
        pred_score.append(float(i[6]))
        # fdrs = i[7].split(';')
        # has_diff = False
        # for fdr in fdrs:
        #     if float(fdr) < fdr_threshold:
        #         has_diff = True
        #         break
        # labels.append(int(has_diff))
        pvalues = [float(j) for j in i[7].split(';')]
        _, merged_p_value = combine_pvalues(pvalues, method='mudholkar_george')
        merged_p_values.append(merged_p_value)
    _, fdrs_merged = fdrcorrection(merged_p_values, alpha=0.05, method='indep', is_sorted=False)
    pred_score = np.array(pred_score)
    labels = np.array(fdrs_merged<fdr_threshold, dtype=int)

    print(f"sum(pred_score < {fdr_threshold}): {sum(pred_score < fdr_threshold)}")
    print(f"positive labels number:negative: {sum(labels)}:{len(labels)-sum(labels)}")
    # make prediction=1 when pred_score < fdr_threshold
    predictions = np.zeros(pred_score.shape)
    predictions[pred_score < fdr_threshold] = 1
    # compute metrics
    acc = accuracy_score(labels, predictions)
    precision = precision_score(labels, predictions)
    recall = recall_score(labels, predictions)
    f1 = f1_score(labels, predictions)
    fbeta = fbeta_score(y_true=labels, y_pred=predictions, beta=2)

    auroc = roc_auc_score(labels, 1-pred_score)
    fpr, tpr, thresholds = roc_curve(labels, 1-pred_score)
    auroc_by_auc = auc(fpr, tpr)

    auprc = average_precision_score(labels, 1-pred_score)
    precision_point, recall_point, thresholds = precision_recall_curve(labels, 1-pred_score, pos_label=1)
    # precision_point[(recall_point==0)] = 1.0
    auprc_by_auc = auc(recall_point, precision_point)
    print(f"accuracy: {acc}, precision: {precision}, recall: {recall}, f1: {f1}, fbeta: {fbeta}, auroc: {auroc}, auprc: {auprc}")

    return fpr, tpr, auroc, precision_point, recall_point, auprc

def plot_roc():
    ROOT_PATH = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    import matplotlib.pyplot as plt
    plt.figure(dpi=300, figsize=(8, 6))
    FONTSIZE = 16
    plt.rc("font", family="Times New Roman")
    params = {"axes.titlesize": FONTSIZE,
                "legend.fontsize": 14,
                "axes.labelsize": FONTSIZE,
                "xtick.labelsize": FONTSIZE,
                "ytick.labelsize": FONTSIZE,
                "figure.titlesize": FONTSIZE,
                "font.size": FONTSIZE}
    plt.rcParams.update(params)
    """
    ['Solarize_Light2', '_classic_test_patch', '_mpl-gallery', '_mpl-gallery-nogrid',
    'bmh', 'classic', 'dark_background', 'fast', 'fivethirtyeight', 'ggplot', 'grayscale',
    'seaborn-v0_8', 'seaborn-v0_8-bright', 'seaborn-v0_8-colorblind', 'seaborn-v0_8-dark',
    'seaborn-v0_8-dark-palette', 'seaborn-v0_8-darkgrid', 'seaborn-v0_8-deep', 'seaborn-v0_8-muted',
    'seaborn-v0_8-notebook', 'seaborn-v0_8-paper', 'seaborn-v0_8-pastel', 'seaborn-v0_8-poster',
    'seaborn-v0_8-talk', 'seaborn-v0_8-ticks', 'seaborn-v0_8-white', 'seaborn-v0_8-whitegrid',
    'tableau-colorblind10']
    """
    # plt.style.use('ggplot')
    # plt.style.use('bmh')
    ax = plt.gca()
    # roc_display1 = RocCurveDisplay.from_predictions(labels, 1-pred_score, name="mine", color="red")
    # roc_display1 = RocCurveDisplay(fpr=fpr, tpr=tpr, roc_auc=auroc)
    # roc_display1.plot(ax=ax, name="mine", color="red")
    # pr_display1 = PrecisionRecallDisplay.from_predictions(labels, 1-pred_score, name="mine", color="red")
    # pr_display1 = PrecisionRecallDisplay(precision=precision_point, recall=recall_point)
    # pr_display1.plot(ax=ax, name="mine", color="red")

    fpr, tpr, auroc, _, _, _ = metrics_compute(os.path.join("my_dci_result", "my_interactions_nonorm_neighbor.csv"), 0.05)
    roc_display1 = RocCurveDisplay(fpr=fpr, tpr=tpr, roc_auc=auroc)
    roc_display1.plot(ax=ax, name="AREDCI", color="#ED1C24")
    fpr, tpr, auroc, _, _, _ = metrics_compute(os.path.join("diffloop_interactions.csv"), 0.05)
    roc_display1 = RocCurveDisplay(fpr=fpr, tpr=tpr, roc_auc=auroc)
    roc_display1.plot(ax=ax, name="diffloop", color="#22B14C")
    fpr, tpr, auroc, _, _, _ = metrics_compute(os.path.join("diffloop_interactions_with_counts_after_myQC.csv"), 0.05)
    roc_display1 = RocCurveDisplay(fpr=fpr, tpr=tpr, roc_auc=auroc)
    roc_display1.plot(ax=ax, name="edgeR_after_myQC", color="#00A2E8")

    plt.plot([-0.05, 1.05], [-0.05, 1.05], color='gray', lw=1, linestyle='--', alpha=0.6)
    ax.set_aspect('equal')
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    # plt.title('Receiver operating characteristic')
    plt.legend(loc='best')
    # plt.tight_layout()
    output_path = os.path.join(ROOT_PATH, "test", "rna_exp")
    plt.savefig(f"{output_path}/roc.pdf", format="pdf", bbox_inches="tight")
    # plt.show()
    plt.close()

def plot_acc_fscore():
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker

    # 假设这是三种方法的accuracy和F1-score
    accuracy = [0.9504319812564065, 0.992963528031373, 0.553704788402401]
    f1_score = [0.0014749262536873158, 0., 0.000491924243666475]

    # 设置柱状图的x轴位置
    x_pos = np.arange(len(accuracy))
    width = 0.35

    # 创建一个图表
    plt.figure(dpi=300, figsize=(8, 6))
    FONTSIZE = 16
    plt.rc("font", family="Times New Roman")
    params = {"axes.titlesize": FONTSIZE,
                "legend.fontsize": FONTSIZE,
                "axes.labelsize": FONTSIZE,
                "xtick.labelsize": FONTSIZE,
                "ytick.labelsize": FONTSIZE,
                "figure.titlesize": FONTSIZE,
                "font.size": FONTSIZE}
    plt.rcParams.update(params)
    # fig, ax1 = plt.subplots()
    ax1 = plt.gca()

    # 绘制accuracy柱状图
    ax1.bar(x_pos - width/2, accuracy, width, label='Accuracy', color="#00A2E8")
    plt.yticks(rotation=90)

    # 创建一个共享相同x轴的第二个轴对象
    ax2 = ax1.twinx()

    # 绘制F1-score柱状图
    ax2.bar(x_pos + width/2, f1_score, width, label='F1-Score', color='#ED1C24')

    # 设置图表标题和轴标签
    # ax1.set_title('Comparison of Methods')
    # ax1.set_xlabel('Methods')
    ax1.set_ylabel('Accuracy', color='#00A2E8')
    ax2.set_ylabel('F1-Score (1E-2)', color='#ED1C24')

    # 设置x轴标签
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(["AREDCI", "diffloop", "edgeR_after_myQC"])

    # 创建一个自定义的刻度格式器
    def custom_formatter(x, pos):
        return f"{float(x/1E-2):.2f}"

    # 应用自定义刻度格式器到y轴
    ax2.yaxis.set_major_formatter(ticker.FuncFormatter(custom_formatter))
    plt.yticks(rotation=90)

    # 设置图例
    # ax1.legend(loc='best')
    # ax2.legend(loc='best')
    # 创建一个包含所有标签的图例
    handles, labels = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(handles + handles2, labels + labels2, loc='best')

    # 显示图表
    output_path = os.path.join(ROOT_PATH, "test", "rna_exp")
    plt.savefig(f"{output_path}/acc_fscore.pdf", format="pdf", bbox_inches="tight")
    # plt.show()
    plt.close()

if __name__ == "__main__":
    ROOT_PATH = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    plot_roc()
    plot_acc_fscore()


