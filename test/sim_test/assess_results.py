
import os
import numpy as np
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, fbeta_score, roc_auc_score, roc_curve, auc, average_precision_score, precision_recall_curve


def compute_metrics(fdrs, labels, fdr_threshold=0.05):
    print(f"sum(pred_score < {fdr_threshold}): {sum(fdrs < fdr_threshold)}")
    print(f"positive labels number:negative: {sum(labels)}:{len(labels)-sum(labels)}")
    # make prediction=1 when fdrs < 0.05
    predictions = np.zeros(fdrs.shape)
    predictions[fdrs < 0.05] = 1
    # compute metrics
    acc = accuracy_score(labels, predictions)
    precision = precision_score(labels, predictions)
    recall = recall_score(labels, predictions)
    f1 = f1_score(labels, predictions)
    fbeta = fbeta_score(y_true=labels, y_pred=predictions, beta=2)

    auroc = roc_auc_score(labels, 1-fdrs)
    fpr, tpr, thresholds = roc_curve(labels, 1-fdrs)
    auroc_by_auc = auc(fpr, tpr)

    auprc = average_precision_score(labels, 1-fdrs)
    precision_point, recall_point, thresholds = precision_recall_curve(labels, 1-fdrs, pos_label=1)
    # precision_point[(recall_point==0)] = 1.0
    auprc_by_auc = auc(recall_point, precision_point)
    print(f"accuracy: {acc}, precision: {precision}, recall: {recall}, f1: {f1}, fbeta: {fbeta}, auroc: {auroc}, auprc: {auprc}")
    return fpr, tpr, auroc_by_auc, precision_point, recall_point, auprc_by_auc


def load_edgeR_sim_results(path):
    # f = open("D:/workspace/dci/data/simulated_data/edgeR_sim_results.csv", "r")
    f = open(path, "r")
    lines = f.readlines()
    f.close()
    pvalues = list()
    fdrs = list()
    for line in lines[1:]:
        line = line.strip().split(",")
        pvalues.append(float(line[4]))
        fdrs.append(float(line[5]))
    return np.array(pvalues), np.array(fdrs)


def load_sim_data(y_path):
    f = open(y_path, "r")
    lines = f.readlines()[1:]
    f.close()
    y = list()
    for line in lines:
        line = [1. if i == "TRUE" else 0. for i in line.strip().split(",")[1:]]
        assert all(line) or (not any(line))
        line = 1 if all(line) else 0
        y.append(line)
    return np.array(y)

def load_my_sim_results(path):
    f = open(path, "r")
    lines = f.readlines()
    f.close()
    pvalues = list()
    fdrs = list()
    for line in lines:
        line = line.strip().split()
        pvalues.append(float(line[0]))
        fdrs.append(float(line[1]))
    return np.array(pvalues), np.array(fdrs)

if __name__ == "__main__":
    ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    SUFFIX = "outlier_0.5"

    # * load simulated data
    SIM_DATA_DIR = os.path.join(ROOT_DIR, "data", "simulated_data")
    SIM_LABELS_DIR = os.path.join(SIM_DATA_DIR, f"sim_DE_labels.csv")
    labels = load_sim_data(SIM_LABELS_DIR)

    # * edgeR results
    # EDGER_RESULT = os.path.join(ROOT_DIR, 'test', "sim_test", f"edgeR_sim_results.{SUFFIX}.csv")
    # pvalues_edgeR, fdrs_edgeR = load_edgeR_sim_results(EDGER_RESULT)

    # * my results
    MY_RESULT = os.path.join(ROOT_DIR, 'test', "sim_test", f"my_sim_result_kde_minus.{SUFFIX}.csv")
    pvalues_kde, fdrs_kde = load_my_sim_results(MY_RESULT)
    # MY_RESULT = os.path.join(ROOT_DIR, 'test', "sim_test", f"my_sim_result_kde_absminus.{SUFFIX}.csv")
    # pvalues_kde, fdrs_kde = load_my_sim_results(MY_RESULT)
    # MY_RESULT = os.path.join(ROOT_DIR, 'test', "sim_test", f"my_sim_result_kde.{SUFFIX}.csv")
    # pvalues_kde, fdrs_kde = load_my_sim_results(MY_RESULT)

    # MY_RESULT = os.path.join(ROOT_DIR, 'test', "sim_test", f"my_sim_result_md_fit_a0.{SUFFIX}.csv")
    # pvalues_md_a0, fdrs_md_a0 = load_my_sim_results(MY_RESULT)

    # MY_RESULT = os.path.join(ROOT_DIR, 'test', "sim_test", f"my_sim_result_md_fit_a1.{SUFFIX}.csv")
    # pvalues_md_a1, fdrs_md_a1 = load_my_sim_results(MY_RESULT)

    # MY_RESULT = os.path.join(ROOT_DIR, 'test', "sim_test", f"my_sim_result_md_fit_a5.{SUFFIX}.csv")
    # pvalues_md_a5, fdrs_md_a5 = load_my_sim_results(MY_RESULT)

    # MY_RESULT = os.path.join(ROOT_DIR, 'test', "sim_test", f"my_sim_result_ssrd_fit.{SUFFIX}.csv")
    # pvalues_ssrd, fdrs_ssrd = load_my_sim_results(MY_RESULT)


    # * compute metrics
    # print("edgeR:")
    # compute_metrics(fdrs_edgeR, labels)
    print("my kde:")
    compute_metrics(fdrs_kde, labels)
    # print("my md_a0:")
    # compute_metrics(fdrs_md_a0, labels)
    # print("my md_a1:")
    # compute_metrics(fdrs_md_a1, labels)
    # print("my md_a5:")
    # compute_metrics(fdrs_md_a5, labels)
    # print("my ssrd:")
    # compute_metrics(fdrs_ssrd, labels)

