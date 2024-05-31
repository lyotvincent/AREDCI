
import os
import numpy as np
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import combine_pvalues

ROOT_PATH = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# print(ROOT_PATH)
edgeR_with_my_QC_PATH = os.path.join(ROOT_PATH, 'test', 'chip_signal', 'diffloop_interactions_with_counts_after_myQC.csv')
mine = os.path.join(ROOT_PATH, 'test', 'chip_signal', 'my_dci_result', 'my_interactions_nonorm_kde.csv')
fdr_threshold = 0.05

edgeR_pred_score = list()
merged_p_values = list()
with open(edgeR_with_my_QC_PATH, 'r') as f:
    edgeR_line = f.readline()
    while edgeR_line:
        fields = edgeR_line.strip().split()
        edgeR_pred_score.append(float(fields[6]))
        pvalues = [float(j) for j in fields[7].split(';')]
        _, merged_p_value = combine_pvalues(pvalues, method='mudholkar_george')
        merged_p_values.append(merged_p_value)
        edgeR_line = f.readline()
_, fdrs_merged = fdrcorrection(merged_p_values, alpha=0.05, method='indep', is_sorted=False)
labels = np.array(fdrs_merged<fdr_threshold, dtype=int)
edgeR_pred_score = np.array(edgeR_pred_score)
edgeR_predictions = np.zeros(edgeR_pred_score.shape, dtype=int)
edgeR_predictions[edgeR_pred_score < fdr_threshold] = 1

my_pred_score = list()
with open(mine, 'r') as f:
    mine_line = f.readline()
    while mine_line:
        fields = mine_line.strip().split()
        my_pred_score.append(float(fields[6]))
        mine_line = f.readline()
my_pred_score = np.array(my_pred_score)
print(my_pred_score[228])
my_predictions = np.zeros(my_pred_score.shape, dtype=int)
my_predictions[my_pred_score < fdr_threshold] = 1
print(my_predictions[228])

print("edgeR_predictions pos:neg : ", np.sum(edgeR_predictions), np.sum(1-edgeR_predictions))
print("my_predictions pos:neg : ", np.sum(my_predictions), np.sum(1-my_predictions))
print("labels pos:neg : ", np.sum(labels), np.sum(1-labels))
print("all pos, pos:neg : ", np.sum(edgeR_predictions & my_predictions & labels))
all_pos = list()
j=0
for i in zip(edgeR_predictions, my_predictions, labels):
    if j == 228:
        print(i, all(i))
    j+=1
    all_pos.append( all(i) )
all_pos = np.array(all_pos)
print("all pos : ", np.sum(all_pos))
print(all_pos[228])

print(np.where(all_pos==1))






