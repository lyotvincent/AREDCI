
import pandas as pd
import numpy as np


def custom_sort(x):
    # sort by int < X < Y
    letter_order = {'X': 1, 'Y': 2}  # Add more if needed
    try:
        return int(x), 0
    except ValueError:
        return float('inf'), letter_order.get(x, 0)

def get_signalValue(cellline_rep, combined_peak):
    signalValue_cellline_rep_4col = list()
    cellline_rep_4col = cellline_rep[["chr", "start", "end", "signalValue"]].values.tolist()
    cellline_rep_4col = [[peak[0], int(peak[1]), int(peak[2]), int(peak[3])] for peak in cellline_rep_4col]
    cellline_rep_4col = sorted(cellline_rep_4col, key=lambda x: (custom_sort(x[0][3:]), x[1], x[2]))
    print(f"{len(cellline_rep_4col)}/{len(combined_peak)}")
    i, j = 0, 0
    temp_signalValue = list()
    while j < len(combined_peak):
        if i < len(cellline_rep_4col):
            if combined_peak[j][0] == cellline_rep_4col[i][0]:
                if combined_peak[j][2] < cellline_rep_4col[i][1]:
                    if temp_signalValue:
                        signalValue_cellline_rep_4col.append(np.mean(temp_signalValue))
                        temp_signalValue = list()
                    else:
                        signalValue_cellline_rep_4col.append(0)
                    j += 1
                elif combined_peak[j][1] <= cellline_rep_4col[i][1] and cellline_rep_4col[i][2] <= combined_peak[j][2]:
                    temp_signalValue.append(int(cellline_rep_4col[i][3]))
                    i += 1
                else:
                    raise ValueError(f"Error: {cellline_rep_4col[i][0]}:{cellline_rep_4col[i][1]}-{cellline_rep_4col[i][2]} in? {combined_peak[j][0]}:{combined_peak[j][1]}-{combined_peak[j][2]}")
            else:
                if temp_signalValue:
                    signalValue_cellline_rep_4col.append(np.mean(temp_signalValue))
                    temp_signalValue = list()
                else:
                    signalValue_cellline_rep_4col.append(0)
                j += 1
        else:
            signalValue_cellline_rep_4col.append(0)
            j += 1
    # i = 0
    # for peak in combined_peak:
    #     print(f"[xx] {cellline_rep_4col[i][0]}:{cellline_rep_4col[i][1]}-{cellline_rep_4col[i][2]} in? {peak[0]}:{peak[1]}-{peak[2]}")
    #     if cellline_rep_4col[i][0] == peak[0] and peak[1] <= cellline_rep_4col[i][1] and cellline_rep_4col[i][2] <= peak[2]:
    #         signalValue_cellline_rep_4col.append(int(cellline_rep_4col[i][3]))
    #         i += 1
    #     else:
    #         signalValue_cellline_rep_4col.append(0)
    #     if 1500000 < peak[1] < 2000000:
    #         exit()
    print(f"{i}/{len(signalValue_cellline_rep_4col)}")
    return signalValue_cellline_rep_4col

def get_combined_peak():
    """
    从10个文件中的第一列chr，第二列start， 第三列end，找出重合的peak，并合并，最后按chr,start,end排序
    """
    # * 先获得所有10个文件中的peak
    k562_rep1 = pd.read_csv("./K562/ENCFF185XRG.bed", sep="\t", header=None, names=["chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"])
    k562_rep2 = pd.read_csv("./K562/ENCFF274YGF.bed", sep="\t", header=None, names=["chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"])
    k562_rep3 = pd.read_csv("./K562/ENCFF628UDH.bed", sep="\t", header=None, names=["chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"])
    k562_rep4 = pd.read_csv("./K562/ENCFF926HBH.bed", sep="\t", header=None, names=["chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"])
    mcf7_rep1 = pd.read_csv("./MCF7/ENCFF438LQM.bed", sep="\t", header=None, names=["chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"])
    mcf7_rep2 = pd.read_csv("./MCF7/ENCFF522NDW.bed", sep="\t", header=None, names=["chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"])
    mcf7_rep3 = pd.read_csv("./MCF7/ENCFF774JMS.bed", sep="\t", header=None, names=["chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"])
    mcf7_rep4 = pd.read_csv("./MCF7/ENCFF835KCG.bed", sep="\t", header=None, names=["chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"])

    # * 取出chr, start, end列；放到一个列表中；并排序
    peak_list = list()
    for df in [k562_rep1, k562_rep2, k562_rep3, k562_rep4, mcf7_rep1, mcf7_rep2, mcf7_rep3, mcf7_rep4]:
        peak_list.extend(df[["chr", "start", "end"]].values.tolist())
    peak_list = [[peak[0], int(peak[1]), int(peak[2])] for peak in peak_list]
    peak_list = sorted(peak_list, key=lambda x: (custom_sort(x[0][3:]), x[1], x[2]))

    # * 合并重合的peak
    combined_peak = list()
    current_peak = peak_list[0]
    for peak in peak_list[1:]:
        if peak[0] == current_peak[0] and peak[1] <= current_peak[2]:
            current_peak[2] = max(current_peak[2], peak[2])
        else:
            combined_peak.append(current_peak)
            current_peak = peak
    combined_peak.append(current_peak)
    # print(len(combined_peak), len(peak_list))
    with open("dnaseseq_k562_mcf7_combined_peak.bed", "w") as f:
        for peak in combined_peak:
            f.write(f"{peak[0]}\t{peak[1]}\t{peak[2]}\n")

    # * 将k562_rep1的signalValue填入combined_peak中
    signalValue_k562_rep1 = get_signalValue(k562_rep1, combined_peak)
    # print(len(signalValue_k562_rep1))
    signalValue_k562_rep2 = get_signalValue(k562_rep2, combined_peak)
    signalValue_k562_rep3 = get_signalValue(k562_rep3, combined_peak)
    signalValue_k562_rep4 = get_signalValue(k562_rep4, combined_peak)
    signalValue_mcf7_rep1 = get_signalValue(mcf7_rep1, combined_peak)
    signalValue_mcf7_rep2 = get_signalValue(mcf7_rep2, combined_peak)
    signalValue_mcf7_rep3 = get_signalValue(mcf7_rep3, combined_peak)
    signalValue_mcf7_rep4 = get_signalValue(mcf7_rep4, combined_peak)

    signalValue_matrix = [signalValue_k562_rep1, signalValue_k562_rep2, signalValue_k562_rep3, signalValue_k562_rep4, signalValue_mcf7_rep1, signalValue_mcf7_rep2, signalValue_mcf7_rep3, signalValue_mcf7_rep4]
    signalValue_matrix = np.array(signalValue_matrix)
    signalValue_matrix = signalValue_matrix.T
    # np.savetxt("chipseq_k562_mcf7_signalValue_matrix.csv", signalValue_matrix, fmt="%d", delimiter=",")
    signalValue_matrix = pd.DataFrame(signalValue_matrix, columns=['K562_REP1', 'K562_REP2', 'K562_REP3', 'K562_REP4', 'MCF7_REP1', 'MCF7_REP2', 'MCF7_REP3', 'MCF7_REP4'])
    signalValue_matrix.to_csv('./dnaseseq_k562_mcf7_signalValue_matrix.csv', sep=',', header=True, index=True)


if __name__ == "__main__":
    get_combined_peak()

