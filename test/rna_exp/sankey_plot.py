
import os, time
from collections import defaultdict
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import combine_pvalues
import numpy as np


class Anchor:
    def __init__(self, anchor_id, chr, start, end):
        self.anchor_id = anchor_id
        self.chr = chr
        self.start = start
        self.end = end

    def __str__(self):
        return f"({self.anchor_id},{self.chr},{self.start},{self.end})"

class Interaction:
    def __init__(self, anchor1, anchor2, pvalue, fdr):
        self.chr = anchor1.chr
        self.anchor1 = anchor1
        self.anchor2 = anchor2
        self.pvalue = pvalue
        self.fdr = fdr
        self.genes = list()
        self.ocrs = list()
        self.ctcfchips = list()

    def __str__(self):
        return f"[{self.chr},{self.anchor1},{self.anchor2},{self.pvalue},{self.fdr}]"

class Gene:
    def __init__(self, gene_id, chr, start, end, pvalue, fdr):
        self.gene_id = gene_id
        self.chr = chr
        self.start = start
        self.end = end
        self.pvalue = pvalue
        self.fdr = fdr

    def __str__(self):
        return f"({self.gene_id},{self.chr},{self.start},{self.end},{self.pvalue},{self.fdr})"

class OCR:
    def __init__(self, ocr_id, chr, start, end, pValue, fdr):
        self.ocr_id = ocr_id
        self.chr = chr
        self.start = start
        self.end = end
        self.pValue = pValue
        self.fdr = fdr

    def __str__(self):
        return f"({self.ocr_id},{self.chr},{self.start},{self.end},{self.pValue},{self.fdr})"

def read_gff3():
    """
    读取ROOT_DIR/gencode.v45.basic.annotation.gff3文件,
    从第三列是gene的行中提取chr(first column), gene_id(9th column, ID=)。
    目前所有基因ID都是ENSG开头的。
    """
    ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    RNA_EXP_DIR = os.path.join(ROOT_DIR, 'data', 'RNA_expression')
    GFF_PATH = os.path.join(RNA_EXP_DIR, 'gencode.v45.basic.annotation.gff3')
    # rna_pvalue_file = os.path.join(RNA_EXP_DIR, 'K562_MCF7_edgeR_RNAcounts_results.csv')
    rna_pvalue_file = os.path.join(RNA_EXP_DIR, 'K562_MCF7_pvalues.csv')
    with open(GFF_PATH, 'r') as f:
        lines = f.readlines()
    gff_info = dict()
    for line in lines:
        if line.startswith('#'):
            continue
        line = line.strip().split('\t')
        if line[2] == 'gene':
            gene_id = line[8].split(';')[0].split('=')[1]
            gff_info[gene_id] = (line[0][3:], int(line[3]), int(line[4])) # chr, start, end

    with open(rna_pvalue_file, 'r') as f:
        lines = f.readlines()
    genes = list()
    for line in lines:
        line = line.strip().split(',')
        gene_id = line[0].strip('"')
        if gene_id in gff_info:
            genes.append(Gene(gene_id, gff_info[gene_id][0], gff_info[gene_id][1], gff_info[gene_id][2], float(line[4]), float(line[5])))
    return genes

def load_k562_mcf7_ctcfs():
    ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    OCR_POS_FILE = os.path.join(ROOT_DIR, 'data', 'CTCF_ChIP_seq', 'chipseq_k562_mcf7_combined_peak.bed')
    # OCR_FDR_FILE = os.path.join(ROOT_DIR, 'data', 'CTCF_ChIP_seq', 'K562_MCF7_edgeR_ChIP_signal_results.csv')
    OCR_FDR_FILE = os.path.join(ROOT_DIR, 'data', 'CTCF_ChIP_seq', 'K562_MCF7_pvalues.csv')
    with open(OCR_POS_FILE, 'r') as f:
        lines = f.readlines()
    with open(OCR_FDR_FILE, 'r') as f:
        lines2 = f.readlines()[1:]
    ocrs_list = list()
    i = 0
    for line, line2 in zip(lines, lines2):
        line = line.strip().split('\t')
        line2 = line2.strip().split(',')
        ocrs_list.append(OCR(i, line[0][3:], int(line[1]), int(line[2]), float(line2[4]), float(line2[5])))
        i += 1
    return ocrs_list

def load_k562_mcf7_ocrs():
    ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    OCR_POS_FILE = os.path.join(ROOT_DIR, 'data', 'DNase_seq', 'dnaseseq_k562_mcf7_combined_peak.bed')
    # OCR_FDR_FILE = os.path.join(ROOT_DIR, 'data', 'CTCF_ChIP_seq', 'K562_MCF7_edgeR_ChIP_signal_results.csv')
    OCR_FDR_FILE = os.path.join(ROOT_DIR, 'data', 'DNase_seq', 'K562_MCF7_pvalues.csv')
    with open(OCR_POS_FILE, 'r') as f:
        lines = f.readlines()
    with open(OCR_FDR_FILE, 'r') as f:
        lines2 = f.readlines()[1:]
    ocrs_list = list()
    i = 0
    for line, line2 in zip(lines, lines2):
        line = line.strip().split('\t')
        line2 = line2.strip().split(',')
        ocrs_list.append(OCR(i, line[0][3:], int(line[1]), int(line[2]), float(line2[4]), float(line2[5])))
        i += 1
    return ocrs_list

def judge_overlap(gene, interaction):
    i_start = min(interaction.anchor1.start, interaction.anchor1.end, interaction.anchor2.start, interaction.anchor2.end)
    i_end = max(interaction.anchor1.start, interaction.anchor1.end, interaction.anchor2.start, interaction.anchor2.end)
    if gene.end < i_start or gene.start > i_end:
        return False
    return True

def judge_overlap2(gene, interaction):
    # if not in two anchors region
    if gene.end < interaction.anchor1.start or \
        (gene.start > interaction.anchor1.end and gene.end < interaction.anchor2.start) or \
        gene.start > interaction.anchor2.end:
        return False
    return True

def load_my_result(my_result_file):
    time1 = time.time()
    ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    TEST_DIR = os.path.join(ROOT_DIR, 'test', 'rna_exp')
    # * Load my result interactions
    with open(os.path.join(TEST_DIR, 'my_dci_result', my_result_file), 'r') as f:
        lines = f.readlines()
    interactions = list()
    for line in lines:
        line = line.strip().split('\t')
        anchor1 = Anchor(0, line[0][3:], int(line[1]), int(line[2]))
        anchor2 = Anchor(0, line[0][3:], int(line[3]), int(line[4]))
        interactions.append(Interaction(anchor1, anchor2, float(line[5]), float(line[6])))
    print(interactions[0])
    print(len(interactions)) # 43710
    print(type(interactions[0].chr))

    # * Load genes from gff3 and rna_pvalue_file
    genes = read_gff3()
    print(genes[0])
    print(len(genes)) # 34155
    print(type(genes[0].chr))
    genes_dict = defaultdict(list)
    for gene in genes:
        genes_dict[gene.chr].append(gene)

    # * Load ctcfs from gff3 and edgeR_rna_file
    ctcfs = load_k562_mcf7_ctcfs()
    print(ctcfs[0])
    print(len(ctcfs)) # 115282
    print(type(ctcfs[0].chr))
    ctcfs_dict = defaultdict(list)
    for ctcf in ctcfs:
        ctcfs_dict[ctcf.chr].append(ctcf)

    # * load ocrs from bed file and edgeR_rna_file
    ocrs = load_k562_mcf7_ocrs()
    print(ocrs[0])
    print(len(ocrs)) # 321256
    print(type(ocrs[0].chr))
    ocrs_dict = defaultdict(list)
    for ocr in ocrs:
        ocrs_dict[ocr.chr].append(ocr)

    # * Add genes to interactions
    i = 0
    for interaction in interactions:
        i += 1
        if i % 1000 == 0:
            print(i)
        for gene in genes_dict[interaction.chr]:
            if judge_overlap(gene, interaction):
                interaction.genes.append(gene)
        for ocr in ocrs_dict[interaction.chr]:
            if judge_overlap(ocr, interaction):
                interaction.ocrs.append(ocr)
        for ctcf in ctcfs_dict[interaction.chr]:
            if judge_overlap2(ctcf, interaction):
                interaction.ctcfchips.append(ctcf)
    interactions_with_ocrs_genes = [interaction for interaction in interactions if interaction.ocrs and interaction.genes and interaction.ctcfchips]
    print(interactions_with_ocrs_genes[0])
    print(len(interactions_with_ocrs_genes)) # 22131

    pred_score = list()
    merged_gene_p_values = list()
    merged_ocr_p_values = list()
    merged_ctcf_p_values = list()
    for i in interactions_with_ocrs_genes:
        pred_score.append(i.fdr)
        gene_pvalues = [gene.pvalue for gene in i.genes]
        ocr_pvalues = [ocr.pValue for ocr in i.ocrs]
        ctcfs_pvalues = [ctcf.pValue for ctcf in i.ctcfchips]
        _, merged_gene_p_value = combine_pvalues(gene_pvalues, method='mudholkar_george')
        _, merged_ocr_p_value = combine_pvalues(ocr_pvalues, method='mudholkar_george')
        _, merged_ctcf_p_value = combine_pvalues(ctcfs_pvalues, method='mudholkar_george')
        merged_gene_p_values.append(merged_gene_p_value)
        merged_ocr_p_values.append(merged_ocr_p_value)
        merged_ctcf_p_values.append(merged_ctcf_p_value)
    merged_gene_p_values = np.array(merged_gene_p_values)
    merged_ocr_p_values = np.array(merged_ocr_p_values)
    merged_ctcf_p_values = np.array(merged_ctcf_p_values)
    merged_gene_p_values[merged_gene_p_values > 1] = 1
    merged_ocr_p_values[merged_ocr_p_values > 1] = 1
    merged_ctcf_p_values[merged_ctcf_p_values > 1] = 1
    _, gene_fdrs_merged = fdrcorrection(merged_gene_p_values, alpha=0.05, method='indep', is_sorted=False)
    _, ocr_fdrs_merged = fdrcorrection(merged_ocr_p_values, alpha=0.05, method='indep', is_sorted=False)
    _, ctcf_fdrs_merged = fdrcorrection(merged_ctcf_p_values, alpha=0.05, method='indep', is_sorted=False)
    pred_score = np.array(pred_score)
    # print(len(fdrs_merged))
    pred_score_logic = pred_score < 0.05
    gene_fdrs_merged_logic = gene_fdrs_merged <= 0.8
    print(len(gene_fdrs_merged_logic), sum(gene_fdrs_merged_logic))
    print(gene_fdrs_merged_logic)
    a = np.where(gene_fdrs_merged_logic)
    print(a)
    print(len(a[0]))
    ocr_fdrs_merged_logic = ocr_fdrs_merged <= 0.05
    print(len(ocr_fdrs_merged_logic), sum(ocr_fdrs_merged_logic))
    print(ocr_fdrs_merged_logic)
    a = np.where(ocr_fdrs_merged_logic)[0]
    a = [i for i in a if 3600<i<3750]
    print("ocr around 3600-3750")
    print(a)
    print(len(a))
    print(interactions_with_ocrs_genes[3604])
    print("=====================================")
    ctcf_fdrs_merged_logic = ctcf_fdrs_merged <= 0.05
    print(len(ctcf_fdrs_merged_logic), sum(ctcf_fdrs_merged_logic))
    print(ctcf_fdrs_merged_logic)
    a = np.where(ctcf_fdrs_merged_logic)[0]
    a = [i for i in a if 3600<i<3750]
    print("ctcf around 3600-3750")
    print(a)
    print(len(a))
    print(interactions_with_ocrs_genes[3613])
    print(interactions_with_ocrs_genes[3614])
    print(interactions_with_ocrs_genes[3615])
    print(interactions_with_ocrs_genes[3616])
    print(interactions_with_ocrs_genes[3623])
    print("=====================================")
    a = np.where(np.logical_and(ocr_fdrs_merged_logic, ctcf_fdrs_merged_logic))
    # a = np.where(np.logical_and(ocr_fdrs_merged_logic, np.logical_and(ctcf_fdrs_merged_logic, gene_fdrs_merged_logic)))
    print(a)
    print(len(a[0]))
    b = [all(i) for i in zip(ocr_fdrs_merged_logic, ctcf_fdrs_merged_logic)]
    print(len(b), sum(b))
    # for i in zip( ocr_fdrs_merged_logic, ctcf_fdrs_merged_logic):
    #     print(i)
    #     print(all(i))
    #     exit()
    # print(np.where(pred_score < 0.05))
    # print(np.where(gene_fdrs_merged < 0.05))
    # intersection
    # print(np.where(np.logical_and(pred_score < 0.05, gene_fdrs_merged < 0.05))) # 3616
    # print(pred_score[3616])
    # print(gene_fdrs_merged[3616])
    # print(interactions_with_ocrs_genes[3616])
    # for g in interactions_with_ocrs_genes[3616].genes:
    #     print(g)

    indices = [3604, 3613, 3616, 3617, 3623, 3625]
    for i in indices:
        print("gene fdr:", gene_fdrs_merged[i])
        print("ocr fdr:", ocr_fdrs_merged[i])
        # print("ocr num:", len(interactions_with_ocrs_genes[i].ocrs))
        print("ctcf fdr:", ctcf_fdrs_merged[i])
        gene_ids = "/".join([g.gene_id for g in interactions_with_ocrs_genes[i].genes])
        print("gene ids:", gene_ids)
# https://cnsknowall.com/#/Home/SankeyBubbleDiagram?pid=20701008

    print("Time:", time.time()-time1)


if __name__ == "__main__":
    load_my_result("my_dci_result_nonorm_neighbor.csv")
