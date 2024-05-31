
import os
from collections import defaultdict

def read_gff3():
    """
    读取ROOT_DIR/gencode.v45.basic.annotation.gff3文件,
    从第三列是gene的行中提取chr(first column), gene_id(9th column, ID=)。
    目前所有基因ID都是ENSG开头的。
    """
    ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    CELL_DIFF_DIR = os.path.join(ROOT_DIR, 'test', 'mouse_neuronal_differentiation')
    GFF_PATH = os.path.join(CELL_DIFF_DIR, 'gencode.vM34.basic.annotation.gff3')
    with open(GFF_PATH, 'r') as f:
        lines = f.readlines()
    # gff_info = dict()
    # for line in lines:
    #     if line.startswith('#'):
    #         continue
    #     line = line.strip().split('\t')
    #     if line[2] == 'gene':
    #         gene_id = line[8].split(';')[0].split('=')[1]
    #         gff_info[gene_id] = (line[0][3:], int(line[3]), int(line[4])) # chr, start, end
    gene_dict = defaultdict(list)
    for line in lines:
        if line.startswith('#'):
            continue
        line = line.strip().split('\t')
        if line[2] == 'gene':
            gene_annotation = line[8].split(';')
            for annotation in gene_annotation:
                if annotation.startswith('gene_name='):
                    gene_name = annotation[10:]
                    break
            gene_dict[line[0]].append((int(line[3]), int(line[4]), gene_name)) # chr, start, end
    return gene_dict

def load_dci_result(file_dir):
    """
    读取file_dir目录下的my_dci_result.csv，第7列是FDR，第1列是chr，第2，3列是anchor1的start和end，第4，5列是anchor2的start和end。
    only keep FDR < 0.05 rows.
    """
    dci_result = list()
    with open(file_dir, 'r') as f:
        lines = f.readlines()
    for line in lines:
        line = line.strip().split()
        if float(line[6]) <= 0.05:
            left = min(int(line[1]), int(line[2]), int(line[3]), int(line[4]))
            right = max(int(line[1]), int(line[2]), int(line[3]), int(line[4]))
            dci_result.append((line[0], left, right))
    return dci_result

def find_gene(file_dir):
    """
    从dci_result中找到gene_dict中的基因。
    """
    gene_dict = read_gff3()
    dci_result = load_dci_result(os.path.join(file_dir, "my_dci_result.csv"))
    result = list()
    for chr, start, end in dci_result:
        for gene in gene_dict[chr]:
            if not ( gene[1] < start or gene[0] > end ):
                result.append(gene[2])
    with open(os.path.join(file_dir, "gene_in_dci.csv"), 'w') as f:
        for gene in result:
            f.write(gene + '\n')
    # return result

if __name__ == '__main__':
    ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    # CELL_DIFF_DIR = os.path.join(ROOT_DIR, 'test', 'mouse_neuronal_differentiation', 'ESC_vs_NSC_result')
    # CELL_DIFF_DIR = os.path.join(ROOT_DIR, 'test', 'mouse_neuronal_differentiation', 'ESC_vs_NPC_result')
    CELL_DIFF_DIR = os.path.join(ROOT_DIR, 'test', 'mouse_neuronal_differentiation', 'NSC_vs_NPC_result')
    find_gene(CELL_DIFF_DIR)
