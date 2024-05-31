import os

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

def read_gff3():
    """
    读取ROOT_DIR/gencode.v45.basic.annotation.gff3文件,
    从第三列是gene的行中提取chr(first column), gene_id(9th column, ID=)。
    目前所有基因ID都是ENSG开头的。
    """
    gff3_file = os.path.join(ROOT_DIR, 'gencode.v45.basic.annotation.gff3')
    with open(gff3_file, 'r') as f:
        lines = f.readlines()
    gene_info = []
    gene_id_set = set()
    for line in lines:
        if line.startswith('#'):
            continue
        line = line.strip().split('\t')
        if line[2] == 'gene':
            gene_id = line[8].split(';')[0].split('=')[1]
            gene_info.append((line[0], gene_id))
            gene_id_set.add(gene_id)
    return gene_info, gene_id_set

def read_gtf():
    """
    读取ROOT_DIR/gencode.v45.basic.annotation.gtf文件,
    从第三列是gene的行中提取chr(first column), gene_id(9th column, gene_id)。
    目前所有基因ID都是ENSG开头的。
    """
    gtf_file = os.path.join(ROOT_DIR, 'gencode.v45.basic.annotation.gtf')
    with open(gtf_file, 'r') as f:
        lines = f.readlines()
    gene_info = []
    gene_id_set = set()
    for line in lines:
        if line.startswith('#'):
            continue
        line = line.strip().split('\t')
        if line[2] == 'gene':
            gene_id = line[8].split(';')[0].split(r'"')[1]
            gene_info.append((line[0], gene_id))
            gene_id_set.add(gene_id)
    # * gene_info: [[chr, gene_id], ...]
    return gene_info, gene_id_set

def read_tsv(file_path):
    """
    读取tsv文件，first column是gene_id，第8列是表达量posterior_mean_count。
    目前所有基因ID都是ENSG开头的。
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()
    gene_exp = []
    gene_id_set = set()
    for line in lines:
        line = line.strip().split('\t')
        gene_id = line[0]
        if not gene_id.startswith('ENSG'):
            continue
        exp = float(line[7])
        gene_exp.append((gene_id, exp))
        gene_id_set.add(gene_id)
    # * gene_info: [[gene_id, exp], ...]
    return gene_exp, gene_id_set


if __name__ == "__main__":
    gene_info, gene_id_set = read_gff3()
    # print(gene_info[:10])
    print(len(gene_info), len(gene_id_set))
    # gene_info2, gene_id_set2 = read_gtf()
    # print(len(gene_info2))
    # print(gene_id_set==gene_id_set2)

    K562_REP1 = os.path.join(ROOT_DIR, "K562", 'ENCFF421TJX.tsv')
    K562_REP2 = os.path.join(ROOT_DIR, "K562", 'ENCFF717PSW.tsv')
    K562_REP3 = os.path.join(ROOT_DIR, "K562", 'ENCFF718ZTE.tsv')
    K562_REP4 = os.path.join(ROOT_DIR, "K562", 'ENCFF928NYA.tsv')
    K562_REP5 = os.path.join(ROOT_DIR, "K562", 'ENCFF986UIO.tsv')

    K562_REP1_exp, K562_REP1_gene_id_set = read_tsv(K562_REP1)
    K562_REP2_exp, K562_REP2_gene_id_set = read_tsv(K562_REP2)
    K562_REP3_exp, K562_REP3_gene_id_set = read_tsv(K562_REP3)
    K562_REP4_exp, K562_REP4_gene_id_set = read_tsv(K562_REP4)
    K562_REP5_exp, K562_REP5_gene_id_set = read_tsv(K562_REP5)
    
    interac = gene_id_set.intersection(K562_REP1_gene_id_set)
    print(f"K562_REP1_exp len: {len(K562_REP1_exp)}, {len(K562_REP1_gene_id_set)}; gene_id interac gene_id_set: {len(interac)}")
    interac = gene_id_set.intersection(K562_REP2_gene_id_set)
    print(f"K562_REP2_exp len: {len(K562_REP2_exp)}; gene_id interac gene_id_set: {len(interac)}")
    interac = gene_id_set.intersection(K562_REP3_gene_id_set)
    print(f"K562_REP3_exp len: {len(K562_REP3_exp)}; gene_id interac gene_id_set: {len(interac)}")
    interac = gene_id_set.intersection(K562_REP4_gene_id_set)
    print(f"K562_REP4_exp len: {len(K562_REP4_exp)}; gene_id interac gene_id_set: {len(interac)}")
    interac = gene_id_set.intersection(K562_REP5_gene_id_set)
    print(f"K562_REP5_exp len: {len(K562_REP5_exp)}; gene_id interac gene_id_set: {len(interac)}")

    MCF7_REP1 = os.path.join(ROOT_DIR, "MCF-7", 'ENCFF721BRA.tsv')
    MCF7_REP1_exp, MCF7_REP1_gene_id_set = read_tsv(MCF7_REP1)
    interac = gene_id_set.intersection(MCF7_REP1_gene_id_set)
    print(f"MCF7_REP1_exp len: {len(MCF7_REP1_exp)}; gene_id interac gene_id_set: {len(interac)}")

    GM12878_REP1 = os.path.join(ROOT_DIR, "GM12878", 'ENCFF084FUG.tsv')
    GM12878_REP2 = os.path.join(ROOT_DIR, "GM12878", 'ENCFF240WBI.tsv')
    GM12878_REP3 = os.path.join(ROOT_DIR, "GM12878", 'ENCFF345SHY.tsv')
    GM12878_REP4 = os.path.join(ROOT_DIR, "GM12878", 'ENCFF362RMV.tsv')
    GM12878_REP5 = os.path.join(ROOT_DIR, "GM12878", 'ENCFF910XWA.tsv')

    GM12878_REP1_exp, GM12878_REP1_gene_id_set = read_tsv(GM12878_REP1)
    GM12878_REP2_exp, GM12878_REP2_gene_id_set = read_tsv(GM12878_REP2)
    GM12878_REP3_exp, GM12878_REP3_gene_id_set = read_tsv(GM12878_REP3)
    GM12878_REP4_exp, GM12878_REP4_gene_id_set = read_tsv(GM12878_REP4)
    GM12878_REP5_exp, GM12878_REP5_gene_id_set = read_tsv(GM12878_REP5)

    interac = gene_id_set.intersection(GM12878_REP1_gene_id_set)
    print(f"GM12878_REP1_exp len: {len(GM12878_REP1_exp)}; gene_id interac gene_id_set: {len(interac)}")
    interac = gene_id_set.intersection(GM12878_REP2_gene_id_set)
    print(f"GM12878_REP2_exp len: {len(GM12878_REP2_exp)}; gene_id interac gene_id_set: {len(interac)}")
    interac = gene_id_set.intersection(GM12878_REP3_gene_id_set)
    print(f"GM12878_REP3_exp len: {len(GM12878_REP3_exp)}; gene_id interac gene_id_set: {len(interac)}")
    interac = gene_id_set.intersection(GM12878_REP4_gene_id_set)
    print(f"GM12878_REP4_exp len: {len(GM12878_REP4_exp)}; gene_id interac gene_id_set: {len(interac)}")
    interac = gene_id_set.intersection(GM12878_REP5_gene_id_set)
    print(f"GM12878_REP5_exp len: {len(GM12878_REP5_exp)}; gene_id interac gene_id_set: {len(interac)}")
    print(K562_REP1_gene_id_set==K562_REP2_gene_id_set==K562_REP3_gene_id_set==K562_REP4_gene_id_set==K562_REP5_gene_id_set==MCF7_REP1_gene_id_set==GM12878_REP1_gene_id_set==GM12878_REP2_gene_id_set==GM12878_REP3_gene_id_set==GM12878_REP4_gene_id_set==GM12878_REP5_gene_id_set)

    import numpy as np
    import pandas as pd
    row_name = [i[0] for i in K562_REP1_exp]
    K562_REP1_exp = np.array([i[1] for i in K562_REP1_exp], dtype=np.float32)
    K562_REP2_exp = np.array([i[1] for i in K562_REP2_exp], dtype=np.float32)
    K562_REP3_exp = np.array([i[1] for i in K562_REP3_exp], dtype=np.float32)
    K562_REP4_exp = np.array([i[1] for i in K562_REP4_exp], dtype=np.float32)
    K562_REP5_exp = np.array([i[1] for i in K562_REP5_exp], dtype=np.float32)
    # 5 columns, each column is a replicate; each row is a gene_id
    K562_exp_matrix = np.stack([K562_REP1_exp, K562_REP2_exp, K562_REP3_exp, K562_REP4_exp, K562_REP5_exp], axis=1)
    print(K562_exp_matrix.shape)
    K562_exp_matrix = pd.DataFrame(K562_exp_matrix, columns=['K562_REP1', 'K562_REP2', 'K562_REP3', 'K562_REP4', 'K562_REP5'], index=row_name)
    # K562_exp_matrix.to_csv(os.path.join(ROOT_DIR, 'K562_exp_matrix.tsv'), sep=',', header=True, index=True)

    GM12878_REP1_exp = np.array([i[1] for i in GM12878_REP1_exp], dtype=np.float32)
    GM12878_REP2_exp = np.array([i[1] for i in GM12878_REP2_exp], dtype=np.float32)
    GM12878_REP3_exp = np.array([i[1] for i in GM12878_REP3_exp], dtype=np.float32)
    GM12878_REP4_exp = np.array([i[1] for i in GM12878_REP4_exp], dtype=np.float32)
    GM12878_REP5_exp = np.array([i[1] for i in GM12878_REP5_exp], dtype=np.float32)
    # 5 columns, each column is a replicate; each row is a gene_id
    GM12878_exp_matrix = np.stack([GM12878_REP1_exp, GM12878_REP2_exp, GM12878_REP3_exp, GM12878_REP4_exp, GM12878_REP5_exp], axis=1)
    print(GM12878_exp_matrix.shape)
    GM12878_exp_matrix = pd.DataFrame(GM12878_exp_matrix, columns=['GM12878_REP1', 'GM12878_REP2', 'GM12878_REP3', 'GM12878_REP4', 'GM12878_REP5'], index=row_name)
    # GM12878_exp_matrix.to_csv(os.path.join(ROOT_DIR, 'GM12878_exp_matrix.tsv'), sep=',', header=True, index=True)

    # combine K562 and GM12878
    exp_matrix = pd.concat([K562_exp_matrix, GM12878_exp_matrix], axis=1)
    exp_matrix.to_csv(os.path.join(ROOT_DIR, 'K562_GM12878_exp_matrix.csv'), sep=',', header=True, index=True)
    # combine K562 and MCF7
    MCF7_REP1_exp = np.array([i[1] for i in MCF7_REP1_exp], dtype=np.float32)
    MCF7_exp_matrix = pd.DataFrame(MCF7_REP1_exp, columns=['MCF7_REP1'], index=row_name)
    exp_matrix = pd.concat([K562_exp_matrix, MCF7_exp_matrix], axis=1)
    exp_matrix.to_csv(os.path.join(ROOT_DIR, 'K562_MCF7_exp_matrix.csv'), sep=',', header=True, index=True)
