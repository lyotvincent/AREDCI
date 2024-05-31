

CHR_LIST = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6",\
            "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",\
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18"\
            "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

GM12878_REP1 = "../GM12878_ChIA-PET2_result/replicate1/GM12878.interactions.MICC"
GM12878_REP2 = "../GM12878_ChIA-PET2_result/replicate2/GM12878.interactions.MICC"
K562_REP1 = "../K562_ChIA-PET2_result/replicate1/K562.interactions.MICC"
K562_REP2 = "../K562_ChIA-PET2_result/replicate2/K562.interactions.MICC"
MCF7_REP1 = "../MCF-7_ChIA-PET2_result/replicate1/MCF-7.interactions.MICC"
MCF7_REP2 = "../MCF-7_ChIA-PET2_result/replicate2/MCF-7.interactions.MICC"



"""
if first line start with "chr", remove first line;
only keep the columns 0, 1, 2, 3, 4, 5, 10, 12
change the extention name from .MICC to .fdr.mango
"""

def parse_micc_to_mango(micc_path, mango_path):
    with open(micc_path, "r") as micc_file:
        with open(mango_path, "w") as mango_file:
            for line in micc_file:
                if line.startswith("chr\tstart\tend"):
                    continue
                else:
                    line = line.strip().split("\t")
                    if line[0] != line[3]: # only keep the interactions within the same chromosome
                        continue
                    if line[0] not in CHR_LIST or line[3] not in CHR_LIST:
                        continue
                    line = "\t".join([line[0], line[1], line[2], line[3], line[4], line[5], line[10], line[12]])
                    mango_file.write(line+"\n")

if __name__ == "__main__":
    parse_micc_to_mango(K562_REP1, "./K562_rep1.interactions.fdr.mango")
    parse_micc_to_mango(K562_REP2, "./K562_rep2.interactions.fdr.mango")
    parse_micc_to_mango(MCF7_REP1, "./MCF7_rep1.interactions.fdr.mango")
    parse_micc_to_mango(MCF7_REP2, "./MCF7_rep2.interactions.fdr.mango")


