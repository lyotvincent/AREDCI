

"""
transform the txt file to MICC format.
complement info for :
GSM1084137_mouse_ESC_RNAPII_ChiaPet_Technical_rep1_mm9.cluster.pet2+.txt -> ESC_rep1.csv
GSM1084138_mouse_ESC_RNAPII_ChiaPet_Technical_rep2_mm9.cluster.pet2+.txt -> ESC_rep2.csv
GSM1084139_mouse_NSC_RNAPII_ChiaPet_Biological_rep1_mm9.cluster.pet2+.txt -> NSC_rep1.csv
GSM1084140_mouse_NSC_RNAPII_ChiaPet_Biological_rep2_mm9.cluster.pet2+.txt -> NSC_rep2.csv
GSM1084141_mouse_NPC_RNAPII_ChiaPet.cluster.pet2+.txt -> NPC_rep1.csv

columns in txt:
chr1    start1  end1    chr2    start2  end2    #PET    PP  FDR

we need:
chr1    start1  end1    chr2    start2  end2    peak1   peak2   depth1  depth2  #PET    PP  FDR

peak1 and peak2 are the index of the (chr, start, end) in the file
depth1 and depth2 are the depth of the peak1 and peak2, unused, just set to 0
"""

import os

def complement_info(input_path, output_path):
    with open(input_path, "r") as f:
        lines = f.readlines()[1:]
    peak_dict = dict()
    for line in lines:
        line = line.strip().split()
        peak1 = (line[0], line[1], line[2])
        peak2 = (line[3], line[4], line[5])
        if peak1 not in peak_dict:
            peak_dict[peak1] = len(peak_dict)+1
        if peak2 not in peak_dict:
            peak_dict[peak2] = len(peak_dict)+1
    with open(output_path, "w") as f:
        for line in lines:
            line = line.strip().split()
            peak1 = (line[0], line[1], line[2])
            peak2 = (line[3], line[4], line[5])
            f.write(f"{line[0]}\t{line[1]}\t{line[2]}\t{line[3]}\t{line[4]}\t{line[5]}\tpeak_{peak_dict[peak1]}\tpeak_{peak_dict[peak2]}\t0\t0\t{line[6]}\t{line[7]}\t{line[8]}\n")

if __name__ == "__main__":
    ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
    complement_info(os.path.join(ROOT_DIR, "GSM1084137_mouse_ESC_RNAPII_ChiaPet_Technical_rep1_mm9.cluster.pet2+.txt"), os.path.join(ROOT_DIR, "ESC_rep1.csv"))
    complement_info(os.path.join(ROOT_DIR, "GSM1084138_mouse_ESC_RNAPII_ChiaPet_Technical_rep2_mm9.cluster.pet2+.txt"), os.path.join(ROOT_DIR, "ESC_rep2.csv"))
    complement_info(os.path.join(ROOT_DIR, "GSM1084139_mouse_NSC_RNAPII_ChiaPet_Biological_rep1_mm9.cluster.pet2+.txt"), os.path.join(ROOT_DIR, "NSC_rep1.csv"))
    complement_info(os.path.join(ROOT_DIR, "GSM1084140_mouse_NSC_RNAPII_ChiaPet_Biological_rep2_mm9.cluster.pet2+.txt"), os.path.join(ROOT_DIR, "NSC_rep2.csv"))
    complement_info(os.path.join(ROOT_DIR, "GSM1084141_mouse_NPC_RNAPII_ChiaPet.cluster.pet2+.txt"), os.path.join(ROOT_DIR, "NPC_rep1.csv"))

