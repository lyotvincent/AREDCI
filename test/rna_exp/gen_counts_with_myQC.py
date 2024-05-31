

import sys, os
ROOT_PATH = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(ROOT_PATH, "src"))
from preprocessing.quality_controller import quality_control
from preprocessing.preprocessing_utils import complement_zero_pet_loops_for_samples
from utils.utils import print_counts

GM12878_REP1 = ROOT_PATH+"/data/GM12878_ChIA-PET2_result/replicate1/GM12878.interactions.MICC"
GM12878_REP2 = ROOT_PATH+"/data/GM12878_ChIA-PET2_result/replicate2/GM12878.interactions.MICC"
K562_REP1 = ROOT_PATH+"/data/K562_ChIA-PET2_result/replicate1/K562.interactions.MICC"
K562_REP2 = ROOT_PATH+"/data/K562_ChIA-PET2_result/replicate2/K562.interactions.MICC"
MCF7_REP1 = ROOT_PATH+"/data/MCF-7_ChIA-PET2_result/replicate1/MCF-7.interactions.MICC"
MCF7_REP2 = ROOT_PATH+"/data/MCF-7_ChIA-PET2_result/replicate2/MCF-7.interactions.MICC"
files = [K562_REP1, K562_REP2, MCF7_REP1, MCF7_REP2]
samples = ["K562_rep1", "K562_rep2", "mcf7_rep1", "mcf7_rep2"]
groups = ["K562", "K562", "mcf7", "mcf7"]

# * default parameters
gap = 500
blacklist_path = ROOT_PATH+"/data/blacklist/hg38-blacklist.v2.bed"
remove_loops_in_blacklist=True
remove_self_ligated_loops=True
fdr_threshold=0.01
intra_only=True
chr_filter=["chrM", "chrX", "chrUn"]
pet_threshold=5
output_path = ROOT_PATH+"/data/my_test_data/"

# * Quality control
data = quality_control(files=files,
                samples=samples,
                groups=groups,
                blacklist_path=blacklist_path,
                gap=gap,
                remove_loops_in_blacklist=remove_loops_in_blacklist,
                remove_self_ligated_loops=remove_self_ligated_loops,
                fdr_threshold=fdr_threshold,
                intra_only=intra_only,
                chr_filter=chr_filter,
                pet_threshold=pet_threshold,
                output_path=output_path)
# data = complement_zero_pet_loops_for_samples(data, func_depth=0)
print_counts(data, os.path.join(ROOT_PATH, "test", "rna_exp", "counts_after_myQC.csv"))

