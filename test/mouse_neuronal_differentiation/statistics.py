
import os
import pandas as pd

ROOT_PATH = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# print(ROOT_PATH)
N_DIFF_PATH = os.path.join(ROOT_PATH, 'test', 'mouse_neuronal_differentiation')

ESC_VS_NSC_PATH = os.path.join(N_DIFF_PATH, 'ESC_vs_NSC_result')
ESC_VS_NPC_PATH = os.path.join(N_DIFF_PATH, 'ESC_vs_NPC_result')
NSC_VS_NPC_PATH = os.path.join(N_DIFF_PATH, 'NSC_vs_NPC_result')

VS_PATH = [ESC_VS_NSC_PATH, ESC_VS_NPC_PATH, NSC_VS_NPC_PATH]
GO_FILE = ["GO_bp_analysis.txt", "GO_mf_analysis.txt"]

def load_bp_and_md(file_path):
    go_dict = dict()
    with open(os.path.join(file_path, GO_FILE[0]), 'r') as f:
        bp_lines = f.readlines()
        bp_lines = bp_lines[12:]
    for l in bp_lines:
        go_term = l.split('\t')[0].strip(")").split(" (")
        term, go_id = go_term[0], go_term[1]
        if go_id == "UNCLASSIFIED":
            continue
        go_dict[go_id] = term
    with open(os.path.join(file_path, GO_FILE[1]), 'r') as f:
        mf_lines = f.readlines()
        mf_lines = mf_lines[12:]
    for l in mf_lines:
        go_term = l.split('\t')[0].strip(")").split(" (")
        term, go_id = go_term[0], go_term[1]
        if go_id == "UNCLASSIFIED":
            continue
        go_dict[go_id] = term
    print(file_path, len(go_dict))
    return go_dict

esc_npc_go_dict = load_bp_and_md(ESC_VS_NPC_PATH)
esc_nsc_go_dict = load_bp_and_md(ESC_VS_NSC_PATH)
nsc_npc_go_dict = load_bp_and_md(NSC_VS_NPC_PATH)

all_go_dict = dict()
all_go_dict.update(esc_npc_go_dict)
all_go_dict.update(esc_nsc_go_dict)
all_go_dict.update(nsc_npc_go_dict)

esc_npc_go_id = set(esc_npc_go_dict.keys())
esc_nsc_go_id = set(esc_nsc_go_dict.keys())
nsc_npc_go_id = set(nsc_npc_go_dict.keys())
print(type(esc_npc_go_id))

go_id_3 = esc_nsc_go_id.intersection(esc_npc_go_id).intersection(nsc_npc_go_id)
print(len(go_id_3))
with open(os.path.join(N_DIFF_PATH, 'intersection_of_the_three.csv'), "w") as f:
    for go_id in go_id_3:
        f.write(go_id+"\t"+ all_go_dict[go_id] + "\n")

go_id_esc_nsc_nsc_npc = esc_nsc_go_id.intersection(nsc_npc_go_id)
print(len(go_id_esc_nsc_nsc_npc))
go_in_esc_npc_nsc_npc = esc_npc_go_id.intersection(nsc_npc_go_id)
print(len(go_in_esc_npc_nsc_npc))
go_in_esc_nsc_esc_npc = esc_nsc_go_id.intersection(esc_npc_go_id)
print(len(go_in_esc_nsc_esc_npc))

# c=go_id_esc_nsc_nsc_npc-go_id_3
# print(c)
c = go_in_esc_nsc_esc_npc - go_id_3
with open(os.path.join(N_DIFF_PATH, 'intersection_of_esc_nsc_esc_npc.csv'), "w") as f:
    for go_id in c:
        f.write(go_id+"\t"+ all_go_dict[go_id] + "\n")

c = esc_npc_go_id - esc_nsc_go_id - nsc_npc_go_id
with open(os.path.join(N_DIFF_PATH, 'intersection_of_esc_npc_only.csv'), "w") as f:
    for go_id in c:
        f.write(go_id+"\t"+ all_go_dict[go_id] + "\n")

c = esc_nsc_go_id - esc_npc_go_id - nsc_npc_go_id
with open(os.path.join(N_DIFF_PATH, 'intersection_of_esc_nsc_only.csv'), "w") as f:
    for go_id in c:
        f.write(go_id+"\t"+ all_go_dict[go_id] + "\n")

c = nsc_npc_go_id - esc_nsc_go_id - esc_npc_go_id
with open(os.path.join(N_DIFF_PATH, 'intersection_of_nsc_npc_only.csv'), "w") as f:
    for go_id in c:
        f.write(go_id+"\t"+ all_go_dict[go_id] + "\n")

