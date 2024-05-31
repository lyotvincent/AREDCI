
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
    bp = []
    mf = []
    with open(os.path.join(file_path, GO_FILE[0]), 'r') as f:
        bp_lines = f.readlines()
        bp_lines = bp_lines[12:]
    for l in bp_lines:
        go_term = l.split('\t')[0].split()[-1].strip("()")
        if go_term == "UNCLASSIFIED":
            continue
        bp.append(go_term)
    # print(bp)
    with open(os.path.join(file_path, GO_FILE[1]), 'r') as f:
        mf_lines = f.readlines()
        mf_lines = mf_lines[12:]
    for l in mf_lines:
        go_term = l.split('\t')[0].split()[-1].strip("()")
        if go_term == "UNCLASSIFIED":
            continue
        mf.append(go_term)
    print(file_path, len(bp), len(mf))
    return bp, mf

bp_esc_nsc, mf_esc_nsc = load_bp_and_md(ESC_VS_NSC_PATH)
bp_esc_npc, mf_esc_npc = load_bp_and_md(ESC_VS_NPC_PATH)
bp_nsc_npc, mf_nsc_npc = load_bp_and_md(NSC_VS_NPC_PATH)
print(len(bp_esc_nsc), len(bp_esc_npc), len(bp_nsc_npc))
print(len(set(bp_esc_nsc)), len(set(bp_esc_npc)), len(set(bp_nsc_npc)))
print(len(set(bp_esc_nsc) & set(bp_esc_npc)))
print(len(set(bp_esc_nsc) & set(bp_nsc_npc)))
print(len(set(bp_esc_npc) & set(bp_nsc_npc)))
print(len(set(bp_esc_nsc) & set(bp_esc_npc) & set(bp_nsc_npc)))

print(len(mf_esc_nsc), len(mf_esc_npc), len(mf_nsc_npc))
print(len(set(mf_esc_nsc)), len(set(mf_esc_npc)), len(set(mf_nsc_npc)))
print(len(set(mf_esc_nsc) & set(mf_esc_npc)))
print(len(set(mf_esc_nsc) & set(mf_nsc_npc)))
print(len(set(mf_esc_npc) & set(mf_nsc_npc)))
print(len(set(mf_esc_nsc) & set(mf_esc_npc) & set(mf_nsc_npc)))


venn_list = list()
for i in list(set(bp_esc_nsc)):
    venn_list.append((i, "ESC_vs_NSC"))
for i in list(set(bp_esc_npc)):
    venn_list.append((i, "ESC_vs_NPC"))
for i in list(set(bp_nsc_npc)):
    venn_list.append((i, "NSC_vs_NPC"))

df1 = pd.DataFrame(venn_list, columns=["GO", "Comparison"])
df1.to_excel(os.path.join(N_DIFF_PATH, "venn_bp_list.xlsx"), index=False)

venn_list = list()
for i in list(set(mf_esc_nsc)):
    venn_list.append((i, "ESC_vs_NSC"))
for i in list(set(mf_esc_npc)):
    venn_list.append((i, "ESC_vs_NPC"))
for i in list(set(mf_nsc_npc)):
    venn_list.append((i, "NSC_vs_NPC"))

df2 = pd.DataFrame(venn_list, columns=["GO", "Comparison"])
df2.to_excel(os.path.join(N_DIFF_PATH, "venn_mf_list.xlsx"), index=False)

df = pd.concat([df1, df2])
# add abundance column, all abundance=1
abundance = [1]*df.shape[0]
df["Abundance"] = abundance
# add color column
# if a GO present three times, it means it is present in all three Comparison/groups(ESC_vs_NSC, ESC_vs_NPC, NSC_vs_NPC)
# set these GO terms to be #000000 color.
# set the GO terms only present in ESC_vs_NPC to be red #ED1C24 color.
# set the GO terms only present in ESC_vs_NSC to be blue #00A2E8 color.
# set the GO terms only present in NSC_vs_NPC to be green #22B14C color.
# set the GO terms present in ESC_vs_NSC and ESC_vs_NPC to be purple #A349A4 color.
# set the GO terms present in ESC_vs_NSC and NSC_vs_NPC to be 青色 #99D9EA.
# set the GO terms present in ESC_vs_NPC and NSC_vs_NPC to be yellow #FFC90E color.
# set the GO terms present in ESC_vs_NSC, ESC_vs_NPC and NSC_vs_NPC to be black #000000 color.
color = []
for i in df["GO"]:
    if df[df["GO"] == i].shape[0] == 3:
        color.append("#000000")
    elif df[df["GO"] == i].shape[0] == 2:
        if "ESC_vs_NSC" in list(df[df["GO"] == i]["Comparison"]) and "ESC_vs_NPC" in list(df[df["GO"] == i]["Comparison"]):
            color.append("#A349A4")
        elif "ESC_vs_NSC" in list(df[df["GO"] == i]["Comparison"]) and "NSC_vs_NPC" in list(df[df["GO"] == i]["Comparison"]):
            color.append("#99D9EA")
        elif "ESC_vs_NPC" in list(df[df["GO"] == i]["Comparison"]) and "NSC_vs_NPC" in list(df[df["GO"] == i]["Comparison"]):
            color.append("#FFC90E")
    elif df[df["GO"] == i].shape[0] == 1:
        if "ESC_vs_NSC" in list(df[df["GO"] == i]["Comparison"]):
            color.append("#00A2E8")
        elif "ESC_vs_NPC" in list(df[df["GO"] == i]["Comparison"]):
            color.append("#ED1C24")
        elif "NSC_vs_NPC" in list(df[df["GO"] == i]["Comparison"]):
            color.append("#22B14C")
df["Color"] = color
# 把Comparison这列换到第三列
df = df[["GO", "Abundance", "Color", "Comparison"]]
# add Condition column
condition = ["ellipse"]*df.shape[0]
df["Condition"] = condition
# set header = "node_name, Abundance, Color, Group, Condition"
df.to_csv(os.path.join(N_DIFF_PATH, "venn_list.csv"), index=False, header=["node_name", "Abundance", "Color", "Group", "Condition"], sep="\t")

