
import os

CUR_PATH = os.path.dirname(os.path.abspath(__file__))

with open(os.path.join(CUR_PATH, "sim_counts.outlier_0.0.csv")) as f:
    lines0 = f.readlines()
with open(os.path.join(CUR_PATH, "sim_counts.outlier_0.1.csv")) as f:
    lines1 = f.readlines()
with open(os.path.join(CUR_PATH, "sim_counts.outlier_0.2.csv")) as f:
    lines2 = f.readlines()
with open(os.path.join(CUR_PATH, "sim_counts.outlier_0.3.csv")) as f:
    lines3 = f.readlines()
with open(os.path.join(CUR_PATH, "sim_counts.outlier_0.4.csv")) as f:
    lines4 = f.readlines()
with open(os.path.join(CUR_PATH, "sim_counts.outlier_0.5.csv")) as f:
    lines5 = f.readlines()

outlier_num1 = 0
outlier_num2 = 0
outlier_num3 = 0
outlier_num4 = 0
outlier_num5 = 0
for line0, line1 in zip(lines0, lines1):
    if line0 != line1:
        outlier_num1 += 1
print(outlier_num1, outlier_num1/len(lines0))

for line0, line2 in zip(lines0, lines2):
    if line0 != line2:
        outlier_num2 += 1
print(outlier_num2, outlier_num2/len(lines0))

for line0, line3 in zip(lines0, lines3):
    if line0 != line3:
        outlier_num3 += 1
print(outlier_num3, outlier_num3/len(lines0))

for line0, line4 in zip(lines0, lines4):
    if line0 != line4:
        outlier_num4 += 1
print(outlier_num4, outlier_num4/len(lines0))

for line0, line5 in zip(lines0, lines5):
    if line0 != line5:
        outlier_num5 += 1
print(outlier_num5, outlier_num5/len(lines0))

