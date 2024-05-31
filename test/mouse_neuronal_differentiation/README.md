# data
【1】 mouse embryonic stem cell E14, Neural stem cell NS5 and neuroshpere cells GSE44067 GEO Accession viewer (nih.gov)
在进化生物学中，细胞系的分支结构通常指的是它们在发育过程中的谱系关系，即它们是如何从原始的干细胞分化而来的。对于您提到的三种细胞系——小鼠胚胎干细胞（E14）、神经干细胞（NS5）和神经球细胞（neurosphere cells），我们可以从它们的起源和分化潜能来理解它们在进化上的分支结构。  
小鼠胚胎干细胞（E14）Embryonic Stem Cells（ESC, 胚胎干细胞）：这些细胞来源于小鼠胚胎的内细胞团（inner cell mass），是多能干细胞（pluripotent stem cells），具有分化成体内所有细胞类型的潜力。在进化上，这些细胞代表了早期胚胎发育阶段的细胞，它们能够通过分化过程产生各种类型的细胞，包括神经干细胞。  
神经干细胞（NS5）Neuronal Stem Cells（NSC, 神经干细胞）：神经干细胞是从胚胎干细胞（如E14）中分化出来的，它们是多能或多潜能干细胞，专门负责生成神经系统的细胞。在进化上，神经干细胞代表了从多能干细胞向特定组织（如神经系统）分化的分支。NS5可能是指特定的神经干细胞系，但具体的进化分支结构需要更详细的信息来确定。  
神经球细胞（neurosphere cells, NPC）：神经球细胞是从神经干细胞中进一步分化而来的，它们通常在体外培养中形成，能够自我更新并分化成多种类型的神经细胞，包括神经元和神经胶质细胞。在进化上，神经球细胞代表了神经干细胞分化为成熟神经细胞的进一步分支。  
这三种细胞系在进化上的分支结构可以概括为：从多能干细胞（E14）分化为神经干细胞（NS5），然后神经干细胞进一步分化为神经球细胞，最终形成成熟的神经细胞。  
| 细胞系 | acc.ver |
| --- | --- |
| ESC_rep1 | SRR768392 |
| ESC_rep2 | SRR768393 |
| ESC_rep2 | SRR768394 |
| NSC_rep1 | SRR768395 |
| NSC_rep1 | SRR768396 |
| NSC_rep2 | SRR768397 |
| NSC_rep2 | SRR768398 |
| NPC_rep1 | SRR768399 |
| NPC_rep1 | SRR768400 |
【2】mouse stem cells, activated B cells, CH12 cell line, and plasmacytomas GSE98119 GEO Accession viewer (nih.gov)
在进化和发育生物学的背景下，细胞系的分支结构通常指的是它们在生物体发育过程中的谱系发展。对于您提到的细胞类型，我们可以从它们的起源和功能来理解它们在进化上的分支结构：  
小鼠干细胞（mESCs）：这些细胞来源于小鼠胚胎的内细胞团（inner cell mass），是多能干细胞，具有分化成体内所有细胞类型的潜力。在进化上，mESCs代表了早期胚胎发育阶段的细胞，它们能够通过分化过程产生各种类型的细胞，包括神经干细胞、血液细胞等。  
激活的B细胞（Activated B cells）：激活的B细胞是免疫系统中的一种细胞，它们在遇到抗原后被激活，能够增殖并分化成产生抗体的浆细胞（plasma cells）。在进化上，B细胞是从造血干细胞（hematopoietic stem cells）分化而来的，它们在适应性免疫系统中起着关键作用。  
CH12细胞系：CH12细胞系是一种特定的B细胞系，通常用于研究B细胞的激活和抗体产生。这些细胞在实验室条件下被用来模拟B细胞的激活过程，但它们在进化上的分支结构与自然状态下的B细胞相同，即从造血干细胞分化而来。  
浆细胞（Plasmacytomas）：浆细胞是高度分化的B细胞，专门负责产生大量的抗体。它们是B细胞在遇到抗原后经过一系列激活和分化过程的最终产物。在进化上，浆细胞代表了B细胞分化的末端阶段，是适应性免疫反应的重要组成部分。

# preprocessing
```python complement_info.py``` transform the txt file to MICC format. generate NPC_rep.csv ESC_rep.csv NSC_rep.csv NPC_rep.csv

# DCI result
run by main.py  

*kde:*
ESC vs NSC: fdrs < 0.05: 2562
ESC vs NPC: fdrs < 0.05: 2708
NSC vs NPC: fdrs < 0.05: 970

*neighbor:*
ESC vs NSC: fdrs < 0.05: 192
ESC vs NPC: fdrs < 0.05: 231
NSC vs NPC: fdrs < 0.05: 230

from https://geneontology.org/ get GO_bp_analysis.txt, GO_mf_analysis.txt

```python process_result.py```提取DCI结果中的基因列表，分别进行GO生物过程和分子功能的富集分析，得到GO_bp_analysis.txt, GO_mf_analysis.txt

plot_venn.py 画venn图的列表，在http://www.ehbio.com/test/venn/#/ 画

statistics.py 统计venn图的各个集合的GO
