# Data records
## Homo sapiens
|cell line|Isogenic replicate|Library|read1 Accession|read2 Accession|Lab|Date added|used in CHASOS|
|--------|--------|--------|--------|--------|--------|--------|--------|
|GM12878|replicate1|ENCLB716IME|ENCFF409UJV|ENCFF810FRZ|Yijun Ruan, JAX|2019-10-02|✓|
|GM12878|replicate2|ENCLB203AEA|ENCFF048QWJ|ENCFF062ISF|Charles Lee|2021-08-25|×|
|K562|replicate1|ENCLB728NYL|ENCFF485GHA|ENCFF116MYA|Yijun Ruan, JAX|2019-09-06|×|
|K562|replicate2|ENCLB587PEE|ENCFF915YCD|ENCFF245TLN|Yijun Ruan, JAX|2019-09-06|✓|
|HCT116|replicate1|ENCLB610CWY|ENCFF695OOH|ENCFF424OGG|Charles Lee|2021-03-18|×|
|HCT116|replicate2|ENCLB557QDS|ENCFF869HLZ|ENCFF933EOS|Charles Lee|2021-03-18|✓|
|MCF-7|replicate1|ENCLB221UOT|ENCFF391NJE|ENCFF465LQF|Yijun Ruan, JAX|2019-11-01|×|
|MCF-7|replicate2|ENCLB978GEE|ENCFF592KNN|ENCFF342KQM|Yijun Ruan, JAX|2019-11-01|✓|

## Mus musculus
|cell line|Isogenic replicate|Library|read1 Accession|read2 Accession|Lab|Date added|used in CHASOS|
|--------|--------|--------|--------|--------|--------|--------|--------|
|Mus musculus strain B6CASTF1/J embryonic stem cell male embryo|replicate1|ENCLB848JLS|ENCFF533ZMZ|ENCFF776ZNO|Charles Lee|2021-08-13|×|
|Mus musculus strain B6CASTF1/J embryonic stem cell male embryo|replicate2|ENCLB844IKU|ENCFF399JRD|ENCFF538TLL|Charles Lee|2021-08-07|×|
|Mus musculus strain C57BL/6J spleen tissue male embryo|replicate1|ENCLB634TDO|ENCFF966NIS|ENCFF273OFU|
Charles Lee|2021-08-15|×|

# ChIA-PET2 command
> conda activate py2

## Homo sapiens
> nohup ChIA-PET2 -g ./hg38/hg38.fa -b ./hg38/hg38.chrom.sizes -f ./raw_data/GM12878/ENCFF409UJV.fastq.gz -r ./raw_data/GM12878/ENCFF810FRZ.fastq.gz -A ACGCGATATCTTATCTGACT -B AGTCAGATAAGATATCGCGT -o GM12878_rep1_ChIA-PET2_result -n GM12878 -m 1 -t 40 -C 0 > nohup.out 2>&1 &

> nohup ChIA-PET2 -g ./hg38/hg38.fa -b ./hg38/hg38.chrom.sizes -f ./raw_data/GM12878/ENCFF048QWJ.fastq.gz -r ./raw_data/GM12878/ENCFF062ISF.fastq.gz -A ACGCGATATCTTATCTGACT -B AGTCAGATAAGATATCGCGT -o GM12878_rep2_ChIA-PET2_result -n GM12878 -m 1 -t 40 -C 0 > nohup.out 2>&1 &

> nohup ChIA-PET2 -g ./hg38/hg38.fa -b ./hg38/hg38.chrom.sizes -f ./raw_data/K562/ENCFF485GHA.fastq.gz -r ./raw_data/K562/ENCFF116MYA.fastq.gz -A ACGCGATATCTTATCTGACT -B AGTCAGATAAGATATCGCGT -o K562_rep1_ChIA-PET2_result -n K562 -m 1 -t 40 -C 0 > nohup.out 2>&1 &

> nohup ChIA-PET2 -g ./hg38/hg38.fa -b ./hg38/hg38.chrom.sizes -f ./raw_data/K562/ENCFF915YCD.fastq.gz -r ./raw_data/K562/ENCFF245TLN.fastq.gz -A ACGCGATATCTTATCTGACT -B AGTCAGATAAGATATCGCGT -o K562_rep2_ChIA-PET2_result -n K562 -m 1 -t 40 -C 0 > nohup.out 2>&1 &

> nohup ChIA-PET2 -g ./hg38/hg38.fa -b ./hg38/hg38.chrom.sizes -f ./raw_data/HCT116/ENCFF695OOH.fastq.gz -r ./raw_data/HCT116/ENCFF424OGG.fastq.gz -A ACGCGATATCTTATCTGACT -B AGTCAGATAAGATATCGCGT -o HCT116_rep1_ChIA-PET2_result -n HCT116 -m 1 -t 40 -C 0 > nohup.out 2>&1 &

> nohup ChIA-PET2 -g ./hg38/hg38.fa -b ./hg38/hg38.chrom.sizes -f ./raw_data/HCT116/ENCFF869HLZ.fastq.gz -r ./raw_data/HCT116/ENCFF933EOS.fastq.gz -A ACGCGATATCTTATCTGACT -B AGTCAGATAAGATATCGCGT -o HCT116_rep2_ChIA-PET2_result -n HCT116 -m 1 -t 40 -C 0 > nohup.out 2>&1 &

> nohup ChIA-PET2 -g ./hg38/hg38.fa -b ./hg38/hg38.chrom.sizes -f ./raw_data/MCF-7/ENCFF391NJE.fastq.gz -r ./raw_data/MCF-7/ENCFF465LQF.fastq.gz -A ACGCGATATCTTATCTGACT -B AGTCAGATAAGATATCGCGT -o MCF-7_rep1_ChIA-PET2_result -n MCF-7 -m 1 -t 40 -C 0 > nohup.out 2>&1 &

> nohup ChIA-PET2 -g ./hg38/hg38.fa -b ./hg38/hg38.chrom.sizes -f ./raw_data/MCF-7/ENCFF592KNN.fastq.gz -r ./raw_data/MCF-7/ENCFF342KQM.fastq.gz -A ACGCGATATCTTATCTGACT -B AGTCAGATAAGATATCGCGT -o MCF-7_rep2_ChIA-PET2_result -n MCF-7 -m 1 -t 40 -C 0 > nohup.out 2>&1 &

## Mus musculus

> nohup ChIA-PET2 -g ./mm39/mm39.fa.gz -b ./mm39/mm39.chrom.sizes -f ./raw_data/ESC/ENCFF533ZMZ.fastq.gz -r ./raw_data/ESC/ENCFF776ZNO.fastq.gz -A ACGCGATATCTTATCTGACT -B AGTCAGATAAGATATCGCGT -o ESC_rep1_ChIA-PET2_result -n ESC -m 1 -t 40 -C 0 > nohup.out 2>&1 &

> nohup ChIA-PET2 -g ./mm39/mm39.fa.gz -b ./mm39/mm39.chrom.sizes -f ./raw_data/ESC/ENCFF399JRD.fastq.gz -r ./raw_data/ESC/ENCFF538TLL.fastq.gz -A ACGCGATATCTTATCTGACT -B AGTCAGATAAGATATCGCGT -o ESC_rep2_ChIA-PET2_result -n ESC -m 1 -t 40 -C 0 > nohup.out 2>&1 &

> nohup ChIA-PET2 -g ./mm39/mm39.fa.gz -b ./mm39/mm39.chrom.sizes -f ./raw_data/spleen/ENCFF966NIS.fastq.gz -r ./raw_data/spleen/ENCFF273OFU.fastq.gz -A ACGCGATATCTTATCTGACT -B AGTCAGATAAGATATCGCGT -o spleen_rep1_ChIA-PET2_result -n spleen -m 1 -t 40 -C 0 > nohup.out 2>&1 &


<style>
  tr:nth-child(4n), tr:nth-child(4n-1) {
    background-color: #f2f2f2;
  }
</style>