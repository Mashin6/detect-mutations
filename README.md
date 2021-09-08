# detect-mutations
Script for detecting mutations in aligned reads (.bam |.sam) e.g. from RNA base conversion  techniques such as TimeLapse, SLAM-seq, TUC-seq

Usage:<br>
`python detect-mutations.py -b <file.bam> [-s <snpFile.txt>] [--reads <PE|SE>] [--mutType <TC[,GA,..]> ][--minDist <int>] [--minQual <int>] [--tracks] [--mutPos]`

Requires: pysam

Options:
```
    -b|--bam    path to input .bam file
    -s|--SNP    path to file with SNPs positions
    --mutType   list of mutations types to analyze comma-separated (defalut: TC)
    --reads     type of reads paird-end (PE) or single end (SE) (default: PE)
    --minDist   smallest allowed distance of mutation for either 5' or 3' end of read
    --minQual   smallest allowed sequencing quality of mutation
    --tracks    make .bedGraph tracks of mutations
    --mutPos    make _cB.csv file with mutation statisctics for every nucleotide
```

Input: 
* .bam file
* list of SNPs for common mutations removal (optional) Format: chrom:position:REF:MUT

Output: 
* *_counts.csv     - mutation statistics for each read or read pair
    * qname  : read (pair) name
    * nA     : number of A in reference sequence covered by the read 
    * nC     : number of C
    * nT     : number of T
    * nG     : number of G,
    * rname  : chromosome
    * FR     : forward or reverse strand
    * sj     : read (pair) contains splice junction 
    * TA     : number of T -> A conversions recorded in read (pair)
    * CA
    * GA
    * NA
    * AT
    * CT
    * GT
    * NT
    * AC
    * TC
    * GC
    * NC
    * AG
    * TG
    * CG
    * NG
    * AN
    * TN
    * CN
    * GN
* *_cB.csv         - mutation statistics for each nucleotide position. Coverage + number of mutations observed
    * rname  : chromosome name
    * gloc   : genomic position of nucleotide
    * trials : number of reads that recorded given nucleotide
    * n      : number of mutations observed at nucleotide
*  *_muts.bedGraph  - browser tracks for viewing mutation counts
    * chromosome name
    * nucleotide start, end
    * count of mutations at position










