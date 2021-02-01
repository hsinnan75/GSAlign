GSAlign: an ultra-fast sequence alignment algorithm for intra-species genome comparison
===================

Developers: Dr. Hsin-Nan Lin and Dr. Wen-Lian Hsu Institute of Information Science, Academia Sinica, Taiwan.

# Introduction

Personal genomics and comparative genomics are two fields that are more and more important in clinical practices and genome researches. Both fields require sequence alignment to discover sequence conservation and structural variation. Though many methods have been developed to handle genome sequence alignment, some are designed for small genome comparison while some are not efficient for large genome comparison. Here, we present GSAlign to handle large genome comparison efficiently. GSAlign includes three unique features: 1) it is the first attempt to use Burrows-Wheeler Transform on genome sequence alignment; 2) it supports parallel computing; 3) it adopts a divide-and-conquer strategy to separate a query sequence into regions that are easy to align and regions that require gapped alignment. With all these features, we demonstrated GSAlign is very efficient and sensitive in finding both the exact matches and differences between two genome sequences and it is much faster than existing state-of-the-art methods. 

# Download

## Conda
Install [Bioconda](https://bioconda.github.io/user/install.html) then type:
```
$ conda install gsalign
```

## Github
  ```
  $ git clone https://github.com/hsinnan75/GSAlign.git
  ```
to download the package of GSAlign.

# Compiling

To compile GSAlign and the index tool, please change to GSAlign's folder and just type 'make' to compile GSAlign and bwt_index. If the compilation or the programs fail, please contact me (arith@iis.sinica.edu.tw), Thanks.

You may run ./run_test.sh to test GSAlign with two E.coli strains.

# Instructions

To index a reference genome, GSAlign requires the target genome file (in fasta format) and the prefix of the index files (including the directory path).

  ```
  $ bin/bwt_index ref_file[ex.ecoli.fa] index_prefix[ex. Ecoli]
  ```
or

  ```
  $ bin/GSAlign index ref_file[ex.ecoli.fa] index_prefix[ex. Ecoli]
  ```

The above command is to index the genome file Ecoli.fa and store the index files begining with ecoli.
If the index files are not mdade beforehand, GSAlign will generate index files istself with the given reference genome sequences.

To align two genome sequences, GSAlign requires two genome files (in fasta format)

  ```
  $ bin/GSAlign -r fa1 -q fa2 -o output
  ```
or with a pre-built index file

  ```
  $ bin/GSAlign -i idx_prefix -q fa2 -o output
  ```

# Datasets

You may download the test datasets at http://bioapp.iis.sinica.edu.tw/~arith/GSAlign/ to test the performance of GSAlign.
You may also use the evaulation tool (Evaluation.cpp) to measure the precision and recall of the resulting VCF files. 
To compile Evaluation.cpp, just type 'g++ Evaluation.cpp -o eva'

# File formats

- Reference and query genome files

    Both the reference genome and query genome files should be in FASTA format.

- Output file

	1. maf/aln file: it shows the pairwise alignments between two sequences (MAF/ALN format).
	2. vcf file: it shows sequence variants between two sequences (VCF format).
	3. ps  file: it shows dotplot between two sequences (gnuplot is required).

# Parameter setting

 ```
-t INT number of threads [8]

-i STR index prefix [BWT based (BWA)]

-r STR reference genome filename [fasta]

-q STR query genome filename [fasta]

-o STR prefix of output files [output]

-sen Sensitive mode [False]

-dp Output Dot-plots [false]. It may not work on MacOS.

-gp STR specify the path of gnuplot, ex: -gp /usr/bin/gnuplot

-fmt INT Set the output format [1]: 1:maf, 2:aln

-idy set the minimal sequnce identity [70]

-one set one on one aligment mode [false] (This option forces the query
sequence is only allowed to be aligned at most one position).

-unique Only output unique alignment [false]

-slen set the minimal seed length [15]. To compare two sequences of less similarity, please use smaller size seed length.

-alen set the minimal alignment length [5000]

-ind set the maximal indel size [25]. The maximal sigle indel size allowed in an alignment.

-clr set the minimal cluster size [250]. A cluster is a group of seeds, and  cluster size is the total length of its seeds. This value can filter out random clusters.

-no_vcf No VCF output [false]. Do not identify sequence variations.

  ```
