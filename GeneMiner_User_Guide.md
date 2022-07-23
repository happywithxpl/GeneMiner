<h1 align="center">GeneMiner Manual</h1>



[TOC]

# 1. Overview

​		GeneMiner is a software for extracting phylogenetic markers from next generation sequencing (NGS) data.With GeneMiner, users can accurately and efficiently obtain a large number of target sequences from a wide range of histological data at a very low cost, such as: all or part of the chloroplast/mitochondrial genome, highly repetitive regions in the nuclear genome (e.g. nrDNA), etc. from shallow whole-genome sequencing data; single-low copy genes from transcriptome sequencing data, etc. GeneMiner improves existing phylogenetic research strategies at the most basic data level, and has significant advantages in reducing experimental costs and expanding phylogenetic marker selection. In addition, GeneMiner can be applied to the extraction of DNA barcodes, customs inspection and quarantine, mining of specific functional genes and other research, which has broad application prospects.

# 2. Download and install

​		GeneMiner is open source under the GPL-3.0 license, which is distributed through the github repository (https://github.com/happywithxpl/GeneMiner/releases). Please be sure to follow our github page to stay up-to-date with the latest code changes. We do not provide any support for previous versions of the code! Version numbers follow the notation x.y.z, where x changes with major code reorganizations, y changes when new features are added, and z changes with bug fixes.

​		GeneMiner is an easy-to-use software written in python, which is is provided for x86-64 systems running GNU/Linux, macOS (version 10.13 or higher) and Windows (64-bit, version 7 or higher).

​		Users on Windows, macOS and Linux computers can run GeneMiner directly from the command line. we also offer a more convenient GUI version for Windows and macOS users.

## 2.1 GeneMiner with GUI

​		We strongly recommend using the GUI version for users who are not familiar with the command line or light use.Download the corresponding version of the packaged GUI from [here](:https://github.com/happywithxpl/GeneMiner/releases) and double click to run it.

![](picture/GeneMiner_GUI-16580658835161.png)

## **2.2 GeneMiner with command line**



### **2.2.1 Cloning the repo  (support)**

Instead of downloading the source distribution as a compressed archive, you could clone the repo and build it as shown below.

```shell
git https://github.com/happywithxpl/GeneMiner.git
cd GeneMiner
python setup.py install --record logName --user
geneminer.py -h
```



### **2.2.2 Source distribution** 

Download the source distribution from a [release](https://github.com/bpp/bpp/releases)  and  install dependencies, use the following commands:

```shell
wget https://github.com/happywithxpl/Geneminer/geneminer-1.0.0-linux-x86_64.tar.gz
tar geneminer_v1.0.0.tar.gz
cd  GeneMiner_v1.0.0
python setup.py install --record logName --user
geneminer.py -h
```



### **2.2.3 Flexible construction**

If both of the above methods fail or you want to have a deeper control of GeneMiner, you can use a more flexible method

- Download the GeneMiner's distribution from  [here](:https://github.com/happywithxpl/GeneMiner/releases).

```shell
wget https://github.com/happywithxpl/GeneMiner/geneminer-1.0.0-linux-x86_64.tar.gz `
tar -zvxf geneminer-1.0.0-linux-x86_64.tar.gz
```

- Use the following commands to make `geneminer.py`  executable.

```shell
cd geneminer-1.0.0
chmod 755 geneminer.py
```

- Install python library `biopython` using pip or conda.

```shell
# install required libs
pip install biopython --user
```

- Add `geneminer.py` to the `$PATH`.    

```shell
echo "export PATH=\$PATH:$(pwd)" >> ~/.bashrc
source ~/.bashrc
geneminer.py -h
```



# 3. Quick start

​		在使用GeneMiner挖掘特定目标序列之前，我们强烈建议您先了解自己的测序数据(来自Illumina、Roche-454、ABI或其他测序平台）的状况，包括测序方式、深度、质量、数据量大小等等，这有助于您对参数的合理选择。

​		Before using GeneMiner to mine specific target sequences, we strongly recommend you to know the status of your sequencing data (from Illumina, Roche-454, ABI or other sequencing platforms), including sequencing method, depth, quality, data volume size, etc., which will help you to make a reasonable choice of parameters.

​		Without any options, GeneMiner takes a reference database and a query sequence file as input and produce phylogenetic markers.We have prepared a simulated dataset of `Arabidopsis thaliana` to help you quickly use GeneMiner.you can download them from https://github.com/happywithxpl/GeneMiner-Test

（1）Mining single target sequence

```shell
 geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta -o out
```

`-1`,`-2`:  the input files with paired-end reads, given in FASTQ format.
`-rtfa`:  reference sequences in fasta format. 
`-o`:        the output directory

注意：

目标序列来源于近缘物种（同一个属或同一个科）的同源基因。

如果一个目标序列有来源于一个或多个物种的同源基因，请把这些同源基因存放在一个fasta格式的文件里。如下：

**NOTE：**Target sequences are derived from homologous genes of closely related species (same genus or family). If a target sequence has homologous genes from one or more species, please store these homologous genes in a fasta format file as a reference file .The file may look like this:

```shell
>species1_GeneA
ATCGATCG
>species2_GeneA
ATCGATCC
>species3_GeneA
ATTGATCC
>species4_GeneA
ATTGATCT
```



(2) Mining multiple target sequences

GeneMiner允许将多个参考文件存放在同一个文件夹以方便批量挖掘目标序列

GeneMiner allows multiple  reference files to be stored in the same folder to facilitate batch mining of target sequences

```shell
#for fasta format
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa fasta_folder -o out
#for GenBank format
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtgb GenBank_folder -o out
```

`-rtfa`:  reference sequences in fasta format

`-rtgb`:   reference sequences in GenBank format



(3) Accuracy Assessment

```shell
geneminer.py -1 data1.fq.gz -2 data2.fq.gz  -rtfa  ref.fasta -bn 100 -o out
```

`-bn` Number of resampling based on base substitution model。

基于核苷酸替换模型，GeneMiner对生成的结果序列反复重采样，从而在不依赖参考序列的情况下对结果进行准确性评估

Based on the base substitution model, GeneMiner iteratively resamples the generated result sequences to evaluate the accuracy of the results without relying on the reference sequences



# 4. Detailed User Guide

## 4.1 Input

### 4.1.1 Reads Files

​		GeneMiner assembles any type of next generation sequencing (NGS) reads, given in the FASTQ format. Giving paired or unpaired reads as input is OK, and keep in mind that GeneMiner will use pairing information. GeneMiner can direclty read compressed files in  `gzip` . Compressed files should end with  `.fq.gz` or .`fastq.gz` .

```
1	@ST-E00600:235:HFCTTALXX:1:1101:4752:1016 1:N:0:GGAGAACA
2	NCCCTATCATTTCTGAGGGGTTACATTCTCATTCTCTCAACATAAATACAGAGACCACTTCACCAAACACACCAACTCTGTCTCTGGATGTTGATATGATAAACACATCAATTCCTGACACTCC
	ATCTTAGTTCTTGGAGAGGCACCACT
3	+
4   #-FFFJJJJJJJJJJJ7JJ7J7JJJJ-J--JJ-JJJ-JJ7JJJJJJJJ-JJJ7J--J-J-JJJ-JJJ--JJJ--J--JJ-JJJJ7--J7--J7JJJJJJ77JJJJ7J-7JJ---J-JJ--JJ7-     J7JJ---JJ-JJ-JJJJJJJ-F>F>F
```



## 4.1.2 Reference files

​		DNA sequences are stored in reference files in FASTA format or GenBank format. Each reference file contains one or more DNA sequence(s),which is the homologous genes from closely related species (same genus or family). In addition, GeneMiner allows multiple reference files to be stored in the same folder to facilitate batch mining of target sequences.

(1)For a reference in FASTA format, it may look like this:

```shell
>Aegopodium_podagraria
TCGAATCCTGTGATAGCAGAACGACCCGCTAACTGGTAAATATATTGGGCAAGCTCATGGGGATTTTATCCCCTGTTGGTGAACCCTTGGTAGGTGGTCACTCCCCGGTTGCCACTGGCC
>Anethum_foeniculum
TCCTCCTTATCATTAGAGGAAGGAGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGG
>Bupleurum_chinense
AACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTCGAATCCTGAATCGAAGAGCGACCCGAGAACATGTTTTAAGACGGGGCCAGCGGTCGTCGGCCTCGGCCTGTCGGCTGCG
>Chamaesium_paradoxum
AAGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTCGATGCCTGCAACAGCAGAATGACCCGTGAACACGTATAAAACATTGGGCTAGTAGATGGGGCGCAAGTTCCCGGA
CATGAACCCCAGGACGGATG
```

**NOTE:** The reference's `filename` will be used as the name of the target sequence in GeneMiner.



(2) For a reference in GenBank format, it may look like this:

```shell
LOCUS       Daucus_carota         155895 bp    DNA     circular UNA 28-FEB-2021
DEFINITION  Daucus carota chloroplast, complete genome.
ACCESSION   urn.local...2i-d2j4m7m
VERSION     urn.local...2i-d2j4m7m
KEYWORDS    .
SOURCE      Daucus carota (carrot)
  ORGANISM  Daucus carota
            Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
            Spermatophyta; Magnoliophyta; eudicotyledons; Gunneridae;
            Pentapetalae; asterids; campanulids; Apiales; Apiaceae; Apioideae;
            Scandiceae; Daucinae; Daucus; Daucus sect. Daucus.
FEATURES             Location/Qualifiers
     source          1..155895
                     /organism="Daucus carota"
                     /organelle="plastid:chloroplast"
                     /mol_type="genomic DNA"
     tRNA            complement(6..77)
                     /gene="trnH-GUG"
                     /product="tRNA-His"
                     /note="anticodon:GUG"
                     /standard_name="trnH-GUG tRNA"
     gene            complement(6..77)
                     /gene="trnH-GUG"
                     /standard_name="trnH-GUG gene"
      ...

   	ORIGIN      
        1 ttgggcgaac gacgggaatt gaacccgcgc gtggtggatt cacaatccac tgccttgatc
       61 cacttggcta catccgcccc gccagttttc ttttatttat ttgcatttca aaggattcct
      121 ttttgatcat tcaaaaatat ttgtttatct aaaaaagtct taaataaata aaaaaggagc
      181 aagaccgcct cttgatagaa caagaaagag gttattgctc cttttttaat atttcaaaaa
      ...
   155821 cgttcactaa aaaaaaaccc ttttgtagca aatcgtttat taagaaaaat tgataacctt
   155881 aacacaaaag cagaa
   //
```

**NOTE:** The `gene` name will be used as the name of the target sequence in GeneMiner.



## **4.2 Explanation of parameters**

Using `geneminer.py -h`, the computer will display all options.

```shell
GeneMiner: a software for extracting phylogenetic markers from next generation sequencing data
Version: 1.0.0
Copyright (C) 2022 Pulin Xie
Please contact <xiepulin@stu.edu.scu.cn> if you have any questions

optional arguments:
  -h, --help          show this help message and exit

Basic option:
  -1                  File with forward paired-end reads (*.fq/.gz/.tar.gz)
  -2                  File with reverse paired-end reads (*.fq/.gz/.tar.gz)
  -s , --single       File with unpaired reads (*.fq/.gz/.tar.gz).
  -o , --out          Specify the result folder 
  -rtfa <file|dir>    References of target sequences, only support fasta format
  -rtgb <file|dir>    References of target sequences, only support GenBank format

Advanced option:
  -k1 , --kmer1       Specify the length of the k-mer to filter reads  [default = 29]
  -k2 , --kmer2       Specify the length of the k-mer to assemble reads   [default = 41]
  -d , --data         Specifies the number of reads to reduce raw data. If you want to use all the data, you can set as 'all' 						  [default = 'all']
  -step_length        Step length of the sliding window on the reads [default = 4]
  -limit_count        limit of kmer count [default=auto]
  -limit_min_ratio    The minimum ratio of contig length to reference average length [default = 1.0]
  -limit_max_ratio    The maximum ratio of contig length to reference average length [default = 2.0]
  -change_seed        Times of changing seed [default = 32]
  -scaffold           Make scaffold
  -max                The maximum length of contigs to be retained [default = 5000]
  -min                The minimum length of contigs to be retained [default = 300]
  -t , --thread       Number of threads [default = 'auto']
  -b , --boundary     The length of the extension along both sides of the target sequence [default = 75]
  -bn , --bootstrap   Number of resampling based on base substitution model
```



### 4.2.1 Basic parameters

- `-1`  and `-2` : The input files with paired-end reads, given in FASTQ format (\*.fq\\\*.fastq\\\*.fq.gz\\\*.fastq.gz).Make sure to keep the correct file extensions

- `-s`  : The input file with unpaired readsgiven in FASTQ format(\*.fq\\\*.fastq\\\*.fq.gz\\\*.fastq.gz).Make sure to keep the correct file extensions.

- `-o , --out` : The output folder.

- `-rtfa`  <file|dir>  

  References of target sequences, only support fasta format.(\*.fasta\\\*.fa\\\*.fas\\\*.fna).Make sure to keep the correct file extensions

  **NOTE:**

  (1) Target sequences are derived from homologous genes of closely related species (same genus or family).

  If a target sequence has homologous genes from one or more species, please store these homologous genes in a fasta format file. As follows：

  ```shell
  >Gene_A species1
  ATCGATCG
  >Gene_A species2
  ATCGATCC
  >Gene_A species3
  ATTGATCC
  ```


​		GeneMiner使用hash的方法存储参考，过多的参考将导致内存溢出。基于系统发育标记的长度特征，我们建议参考序列长度在300bp~5000bp，总条目不		超过100000。

​        (2) GeneMiner uses a hash method to store references, and too many references will lead to memory overflow. Based on the length     				   		characteristics of phylogenetic markers, we suggest that the reference sequence length is 300~5000bp and the total entries do not exceed 100000.

- `-rtgb`  <file|dir>  		

  References of target sequences, only support GenBank format(\*.gb).Make sure to keep the correct file extensions.

  **NOTE:**

  GeneMiner不会组装整个质体基因组，仅仅是挖掘GenBank文件中的Gene部分。用户还可以编辑genbank文件，从而指定特定的目标序列

  GeneMiner does not assemble the entire plastid genome, but only mines the `gene` portion (feature=="gene") of the GenBank file.

  Users can also edit the genbank file to specify specific target sequences.

  



### **4.2.2 Advanced parameters**

- `-k1 , --kmer1`  :用于过滤reads的k-mer长度.该参数用于筛选与参考序列高度匹配的reads。我们将参考序列拆分为K-mers并存储为hash表.随后将原始数据中的每条reads，拆分成(l-k1+1)条K-mer，其中，l为测序的读长，k1的大小与参考序列的设定一致。如果某条read的子序列在参考序列哈希表中出现，则把该条read保留并分配给特定目标序列的过滤数据集中。k1的取值与提供原始测序数据的物种和提供参考序列的物种之间的亲缘关系有关.  两者亲缘关系越近，k1的能够接受的取值上限越大；两者亲缘关系越远，k1的能够接受的取值上限越小。较大的k1有利于提高拼接的准确度

- `-k2 , --kmer2`:用于组装reads的k-mer长度.k2 是de Bruijn图中节点的长度, GeneMiner 根据节点之间 k-1 mer的overlap 将k-mer组装为contig.一般来说，较小的k-mer有利于序列在较低覆盖度时的延申，但可能引入错误k-mers；较大的k2能有效处理反向重复，但可能提前终止序列的延申。k2 默认大小为31，k2的取值范围在17~127之间	

  

- `-d , --data` : reads的数目。GeneMiner 允许使用原始测序数据中的一部分挖掘目标序列，特别是具有中-高拷贝数的目标序列。如果您需要将全部的原始测序数据作为输入，您可以将-n设置为‘all’ .-d 默认设定为all

Specify the number of reads to reduce raw data.GeneMiner allows mining target sequences, especially those with medium-high copy numbers, using a portion of the raw sequencing data. If you need to use the entire raw sequencing data as input, you can set -n as `all` .default = `all`

- `-step_length`指定reads过滤过程中滑动窗口的长度,默认值为4



Step length of the sliding window on the reads. With a sequence is ATCGAATTCA, when `step_length` is 1 and `kmer_size` is 5,  we can get  ATCGA, TCGAA, CGAAT, GAATT, AATTC and ATTCA; when `step_length`is 2, we will get ATCGA, GAATT, AATTCA. This parameter can be used to reduce the program runtime when the dataset is large with sufficient coverage

- `-limit_count`:	k-mer频次最低阈值（limit）。在组装过程中，GeneMiner会将过滤后reads拆分为长度为`k1`的子序列(k-mer)，并统计这些k-mers出现的次数.频次低于最低阈值（limit_count）的k-mer将被剔除用户既可以自行设置`-limit_count` ，也可以选择`auto`模式.在`auto`模式下GeneMiner会基于目标序列的k-mers频数分布自行为每一个目标序列分配合理的`-limit_count`. 默认值=`auto`

limit of k-mer count. This parameter is used to remove erroneous, low quality k-mers.During the assembly process, GeneMiner will split the filtered reads into subsequences (k-mers) of length `k1` and count the number of occurrences of these k-mers. The k-mers whose frequency is below the `limit_count`) will be removed.Users can either set a uniform `-limit_count` or choose the `-auto` mode. In `auto` mode, GeneMiner will assign a reasonable `-limit_count` to each target sequence based on the k-mers frequency distribution . default=`auto`





- `-limit_min_ratio`:恢复的目标序列相较于参考序列平均长度的最小比值。在组装过程中，GeneMiner会舍弃比值小于`limit_min_ratio`的序列



The minimum ratio of the recovered target sequence compared to the average length of the reference sequence. During the assembly process, GeneMiner discards sequences with ratios smaller than `limit_min_ratio`





- `-limit_max_ratio`:恢复的目标序列相较于参考序列平均长度的最大比值。在组装过程中，GeneMiner会舍弃比值大于`limit_min_ratio`的序列



The maximum ratio of the recovered target sequence compared to the average length of the reference sequence. During the assembly process, GeneMiner discards sequences with ratios larger than `limit_max_ratio`



- `-change_seed`:Times of changing seed. GeneMiner会自动更换候选种子以达到最优组装效果.default=32

  GeneMiner automatically replaces candidate seeds to achieve optimal assembly





- `-max` :指定GeneMiner挖掘出的基因的最大长度，默认为5000bp.

The maximum length of contigs to be retained.default=5000 bp

- `-min`: 指定GeneMiner挖掘出的基因的最小长度，默认为300bp.

The minimum length of contigs to be retained



- `-t , --thread`  :指定线程数量,如果不指定的话，GeneMiner会根据计算机性能自动选择合适的线程数量.默认值为"auto"

Specify the number of threads, if not specified, GeneMiner will automatically select the appropriate number of threads based on computer performance. default=`auto`



- `-b , --boundary`  指定软边界的长度。沿着目标序列两侧延申的长度当沿着恢复的目标序列的两侧延申的时候，随着延申长度的增加，碱基的准确率逐渐降低.然而，准确率的降低不是断崖式的，而是在某一定长度的缓冲区内逐渐下滑。我们将保留缓冲区并将这段缓冲区称软边界。软边界取值范围在0~200之间，我们推荐大小为 05 * reads length。该参数与-rtgb参数同时使用时才生效.



The length of the extension along both sides of the target sequence (length of `soft boundary`).

When extending along both sides of the recovered target sequence, the accuracy of the bases gradually decreases as the extension length increases.

However, the decrease in accuracy is not precipitous, but a gradual decline within a buffer of some certain length. We will keep the buffer and call this buffer  `soft boundary`.

`-b` takes a range of values from 0 to 200 bp,default =75 bp. Recommended length is 0.5 * reads length

NOTE: `-b` only takes effect when the user uses the `-rtgb` .



- `-bn , --bootstrap`: 重采样的次数基于核苷酸替换模型。GeneMiner从Bootstrap借鉴了思想，然后创新性的提出了基于核苷酸替换模型的迭代校验方法，可以有效的评估结果的准确性[默认值=100]

​	Number of resampling. GeneMiner borrowed ideas from Bootstrap, and then innovatively proposed an iterative verification method based on base substitution model, which can effectively evaluate the accuracy of the results





- `-k1 , --kmer1`: Length of k-mer for filtering reads. This parameter is used to filter the reads that are highly matched with the reference sequence.We split the reference sequence into k-mers and store them as hash tables. Each reads in the original data is subsequently split into (`L`-`k1`+1) k-mers, where `L` is the read length and `k1` is of the same size as the setting of the reference sequence. If a subsequence of a read appears in the reference sequence hash table, the read is retained and assigned to a filtered dataset of a specific target sequence.The value of `k1` is related to the relationship between the species providing the original sequencing data and the species providing the reference sequence.  The closer the relationship between the two species, the larger the acceptable upper limit of `k1`; the more distant the relationship between the two species, the smaller the acceptable upper limit of `k1`. A larger `k1` is beneficial to improve the accuracy of the assembly.  `k1` takes a range of values from 17 to 127 bp,default=29.

- `-k2 , --kmer2`: Length of k-mer for assembling reads. `k2` is the length of nodes in the de Bruijn graph, and GeneMiner assembles k-mer into contigs based on the overlap of k-1 mer between nodes. In general, smaller k-mer facilitates sequence extension at lower coverage, but may introduce erroneous k-mers; larger k-mer handles reverse repetition efficiently, but may terminate sequence extension early.For Illumina reads(150bp) with sufficient coverage (> 40x), we have good results with k = 41. `k2` takes a range of values from 17 to 127 bp,default=41

- `-d , --data`: Specify the number of reads to reduce raw data.GeneMiner allows mining target sequences, especially those with medium-high copy numbers, using a portion of the raw sequencing data. If you need to use the entire raw sequencing data as input, you can set `-d` as `all` .default = `all`

- `-step_length`: Step length of the sliding window on the reads. With a sequence is ATCGAATTCA, when `step_length` is 1 and `kmer_size` is 5,  we can get  ATCGA, TCGAA, CGAAT, GAATT, AATTC and ATTCA; when `step_length`is 2, we will get ATCGA, GAATT, AATTCA. This parameter can be used to reduce the program runtime when the dataset is large with sufficient coverage



- `-limit_count`: limit of k-mer count. This parameter is used to remove erroneous, low quality k-mers.During the assembly process, GeneMiner will split the filtered reads into subsequences (k-mers) of length `k1` and count the number of occurrences of these k-mers. The k-mers whose frequency is below the `limit_count`) will be removed.Users can either set a uniform `-limit_count` or choose the `-auto` mode. In `auto` mode, GeneMiner will assign a reasonable `-limit_count` to each target sequence based on the k-mers frequency distribution .default=`auto`

- `-limit_min_ratio`: The minimum ratio of the recovered target sequence compared to the average length of the reference sequence. During the assembly process, GeneMiner discards sequences with ratios smaller than `limit_min_ratio`

- `-limit_max_ratio`: The maximum ratio of the recovered target sequence compared to the average length of the reference sequence. During the assembly process, GeneMiner discards sequences with ratios larger than `limit_max_ratio`

- `-change_seed`:Times of changing seed. GeneMiner automatically replaces candidate seeds to achieve optimal assembly. default=32



- `-max` :The maximum length of contigs to be retained.default=5000 bp

- `-min`: The minimum length of contigs to be retained.default=300 bp

- `-t , --thread`: Specify the number of threads, if not specified, GeneMiner will automatically select the appropriate number of threads based on computer performance. default=`auto`

- `-b , --boundary` : The length of the extension along both sides of the target sequence (length of `soft boundary`).When extending along both sides of the recovered target sequence, the accuracy of the bases gradually decreases as the extension length increases.However, the decrease in accuracy is not precipitous, but a gradual decline within a buffer of some certain length. We will keep the buffer and call this buffer  `soft boundary`. Recommended length is 0.5 * reads length.`-b` takes a range of values from 0 to 200 bp,default =75 bp

  **NOTE**: `-b` only takes effect when the user uses the `-rtgb` .

- `-bn , --bootstrap`: Number of resampling .Number of resampling. GeneMiner borrowed ideas from Bootstrap, and then innovatively proposed an iterative verification method based on base substitution model, which can effectively evaluate the accuracy of the results.default=100







## 4.3 Output

The output directory contains the`reference_database`, `filtered_out`, `assembled_out` , `GM_results` , `bootstrap_out` and `results.csv`

### 4.3.1 reference_database

`reference_database`  \<folder\>: 用于存放参考序列数据集. GeneMiner对用户提供的近缘类群的参考序列进行预处理，并将这些目标序列写入fasta格式的文件。Used to store reference sequence datasets. GeneMiner preprocesses the reference sequences of the close taxa provided by the user and writes these target sequences to fasta-format files.

### 4.3.2 filtered_out

`filtered_out`  \<folder\>: 用于存放过滤数据集. GeneMiner基于参考数据集检索原始测序数据，将与参考序列高度匹配的reads分配给特定目标序列的过滤数据集。如果某个特定目标序列的过滤数据集文件大小超过10MB，GeneMiner会将它保存在big_reads中，并对其进行重新过滤.

Used to store filtered datasets. GeneMiner search the raw sequencing data based on the reference dataset and assign the reads that are highly matched with the reference sequence to the filtered dataset of a specific target sequence (`specific filtered dataset`). If the file size of `specific filtered dataset` exceeds 10MB, GeneMiner will save it in `big_reads` and refilter.

### 4.3.3 assembled_out

`assembled_out`  \<folder\>: 用于存放组装的contigs。GeneMiner利用specific filtered dataset组装contigs.根据组装完成度的不同，assembled_out下还可以细分为`short_contig`和`contig`将恢复的目标序列同参考序列的平均长度相比，其比值大于`-limit_min_ratio`的contig 将放入contig folder, 比值小于`-limit_min_ratio`的contig 将放入short_contig folder.

Used to store assembled contigs.GeneMiner assembles contigs using `specific filtered datasets`. Depending on the completion of the assembly, it can also be subdivided into `short_contig` and `contig` under `assembled_out`.Compared the recovered target sequence with the average length of the reference sequence, the contig with a ratio greater than `-limit_min_ratio` will be put into the `contig` folder, and the contig with a ratio less than `-limit_min_ratio` will be put into the `short_contig` folder.

### 4.3.4 GM_results 

`GM_results` \<folder\>:用于存放GeneMiner挖掘出的所有的目标序列. GM_results是最重要的文件夹之一

used to store all the target sequences mined by GeneMiner.`GM_results` is one of the most important folders

### 4.3.5 bootstrap_out

`bootstrap_out` \<folder\>:

GeneMiner对GM_results基于核苷酸替换模型对GM_results重采样，并进行迭代校验。bootstrap_out用于存放评估结果。bootstrap_out文件夹下包括：

Used to store the evaluation records.GeneMiner resamples `GM_results` based on the base substitution model and performs iterative checks on `GM_results`. The `bootstrap_out` folder contains `reference_database`,`filtered_out`,`assembled_out`,`GM_results`,`high_quality_results `and `bootstrap.csv`

**NOTE**：`bootstrap_out` is only generated when the ­ `-bn`/`--bootstrap` option is used

(1)`reference_database` \<folder\>:  

​		存放变异的参考序列数据集。首次获取目标序列后，与参考序列进行比对，得到变异率v。基于由参考序列的碱基组成建立碱基替代模型，将目标序列的以变异率v和碱基替代模型进行重采样，获取变异率同样为v的ref~1~, ref~2~.....ref~n~ 

Used to store the mutated reference sequence dataset (MRSD). After the first acquisition of the target sequence, it is compared with the reference sequence to obtain the variation rate v. Based on the base substitution model built from the base composition of the reference sequence, the target sequence is resampled with the variation rate v to obtain ref~1~, ref~2~ ..... ref~n~  



(2)`filtered_out` \<folder\>:  

​		将变异的参考序列数据集作为新的参考序列数据集，使用GeneMiner中相同的步骤，获得新的过滤数据集.

Use the  MRSD as the new reference sequence dataset and use the same steps in GeneMiner to obtain the new filtered dataset

(3)`assembled_out` \<folder\>:

​		GeneMiner利用新的过滤数据集完成组装，并存放新的组装结果

​		GeneMiner completes the assembly using the new filtered dataset and stores the new assembly results

(4)`GM_results` \<folder\>:

​		用于存放使用基于核苷酸替换模型的迭代校验方法挖掘出的所有目标序列

Used to store all target sequences mined by GeneMiner.  Here, GeneMine used an iterative verification method based on base substitution model.



(5)`high_quality_results` \<folder\>:

​	GeneMiner 比较 GM_results 和bootstrap_out/GM_results 中的目标序列，并对每一个目标序列给出得分.score=identity*100, identity代表前后生成的两条目标序列之间的一致度，得分大于99的目标序列将被存放在high_quality_results文件夹中



GeneMiner compares the target sequences in `GM_results` and `bootstrap_out/GM_results`, and gives a score for each target sequence.

`score=identity*100` The target sequences with score higher than 99 will be stored in the high_quality_results folder

(6)`bootstrap.csv` \<file\>:

​		Assessment Record Sheet

### 4.3.6 results.csv 

`results.csv`  \<file\>:

​		Record various information about the target sequences

​		This file consists of comma-separated columns containing various information on each target  sequence found. The file can be easily imported into programs such as Excel. The contents of the columns (from left to right) are explained in this table:

| **Column**          | **Description**                                              |
| ------------------- | ------------------------------------------------------------ |
| gene                | Target sequence's name                                       |
| k1                  | Length of kmer for filtering reads                           |
| re_k1               | Length of kmer for filtering reads after re-filtering        |
| richness            | Approximate sequencing depth based on the target sequence's filtered dataset |
| limit               | Limit of k-mer count                                         |
| seed                | GeneMiner在组装contigs的时候使用的起始序列 The starting sequence used by GeneMiner in assembling contigs |
| k2                  | Length of kmer for assembling reads                          |
| ref_length          | Average length of the reference sequences                    |
| short_contig_length | Length of the short contig                                   |
| contig_length       | Length of the contig                                         |
| scaffold_length     | Length of the scaffold                                       |
| bootstrap_number    | Number of resampling                                         |
| score               | Target sequence assessment score                             |

NOTE: `bootstrap_number` and `score` are  only printed when the ­ `-bn`/`--bootstrap`  option is used

​		



## 5.4 example

### 5.4.1 Mining chloroplast genes from skimming whole genome sequencing (WGS).

当有较为充足数据量和近缘的参考序列时，本软件几乎能从浅层基因组数据中恢复所有的叶绿体基因

GeneMiner can recover almost all chloroplast (cp) genes from shallow genomic data when there is a sufficient coverage (>10x) and closely related reference sequences.

```shell
#Mining single cp gene : matk
geneminer.py -1 skimming_data1.fq.gz -2 skimming_data2.fq.gz -rtfa matK.fasta -t 4 -o matK_out
#Mining multiple cp genes
geneminer.py -1 skimming_data1.fq.gz -2 skimming_data2.fq.gz -rtfa Ref_cp_fasta -t 4 -o Ref_cp_fasta_out
#Mining all cp genes
geneminer.py -1 skimming_data1.fq.gz -2 skimming_data2.fq.gz -rtgb chloro.gb -t auto -b 0 -min 300 -max 5000 -o chloro_gb_out 
```

### 5.4.2 Mining mitochondrial genes from skimming WGS：

​		线粒体基因恢复难度往往会比叶绿体基因大，其一般原因在于：（1）线粒体基因变异较大（2）原始数据未包含太多线粒体基因。针对这种情况，您可以增大原始数据大小以及选择更近源的参考序列.

Mitochondrial (mito)   genes are often more difficult to recover than chloroplast genes, generally because (1) mitochondrial genes are more variable (2) the raw data does not contain too many mitochondrial genes. In this case, you can increase the size of the raw data and choose more closely related reference sequences.

```shell
#Mining mito genes
geneminer.py -1 skimming_data1.fq.gz -2 skimming_data2.fq.gz -rtgb mito.gb -b 0 -min 300 -max 5000 -o mito_gb_out 
```

### 5.4.3 Mining nuclear genes from skimming WGS：

GeneMiner可以挖掘核基因组中的高度重复区 (如nrDNA)

GeneMiner can mine highly repetitive regions in the nuclear genome (e.g. nrDNA)

```shell
#mining ITS
geneminer.py -1 skimming_data1.fq.gz -2 skimming_data2.fq.gz -t 4 -rtfa ITS.fasta -o ITS_out
```



### 5.4.4 Mining Angiosperm 353 genes from transcriptome data

被子植物353基因是近年来在被子植物系统发育研究中被广泛采用的分子标记，由353条基因组成

The angiosperm 353 genes are molecular markers that have been widely adopted in recent years in angiosperm phylogenetic studies and consist of 353 genes

```shell
#Mining Angiosperm 353 genes
geneminer.py -1 Arabidopsis_thaliana_sim_353_data1.fq.gz -2 Arabidopsis_thaliana_sim_353_data1.fq.gz -rtfa Ref_353 -k1 29
-k2 41 -t 4 -o Angiosperm_353_out4
```





# 6. Methods

针对短片段的系统发育标记提取，GeneMiner结合了有参过滤和无参拼接的优势.

(1)首先基于近缘的参考序列采用K-mer过滤的方法获取可靠的reads;

(2)然后采用a weighted seed selection and extension algorithm based on De Brujin Graph方法通过组装reads的小目标区域恢复目标序列;

(3)最后使用基于核苷酸替代模型的迭代方法对结果进行校验.



​		For phylogenetic marker extraction of short sequences, GeneMiner combines the advantages of reference-based filtering and reference-free assembly. The core process of GeneMiner is divided into three main steps:

​		(1)Firstly, reliable reads were obtained based on the reference sequences of closely related taxa using a K-mer filtering method;

​		(2) then a weighted seed selection and extension algorithm based on De Brujin Graph was used to recover the target sequences by assembling small target regions of the reads;

​		(3)Finally, the results are verified using an iterative method based on the base substitution model.



![](picture/geneminer_workflow.svg)





## 6.1 Filtering reads

​		GeneMiner filters reads that are highly similar to the target sequences based on the reference sequences of closely related taxa. The process consists of two stages: building the hash table and filtering reads. In building the hash table, we split the reference sequences into K-mers (k-mer is the iterative division of the reads into sequences containing K bases) and record the position information, the number of occurrences and the corresponding labels of these subsequences, finally store them as a hash table. The total number of k-mer is![img](picture/clip_image002.png), where L~i~ is the length of the i-th reference sequence, n is the total number of reference sequences, and k~f~  is set by the user.

## 6.2 Assembly

​		Geneminer developed a weighted seed selection and extension algorithm based on de Brujin Graph, the general flow is as follows：

​		(Ⅰ) Make kmers. Split the filtered reads into k-mers and build a k-mer set recorded as T. The ka used here may be different from the kf used in the filtering.

​		 (II) Remove low quality kmers. The program will automatically fit the k-mers frequency distribution of the set T, and assign a minimum threshold of k-mer frequency (limit) to each target sequence based on the bottom of the first L-peak, and k-mer with frequency below the limit will be removed.

​		 (III) Choose seed. Split the reference sequences into k-mers using ka and build a k-mer set recorded as R. The k-mer that occurs with high frequency will be considered as a presumed conservative region and kept as a candidate seed if the k-mer also occurs in the set T. GeneMiner automatically replaces candidate seeds to achieve optimal assembly. 

​		(Ⅳ)Seed extend. Take each k-mer as a node and assign a weighted score according to its frequency and position in the reference sequences. The weight score is![img](picture/clip_image002-16585694867633.png), where count represents the frequency of k-mer in the set T, Pos~1~ represents the current assembly position of k-mer, and Pos~2~ represents the average position of k-mer in the set R. Using the candidate seeds as starting nodes, search the set T to find nodes with k-1 overlapping bases as adjacent nodes, where the nodes with high weight scores will be prioritized. Connect the neighboring nodes until no neighboring nodes can be found, and finally the overlap cluster (contig) with the highest cumulative weighted score will be the output.



## 6.3 Verifying results 

​		GeneMiner innovatively proposes an iterative verification method based on base substitution model, which can effectively assess the accuracy of the results. 

​		(I) After the first acquisition of the recovered target sequences, they were compared with the reference sequence to obtain the variation rate (v). 

​		(II) Combine the base substitution model built from the base composition of the reference sequences to resample the recovered target sequence with v to obtain ref~1~, ref~2~ .... ref~n~ 

​		(III) Use ref~1~, ref~2~ .... ref~n~ as the new reference sequences to re-run the whole process and get the new recovered target sequences

​		(IV) The identity of the two results will be used to evaluate the accuracy



# 7. Get help

Please check [Geneminer's homepage](https://github.com/happywithxpl/GeneMiner) first. If your question is running specific,please do not be surprised and report it to us. We usually have quick response to bugs.

- Find Questions & Answers at [GeneMiner Discussions](https://github.com/happywithxpl/GeneMiner/discussions/categories/q-a): **Recommended**

- Report Bugs & Issues at [GetOrganelle Issues](https://github.com/Kinggerm/GetOrganelle/issues):

  Please avoid repetitive or irrelevant issues

- QQ group (ID: 78266311): Mainly for mutual help, responses are likely to be not timely

**Do NOT** directly write to us with your questions, instead please post the questions **publicly**, using above platforms  Our emails (xiepulin@scu.edu.cn, 1791173948@qq.com) are only for receiving public question alert and private data (if applied) associated with those public questions. When you send your private data to us, enclose the email with a link where you posted the question. Our only reply emails will be a receiving confirmation, while our answers will be posted in a public place.



# 8. Citation

When you use GeneMiner please cite:

**GeneMiner : a software for extracting phylogenetic markers from next generation sequencing data**

