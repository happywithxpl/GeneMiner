# GeneMiner
​		Thank you for choosing miner. This document will help you learn how to use GeneMiner. Please note that GeneMiner is still in the testing stage, so this document may be updated in a future version upgrade.
If you have any questions, you can contact xiepulin@stu.scu.edu.cn  or 15162176893@163.com

[TOC]



# 1. About GeneMiner

​	 GeneMiner (gene miner) is a software used to mine genes from next-generation sequencing (NGS) data. It can obtain high-quality specific target genes from very low-quality and deep source data. For example, accurately extract all or part of chloroplast / mitochondrial genes and highly repetitive regions in nuclear genome (such as nrDNA) from skimming genome sequencing data; A large number of single copy phylogenetic markers were extracted from transcriptome sequencing data; Obtaining environmental response genes of specific microorganisms from macrogenomic data. It can be widely used in phylogenetic and evolutionary research, customs inspection and quarantine, mining of specific functional genes, etc. it has significant advantages in reducing sequencing cost and expanding gene selection. GeneMiner's result is very accurate. In the real experimental verification, it reaches or is close to the result of next-generation sequencing. Even if the similarity between the reference sequence and the target sequence is less than 90%, it can still rely on the gradient approximation algorithm to obtain 100% accurate results. The software innovatively proposes a verification method based on self-developed detection, which can repeatedly verify the target sequence without relying on the reference sequence, and output a more reliable consistent sequence. Based on a large number of optimization at the algorithm level, GeneMiner has excellent computing speed and memory consumption. It supports multi-threaded parallelism and makes full use of computer resources. It can be deployed on high-performance computing clusters or ordinary personal computers. GeneMiner is very user-friendly and supports various mainstream operating systems of windows, MAC and Linux. Users can choose command-line interface or graphical interface.

# 2. Downloading GeneMiner
​		GeneMiner is open source under MIT license. It is distributed through the github Repository: https://github.com/sculab/GeneMiner , you can download the latest version at any time. Be sure to keep an eye on github to keep your code up to date. We do not provide any support for previous versions of code! The version number follows the symbol x.y.z, where x changes with major code reorganization, y changes when new features are added, and changes with bug fixes

# 3. Installing GeneMiner

## 3.1 For Linux users

Automatic installation (recommended)

```shell
tar -zxvf geneminer.tar.gz  # decompression
cd geneminer
python setup.py   #according to the script prompt, automatically install the dependency and write GeneMiner into the environment variable
```

Manual installation (used when automatic installation encounters problems)

```shell
tar -zxvf geneminer.tar.gz # decompression
cd geneminer
# Manually install the required libraries
pip3 install  biopython --user
pip3 install  pandas --user
pip3 install  tqdm --user
pip3 install  openpyxl --user
pip3 install  pysimplegui --user
#It can also be installed in batches according to the dependent file provided by the software
pip3 install -r requirements.txt --user
#Configure environment variables
echo "export PATH=\$PATH:$(pwd)" >> ~/.bashrc
source ~/.bashrc
```

Close and restart the terminal to check whether the configuration is successful

```shell
geneminer -h # command line interface
```

## 3.2 For Mac users

Download the corresponding version of the packaged graphical interface app and double-click it to run. (recommended)
To manually configure the command line and GUI versions, use the following commands for automatic installation:

```shell
tar -zxvf geneminer.tar.gz # decompression
cd geneminer 
python3 setup.py #according to the script prompt, automatically install the dependency and write GeneMiner into the environment variable
```

Install or use the following command manually:

```shell
tar -zxvf geneminer.tar.gz # decompression
cd geneminer
# Manually install the required libraries
pip3 install  biopython --user
pip3 install  pandas --user
pip3 install  tqdm --user
pip3 install  openpyxl --user
pip3 install  pysimplegui --user
#It can also be installed in batches according to the dependent file provided by the software
pip3 install -r requirements.txt --user

#Configure environment variables
#For MacOS Catalina (10.15) and later systems:
echo "export PATH=\$PATH:$(pwd)" >> ~/.zshrc
source ~/.zshrc
#For systems before MacOS Catalina (10.15):
echo "export PATH=\$PATH:$(pwd)" >> ~/.bash_profile
source ~/.bash_profile
```

关闭并重启终端，检测是否配置成功

```shell
geneminer -h # command line interface
```

## 3.3 For Windows users

### 3.3.1 Install WSL for windows

​		Because geneminer needs to use WSL in windows, it only supports windows 10 and later operating systems. WSL is not the default function of the system. You need to install WSL in the system first. For windows 10 version 2004 and later (build 19041 and later) or windows 11, you can use the following steps to install WSL. During the installation process, you need to connect to the Internet:

- Run PowerShell as administrator: find windows PowerShell in the start menu, right-click and select run as administrator.
- In the command line window that opens, enter:

```shell
wsl --install
```

- After the installation is completed, run WSL in the command line window to confirm that the installation is successful.

```shell
wsl
```

​		The first time you start a newly installed Linux distribution, a console window opens asking you to wait for the files to be extracted and stored on your computer. All future start-up times should be less than one second.

​		For the old version of windows 10, it is recommended that you upgrade to the latest version or use the manual installation method of the old version. For details, please refer to the technical documents of Microsoft:

- https://docs.microsoft.com/zh-cn/windows/wsl/install-manual

### 3.3.2 Install packaged GUI applications (recommended)

Download the corresponding version of the graphical interface application of geneminer and double-click to run it. Geneminer Windows version maual, see xxxxxxxxxxxxxxx for details

### 3.3.3 Configuring command line interface applications

​	To configure geneminer in the command line interface, you need to install Python version 3.6 or above in the system. You can also install Anaconda or miniconda. For the specific installation method, please refer to the official technical documents of Python and its different distributions. You can refer to the following steps to configure geneminier:
Automatic installation

- Download the windows installation package of geneminer and double-click to install automatically.
  If the automatic installation encounters problems, you can perform a manual installation:


Manual installation

- Unzip: unzip the downloaded windows installation package to a specified folder, such as geneminer.
- Open the command line and manually install the required libraries:

```shell
pip3 install  biopython --user
pip3 install  pandas --user
pip3 install  tqdm --user
pip3 install  openpyxl --user
pip3 install  pysimplegui --user
#It can also be installed in batches according to the dependent file provided by the software
pip3 install -r requirements.txt --user
```

- 将geneminer文件夹加入用户环境变量path中。
- 关闭并重启终端，检测是否配置成功

```shell
geneminer -h # command line interface
```

# 4. Quick  start

​	Before mining your sequencing data, we suggest you first understand your sequencing data, including sequencing method, depth, quality, data volume, etc. Our software is mainly applicable to the second-generation sequencing data returned by Illumina platform. Verified by a large number of test data, even for low sequencing depth (below 10x), GeneMiner can mine single copy nuclear gene, chloroplast gene and mito gene, and can also be used to mine single copy genes from transcriptome data.
​	 The data used for hand training is stored in geneminer/example/

(1) Mining single target gene

```shell
 geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta 
```

data1.fq and data2.fq is the sequencing data returned from second-generation sequencing, and ref.fasta is the homologous gene of related species (same genus or same family). In this case,
Note: each FASTA file can only store the same gene (it can be from different species) . As shown in the following example:

```
>GeneA_species1
ATCGATCG
>GeneA_species2
ATCGATCC
>GeneA_species3
ATTGATCC
```

(2) Mining multiple target genes
Multiple files in FASTA format can be placed in one folder, and genes can be mined in batches

```shell
 geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta
```

（3）Batch mining of different types of genes

```shell
 geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rcp cp.gb -rmito mito.gb
```

cp.gb is the chloroplast reference genome in GenBank format, Mito GB is the mitochondrial reference genome in GenBank format, - rcp and - rmito are used to specify the reference genome type


（4）Evaluate reliability

```shell
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta  -bn 10
```

-bn : number of self-development tests. The method based on bootstrap idea is used to evaluate the accuracy of the results. The larger the self expanding detection value, the more time-consuming the calculation

# 5.Detailed usage guide

## 5.1 Parameter interpretation

With geneminer-h, the computer displays all the options, and I'll explain how to use them and some tips below

```shell
geneminer -h
usage: GeneMiner <-1 -2|-s|-12>  <-rn|rcp|-rmito|rt>  [options]

GeneMiner: a software for extracting phylogenetic markers from skimming genome
Version: 1.0
Copyright (C) 2021 Pu-lin Xie
Please contact <xiepulin@stu.edu.scu.cn>, if you have any bugs or questions

optional arguments:
  -h, --help            show this help message and exit

Basic option:
  -1                    one end of paired-end reads,support fastq/fastq.gz/fastq.bz2 
  -2                    another end of paired-end reads,support fastq/fastq.gz/fastq.bz2
  -12                   interlaced forward and reverse paired-end reads,support fastq/fastq.gz/fastq.bz2
  -s , --single         single-read,support fastq/fastq.gz/fastq.bz2
  -o , --out            Specify the result folder [default='auto']
  -rcp <file|dir>       reference of chloroplast genome,only support GenBank-format
  -rmito <file|dir>     reference of mitochondrial genome,only support GenBank-format
  -rtfa <file|dir>      References of target genes, only support fasta format
  -rtgb <file|dir>      References of target genes, only support GenBank format
  
Advanced option:
  -n , --number         The number of rows of raw data from skimming genomes,default=1000000
  -k , --kmer           size of a kmer  [default =31]
  -max                  gene maximum length  [default =5000]
  -min                  gene minimum length  [default =300]
  -t , --thread         Specify the number of threads you want to run [default='auto']
  -b , --boundary       extend the length to both sides of the gene while extracting                            					 genes from  Genebank file [default=75]
  -sf                   Select the reference sequences to reduce the computation.
                        s1: do nothing;
                        s2,3,4: only use the reference sequence with the shortest/median/longest length;
                        s5: remove sequences with abnormal length.[default = 's1']

Gradient approximation option:
  -in , --iterative_number
                        Specify the number of iterative loop to gradually approximate the best results

Bootstrap option:
  -bn , --bootstrap_number 
                        Specify the bootstrap number.Evaluate the results based on the bootstrap method
```

### 5.1.1 Basic parameters

```shell
-1				
One end of paired-end reads,support fastq/fastq.gz/fastq.bz2 format.Be sure to keep the correct file extension.
example: geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa  ref.fasta
    
-2   			
Another end of paired-end reads,support fastq/fastq.gz/fastq.bz2.Be sure to keep the correct file extension.
example: geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta

-12  			
Interlaced forward and reverse paired-end reads,support fastq/fastq.gz/fastq.bz2.Be sure to keep the correct file extension.
example: geneminer.py -12 data.fq.gz -rtfa  ref.fasta

-s   			
Single-read,support fastq/fastq.gz/fastq.bz2.Be sure to keep the correct file extension.
example: geneminer.py -s data1.fq.gz -rtfa  ref.fasta

-rcp 	<file|dir>	
Reference of chloroplast genome,only support GenBank-format. 
You can input a GenBank file, which can contain either one species or multiple species; You can also input multiple Genebank files into a single folder.It is worth noting that :(I) users can modify GeneBank files to keep only parts of the genes they are interested in. (ⅱ) When selecting the reference genome, try to select the closest reference genome
example: geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rcp cp.gb

-rmito    <file|dir>	
Reference of mitochondrial genome,only support GenBank-format.The specific usage is the same as -rcp
example: geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rmito mito.gb

-rtfa     <file|dir>  		
References of target genes, only support fasta format. It can be used to find genes of interest
example: geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta
A fasta-format file can only store the data of the same gene, as shown in the following example:
>species_a ITS
AGCTAGCT
>species_b ITS 
AGCTAGCC
>species_c ITS
AGCTAGCA
t4 species_d ITS
AGCTAGAA

-rtgb     <file|dir>  		
References of target genes, only support GenBank format.It can be used to find genes of interest.
-rtgb includes -rcp and -rmito. The reason why -rcp and -mito are independent is for the convenience of users and the extension of software functions in the later period
example: geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtgb gb_folder

-o , --out            
Specify the output folder. If you do not specify an output folder, geneminer.py will automatically use 'GM+ timestamp 'as the output folder name. But Windows and MAC GUI versions must specify the output folder.
example: geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta -o geneminer_out

-sf   
Select the reference sequences to reduce the computation.[default = 's1']
When the number of reference sequences is too large or the difference between them is too large, it not only increases the calculation time but also affects the accuracy of the results. GeneMiner currently provides five strategies for filtering reference sequences.
strategy 1(s1):do nothing;  
strategy 2(s2):only use the reference sequence with the shortest length
strategy 3(s3):only use the reference sequence with the median length
strategy 4(s4):only use the reference sequence with the longest length
strategy 5(s5): remove sequences with abnormal length
example: geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta -sf s5
```

### 5.1.2 Advanced parameters

```shell
-n , --number     
Enter the number of rows for the raw data amount. The default value is 10 million rows. If you need to enter all raw data, you can set -n as all.
After testing, good results can be obtained by selecting only a part of the original data(100w~1000w), while greatly reducing the running time of the software
For the second-generation sequencing data with a read length of 150bp, the data volume of 1000W lines is about 800~1000MB.
If you are interested in the number of rows of raw data, you can view your data using the following command
zcat your_data.fq.gz|wc -l 
example:
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta -n 2000000
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta -n all

-k , --kmer   
Specifies the length of k-mer, which is the length of the node in the de Bruijn diagram. The value of kmer depends heavily on the data set.
The default value is 31. The value of kmer ranges from 15 to 127
example:
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta -n 2000000 -k 43

-max 
Specifies the maximum length of the mined gene. Default is 5000bp.
example:
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rcp cp.gb -n 2000000 -k 43 -max 3000

-min    
Specifies the minimum length of the mined gene. Default is 300bp.
example:
geneminer.py  -1 data1.fq.gz  -2 data2.fq.gz -rcp mito.gb -n 2000000 -k 43 -max 3000     -min 200

-t , --thread  
example:
Specify the number of threads. If not, the software automatically selects the appropriate number of threads based on computer performance.
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rcp cp.gb -n 200000 -k 43 -max 3000 -min 200 -t 8

-b , --boundary   
Specifies the length of the soft boundary.
As the lines were traced along the sides of the excavated gene, the accuracy of the bases became less and less accurate as the lines increased. The descent is not precipitous, but gradual within a buffer zone of a certain length. We will preserve a buffer and call it a soft boundary. The recommended size is 0.5*reads‘ length, and the soft boundary ranges from 0 to 200. This parameter can be used together with -rcp,-rmito, and -rtgb
example:
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rcp cp.gb -n 2000000 -k 43 -max 3000 -min 500 -b 75
```



### 5.1.3 Gradient approximation parameters (optional)

```shell
-in , ----iterative_number
Specify the number of iterations. As the number of iterations increases, the result gradually approaches the optimal answer.
The basic principle is: The results of GeneMiner are taken as the reference sequence of the next input. If the results of the next output are highly consistent with the results of the last output, it will stop. Otherwise, repeat the process. The maximum number of repetitions is the number specified by the user. [Default =2]
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta  -cn 2
```

### 5.1.4 **Method based on bootstrap** (optional)

```shell
-bn,--bootstrap_number  
The verification method based on bootstrap can evaluate the accuracy of the result and verify the target sequence repeatedly without relying on the reference sequence
[Default =10]
example:
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta  -bn 20
```



## 5.2 Interpretation of results

```shell
GeneMiner will produce a large number of files, which are divided into three main parts:
# Part 1
Filtered sequencing data: "data1.fq" and "data2.fq"

# Part 2
Results: Excel files ("results_information.xlsx") and logs ("log.txt") containing various statistics. 

"results_information.xlsx":

| Column name			    | Example    |  Explanation                                                    
| ------------------------  | ---------- | ----------------------------------------------------------- 	
| nuclear_gene_name         | ITS        | gene name                                                       
| filtered_reads_number     | 300        | The number of reads after filtration                             
| richness                  | 75.88      | richness,The average depth of genes. richness =(n*L1)/L2 
                                           N: average length of reads,L2: The length of the gene             
| graph_construction        | 30.766     | graph construction step(minia)                               
| assembled_percentage      | 0.186      | Percentage of reads used for assembly                             
| assembled_max_length      | 1181       | Maximum length of sequence after assembly
| result_max_length         | 1181       | The maximum length of the aligned-cut(trimmed) sequence  
| identity_trimmed_sequence | 85.83%     | The identity bewteen the reference sequence and the trimmed 	    
                                           sequence
| coverage_trimmed_sequence | 98.52%     | The coverage of the reference sequence over the trimmed sequence
| Gene_extraction           | successful | Whether the gene was extracted successfully                      
| Gene_aligned_cut          | successful | Whether the gene can be aligned and cut successfully            


# Part 3
Result summary files, different types of reference sequences will generate result summary files with different names.-rtfa,-rtgb, -rcp, and -rmito correspond to "tfa_genes", "tgb_genes", "cp_genes", and "mito_genes" respectively

The results summary file can also be subdivided into:
(1) Reference sequence folder  ("reference_database")
	GeneMiner preprocessed the reference sequence provided by the user and wrote each gene into a fasta-format file separately.
(2) filtered folder ("filtered_out")
	Stores reads filtered by the Filter script in fastq format
(3) Assembled folder (" assembled_out ")
	Store contigs and untigs generated by minia  assembly
(4) Results folder ("GM_results")
	Put together all the genes that have been mined.
	First, if the target gene is not mined, no result file will be generated. If the target gene is mined, it is retained as the original result ("xxx_raw.fasta") without any processing. XXX represents a gene, and the following statement is the same.
	Then, if the target gene is excavated but the gene does not meet various subsequent screening conditions, GeneMiner will select the best result from the original result according to the reference sequence and store it as the best original result ("xxx_raw_best.fasta"). If the target gene is mined and the gene meets various subsequent screening conditions, a) when only one sequence is obtained, a unique result ("xxx.fasta") will be generated. B) Generate candidate results ("xxx_options.fasta") when the result has multiple sequences.
	Finally, GeneMiner trimmed the sequence from the unique result (“xxx_fasta”) or  candidate results ("xxx_options.fasta") with the references, and selected the best one  as the trimmed result ("xxx_trimmed.fasta"). It is worth noting that not all result sequences can be align and cut, and users can also align and cut depend on their specific needs.


"xxx_raw. Fasta" 		: raw result
"xxx_raw_best.fasta "	: best raw result
"xxx. fasta" 			: unique result (recommended)
"xxx_options.fasta" 	: candidate result
"xxx_trimmed.aasta" 	: Trimmed results (recommended)

(5) bootstrap folder (" bootstrap_out ", optional)
	Store results generated by method based on bootstrap, including the reference sequence library (" reference_Ddatabase "), filtering results (" filtered_out "), assembly results (" assembled_out "), and final mining results ("GM_results")
	It is worth noting that only the trimmed sequences can be verified using the method based on bootstrap.
(6) Gradient approximation folder ("interative_out", optional)
	Store the results of each iteration

example:
geneminer.py -1 data1.fq  -2 data2.fq  -rtfa ref.fasta -bn 5 -cn 2 -o GeneMiner_out
tree -L 4 GeneMiner_result/   #View the resulting file directory structure
```

## 5.3 example

### 5.3.1 Chloroplast gene extraction:

​		When there is a relatively sufficient amount of data and a close reference sequence, the software can almost extract all chloroplast genes from the skimming genome data at the same time in the phylogenetic research, which This provides another way to solve the problem that chloroplast assembly results are not cyclic.

```shell
#example:
geneminer.py -1 data1.fq.gz -2 data2.fq.gz  -rmito ref_cp.gb -b 0 -max  3000    -min 300 -o results
```

### 5.3.2 Extraction of mitochondrial genes:

​		Mitochondrial genes are often much more difficult to excavate than chloroplast genes for the general reason that the original data itself does not contain many mitochondrial genes. In this case, you can appropriately increase the size of the raw data and choose a reference sequence closer to the source.

```shell
#example
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rmito mito.gb  -n 15000000
```

### 5.3.3 Nuclear gene Extraction:

​		Even in low sequencing depth data, GeneMiner can mine medium and high copy number genes, such as rDNAs.

```shell
#example
#ITS
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rn ITS_ref.fasta -o ITS_out
#18S
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rn 18S_ref.fasta -o 18S_out
#ETS
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rn ETS_ref.fasta -o ETS_out
#26S
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rn 26S_ref.fasta -o 26S_out
```

# 

# 6.Methods

The GeneMiner core process is divided into three steps:

![流程图改改改](https://gitee.com/xiepulin/picgo_xpl/raw/master/GeneMiner_picture/流程图改改改-16453437397731.svg)



## 6.1 Data Filtering

​		Different from the traditional method of first assembling sequences and then mapping to the reference genome, GeneMiner chose to map the sequence to the reference genome and then assemble the sequence. First, the FASTQ format of the next generation data as the original input data, according to the actual situation, select a suitable size of data volume. Then Ukkonen's algorithm was used to construct a suffix tree from the user-provided reference genome. Subsequently, reads of length L in the original data were split into (L-K +1) sub-sequences of length K, where the size of L depends on the sequencing method and the size of K depends on the kmer set by the user. Finally, if a kmer of an reads is detected in the suffix tree, that reads is retained.

​		Due to the strategy of filtering before assembling, the amount of reads data to be filtered is quite large, while the amount of reads data to be assembled is significantly reduced, which means that assembling is no longer the bottleneck of calculation, and filtering is the core step to be broken through in this method. Fortunately, we reduced the algorithmic complexity of the filtering step to O (n), greatly reducing the time cost.

## 6.2 Assembling and verifying

### 6.2.1 Assembling reads

The filtered reads were assembled into contigs, and minia3 was selected as the sequence assembling software. The main reason is that, different from other sequence assembling software, such as Velvet, Spades, Soapdenovo, Minia3 has high accuracy, fast speed, smaller memory consumption, and is more suitable for short sequence assembling.

### 6.2.2 Data verification

​		Minia3 may generate multiple contigs, and these sequences are not all the ones we need, so GeneMiner performs further processing on these contigs.

**(1) Length correction**

​		Contigs that are too long or too short are filtered based on the average length of the reference sequence.

 **(2) Direction correction**

​		GeneMiner uses makeblastdb and blastn tools in BLAST+ suite to retain contigs in the same direction as the reference sequence and reverse complement contigs in the different directions.

**(3) Align and cut contigs**

​		In order to facilitate users to directly use the results of GeneMiner, we cut contigs and GeneMiner processed contigs as follows: (1) A contig was taken out and denoted as sequence A. (2) A reference sequence was taken out and denoted as sequence B. (3) Local comparison was made between sequence A and sequence B to record the starting and ending sites of high-score fragments (HSP, High Scoring Pair). (3) If the high-scoring fragment is close to the length of the reference gene sequence, write gene_trimmed. fasta, otherwise no gene_trimmed. fasta file is generated

**(4) Get the best result** 

​		If there are more than one contig in xxx_trimmed.fasta,these contigs will be  pairwise compared with reference sequences. The contig with the highest consistency and the longest length is reserved as the best result.

### **6.2.3** Gradient approximation

​		GeneMiner optimizes the results through continuous iteration, which we call gradient approximation.

​		Here, we define the complete process of the user entering the reference sequence into GeneMiner to generate the best result as an iteration. In the next iteration, the reference sequence is replaced by the results generated in the previous iteration, and the previous data stitching and verification process is repeated. This process is repeated for a predefined number of iterations, or until no new sequences are found. As the number of iterations increases, more reads can be matched and assembled, so the generated sequence usually grows with each iteration. The iterative method can not only extend the sequence length and retain more phylogenetic information, but also make up for some defects of the reference sequence, especially when the reference sequence is relatively distant or there are holes in the reference sequence. Through the iterative method, it is possible for us to get the optimal solution gradually

## 6.3 Method based on bootstrap

​		GeneMiner innovatively proposed a verification method based on self-propagation detection, which can evaluate the accuracy of the results without relying on the reference sequence, and verify the target sequence repeatedly, so as to output a more reliable consistent sequence.

(ⅰ) After the target gene was obtained for the first time, it was compared with the reference sequence to obtain the mutation rate V. In order to prevent inconvenient calculation of the difference degree due to multiple contigs or reference sequences in the results, GeneMiner selected the contig with the highest contig fit and reference sequence as the standard to calculate the difference degree according to the consistency and coverage between the sequences.

(ⅱ) Each locus of the target gene was randomly resamaged at the variation rate V to obtain ref1 and ref2..... with the same variation rate V refn

(ⅲ) use ref1, ref2..... refn is used as a reference sequence to re-run the whole process and obtain new target genes target1 and target2.... targetn

(Ⅳ) on target1 target2... targetn obtains the consistency sequence, and the consistency ratio of the site is the approval rate of the site

# 7. Frequently asked questions

**Q: ****How do I verify the reliability of the results?**

**A**: Our software sets a method based on bootstrap to check the results. If you have questions about your results, you can use the -bn option and set the self-expansion times. However, self-expansion verification requires repeated calls to GeneMiner internal scripts. At the same time,next generation sequencing data on NCBI or real first-generation data can be used as supplementary verification

**Q: What if I don't get the results?**

**A**: The main reasons affecting the success rate of gene mining are as follows: First, the quality of original data. Clean_data works better than raw_data. Secondly, the amount of data is too small, which may lead to insufficient gene abundance and cannot be assembled. A large amount of data may lead to a variety of splicing situations. You are advised to set -n to 100w to 1000w. Third, the choice of reference sequence. Better results can be obtained by selecting sequences of related genera or different species of the same genus as reference sequences. At the same time, when there are multiple sequences as reference sequences, the difference between reference sequences should not be too large.

**Q: What is the appropriate value for “-n“**

**A:** Generally speaking, the default parameter n=10000000, that is, 2500000 reads, can well meet the needs. However, for single-copy or low-copy sequences, the amount of data can be appropriately increased. We do not recommend entering all the data, because not only will it slow down a lot, but it will not get good results due to the increasing number of concatenation cases.

# 8.Citation

GeneMiner : a software for extracting phylogenetic markers from next generation sequencing data

