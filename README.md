
# Overview

GeneMiner is a software for extracting phylogenetic markers from next generation sequencing (NGS) data. With GeneMiner, users can accurately and efficiently obtain a large number of target sequences from a wide range of histological data at a very low cost, such as: all or part of the chloroplast/mitochondrial genome, highly repetitive regions in the nuclear genome (e.g. nrDNA), etc. from shallow whole-genome sequencing data; single-low copy genes from transcriptome sequencing data, etc. GeneMiner improves existing phylogenetic research strategies at the most basic data level, and has significant advantages in reducing experimental costs and expanding phylogenetic marker selection. In addition, GeneMiner can be applied to the extraction of DNA barcodes, customs inspection and quarantine, mining of specific functional genes and other research, which has broad application prospects.



# Dependencies

[Python](https://www.python.org/downloads/) 3.6 or later, along with the Python libraries

- [biopython](http://biopython.org/wiki/Main_Page) 1.79 or later



# Download and install

GeneMiner is an easy-to-use software written in python, which is is provided for x86-64 systems running GNU/Linux, macOS (version 10.13 or higher) and Windows (64-bit, version 7 or higher).

Users on Windows, macOS and Linux computers can run GeneMiner directly from the command line. we also offer a more convenient GUI version for Windows and macOS users.

## GeneMiner with GUI

We strongly recommend using the GUI version for users who are not familiar with the command line or light use.Download the corresponding version of the packaged GUI from [here](:https://github.com/happywithxpl/GeneMiner/releases) and double click to run it.

![图片](https://github.com/happywithxpl/GeneMiner-Test/blob/main/GeneMiner_GUI.png)

## **GeneMiner with command line**

- option1  **Cloning the repo**
- option2  **Source distribution** 
- option3  **Flexible construction**



**Cloning the repo**  (support)

Instead of downloading the source distribution as a compressed archive, you could clone the repo and build it as shown below.

```shell
git https://github.com/happywithxpl/GeneMiner.git
cd GeneMiner
python setup.py install --record logName --user
```



**Source distribution** 

To download the source distribution from a [release](https://github.com/bpp/bpp/releases) and  install dependencies, use the following commands:

```shell
wget https://github.com/happywithxpl/Geneminer/geneminer-1.0.0-linux-x86_64.tar.gz
tar geneminer_v1.0.0.tar.gz
cd  GeneMiner_v1.0.0
python setup.py install --record logName --user
geneminer.py -h
```



**Flexible construction**

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



# Running GeneMiner

Without any options, GeneMiner takes a reference database and a query sequence file as input and produce phylogenetic markers.We have prepared a simulated dataset of `Arabidopsis thaliana` to help you quickly use GeneMiner.you can download them from [here](https://github.com/happywithxpl/GeneMiner-Test).



## GeneMiner-cmd

- Mining single target sequence from skimming whole genome sequencing (WGS).

```shell
geneminer.py -1 skimming_data1.fq.gz -2 skimming_data2.fq.gz -rtfa ITS.fasta -o ITS_out
```

`-1`,`-2`:  the input files with paired-end reads, given in FASTQ format.

`-rtfa`:  reference sequences in fasta format. 

`-o`:        the output directory



- Mining multiple target sequences  using from skimming WGS

```shell
#for fasta format
geneminer.py -1 skimming_data1.fq.gz -2 skimming_data2.fq.gz -rtfa Ref_cp -o cp_out1 

#for GenBank format
geneminer.py -1 skimming_data1.fq.gz -2 skimming_data2.fq.gz -rtgb cp.gb -o cp_out2 

```

`-rtfa`:  reference sequences in fasta format

`-rtgb`:   reference sequences in GenBank format



- Mining multiple target sequences from skimming WGS and evaluating the accuracy of the results.

```shell
geneminer.py -1 skimming_data1.fq.gz  -2 skimming_data2.fq.gz -rtfa Ref_cp -min 300 -max 5000 -limit_count 3 -t 4 -bn 20 -o cp_verify_out3
```

`min`: The maximum length of contigs

`max`: The minimum length of contigs 

`limit_count`: The minimum threshold of kmer count is used to erroneous, low-abundance K-mers

`-t`:   The number of threads

`-bn`:   Specify the bootstrap number. Evaluate the results based on the bootstrap method



- Mining Angiosperm 353 genes from  transcriptome data

```shell
geneminer.py -1 Arabidopsis_thaliana_sim_353_data1.fq.gz -2 Arabidopsis_thaliana_sim_353_data1.fq.gz -rtfa Ref_353 -k1 29 -k2 41 -t 4  -o Angiosperm_353_out4
```

`k1`:  Specify the length of the k-mer to filter reads

`k2`:  Specify the length of the k-mer to assemble reads



## GeneMiner-GUI

GeneMiner-GUI is easy to use and suitable for lightweight development. Considering the memory and time overhead, we recommend you to use GeneMiner-cmd when you need to run large scale data

For GeneMiner-GUI, you must set `Data1`, `Data2` or `Single reads` under the `Data` module、`Ref.(fasta)` or `Ref.(gb)` under the `Reference` module and  `Output Folder` under the `Outout` module as shown in the image below

![图片](https://github.com/happywithxpl/GeneMiner-Test/blob/main/run_GeneMiner_GUI.png)



## View results

Now, you can view the results of GeneMiner in the output directory.

The output directory contains  `reference_database`、`filtered_out <folder>` 、`assembled_out`、 `GM_results`、`bootstrap_out `  and `results.csv` .

`reference_database <folder>  `: Used to store reference sequence datasets

`filtered_out <folder>`：Used to store filtered reads

`assembled_out <folder>`: Used to store assembled contigs 

`GM_results <folder>`:Used to store all the target sequences mined by GeneMiner

`bootstrap_out <folder> ` :Used to store assessment results, If you have used the `-bn` parameter

`results.csv <file>` :Record various information about the target sequences



# User manual

A more complete manual is here: https://github.com/happywithxpl/GeneMiner/blob/main/GeneMiner_User_Guide.md

# Contact

Please check [Geneminer's homepage](https://github.com/happywithxpl/GeneMiner) first. If your question is running specific, please do not be surprised and report it to us. We usually have quick response to bugs.

- Find Questions & Answers at [GeneMiner Discussions](https://github.com/happywithxpl/GeneMiner/discussions/categories/q-a): **Recommended**
- Report Bugs & Issues at [GeneMiner Issues](https://github.com/happywithxpl/GeneMiner/issues): Please avoid repetitive or irrelevant issues
- QQ group (ID: 78266311): Mainly for mutual help, responses are likely to be not timely

**DO NOT** directly write to us with your questions, instead please post the questions **publicly**, using above platforms  Our emails (xiepulin@scu.edu.cn, 1791173948@qq.com) are only for receiving public question alert and private data (if applied) associated with those public questions. When you send your private data to us, enclose the email with a link where you posted the question. Our only reply emails will be a receiving confirmation, while our answers will be posted in a public place.



# Citation

When you use GeneMiner please cite:

**GeneMiner : a software for extracting phylogenetic markers from next generation sequencing data**

