
# Overview

GeneMiner is a software for extracting phylogenetic markers from next-generation sequencing (NGS) data. (i) With GeneMiner, users can accurately and efficiently obtain large numbers of target gene fragments from NGS data at a meager cost. For example, extraction of all or part of mitochondrial/chloroplast genes and highly repetitive regions (e.g., nrDNA) in the nuclear genome from genome skimming data. Single-low copy genes can also be extracted from transcriptome sequencing data, etc.GeneMiner broadens the choice of phylogenetic markers from the most basic data level. (ii) We propose a novel verification method based on the base substitution model and repetitive resampling, which can statistically evaluate the impact of the reference on the assembly. (iii) GeneMiner provides a cross-platform graphical interface and can be easily docked to downstream phylogenetic analysis processes. In addition, GeneMiner can be applied to the research of DNA barcode extraction, customs quarantine, specific functional gene exploration, and other research, which has broad application prospects.



# Dependencies

[Python](https://www.python.org/downloads/) 3.6 or later, along with the Python libraries.

- [biopython](http://biopython.org/wiki/Main_Page) 1.79 or later



# Download and install

GeneMiner is an easy-to-use software written in python3, which is provided for x86-64 systems running GNU/Linux, macOS (version 10.13 or higher), and Windows (64-bit, version 7 or higher).

Users on Windows, macOS, and Linux can run GeneMiner directly from the command line. We also offer a more convenient GUI version for Windows and macOS users.

## GeneMiner with GUI

For individuals who are not accustomed to utilizing the command line or for light use, we strongly advise using the GUI version. Download the corresponding version of the packaged GUI from [here](:https://github.com/happywithxpl/GeneMiner/releases) and double-click to run it.

![图片](https://github.com/happywithxpl/GeneMiner-Test/blob/main/GeneMiner_GUI.png)

## **GeneMiner with command line**

- option1  **Cloning the repository**
- option2  **Source distribution** 
- option3  **Flexible construction**



**Cloning the repository**  (support)

Clone GeneMiner's repository directly and build it as below:

```shell
git clone https://github.com/happywithxpl/GeneMiner.git
cd GeneMiner
python setup.py install --record logName --user #Add 'geneminer.py' to the '$PATH' 
geneminer.py -h
```

If you want to  uninstall GeneMiner and remove `geneminer.py` from the `$PATH`, you can:

```
cat logName | xargs rm -rf  
```



**Source distribution** 

 Download the source distribution from the [release](https://github.com/bpp/bpp/releases) and  install dependencies:

```shell
wget https://github.com/happywithxpl/Geneminer/geneminer-1.0.0-linux-x86_64.tar.gz
tar geneminer_v1.0.0.tar.gz
cd  GeneMiner_v1.0.0
python setup.py install --record logName --user
geneminer.py -h
```

If you want to  uninstall GeneMiner and remove `geneminer.py` from the `$PATH`, you can:

```
cat logName | xargs rm -rf  
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

- Add `geneminer.py` to the `$PATH`. The following is an example for linux users:

```shell
echo "export PATH=\$PATH:$(pwd)" >> ~/.bashrc  
source ~/.bashrc
geneminer.py -h
```



# Running GeneMiner

GeneMiner takes the reference and fastq format sequencing files as input and the recovered phylogenetic markers as output. We have prepared a simulated dataset of `Arabidopsis thaliana` to help you quickly grasp the main usage of GeneMiner. You can download them from [here](https://github.com/happywithxpl/GeneMiner-Test).



## GeneMiner-cmd

- Mining single target gene from genome skimming.

```shell
geneminer.py -1 skimming_data1.fq.gz -2 skimming_data2.fq.gz -rtfa ITS.fasta -o ITS_out
```

`-1`,`-2`:  the input files with paired-end reads, given in FASTQ format

`-rtfa`:  reference sequences in fasta format

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



- Mining multiple target genes from genome skimming and evaluating the accuracy of the assembly results.

```shell
geneminer.py -1 skimming_data1.fq.gz  -2 skimming_data2.fq.gz -rtfa Ref_cp -min 300 -max 5000 -limit_count 3 -t 4 -bn 20 -o cp_verify_out3
```

`min`: The minimum length of recovered target genes to be retained [default = 0]

`max`: The maximum length of recovered target genes to be retained [default = 5000]

`limit_count`: Minimum threshold for the k-mer count to remove erroneous and low richness k-mers  [default=auto]

`-t`:   The number of threads

`-bn`:   Specify the bootstrap number. Evaluate the results based on repetitive resampling and base substitution model.



- Mining Angiosperms353 genes from  transcriptome data

```shell
geneminer.py -1 Arabidopsis_thaliana_sim_353_data1.fq.gz -2 Arabidopsis_thaliana_sim_353_data1.fq.gz -rtfa Ref_353 -k1 29 -k2 41 -t 4  -o Angiosperm_353_out4
```

`k1`:  Length of kmer for filtering reads [default = 29]

`k2`:  Length of kmer for assembling reads [default = 41]



## GeneMiner-GUI

GeneMiner-GUI is straightforward, simple to use, and ideal for lightweight users. Considering the memory limit and the excessive time overhead, we recommend you utilize GeneMiner-cmd when you run large-scale data

For GeneMiner-GUI, you must set `Data1`, `Data2` or `Single reads` under the `Data` module、`Ref.(fasta)` or `Ref.(gb)` under the `Reference` module and  `Output Folder` under the `Outout` module as below:

![图片](https://github.com/happywithxpl/GeneMiner-Test/blob/main/run_GeneMiner_GUI.png)



## View results

Now, you can view the results of GeneMiner in the output directory.  

The output directory contains  `reference_database`、`filtered_out `、`assembled_out`、 `GM_results`、`bootstrap_out `  and `results.csv` 

`reference_database <folder>  `: Used to store reference sequences

`filtered_out <folder>`：Used to store filtered reads

`assembled_out <folder>`: Used to store assembled contigs 

`GM_results <folder>`:Used to store all the target genes mined by GeneMiner

`bootstrap_out <folder> ` :Used to store assessment of assembly results, If you have used the `-bn` parameter

`results.csv <file>` :Record various information about the recovered target genes



# User manual

A more complete manual is here: https://github.com/happywithxpl/GeneMiner/blob/main/GeneMiner_User_Guide.md

# Contact

Please check [Geneminer's homepage](https://github.com/happywithxpl/GeneMiner) first. If something is wrong when running GeneMiner, please do not be surprised and report it to us. We usually have quick response to bugs.

- Find Questions & Answers at [GeneMiner Discussions](https://github.com/happywithxpl/GeneMiner/discussions/categories/q-a): **Recommended**
- Report Bugs & Issues at [GeneMiner Issues](https://github.com/happywithxpl/GeneMiner/issues): Please avoid repetitive or irrelevant issues
- QQ group (ID: 78266311): Mainly for mutual help, responses are likely to be not timely

**DO NOT** directly write to us with your questions, instead please post the questions **publicly**, using above platforms  Our emails (xiepulin@scu.edu.cn, 1791173948@qq.com) are only for receiving public question alert and private data (if applied) associated with those public questions. When you send your private data to us, enclose the email with a link where you posted the question. Our only reply emails will be a receiving confirmation, while our answers will be posted in a public place.



# Citation

When you use GeneMiner please cite:

**GeneMiner : a software for extracting phylogenetic markers from next-generation sequencing data**



