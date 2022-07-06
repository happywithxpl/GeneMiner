

# GeneMiner

​	感谢您选择GeneMiner(GM)。本文档将帮助您学会如何使用GeneMiner。请注意，GeneMiner仍处于测试阶段，因此本文档可能会在未来的版本升级中更新。

​	如果您有任何问题，可以联系**xiepulin@stu.scu.edu.cn**  或者 15162176893@163.com

[TOC]

# 1. 简介

​		GeneMiner(Gene-Miner，基因矿工)是一款用于从二代测序(NGS)数据中挖掘基因的软件，能够从极低质量和深度的源数据中获取高质量的特定目标基因。例如从浅层基因组测序数据中准确提取叶绿体/线粒体的全部或者部分基因、核基因组中高度重复区 (如*nr*DNA)等；从转录组测序数据中提取大量的单拷贝系统发育标记；从宏基因组数据中获取特定微生物的环境响应基因等。可广泛应用于系统发育与进化研究、海关检验检疫、特定功能基因的挖掘等，在降低测序成本，扩大基因选择方面具有显著优势。GeneMiner的结果十分准确，在真实的实验验证中达到或接近了一代测序的结果，即便参考序列与目标序列相似度低于90%，依然能够依靠梯度逼近算法获得100%准确的结果。软件创新性的提出了基于自展检测的校验方法，可以在不依赖参考序列的情况下对目标序列进行重复验证，输出更加可靠的一致性序列。基于算法层面的大量优化，GeneMiner的运算速度和内存消耗都非常优越，支持多线程并行，充分调用计算机资源，既能够部署在高性能运算集群上，也可以在普通个人电脑上运行。GeneMiner对用户非常友好，支持Windows、Mac和Linux各种主流的操作系统，用户可以选择命令行界面或图形界面进行使用。

# 2.下载

GeneMiner 在MIT License下是开源的。它通过github存储库分发:https://github.com/sculab/GeneMiner，你可以随时下载最新的的版本。请务必关注github，以保持最新的代码更改。我们对以前版本的代码不提供任何支持!版本号遵循符号x.y.z，其中x随着主要的代码重组而改变，y当添加新特性时发生更改，z伴随着bug修复而发生更改

# 3. 安装

## 3.1 For Linux users

自动安装（推荐）

```shell
tar -zxvf geneminer.tar.gz # 解压缩
cd geneminer
python install.py   #根据脚本提示完成，自动安装依赖并将GeneMiner写入环境变量
```

手动安装（自动安装遇到问题时使用）

```shell
tar -zxvf geneminer.tar.gz # 解压缩
cd geneminer
#也可以根据软件提供的依赖文件安装
pip3 install -r requirements.txt --user

# 配置环境变量
echo "export PATH=\$PATH:$(pwd)" >> ~/.bashrc
source ~/.bashrc
```

关闭并重启终端，检测是否配置成功

```shell
geneminer.py -h # 命令行界面
```

## 3.2 For Mac users

下载对应版本的打包好的图形界面app，双击即可运行。（推荐）
要手动配置命令行和图形界面版本，使用如下命令进行自动安装：

```shell
tar -zxvf geneminer.tar.gz # 解压缩
cd geneminer 
python3 install.py #脚本会自动安装依赖，并将GeneMiner写入环境变量
```

或者使用如下命令手动安装：

```shell
tar -zxvf geneminer.tar.gz # 解压缩
cd geneminer
pip3 install -r requirements.txt --user
# 配置环境变量
# 对于macOS Catalina (10.15) 及其之后的系统：
echo "export PATH=\$PATH:$(pwd)" >> ~/.zshrc
source ~/.zshrc
# 对于macOS Catalina (10.15) 之前的系统：
echo "export PATH=\$PATH:$(pwd)" >> ~/.bash_profile
source ~/.bash_profile
```

关闭并重启终端，检测是否配置成功

```shell
geneminer.py -h # 命令行界面
```

## 3.3 For Windows users

 安装打包好的图形界面应用程序（推荐）

下载GeneMiner对应版本的图形界面应用程序，双击即可运行。GeneMiner windows版maual详见：xxxxxxxxxxxxxxxxxx

# 4.快速入门

​		在对您的测序数据进行挖掘之前，我们建议您先了解自己的测序数据，包括测序方式、深度、质量、数据量大小等等。我们的软件主要适用于Illumina平台所返回的二代测序数据。经大量的测试数据验证，即便对于较低的测序深度(10x以下)，GeneMiner也能从中挖掘到单拷贝核基因(nuclear gene)、叶绿体基因(cp gene)以及线粒体基因(mito gene)等，也可以用于从转录组数据中挖掘单-低拷贝基因，如被子植物353基因。

​		用于练手的数据存放在 GeneMiner/example/

（1）挖掘单个目标基因

```shell
 geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta -o out
```

其中data1.fq.gz和data2.fq.gz是二代测序返回的测序数据，ref.fasta是近缘物种（同属或者同科）的同源基因。
**注意：每一个fasta文件中只能存放同一种基因（可以是不同物种的）**。如以下样例：

```
>GeneA species1
ATCGATCG
>GeneA species2
ATCGATCC
>GeneA species3
ATTGATCC
```

（2）挖掘多个目标基因

可以将多个fasta格式的文件放在同一个文件夹中，便能够批量挖掘基因

```shell
#for fasta format
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa fasta_folder -o out
#for GenBank format
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtgb fasta_folder -o out
```

-rtfa : 指定fasta格式的参考序列

-rtgb ：指定GenBank格式的参考序列

（3）评估可靠性

```shell
geneminer.py -1 data1.fq.gz -2 data2.fq.gz  -rtfa  ref.fasta -bn 100 -o out
```

-bn 自展检测的次数。采用基于Bootstrap思想的方法评估结果的准确性。在不依赖参考序列的情况下对目标序列进行重复验证



# 5.详细使用指南

## 5.1参数解读

使用geneminer -h ，计算机将显示所有选项，后面我将详细介绍它们的使用方法和一些小窍门

```shell
GeneMiner: a software for extracting phylogenetic markers from next generation sequencing data
Version: 1.0.0
Copyright (C) 2022 Pulin Xie
Please contact <xiepulin@stu.edu.scu.cn> if you have any questions

optional arguments:
  -h, --help          show this help message and exit

Basic option:
  -1                  One end of the paired-end reads, support fastq format
  -2                  Another end of the paired-end reads, support fastq format
  -s , --single       Single reads, support fastq format
  -o , --out          Specify the result folder
  -rtfa <file|dir>    References of target genes, only support fasta format
  -rtgb <file|dir>    References of target genes, only support GenBank format

Advanced option:
  -k1 , --k-mer1       Specify the size of the wordsize to filter reads  [default = 29]
  -k2 , --k-mer2       Specify the size of the k-mer to assemble reads  [default = 31]
  -d , --data_size    Specifies the number of reads to reduce raw data. If you want to                         use all the data, you can set as 'all' [default = 'all']
  -step_length        the length of the sliding window on the reads [default = 4]
  -limit_count        limit of k-mer count [default=auto]
  -limit_min_length   limit of contig length
  -limit_max_length   limit of contig length
  -change_seed        times of changing seed [default = 32]
  -scaffold           make scaffold
  -max                The maximum length of contigs [default = 5000]
  -min                The minimum length of contigs [default = 300]
  -t , --thread       The number of threads [default = 'auto']
  -b , --boundary     Extend the length to both sides of the gene while extracting genes 					   from Genbank file [default = 75]
  -bn , --bootstrap   Specify the bootstrap number. Evaluate the results based on the 			              bootstrap method

```



### 5.1.1 基础参数

```shell
-1				
双末端测序数据其中一端的数据，支持fastq/fastq.gz格式。务必保留正确的文件扩展名。
example: geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa  ref.fasta -o out

-2   			
双末端测序数据中另一端的数据，支持fastq/fastq.gz/fastq.bz2格式。务必保留正确的文件扩展名。
example: geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta -o out

-s   			
单端测序数据，支持fastq/fastq.gz格式。务必保留正确的文件扩展名。
example: geneminer.py -s data1.fq.gz -rtfa  ref.fasta -o out

-o , --out            
指定输出文件夹。
example: geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta -o out

-rtfa     <file|dir>  
目标基因参考序列，仅支持fasta格式。可用于寻找自己感兴趣的基因。
example: geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta -o out
同一个fasta格式的文件里，只能存放同一个基因的数据，例子如下：
>species_a ITS
AGCTAGCT
>species_b ITS 
AGCTAGCC
>species_c ITS
AGCTAGCA
>species_d ITS
AGCTAGAA

-rtgb     <file|dir>  		
目标基因参考序列，仅支持GenBank格式。可用于寻找自己感兴趣的基因。
example: geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtgb gb_folder -o out
```



### **5.1.2 高级参数**

```shell
-k1 , --k-mer1  
	指定wordsize的大小。该参数用于筛选与参考序列高度匹配的reads。
	geneminer将长度为L的read拆分为(L-k1+1)条长度为k1的子序列(记为word),如果其中某一个word能与参考序列高度匹配，则将这条read保留。k1的大小设定取决于待挖掘目标基因的物种与提供参考序列的物种之间亲缘关系，亲缘关系越近，待挖掘基因与参考序列相似度越大，k1取值越大。
	k1 默认大小为29，k-mer的取值范围在17~127之间。
-k2 , --k-mer2   
	指定k-mer长度，k-mer是de Bruijn图中节点的长度,根据k-1 mer的overlap可将k-mer组装为contig.一般来说，较小的k-mer有利于序列在较低覆盖度时的延申，但可能引入错误k-mer；较大的k-mer能有效处理重复k-mer对拼接的影响，但可能导致序列延申过程中提前终止。
	k2 默认大小为31，k-mer的取值范围在17~127之间	
example:
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta -n 2000000 -k1 17 -k2 31 -o out


-d , --data     
	指定reads的数目。GeneMiner允许使用原始数据的一部分用于挖掘目标基因，特别是具有中-高拷贝数的目标基因。如果您需要将全部的原始数据作为输入，您可以将-n设置为all
	-n 默认设定为all
example:
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta -n 2000000 -o out
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta -n all -o out


-step_length
	指定reads过滤过程中滑动窗口的长度,默认值为4
example:
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta -step_length 4 -o out


-limit_count
	k-mer频次最低阈值（limit）。
	在序列拼接的过程中，GeneMiner会将过滤后reads拆分为长度为k1的子序列(k-mer)，并统计这些k-mer出现的次数.低频的k-mer可能由测序误差导致，具有较低的置信度。
	在"auto"模式下，GeneMiner会自动估计目标loci相关的k-mers的频数分布，基于第一个L型峰的底部为每一个目标loci分配一个k-mer频次最低阈值（limit），频次低于最低阈值（limit）的k-mer将被剔除
	limit值的最小值为2，默认值为"auto"
example:
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta -limit_count auto -o out
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta -limit_count 3 -o out


-limit_min_length
	目标序列占参考序列平均长度的最小比率,剔除拼接过程中过短的contigs.
	limit_min_length=目标序列长度/参考序列平均长度,默认值=1
-limit_max_length
	目标序列占参考序列平均长度的最大比率,剔除拼接过程中过长的contigs.
	limit_max_length=目标序列长度/参考序列平均长度,默认值=2
example:
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta -limit_min_length 0.8 -limit_max_length 1.5 -o out


-max 
	指定GeneMiner挖掘出的基因的最大长度，默认为5000bp.
-min    
	指定GeneMiner挖掘出的基因的最小长度，默认为300bp.
example:
geneminer.py  -1 data1.fq.gz  -2 data2.fq.gz -rtfa ref.fasta  -max 3000  -min 500 -o out


-t , --thread  
	指定线程数量,如果不指定的话，软件会根据计算机性能自动选择合适的线程数量.
	默认值为"auto"
example:
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta  -t 8 -o out

-b , --boundary   
	指定软边界的长度。
	当沿着挖掘出的基因的两侧延申的时候，随着延申长度的增加，碱基的准确率越来越低。但这种下降不是断崖式的，而是在某一定长度的缓冲区内逐渐下滑。我们将保留这一段缓冲区，并将这段缓冲区称为软边界。推荐大小为0.5*reads的长度，软边界取值范围在0~200之间。该参数与-rtgb参数同时使用时才生效.
example:
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtgb ref.gb -b 75 -o out

-bn , --bootstrap
	指定迭代的次数
	GeneMiner创新性的提出了基于碱基替代模型和迭代的校验方法，可以有效的评估获取结果的准确性
	[默认值=100]
example:
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rtfa ref.fasta  -bn 100 -o out
```



## 5.2 结果解读

根据用户不同的参数选择，GeneMiner将产生大量的文件。

**reference_database**  \<folder\>: 	

​		 存放参考序列数据集. GeneMiner对用户提供的近缘类群的参考序列进行预处理，并将目标基因写入fasta格式的文件。



**filtered_out**  \<folder\>: 

​		  存放过滤数据集. GeneMiner基于参考序列数据集对目标片段的相关reads进行筛选，将与参考序列高度匹配的reads分配给特定目标片段的过滤数据集中。如果某个特定目标片段的过滤数据集文件大小超过10MB，GeneMiner会将它保存在big_reads中，并对其进行重新过滤.



**assembled_out**  \<folder\>:

​		 GeneMiner利用过滤数据集完成从头组装，并存放组装结果。根据组装完成度的不同，assembled_out下还可以细分为short_contig,contig,scaffold.

(1) short_contig \<folder\>  :

​		存放较短的contigs.  rate=结果序列长度/参考序列平均长度，如果rate<limit_min_length (see 高级参数)，则把该条结果序列放入short_contig 文件夹

(2) contig \<folder\>: 

​		 存放contigs.  rate=结果序列长度/参考序列平均长度，如果rate>=limit_min_length (see 高级参数)，则把该条结果序列放入contig文件夹

(3) scaffold \<folder\>  :

​		 存放scaffold.

​		NOTE：This folder is only generated when the ­ -scaffold  option is used



**GM_results** \<folder\>:

​		GM_results是最重要的文件夹，该文件夹下存放着GeneMiner挖掘出的所有可用于系统发育研究的目标片段。



**bootstrap_out** \<folder\>:

​		GeneMiner基于碱基替代模型和迭代的方法对挖掘出的目标片段进行校验，并保存评估结果。bootstrap_out文件夹下包括：

(1)reference_database \<folder\>:  

​		存放变异的参考序列数据集。首次获取目标序列后，与参考序列进行比对，得到变异率v。基于参考序列的碱基组成建立碱基替代模型，将目标序列的以变异率v和碱基替代模型进行重采样，获取变异率同样为v的ref~1~, ref~2~.....ref~n~ 

(2)filtered_out \<folder\>:  

​		将变异的参考序列数据集作为新的参考序列数据集，调用相同的流程，获得新的过滤数据集

(3)assembled_out \<folder\>:

​		GeneMiner利用新的过滤数据集完成从头组装，并存放新的组装结果

(4)GM_results \<folder\>:

​		该文件夹下存放着GeneMiner挖掘出的所有可用于系统发育研究的新的目标片段。

(5)high_quality_results \<folder\>:

​		GeneMiner对两次生成的目标片段进行比对，并对目标片段赋予一个评估得分score=identity*100，identity 代表两次生成的目标片段之间的一致度。得分大于99的目标片段将被存放在high_quality_results文件夹中

(6)bootstrap.csv \<folder\>:

​		校验信息记录表

NOTE：This folder is only generated when the ­ -bn/--bootstrap option is used

**results.csv**  \<file\>:

​		This file consists of comma-separated columns containing various information on each target  

sequence found. The file can be easily imported into programs such as Excel. The contents of 

the columns (from left to right) are explained in this table:

| **Column**          | **Description**                                            |
| ------------------- | ---------------------------------------------------------- |
| gene                | 基因名                                                     |
| k1                  | wordsize的大小 (see Advanced parameters：-k1).             |
| re_k1               | 重过滤后wordsize的大小                                     |
| richness            | 基于过滤数据集计算的目标序列大致的测序深度                 |
| limit               | k-mer频次最低阈值  (see Advanced parameters：-limit_count) |
| seed                | GeneMiner组装过程中最终使用的种子片段                      |
| k2                  | k-mer的大小(see Advanced parameters：-k2).                 |
| ref_length          | 参考序列长度（中位数）                                     |
| short_contig_length | 短目标contig的长度                                         |
| contig_length       | 目标contig的长度                                           |
| scaffold_length     | 目标scaffold的长度                                         |
| bootstrap_number    | 迭代次数 (see Advanced parameters：-bn).                   |
| score               | 目标序列评估得分                                           |

NOTE: bootstrap_number and score are  only printed when the ­ -bn/--bootstrap  option is used

​		













## 5.3 例子

### 5.3.1 提取叶绿体基因：

​		当有较为充足的数据量，近缘的参考序列时，本软件几乎能从浅层基因组数据中提取所有的叶绿体基因同时在系统发育研究中.这为叶绿体组装不成环提供了另一种解题思路。

```shell
#example
geneminer.py -1 data1.fq.gz -2 data2.fq.gz  -rmito ref_cp.gb -b 0 -max 3000 -min 300 -o out
```

### 5.3.2 提取线粒体基因：

​		线粒体基因挖掘难度往往会比叶绿体基因大很多，其一般原因在于：原始数据本身未包含太多线粒体基因，线粒体基因变异较大。针对这种情况，您可以适当增加原始数据量大小以及选择更近源的参考序列.

```shell
#example
geneminer.py -1 data1.fq.gz -2 data2.fq.gz -rmito mito.gb  -n 15000000 -o out
```

### 5.3.3 提取核基因：

​        在测序深度较低的测序数据中，GeneMiner也能挖掘出中高拷贝数基因，如rDNAs

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

# 6.方法

​		Geneminer结合了有参拼接和无参拼接的优势，针对短片段的系统发育标记提取，首先基于近缘的参考序列，采用K-mer过滤的方法获取可靠的reads，再采用a weighted seed extension algorithm based on de Brujin Grap方法通过装配短读基因组数据集的小目标区域的方法挖掘目标片段，最后基于碱基替代模型的迭代方法对结果进行校验

GeneMiner 核心流程主要分为三个步骤：基于权重模型和自动种子选取的assembly,基于权重模型和自动种子选取的assembly,基于碱基替代模型和迭代的校验



图片





## **6.1数据过滤：**

​		在这一步骤中，GeneMiner基于近缘类群的参考列对目标片段的相关reads进行筛选。该过程由两阶段组成：建立哈希表和过滤reads。（此处可以补充一句K-mer的定义）在建立哈希表时，我们将参考序列拆分为K-mers的同时记录这些子序列的位置信息、出现的次数和对应的参考序列名，并存储为hash表，K-mer的总数量为![img](C:/Users/xiepulin/AppData/Roaming/Typora/draftsRecover/picture/clip_image002.gif),其中Li为第i条参考序列的长度，参考序列的总数量为n，kf由用户设定。在过滤reads时，原始输入为fastq格式的下一代测序(来自Illumina、Roche-454、ABI或其他测序平台）文件，对于数据中的每条reads，也拆分成![img](C:/Users/xiepulin/AppData/Roaming/Typora/draftsRecover/picture/clip_image004.gif)条K-mer，其中，l为测序的读长，kf的大小与参考序列的设定一致。如果某条read的子序列在参考序列哈希表中出现，则把该条read保留并分配给特定目标片段的过滤数据集中



# 6.2 拼接

​		Geneminer developed a weighted seed extension algorithm based on de Brujin Graph。其大致流程如下：

(Ⅰ)make kmers，把过滤出的reads拆分为k-mers，建立k-mer集合T，此处采用的k与过滤时采用的kf可能不同，计为ka。

(Ⅱ)remove low quality kmers，程序会自动估计目标loci相关的k-mers的频数分布，基于第一个L型峰的底部为每一个目标loci分配一个k-mer频次最低阈值（limit），频次低于最低阈值（limit）的k-mer将被剔除。

(Ⅲ)choose seed，将参考序列使用ka建立k-mer集合R，其中高频出现的k-mer会被视为假定的保守区域，如果该条k-mer也在集合T中出现时，则保留为候选种子。GeneMiner会自动更换候选种子以达到最优拼接效果。

(Ⅳ)seed extend，将每个k-mer作为一个节点，并根据其出现的频次和在参考序列中的位置信息赋予一个权重分：score=log(x1+1) * x2 ,x1代表kmer在过滤后reads中出现的次数，x2代表kmer在参考序列中出现的次数，得分高的候选种子会优先用于拼接。以候选种子作为起始节点，以与其存在k-1碱基的重叠的节点作为相邻节点，其中权重得分高的节点将被优先考虑，并重复加入新的节点的过程，直至无法找到满足要求的新节点。累计权重得分最高的重叠群（contig）将作为输出结果。



## 6.3自展检测

​	GeneMiner创新性的提出了基于碱基替代模型和迭代的校验方法，可以有效的评估获取结果的准确性。

(Ⅰ)首次获取目标序列后，与参考序列进行比对，得到变异率v。

(Ⅱ)基于参考序列的碱基组成建立碱基替代模型，将目标序列的以变异率v和碱基替代模型进行重采样，获取变异率同样为v的ref~1~, ref~2~.....ref~n~ 

(Ⅲ)使用ref~1~, ref~2.~....ref~n~作为新的参考序列重新运行整个流程，得到新的目标序列

(Ⅳ)两次结果的一致度将用于评估其准确性。





# 7. 常见问题

## 7.1 如何验证结果的可靠性？

​		我们的软件设置了基于自展检测的方法对结果进行检验，如果您对您的结果产生了质疑，可以使用-bn选项，并设置自展次数，但自展检验需要反复的调用GeneMiner内部脚本，因此，建议您设置一个合适的值，以免耗费过多时间。同时也可以利用NCBI上的二代测序数据或者真实的一代数据作为补充验证

## 7.2 没有得到结果怎么办？

​		影响基因挖掘成功率的的主要原因大概有如下方面：第一，原始数据质量的好坏。相较于raw_data，clean_data的效果更好。第二，数据量的大小，数据量太小可能导致基因丰度不够无法拼接；数据量太大可能导致多种拼接情况的出现。我们推荐将-n设置为100w~1000w。第三，参考序列的选择。选择近缘属或者同属不同种的序列作为参考序列，往往能够得到较好的结果。同时，当存在多条序列作为参考序列的时候，参考序列之间的差异度不要太大。

## 7.3  -n 设置多少合适？

​		一般来说使用默认参数n=10000000，即2500000条reads即可以很好的满足需要。但是对于单拷贝或者低拷贝序列，可以适当的增加数据量。我们不建议将全部数据输入，因为不仅会大大降低速度，而且由于拼接出的情形增多，反而得不到好的结果。

# 8.引用

当你使用GeneMiner的时候请引用：

GeneMiner : a software for extracting phylogenetic markers from next generation sequencing data

