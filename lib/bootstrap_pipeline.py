#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2021/9/14 16:19
# @Author  : xiepulin
# @File    : bootstrap_pipeline.py
# @Software: PyCharm
import shutil
import time
import sys
import os
from concurrent import futures
import random
from  collections import  defaultdict
from lib.basic import get_basename,get_fasta_file,get_file_list,is_exist,mylog,cutting_line,get_files,get_identity,get_identity_and_mutate_model
from lib.my_filter import  my_filter_main
from lib.my_assemble import  my_assemble_main

################################################
#################################################


'''
将参考序列拆分为kmer 
数据结构 {kmer:[gene1,gene2]}
'''
def get_bootstrap_hashdict(reference, merSize):
    kmer_dict = defaultdict(list)
    infile = open(reference, 'r', encoding='utf-8', errors='ignore')
    name=""
    seq=""
    my_list=[]
    while True:
        line = infile.readline()
        line = line.strip()
        if (line.startswith('>') or not line) and name:  # 保证最后一条序列能能在保存后退出
            temp = {}
            temp = {name: seq}
            my_list.append(temp)
        if line.startswith('>'):
            name = line[1:]
            seq = ''
        else:
            seq += line
        if not line:
            break

    for i in my_list:
        gene_name=list(i.keys())[0]
        refseq=list(i.values())[0]
        for j in range(0,len(refseq)-merSize+1):
            temp_list, kmer = [], refseq[j:j + merSize]
            # print(kmer,j)
            if kmer in kmer_dict:
                if gene_name not in kmer_dict[kmer]:
                    kmer_dict[kmer].append(gene_name)
            else:
                kmer_dict[kmer]=[gene_name]
    return  kmer_dict

#根据名字获得对应序列
def get_seq_from_name(file,name_list):
    infile = open(file, 'r', encoding='utf-8', errors='ignore')
    name = ""
    seq = ""
    my_list=[]

    while True:
        line = infile.readline()
        line = line.strip()
        if (line.startswith('>') or not line) and (name in name_list):  # 保证最后一条序列能能在保存后退出,保证第一条是有效的
            temp = {}
            temp = {name: seq}
            my_list.append(temp)
        if line.startswith('>'):
            name = line[1:]
            seq = ''
        else:
            seq += line
        if not line:
            break
    return  my_list  #[{gene:seq}]

'''
从fasta序列中获取seq，可指定条目
'''
def get_seq(fasta_file, max_seq_number=100, seq_count_limit=False):
    infile = open(fasta_file, 'r', encoding='utf-8', errors='ignore')
    seq, name = "", ""
    my_list = []
    seq_number = 0
    while True:
        line = infile.readline()
        line = line.strip()

        line_back=line                #防止一整行都是"N"导致程序中断，从而保证退出条件仅为 "文件读完毕"
        line = line.replace("N", "")  # 特定对于scaffold

        if (line.startswith('>') or (not line and not line_back)  ) and name:  # 保证最后一条序列能能在保存后退出
            temp = {}
            temp = {name: seq}
            my_list.append(temp)
            seq_number += 1
        if line.startswith('>'):
            name = line[1:]  # 包含了description
            seq = ''
        else:
            seq += line

        if not line and not line_back: #唯一退出条件 读完文件，而非N导致
            break
        if seq_count_limit and max_seq_number:
            if seq_number >= max_seq_number:
                break

    infile.close()
    return my_list

def mutate(dna,nuc_model):
    dna=bytearray(dna,"utf-8")
    for index in range(len(dna)):
        rate_list=nuc_model[chr(dna[index])]
        rand_num=random.random()
        for j in range(4):
            if rand_num <rate_list[j]:
                dna[index]=ord(("A","T","C","G")[j])
                break
    return  dna.decode("ascii")
def bootstrap_mutate(file,mutate_model,output,bootstrap_number):
    file_name=get_basename(file)
    # 随机替换的对应表
    seq=get_seq(file)
    sequence=list(seq[0].values())[0]
    gene_name=list( seq[0].keys())[0]
    #将自展基因写在用一个文件下
    for i in range(1,bootstrap_number+1):
        bootstrap_out=file_name+".fasta"
        bootstrap_out_path=os.path.join(output,bootstrap_out)
        seq_var = mutate(sequence, mutate_model)
        with open(bootstrap_out_path,"a") as f:
            f.write('>'+gene_name+"_bootstrap_"+str(i)+"\n")
            f.write(seq_var+"\n")



def split_mutate_sequences(input_file,out_dir):
    '''
    :param input_file: 合并后的变异序列
    :param out_dir: 拆分后的序列，输出在一个文件夹下，该文件夹需要提前创建
    :return:
    '''
    infile = open(input_file, 'r', encoding='utf-8', errors='ignore')
    name=""
    seq=""
    my_list=[]
    while True:
        line = infile.readline()
        line = line.strip()
        if (line.startswith('>') or not line) and name:  # 保证最后一条序列能能在保存后退出
            temp = {}
            temp = {name: seq}
            my_list.append(temp)
        if line.startswith('>'):
            name = line[1:]
            seq = ''
        else:
            seq += line
        if not line:
            break
    infile.close()

    if my_list==[]:
        return 0
    number=1
    for i in my_list:
        header= list(i.keys())[0]
        sequence=list(i.values())[0]
        temp=[">"+header+"\n"+sequence]

        name="bootstrap"+"_"+str(number)
        path=os.path.join(out_dir,name+".fasta")
        with open(path,"w") as f:
            f.writelines(temp)
        number+=1

'''
合并fasta文件
'''
def combine_fasta(input_file_list,out_put):
    '''
    :param input_file_list: 输入fasta文件列表
    :param out_put:  输出，如果该文件在某一文件夹下，需要提前创建
    :return:
    '''
    records=[]
    for i in input_file_list:
        infile = open(i, 'r', encoding='utf-8', errors='ignore')
        seq, name = "", ""
        my_list = []
        while True:
            line = infile.readline()
            line = line.strip()
            line = line.replace("N", "")  # 特定对于scaffold
            if (line.startswith('>') or not line) and name:  # 保证最后一条序列能能在保存后退出
                temp = {}
                temp = {name: seq}
                my_list.append(temp)
            if line.startswith('>'):
                name = line[1:]
                seq = ''
            else:
                seq += line
            if not line:
                break
        infile.close()

        if my_list == []:
            continue
        for i in my_list:
            header = list(i.keys())[0]
            sequence = list(i.values())[0]
            temp = [">" + header + "\n" + sequence+"\n"]
            records.extend(temp)
    if records==[]:
        return 0
    with open(out_put,"w") as f:
        for i in records:
            f.writelines(i)

'''
获得geneminer my_assmble的结果序列
'''
def get_geneminer_assembled_result(input_file):
    '''
    :param input_file: assembled_out 路径
    :return: 结果路径列表
    '''
    path1 = os.path.join(input_file, "short_contig")
    path2 = os.path.join(input_file, "contig")
    path3 = os.path.join(input_file, "scaffold")
    path_list = [path1, path2, path3]  # short contig scaffold

    results = []  # 冗余
    assembled_results_list = []
    gene_name_list = []
    for i in path_list:
        if is_exist(i):
            fasta_file = get_fasta_file(i)
            results.extend(fasta_file)
    if results == []:
        return []

    for i in results:  #scaffold质量较差，同时存在contig和scaffold时舍弃
        name = get_basename(i)
        if name not in gene_name_list:
            gene_name_list.append(name)
            assembled_results_list.append(i)
    return  assembled_results_list



'''
自展检测  每次检测一个基因
'''
##################################################
#################################################

class BootstrapPipeLine():
    def __init__(self,configuration_information):
        self.configuration_information=configuration_information
        self.data1=configuration_information["data1"]
        self.data2=configuration_information["data2"]
        self.single=configuration_information["single"]
        self.out_dir = configuration_information["out_dir"]

        self.k1=configuration_information["k1"]
        self.k2=configuration_information["k2"]
        self.step_length=configuration_information["step_length"]
        self.data_size = configuration_information["data_size"]
        self.limit_count = configuration_information["limit_count"]
        self.limit_min_length = configuration_information["limit_min_length"]
        self.limit_max_length = configuration_information["limit_max_length"]
        self.change_seed = configuration_information["change_seed"]
        self.scaffold_or_not = configuration_information["scaffold_or_not"]
        self.max_length=configuration_information["max_length"]
        self.min_length=configuration_information["min_length"]
        self.thread_number=configuration_information["thread_number"]
        self.bootstrap_number=configuration_information["bootstrap_number"]
        self.quiet=configuration_information["quiet"]

        self.reference_database=configuration_information["reference_database"]
        self.filtered_out=configuration_information["filtered_out"]
        self.assembled_out=configuration_information["assembled_out"]
        self.GM_results=configuration_information["GM_results"]
        self.boostrap_out=configuration_information["bootstrap_out"]
        self.results_log=configuration_information["results_log"]
        self.my_software_name=configuration_information["my_software_name"]



        '''
        路径
        '''
        #软件路径
        self.filter_path = configuration_information["filter_path"]
        self.assemble_path = configuration_information["assemble_path"]
        #一级路径
        self.reference_database_path= os.path.join(self.out_dir,self.reference_database)
        self.filtered_out_path=os.path.join(self.out_dir,self.filtered_out)
        self.assembled_out_path=os.path.join(self.out_dir,self.assembled_out)
        self.GM_results_path=os.path.join(self.out_dir,self.GM_results)
        self.boostrap_out_path=os.path.join(self.out_dir,self.boostrap_out)


        #二级路径
        self.boostrap_out_reference_database_path=os.path.join(self.out_dir,self.boostrap_out,self.reference_database) #out_dir\bootstrap_out\reference_database
        self.boostrap_out_filtered_out_path=os.path.join(self.out_dir,self.boostrap_out,self.filtered_out)
        self.boostrap_out_assembled_out_path = os.path.join(self.out_dir, self.boostrap_out, self.assembled_out)
        self.boostrap_out_GM_results_path = os.path.join(self.out_dir, self.boostrap_out, self.GM_results)

        #新方法补充
        #(1)参考
        self.reference_database_split="reference_database_split"  #将合并的参考序列拆开，作为重过滤参考
        #(2)过滤
        self.filtered_out_split ="filtered_out_split"             #重过滤
        self.assembled_out_combined ="assembled_out_combined"    #合并重拼接结果，作为求一致序列的reads序列
        self.assembled_out_split="assembled_out_split"            #归纳重拼接结果，方便计算bootstrap score （均值）
        #(3)拼接
        self.assembled_out_consensus="assembled_out_consensus"    #用基于kmer的方法获得一致序列
        #(4)高置信度结果位置
        self.high_quality_results="high_quality_results"




    def get_mutated_sequence(self,gm_result_path, ref_path):
        kmer=self.k2   #拼接的kmer
        output=self.boostrap_out_reference_database_path
        bootstrap_number=self.bootstrap_number

        # hash 建库
        ref_kmer_dict = get_bootstrap_hashdict(ref_path, kmer)
        gm_kmer_dict = get_bootstrap_hashdict(gm_result_path, kmer)
        kmer_count = defaultdict(int)
        for i in gm_kmer_dict:
            if i in ref_kmer_dict:
                for z in ref_kmer_dict[i]:
                    kmer_count[z] += 1

        # 获得中位数kmercount 参考序列作为平均变异度的计算
        sorted_list = sorted(kmer_count.items(), key=lambda x: x[1], reverse=True)
        name_list = []
        length = len(sorted_list)  # 144
        name_list.append(sorted_list[int(length/2)][0])  # 取kmercount中位数

        if not kmer_count:
            limit_kmer = 17
            kmer = limit_kmer
            ref_kmer_dict_limit = get_bootstrap_hashdict(ref_path, kmer)
            gm_kmer_dict_limit = get_bootstrap_hashdict(gm_result_path, kmer)
            for i in gm_kmer_dict_limit:
                if i in ref_kmer_dict_limit:
                    for z in ref_kmer_dict_limit[i]:
                        kmer_count[z] += 1
        if not kmer_count:
            return 0
        # 获得一致度，变异度
        ref_seq = get_seq_from_name(ref_path, name_list)
        ref_seq = list(ref_seq[0].values())[0]
        gm_seq = get_seq(gm_result_path)
        gm_seq = list(gm_seq[0].values())[0]


        identity,mutate_model = get_identity_and_mutate_model(gm_seq,ref_seq)
        #加载核苷酸变异模型
        bootstrap_mutate(gm_result_path, mutate_model, output, bootstrap_number)

    # 如果GM_results不存在，就没有后续了
    def get_mutated_sequence_parallel(self):
        GM_results_path=self.GM_results_path
        reference_database_path = self.reference_database_path
        thread_number=self.thread_number
        bootstrap_out_path=self.boostrap_out_path
        bootstrap_out_reference_database_path=self.boostrap_out_reference_database_path
        files=get_fasta_file(GM_results_path)
        if files==[]:
            return 0
        if not os.path.isdir(bootstrap_out_path):
            os.mkdir( bootstrap_out_path)
        if not os.path.isdir(bootstrap_out_reference_database_path):
            os.mkdir(bootstrap_out_reference_database_path)
        task_pool = []
        results=[]
        executor = futures.ProcessPoolExecutor(max_workers=thread_number)  #24s
        # executor = futures.ThreadPoolExecutor(max_workers=thread_number) #40s
        for i in files:
            name=get_basename(i)
            ref_path=os.path.join(reference_database_path,name+".fasta")
            task_pool.append(executor.submit(self.get_mutated_sequence,i,ref_path))
        total = len(task_pool)
        number = 1
        for i in task_pool:
            if number < total:
                sys.stdout.write('\r' + "{0:}:{1:>4}/{2}".format("Preparing bootstrap data",number, total))
                sys.stdout.flush()
            else:
                sys.stdout.write('\r' + "{0:}:{1:>4}/{2}".format("Preparing bootstrap data", number, total) + "\n")
                sys.stdout.flush()
            results.append(i.result())
            number = number + 1

        executor.shutdown()
        return  1


    def bootstrap_filter(self):
        data1 = self.data1
        data2 = self.data2
        single = self.single
        boostrap_out_reference_database_path = self.boostrap_out_reference_database_path
        boostrap_out_filtered_out_path=self.boostrap_out_filtered_out_path
        k1 = self.k1
        step_length = self.step_length
        thread_number = self.thread_number
        data_size = self.data_size
        quiet=False  #第一次重过滤允许打印
        filter_path = self.filter_path

        files = get_files(boostrap_out_reference_database_path)
        if files == []:
            return 0

        #内部调用
        filter_configuration_information = {
            "data1": data1, "data2": data2, "single": single,
            "thread_number": thread_number, "k1": k1,
            "out_dir": boostrap_out_filtered_out_path, "step_length": step_length,
            "reference": boostrap_out_reference_database_path, "data_size": data_size,
            "quiet":quiet
        }
        my_filter_main(filter_configuration_information)


    def bootstrap_assemble_parallel(self):
        boostrap_out_reference_database_path = self.boostrap_out_reference_database_path
        boostrap_out_filtered_out_path = self.boostrap_out_filtered_out_path
        boostrap_out_assembled_out_path = self.boostrap_out_assembled_out_path
        thread_number = self.thread_number
        limit_count = self.limit_count
        limit_min_length = self.limit_min_length
        limit_max_length = self.limit_max_length
        scaffold_or_not = False
        change_seed = self.change_seed
        k2 = self.k2
        out_dir=self.out_dir
        quiet=False #拼接 合并拼接可以打印
        assemble_path = self.assemble_path  # assemble.py
        files = get_file_list(boostrap_out_filtered_out_path)
        if files == []:
            return 0
        #内部调用
        assemble_configuration_information = {
            "thread_number": thread_number, "k2": k2, "assembled_out_path": boostrap_out_assembled_out_path,
            "reference": boostrap_out_reference_database_path,
            "input": boostrap_out_filtered_out_path,
            "limit_count": limit_count,
            "limit_min_length": limit_min_length, "limit_max_length": limit_max_length,
            "change_seed": change_seed,
            "scaffold_or_not": scaffold_or_not,
            "out_dir": out_dir,
            "quiet":quiet
        }
        my_assemble_main(assemble_configuration_information)


    def bootstrap_get_results_contig(self):
        path1 = self.boostrap_out_assembled_out_path
        path2 = os.path.join(path1, "short_contig")
        path3 = os.path.join(path1, "contig")
        path4 = os.path.join(path1, "scaffold")

        boostrap_out_GM_results_path = self.boostrap_out_GM_results_path
        if not os.path.isdir(boostrap_out_GM_results_path):
            os.mkdir(boostrap_out_GM_results_path)
        path_list = [path2, path3, path4]  # short contig scaffold

        results = []  # 冗余
        GM_results_list = []
        gene_name_list = []
        for i in path_list:
            if is_exist(i):
                fasta_file = get_fasta_file(i)
                results.extend(fasta_file)
        if results == []:
            return 0
        for i in results:
            name = get_basename(i)
            if name not in gene_name_list:
                gene_name_list.append(name)
                GM_results_list.append(i)

        for i in GM_results_list:
            new_path = os.path.join(boostrap_out_GM_results_path, get_basename(i) + ".fasta")
            shutil.copy(i, new_path)






    #默认剔除N后计算一致度
    def get_bootstrap_information(self,gm_result, bootstrap_result, log_path, bootstrap_number, gene_name):
        gene = gene_name
        bootstrap_number = bootstrap_number
        if is_exist(gm_result) and is_exist(bootstrap_result):
            seq1 = get_seq(gm_result, max_seq_number=1, seq_count_limit=True)
            gm_result_seq = list(seq1[0].values())[0]
            seq2 = get_seq(bootstrap_result, max_seq_number=1, seq_count_limit=True)
            bootstrap_result_seq = list(seq2[0].values())[0]
            score = format(get_identity(gm_result_seq, bootstrap_result_seq)*100,".2f") #format控制位数比round好
        else:
            score = "failed"
        sth = [gene, bootstrap_number, str(score)]
        mylog(log_path, sth)

        return (gene_name,score)

    def get_bootstrap_information_paralle(self):
        GM_results_path=self.GM_results_path
        thread_number=self.thread_number
        boostrap_out_path=self.boostrap_out_path
        boostrap_out_GM_results_path=self.boostrap_out_GM_results_path
        bootstrap_csv_path=os.path.join(boostrap_out_path,"bootstrap.csv")
        bootstrap_number=self.bootstrap_number
        high_quality_results=self.high_quality_results

        files=get_fasta_file(GM_results_path)
        if files==[]:
            return 0
        task_pool=[]
        results=[]
        executor=futures.ProcessPoolExecutor(thread_number)

        header=["gene","bootstrap_number","score"]
        mylog(bootstrap_csv_path,header)
        for i in files:
            name=get_basename(i)
            bootstrap_result=os.path.join(boostrap_out_GM_results_path,name+".fasta")
            task_pool.append(executor.submit(self.get_bootstrap_information,i,bootstrap_result,bootstrap_csv_path,bootstrap_number,name))
        for i in task_pool:
            results.append(i.result())
        executor.shutdown()


        #分拣高质量序列
        high_quality_results_path=os.path.join(boostrap_out_path,high_quality_results)
        if not os.path.isdir(high_quality_results_path):
            os.makedirs(high_quality_results_path)
        bootstrap_score_threshold_value=99.5
        if results:
            for i in results:
                gene_name=i[0]
                if i[1]=="failed":
                    continue
                else:
                    score=float(i[1])
                    if score>=bootstrap_score_threshold_value:
                        old_path=os.path.join(boostrap_out_GM_results_path,gene_name+".fasta")
                        new_path=os.path.join(high_quality_results_path,gene_name+".fasta")
                        shutil.copy(old_path,new_path)






def my_bootstrap_pipeline_main(configuration_information):
    t1=time.time()
    print("")
    cutting_line(" Bootstrap ")
    print("Using GeneMiner...")
    my_bootstrap_pipeline = BootstrapPipeLine(configuration_information)
    flag=my_bootstrap_pipeline.get_mutated_sequence_parallel()
    if flag:
        my_bootstrap_pipeline.bootstrap_filter()
        my_bootstrap_pipeline.bootstrap_assemble_parallel()
        my_bootstrap_pipeline.bootstrap_get_results_contig()
        my_bootstrap_pipeline.get_bootstrap_information_paralle()
        t2=time.time()
        print(" " * 50, flush=True, end='\r')  # 之前的信息比较长，冲刷掉
        used_temp_time = format((t2 - t1), ".2f")
        print("Bootstrap total time used: {}s".format(used_temp_time))
    else:
        t2 = time.time()
        print(" " * 50, flush=True, end='\r')  # 之前的信息比较长，冲刷掉
        used_temp_time = format((t2 - t1), ".2f")
        print("Bootstrap Failed: {}s".format(used_temp_time))




if __name__ == '__main__':
    #较慢的测试
    # data1=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeeee10 重构bootstrap\example\data1.fq"
    # data2 =r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeeee10 重构bootstrap\example\data2.fq"
    # single=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeeee10 重构bootstrap\example\data1.fq"


    #快速测试
    data1 = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeeee10 重构bootstrap\example\data1_100w.fq"
    data2 = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeeee10 重构bootstrap\example\data1_100w.fq"
    single = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeeee10 重构bootstrap\example\data1_100w.fq"
    out_dir = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeeee10 重构bootstrap\example\cp_out"
    rtfa = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeeee10 重构bootstrap\example\cp_gene"
    rtgb = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeee9 重构filter\example\ref_gb\chuanxiong.gb"


    #353大型测试
    # data1 = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeeee10 重构bootstrap\example\A_thaliana_s_2g_single.fastq"
    # data2 = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeeee10 重构bootstrap\example\A_thaliana_s_2g_single.fastq"
    # single = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeeee10 重构bootstrap\example\A_thaliana_s_2g_single.fastq"


    # out_dir = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeeee10 重构bootstrap\example\geneminer_2g_t20_k2931_bn100_004"
    # rtfa = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeeee10 重构bootstrap\example\ref_Brassicaceae"
    # rtgb = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeee9 重构filter\example\ref_gb\chuanxiong.gb"


    k1=29
    k2=31
    data_size='all'
    step_length=4
    limit_count="auto"
    limit_min_length=0.5
    limit_max_length=2
    change_seed=32
    scaffold_or_not=True
    max_length = 50000
    min_length = 0
    thread_number=4
    soft_boundary = 0
    bootstrap_information=[True,10]
    bootstrap=bootstrap_information[0]
    bootstrap_number=bootstrap_information[1]
    quiet=True

    reference_database = "reference_database"
    filtered_out = "filtered_out"
    assembled_out = "assembled_out"
    bootstrap_out = "bootstrap_out"
    GM_results = "GM_results"
    results_log = "results.log"

    filter_path=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeee9 重构filter\lib\my_filter.py"
    assemble_path=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeee9 重构filter\lib\my_assemble.py"
    muscle_path=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeee9 重构filter\lib\muscle3"

    #其他信息
    my_software_name = "GM"
    configuration_information = {"out_dir": out_dir,
                                 "data1": data1, "data2": data2, "single": single,
                                 "rtfa": rtfa, "rtgb": rtgb,
                                 "k1": k1, "k2": k2, "thread_number": thread_number,
                                 "step_length": step_length,
                                 "limit_count": limit_count,
                                 "limit_min_length": limit_min_length,
                                 "limit_max_length": limit_max_length,
                                 "change_seed": change_seed,
                                 "scaffold_or_not":scaffold_or_not,
                                 "max_length": max_length, "min_length": min_length,
                                 "soft_boundary": soft_boundary, "data_size": data_size,
                                 "bootstrap": bootstrap_information[0], "bootstrap_number": bootstrap_information[1],
                                 "reference_database": reference_database,
                                 "filtered_out": filtered_out, "assembled_out": assembled_out,
                                 "bootstrap_out": bootstrap_out,
                                 "GM_results": GM_results,
                                 "results_log": results_log,
                                 "my_software_name": my_software_name,
                                 "filter_path": filter_path, "assemble_path": assemble_path, "muscle_path": muscle_path,
                                 "quiet":quiet
                                 }

    # my_bootstrap_pipeline=BootstrapPipeLine(configuration_information)
    # my_bootstrap_pipeline.get_mutated_sequence_parallel()
    # my_bootstrap_pipeline.bootstrap_filter()
    # my_bootstrap_pipeline.bootstrap_assemble_parallel()
    # my_bootstrap_pipeline.bootstrap_get_results_contig()
    # my_bootstrap_pipeline.get_bootstrap_information_paralle()
    my_bootstrap_pipeline_main(configuration_information)





