#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2021/10/20 16:48
# @Author  : xiepulin
# @File    : core_pipeline.py
# @Software: PyCharm
import os
import shutil
from basic import get_file_list,is_exist,get_fasta_file,get_basename
from my_filter import  my_filter_main
from my_assemble import  my_assemble_main
class CorePipeLine():
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
        self.limit_count=configuration_information["limit_count"]
        self.limit_min_length=configuration_information["limit_min_length"]
        self.limit_max_length = configuration_information["limit_max_length"]
        self.change_seed=configuration_information["change_seed"]
        self.scaffold_or_not=configuration_information["scaffold_or_not"]
        self.max_length=configuration_information["max_length"]
        self.min_length=configuration_information["min_length"]
        self.thread_number=configuration_information["thread_number"]
        self.quiet=configuration_information["quiet"]


        self.reference_database=configuration_information["reference_database"]
        self.filtered_out=configuration_information["filtered_out"]
        self.assembled_out=configuration_information["assembled_out"]
        self.GM_results=configuration_information["GM_results"]
        self.results_log=configuration_information["results_log"]
        self.my_software_name=configuration_information["my_software_name"]


        #路径
        self.filter_path = configuration_information["filter_path"]
        self.assemble_path = configuration_information["assemble_path"]
        self.muscle_path = configuration_information["muscle_path"]
        self.reference_database_path= os.path.join(self.out_dir,self.reference_database)
        self.filtered_out_path=os.path.join(self.out_dir,self.filtered_out)
        self.assembled_out_path=os.path.join(self.out_dir,self.assembled_out)
        self.GM_results_path=os.path.join(self.out_dir,self.GM_results)


    def filter_pipeline(self):
        data1=self.data1
        data2=self.data2
        single=self.single
        reference_database_path = self.reference_database_path
        filtered_out_path = self.filtered_out_path
        k1=self.k1
        step_length=self.step_length
        thread_number = self.thread_number
        data_size=self.data_size
        filter_path=self.filter_path
        quiet=self.quiet

        filter_configuration_information={
            "data1":data1,"data2":data2,"single":single,
            "thread_number":thread_number,"k1":k1,
            "out_dir":filtered_out_path,"step_length":step_length,
            "reference":reference_database_path,"data_size":data_size,
            "quiet":quiet
        }
        my_filter_main(filter_configuration_information)


    def assemble_pipeline(self):
        reference_database_path = self.reference_database_path
        filtered_out_path = self.filtered_out_path
        assembled_out_path = self.assembled_out_path
        thread_number = self.thread_number
        limit_count = self.limit_count
        limit_min_length = self.limit_min_length
        limit_max_length = self.limit_max_length
        scaffold_or_not = self.scaffold_or_not
        change_seed = self.change_seed
        k2 = self.k2
        out_dir = self.out_dir
        quiet=self.quiet
        assemble_path = self.assemble_path  # assemble.py


        files = get_file_list(filtered_out_path)
        if files == []:
            return 0

        assemble_configuration_information = {
            "thread_number": thread_number, "k2": k2, "assembled_out_path": assembled_out_path,
            "reference": reference_database_path,
            "input": filtered_out_path,
            "limit_count": limit_count,
            "limit_min_length": limit_min_length, "limit_max_length": limit_max_length,
            "change_seed": change_seed,
            "scaffold_or_not": scaffold_or_not,
            "out_dir": out_dir,
            "quiet":quiet
        }
        my_assemble_main(assemble_configuration_information)



    def check_contig_pipeline(self):
        pass


    def get_results_contig(self):
        path1=self.assembled_out_path
        path2=os.path.join(path1,"short_contig")
        path3=os.path.join(path1,"contig")
        path4=os.path.join(path1,"scaffold")

        GM_results_path= self.GM_results_path
        if not os.path.isdir(GM_results_path):
            os.mkdir(GM_results_path)
        path_list=[path2,path3,path4] #short contig scaffold
        results=[] #冗余
        GM_results_list=[]
        gene_name_list=[]
        for i in path_list:
            if is_exist(i):
                fasta_file=get_fasta_file(i)
                results.extend(fasta_file)

        if results==[]:
            return 0
        #contig和scaffold同时存在时  去除scaffold
        for i in results:
            name=get_basename(i)
            if name not in gene_name_list:
                gene_name_list.append(name)
                GM_results_list.append(i)

        for i in GM_results_list:
            new_path=os.path.join(GM_results_path,get_basename(i)+".fasta")
            shutil.copy(i,new_path)


if __name__ == '__main__':
    data1 = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeeee10 重构bootstrap\example\data1_100w.fq"
    data2 = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeeee10 重构bootstrap\example\data1_100w.fq"
    single = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeeee10 重构bootstrap\example\data1_100w.fq"

    out_dir = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeeee10 重构bootstrap\example\cp_out"
    rtfa = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeeee10 重构bootstrap\example\cp_gene"
    rtgb = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeee9 重构filter\example\ref_gb\chuanxiong.gb"
    # rtgb = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeee9 重构filter\example\ref_gb"

    k1=17
    k2=31
    data_size='all'
    step_length=4
    limit_count=-1
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
    quiet=False



    reference_database = "reference_database"
    filtered_out = "filtered_out"
    assembled_out = "assembled_out"
    bootstrap_out = "bootstrap_out"
    GM_results = "GM_results"
    results_log = "results.log"
    bootstrap_data_set = "bootstrap_data_set.fasta"
    bootstrap_concensus = "bootstrap_concensus.fasta"


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
                                 "bootstrap_data_set": bootstrap_data_set,
                                 "bootstrap_concensus": bootstrap_concensus,
                                 "my_software_name": my_software_name,
                                 "filter_path": filter_path, "assemble_path": assemble_path, "muscle_path": muscle_path,
                                 "quiet":quiet

                                 }
    # print(configuration_information)
    my_core_pipeline=CorePipeLine(configuration_information)
    my_core_pipeline.filter_pipeline()
    my_core_pipeline.assemble_pipeline()
    my_core_pipeline.get_results_contig()

