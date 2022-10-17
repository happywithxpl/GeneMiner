#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2021/10/28 11:39
# @Author  : xiepulin
# @File    : pack_results.py
# @Software: PyCharm


import csv
from collections import defaultdict
import os
from lib.basic import mylog, is_exist,get_fasta_file
# from basic import mylog, is_exist,get_fasta_file


def parse_csv(path):
    my_info = []
    with open(path, "r", newline='') as f:
        reader = csv.reader(f)
        number = 0
        for row in reader:
            temp = {}
            if number == 0:
                number += 1
                continue
            number += 1
            temp[row[0]] = row[1:]
            my_info.append(temp)

    return my_info


def combine_csv_info(csv_info_list, bn=False):
    my_dict = defaultdict(list)
    my_list = []
    csv_header_length_bn = 13  # ["gene", "k1", "re_k1", "richness", "limit", "seed", "k2", "ref_length", "short_contig_length","contig_length", "scaffold_length", "bootstrap_number", "score"]
    csv_header_length = 11

    for i in csv_info_list:
        for j in i:
            key = list(j.keys())[0]
            value = list(j.values())[0]

            if key not in my_dict:
                my_dict[key] = value
            else:
                my_dict[key].extend(value)

    for key, value in my_dict.items():
        temp = [key]
        temp.extend(value)
        length = len(temp)
        if bn:
            if length < csv_header_length_bn:
                temp.extend(["None"] * (csv_header_length_bn - length))
        else:
            if length < csv_header_length:
                temp.extend(["None"] * (csv_header_length - length))

        my_list.append(temp)

    return my_list


class PackResultsPipeline():
    def __init__(self, configuration_information):
        self.configuration_information = configuration_information
        self.out_dir = configuration_information["out_dir"]
        self.filtered_out = configuration_information["filtered_out"]
        self.assembled_out = configuration_information["assembled_out"]
        self.boostrap_out = configuration_information["bootstrap_out"]
        self.GM_results = configuration_information["GM_results"]

        '''
        路径
        '''
        # 二级路径
        self.filter_csv_path = os.path.join(self.out_dir, self.filtered_out, "filter.csv")
        self.assemble_csv_path = os.path.join(self.out_dir, self.assembled_out, "assemble.csv")
        self.bootstrap_csv_path = os.path.join(self.out_dir, self.boostrap_out, "bootstrap.csv")
        self.results_csv_path = os.path.join(self.out_dir, "results.csv")
        self.GM_results_path=os.path.join(self.out_dir,self.GM_results)

    def get_csv_all_info(self):
        path1 = self.filter_csv_path
        path2 = self.assemble_csv_path
        path3 = self.bootstrap_csv_path
        path4 = self.results_csv_path
        csv_info_list = []
        for i in [path1, path2, path3]:
            if is_exist(i):
                csv_info_list.append(parse_csv(i))
        if csv_info_list != []:
            if is_exist(path3):
                my_list = combine_csv_info(csv_info_list, bn=True)
                header = ["gene", "k1", "re_k1", "richness", "limit", "seed", "k2", "ref_length", "short_contig_length",
                          "contig_length", "scaffold_length", "bootstrap_number", "score"]
                mylog(path4, header)
                mylog(path4, my_list, row=False)

            else:
                my_list = combine_csv_info(csv_info_list, bn=False)
                header = ["gene", "k1", "re_k1", "richness", "limit", "seed", "k2", "ref_length", "short_contig_length",
                          "contig_length", "scaffold_length"]
                mylog(path4, header)
                mylog(path4, my_list, row=False)
        else:
            pass

    def print_recovered_info(self):
        path=self.GM_results_path
        recovered_genes=0
        if is_exist(path):
            files=get_fasta_file(path)
            recovered_genes=len(files)
        else:
            pass
        return  recovered_genes




def my_pack_results_pipeline_main(configuration_information):
    my_pack_results_pipeline = PackResultsPipeline(configuration_information)
    my_pack_results_pipeline.get_csv_all_info() #get csv file
    recovered_genes=my_pack_results_pipeline.print_recovered_info()#get recovered genes's number
    return recovered_genes



if __name__ == '__main__':
    data1 = r"E:\Computer\python\GeneMiner\eeeeeeeee10 重构bootstrap\example\data1_500w.fq"
    data2 = r"E:\Computer\python\GeneMiner\eeeeeeeee10 重构bootstrap\example\data2_500w.fq"
    single = r"E:\Computer\python\GeneMiner\eeeeeeeee10 重构bootstrap\example\data1_500w.fq"

    '''
    for bootstrap
    '''
    # out_dir = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeee9 重构filter\example\demo"
    # rtfa = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeee9 重构filter\example\cp_gene"
    # rtgb = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeee9 重构filter\example\ref_gb\chuanxiong.gb"

    # out_dir = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeee9 重构filter\example\matk_bootstrap"
    # rtfa = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeee9 重构filter\example\matk_ref"
    # rtgb = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeee9 重构filter\example\ref_gb\chuanxiong.gb"


    '''
    for pack_results
    '''
    out_dir = r"E:\Computer\python\GeneMiner\eeeeeeeee10 重构bootstrap\example\cp_out"
    rtfa = r"E:\Computer\python\GeneMiner\eeeeeeeee10 重构bootstrap\example\cp_gene"
    rtgb = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeee9 重构filter\example\ref_gb\chuanxiong.gb"

    k1 = 29
    k2 = 41
    data_size = 'all'
    step_length = 4
    limit_count = "auto"
    limit_min_length = 0.5
    limit_max_length = 2
    change_seed = 32
    scaffold_or_not = True
    max_length = 50000
    min_length = 0
    thread_number = 4
    soft_boundary = 0
    bootstrap_information = [True, 10]
    bootstrap = bootstrap_information[0]
    bootstrap_number = bootstrap_information[1]

    reference_database = "reference_database"
    filtered_out = "filtered_out"
    assembled_out = "assembled_out"
    bootstrap_out = "bootstrap_out"
    GM_results = "GM_results"
    results_log = "results.log"

    filter_path = r"E:\Computer\python\GeneMiner\eeeeeeeee10 重构bootstrap\lib\my_filter.py"
    assemble_path = r"E:\Computer\python\GeneMiner\eeeeeeeee10 重构bootstrap\lib\my_assemble.py"

    # 其他信息
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
                                 "scaffold_or_not": scaffold_or_not,
                                 "max_length": max_length, "min_length": min_length,
                                 "soft_boundary": soft_boundary, "data_size": data_size,
                                 "bootstrap": bootstrap_information[0], "bootstrap_number": bootstrap_information[1],
                                 "reference_database": reference_database,
                                 "filtered_out": filtered_out, "assembled_out": assembled_out,
                                 "bootstrap_out": bootstrap_out,
                                 "GM_results": GM_results,
                                 "results_log": results_log,
                                 "my_software_name": my_software_name,
                                 "filter_path": filter_path, "assemble_path": assemble_path
                                 }


    # my_pack_results_pipeline=PackResultsPipeline(configuration_information)
    # my_pack_results_pipeline.get_csv_all_info()
    # recovered_genes = my_pack_results_pipeline.print_recovered_info()
    recovered_genes=my_pack_results_pipeline_main(configuration_information)
    print(recovered_genes)



