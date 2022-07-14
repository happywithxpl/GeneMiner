#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2022/5/8 10:05
# @Author  : xiepulin
# @File    : geneminer.py
# @Software: PyCharm

import argparse
import sys
import subprocess
import datetime
import os
import logging
import re
import multiprocessing
import signal
import threading
import time
import math

'''
导入第三方库（非标准库）
'''
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import pairwise2

from concurrent import futures

# import PySimpleGUI as sg
my_version = 'Version 1.0b build 20220122'
my_cite = 'Cite: https://github.com/happywithxpl/GeneMiner'
cur_path = os.path.realpath(sys.argv[0])  # 脚本当前路径
father_path = os.path.dirname(cur_path)  # 脚本的父目录
sys.path.append(os.path.join(father_path, "lib"))

from lib.global_var import get_init,set_value,get_value
from lib.basic import *
from lib.verify_parameters import *
from lib.build_reference_database import *
from lib.core_pipeline import *
from lib.bootstrap_pipeline import *
from lib.pack_results import *

get_init()  # 在basic 中已经申明过了

def main(args):
    t1=time.time()
    set_value("my_gui_flag", 1)  # 用于判定GUI是否处于运行状态，1代表运行，0代表没有运行
    '''
    从程序外获得的参数信息
    '''
    data1 = args.data1
    data2 = args.data2
    single = args.single
    target_reference_fa = args.target_reference_fa
    target_reference_gb = args.target_reference_gb
    out_dir = args.out
    thread_number = args.thread
    k1 = args.kmer1
    k2 = args.kmer2
    step_length = args.step_length
    limit_count = args.limit_count
    limit_min_length = args.limit_min_length
    limit_max_length = args.limit_max_length
    scaffold_or_not = args.scaffold

    change_seed = args.change_seed
    soft_boundary = args.soft_boundary
    max_length = args.max
    min_length = args.min
    data_size = args.data_size
    bootstrap_number = args.bootstrap_number
    quiet=False


    data1 = get_absolute(data1)
    data2 = get_absolute(data2)
    single = get_absolute(single)
    out_dir = get_absolute(out_dir)
    target_reference_fa = get_absolute(target_reference_fa)
    target_reference_gb = get_absolute(target_reference_gb)

    # 校验信息
    check_python_version()
    out_dir = check_out_dir(out_dir)  # 输出文件夹检测
    set_value("out_dir", out_dir)

    check_input(data1, data2, single)
    check_reference(target_reference_fa, target_reference_gb)
    thread_number = check_threads_number(thread_number)  # 确定线程数量
    k1 = check_k1(k1)  # 检测k1 filter
    k2 = check_k1(k2)  # 检测k1 filter
    step_length = check_step_length(step_length)  # 步长
    limit_count = check_limit_count(limit_count)  # 限定reads中最低kmercount
    limit_min_length, limit_max_length = check_limit_length(limit_min_length,
                                                            limit_max_length)  # 限定组装contig最短百分比长度，较于ref
    change_seed = check_change_seed(change_seed)  # 种子更换策略最大次数
    scaffold_or_not = check_scaffold(scaffold_or_not)  # 是否做scaffold
    soft_boundary = check_soft_boundary(soft_boundary)  # 检测软边界
    check_max_min_length(max_length, min_length)  # 检测基因最大最小长度
    data_size = check_datasize(data_size)  # 检测数据量 大小
    bootstrap_information = check_bootstrap_parameter(bootstrap_number)  # 自展次数

    # 脚本信息
    cur_path = os.path.realpath(sys.argv[0])  # 脚本当前路径
    cur_path = os.path.dirname(cur_path)  # 脚本的父目录,father_path 覆盖
    filter_path = os.path.join(cur_path, "lib", "my_filter.py")
    assemble_path = os.path.join(cur_path, "lib", "my_assemble.py")
    muscle_path = os.path.join(cur_path, "lib", "muscle3")  # muscle3
    # 文件夹目录
    reference_database = "reference_database"
    filtered_out = "filtered_out"
    assembled_out = "assembled_out"
    bootstrap_out = "bootstrap_out"
    GM_results = "GM_results"

    # 文件目录
    results_log = "results.log"  # bootstarp  将所有Bootstrap结果写在一起的fasta
    bootstrap_data_set = "bootstrap_data_set.fasta"
    bootstrap_concensus = "bootstrap_concensus.fasta"  # 导出共识序列

    # 其他信息
    my_software_name = "GM"
    system = get_platform()

    printinfo = {"project name": out_dir,
                 "data1": data1, "data2": data2, "single": single,
                 "reference (fa)": target_reference_fa, "reference (gb)": target_reference_gb,
                 "k1": k1, "k2": k2, "threads": thread_number,
                 "step length": step_length,
                 "limit count": limit_count,
                 "limit min lenght": limit_min_length,
                 "limit max lenght": limit_max_length,
                 "change seed": change_seed,
                 "max length": max_length, "min length": min_length,
                 "soft boundary": soft_boundary, "data size": data_size,
                 "bootstrap": bootstrap_information[0],
                 "bootstrap number": bootstrap_information[1],
                 }
    configuration_information = {"out_dir": out_dir,
                                 "data1": data1, "data2": data2, "single": single,
                                 "rtfa": target_reference_fa, "rtgb": target_reference_gb,
                                 "k1": k1, "k2": k2, "thread_number": thread_number,
                                 "step_length": step_length,
                                 "limit_count": limit_count,
                                 "limit_min_length": limit_min_length,
                                 "limit_max_length": limit_max_length,
                                 "scaffold_or_not": scaffold_or_not,
                                 "change_seed": change_seed,
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
                                 "system": system,
                                 "filter_path": filter_path, "assemble_path": assemble_path, "muscle_path": muscle_path,
                                 "quiet":quiet

                                 }

    print_parameter_information(printinfo)

    # 构建参考序列基因库
    my_bulid_reference_database_pipeline(configuration_information)
    # 核心流程 过滤 拼接 校验
    my_core_pipeline = CorePipeLine(configuration_information)
    my_core_pipeline.filter_pipeline()
    my_core_pipeline.assemble_pipeline()
    my_core_pipeline.get_results_contig()

    # 自展检测
    if bootstrap_number:
        my_bootstrap_pipeline_main(configuration_information)

    my_pack_results_pipeline_main(configuration_information)

    t2=time.time()

    whole_time=format(t2-t1,"6.2f")

    print("whole time: {}s".format(whole_time))





if __name__ == "__main__":
    # signal.signal(signal.SIGINT, signal_handler)  # 检测函数终止退出（ctrl+c） #必须放在主程序中
    # multiprocessing.freeze_support()  # windows上Pyinstaller打包多进程程序需要添加特殊指令
    # set_value("my_gui_flag", 0)  # 用于判定脚本是否跑完，还可以防止run双击覆盖事件
    parser = argparse.ArgumentParser(usage="%(prog)s <-1 -2|-s>  <-rtfa|rtgb>  <-o>  [options]",
                                     description="GeneMiner: a software for extracting phylogenetic markers from next generation sequencing data\n"
                                                 "Version: 1.0.0\n"
                                                 "Copyright (C) 2022 Pulin Xie\n"
                                                 "Please contact <xiepulin@stu.edu.scu.cn> if you have any questions",

                                     formatter_class=argparse.RawTextHelpFormatter,
                                     # help信息中会自动取消掉换行符和空格，argparse.RawTextHelpFormatter
                                     # 可以将help信息恢复为原始的文本信

                                     # epilog="xiepulin", #参数说明之后，显示程序的其他说明
                                     # add_help=False   #禁用帮助信息
                                     )
    # 原始数据输入部分
    basic_option_group = parser.add_argument_group(title="Basic option")
    basic_option_group.add_argument("-1", dest="data1",
                                    help="One end of the paired-end reads, support fastq format", metavar="")
    basic_option_group.add_argument("-2", dest="data2",
                                    help="Another end of the paired-end reads, support fastq format", metavar="")

    basic_option_group.add_argument("-s", "--single", dest="single",
                                    help="Single reads, support fastq format", metavar="")
    basic_option_group.add_argument("-o", "--out", dest="out", help="Specify the result folder ",
                                    metavar="", required=True)

    basic_option_group.add_argument("-rtfa", dest="target_reference_fa",
                                    help="References of target genes, only support fasta format", metavar="<file|dir>")
    basic_option_group.add_argument("-rtgb", dest="target_reference_gb",
                                    help="References of target genes, only support GenBank format",
                                    metavar="<file|dir>")

    # 高级参数部分
    advanced_option_group = parser.add_argument_group(title="Advanced option")

    advanced_option_group.add_argument("-k1", "--kmer1", dest="kmer1",
                                       help="Specify the size of the wordsize to filter reads  [default = 29]",
                                       default=29,
                                       type=int, metavar="")
    advanced_option_group.add_argument("-k2", "--kmer2", dest="kmer2",
                                       help="Specify the size of the kmer to assemble reads  [default = 31]",
                                       default=31,
                                       type=int, metavar="")

    advanced_option_group.add_argument("-d", "--data", dest="data_size",
                                       help="Specifies the number of reads to reduce raw data. If you want to use all the data, you can set as 'all' [default = 'all']",
                                       default='all', metavar="")

    advanced_option_group.add_argument("-step_length", metavar="", dest="step_length", type=int,
                                       help="the length of the sliding window on the reads [default = 4]", default=4)
    advanced_option_group.add_argument('-limit_count', metavar='', dest='limit_count',
                                       help='''limit of kmer count [default=auto]''', required=False,
                                       default='auto')
    advanced_option_group.add_argument('-limit_min_length', metavar='', dest='limit_min_length', type=float,
                                       help='''limit of contig length''',
                                       required=False, default=1)
    advanced_option_group.add_argument('-limit_max_length', metavar='', dest='limit_max_length', type=float,
                                       help='''limit of contig length''',
                                       required=False, default=2)

    advanced_option_group.add_argument("-change_seed", metavar="", dest="change_seed", type=int,
                                       help='''times of changing seed [default = 32]''', required=False,
                                       default=32)
    advanced_option_group.add_argument('-scaffold', metavar="", dest="scaffold", type=str, help='''make scaffold''',
                                       default=False)

    advanced_option_group.add_argument("-max", dest="max", help="The maximum length of contigs [default = 5000]",
                                       default=5000,
                                       type=int, metavar="")
    advanced_option_group.add_argument("-min", dest="min", help="The minimum length of contigs [default = 300]",
                                       default=300,
                                       type=int, metavar="")
    advanced_option_group.add_argument("-t", "--thread",
                                       help="The number of threads [default = 'auto']",
                                       default="auto", metavar="")
    advanced_option_group.add_argument("-b", "--boundary", dest="soft_boundary",
                                       help="Extend the length to both sides of the gene while extracting genes from Genbank file [default = 75]",
                                       default=75, type=int, metavar="")

    advanced_option_group.add_argument("-bn", "--bootstrap", dest="bootstrap_number", type=int,
                                       help="Specify the bootstrap number. Evaluate the results based on the bootstrap method",
                                       metavar="")
    args = parser.parse_args()

    main(args)
    # geneminer_GUI()



