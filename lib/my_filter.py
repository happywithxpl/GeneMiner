#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2022/3/6 15:58
# @Author  : xiepulin
# @File    : my_filter.py
# @Software: PyCharm

import sys
import argparse
from Bio import SeqIO
import gzip
import os
import re
import time
from collections import defaultdict
import platform
import gc
import multiprocessing
from multiprocessing import Process
import shutil
from concurrent.futures import ProcessPoolExecutor
import csv


def get_absolute(path):
    if path == None:
        return None
    else:
        if os.path.isabs(path):
            return path
        else:
            path = os.path.abspath(path)
            return path


def get_file_physical_size(path):
    size = os.path.getsize(path)
    return size


def get_platform():
    current_os = platform.system().lower()  # linux windows darwin
    return current_os


'''
返回参考序列列表   参考基因组可能是文件，也可能是文件夹
'''


def get_file_list(path):
    file_list = []
    if os.path.isdir(path):
        files = get_files(path)
        for i in files:
            if os.path.getsize(i) == 0:
                pass
            else:
                file_list.append(i)
    elif os.path.isfile(path):
        size = os.path.getsize(path)
        if size == 0:
            pass
        else:
            file_list.append(path)
    else:
        pass
    return file_list


'''
获取目录下所有文件
'''


def get_files(ref):
    file_path_list = []
    for root, dirs, files in os.walk(ref):
        for file in files:
            file_path = os.path.join(root, file)
            file_path_list.append(file_path)
    return file_path_list


'''
获得目标序列文件大小
'''


def get_size(path):
    size = 0
    files = get_file_list(path)
    for i in files:
        size += os.path.getsize(i)

    size = size
    return size


'''
文件是否真实存在
'''


def is_exist(file):
    if os.path.isfile(file):
        if os.path.getsize(file) > 0:
            flag = 1
        else:
            flag = 0
    elif os.path.isdir(file):
        files = get_files(file)
        if files == []:
            flag = 0
        else:
            flag = 1
            for i in files:
                if os.path.getsize(i) > 0:
                    continue
                else:
                    flag = 0
                    break
    else:
        flag = 0

    return flag


'''
获取fasta文件的文件名 即基因名
'''


def get_basename(file):
    if is_exist(file):
        basename = os.path.basename(file)

        if ".fasta" in basename:
            basename = basename.split(".fasta")[0]
        elif ".fas" in basename:
            basename = basename.split(".fas")[0]
        elif ".fa" in basename:
            basename = basename.split(".fa")[0]
        else:
            basename = basename
        return basename


# 记录日志信息
def mylog(file_path, sth, row=True):
    if row == True:
        with open(file_path, "a",
                  newline='') as f:  # 不同操作系统换行符不一样 linux:\n  windows:\r\n  mac \r newline参数可以统一换行符  否则中间会出现默认空白行
            writer = csv.writer(f)
            writer.writerow(sth)
    else:
        with open(file_path, "a",
                  newline='') as f:  # 不同操作系统换行符不一样 linux:\n  windows:\r\n  mac \r newline参数可以统一换行符  否则中间会出现默认空白行
            writer = csv.writer(f)
            writer.writerows(sth)


# 分割线 文字内容不得大于分割线
def cutting_line(message=""):
    element = ["="]
    cutting_line_number = 60
    element_all = element * cutting_line_number
    length = len(message)
    start = int((cutting_line_number - length) / 2)

    for i in range(len(message)):
        element_all[i + start] = message[i]
    element_all = "".join(element_all)
    print(element_all)
    return element_all

'''
控制信息的打印状况与保存
'''
def My_Log_Recorder(message,keep_quiet=True,path="",stage="",printout_flag=True,stored_flag=False,make_a_newline=True,use_cutting_line=False):
    '''
    :param message: 需要打印/记录的信息
    :param path:  log路径
    :param stage:  当前阶段（注释信息）
    :param printout_flag: 是否打印
    :param stored_flag:  是否保存记录
    :param make_a_newline:  是否换行
    :param use_cutting_line:  使用函数，打印一些标志性的东西
    :return:
    '''

    if keep_quiet==True :
        return 0
    if printout_flag==True:
        if make_a_newline==True:
            print(message,flush=True)
        else:
            print(message,flush=True,end='\r')
    if stored_flag==True and path:
        my_time=time.localtime()
        time_header=str(my_time.tm_year)+"-"+str(my_time.tm_mon)+"-"+str(my_time.tm_mday)+"-"+str(my_time.tm_hour)+"-"+str(my_time.tm_min)+"-"+str(my_time.tm_sec)
        stored_message=time_header+" "+stage+":"+message+"\n"
        with open(path,"a") as f:
            f.write(stored_message)
    if use_cutting_line:
        cutting_line(message)


# 获取反向互补
def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return str("".join(complement.get(base, base) for base in reversed(seq)))


# 简化反向互补
def reverse_complement_limit(seq):
    return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


# 创建哈西字典,将reference存入字典,并记录位置信息
def gethashdict(reference, merSize, get_rc=False, pos=False, print_info=False,keep_quiet=False):
    files = get_file_list(reference)
    if files == []:
        return 0

    gene_number = 0
    kmer_dict = defaultdict(list)
    for file in files:
        if os.path.isfile(file) == False:
            continue
        infile = open(file, 'r', encoding='utf-8', errors='ignore')
        file_name = get_basename(file)
        infile.readline()
        for line in infile:
            temp_str = []
            while line and line[0] != '>':
                temp_str.append(line)
                line = infile.readline()
            gene_number += 1
            refseq = ''.join(filter(str.isalnum, ''.join(temp_str).upper()))  # 将一个文件下所有seq连在一起
            for j in range(0, len(refseq) - merSize + 1):
                temp_list, kmer = [], refseq[j:j + merSize]
                if kmer in kmer_dict:
                    temp_list = kmer_dict[kmer]
                    temp_list[0] += 1
                    if pos:
                        temp_list[1] += (j + 1) / len(refseq)  # (j+1)/len(refseq) 记录位置信息，如果反复出现，我们可以理解为在求一个位置的均值 （估计值）
                    if file_name not in temp_list:
                        temp_list.append(file_name)
                else:
                    temp_list = [1, (j + 1) / len(refseq), file_name] if pos else [1, 0, file_name]  # 记录位置信息
                    kmer_dict[sys.intern(kmer)] = temp_list
            if get_rc:
                try:
                    refseq = reverse_complement_limit(refseq)
                except:
                    pass
                for j in range(0, len(refseq) - merSize + 1):
                    temp_list, kmer = [], refseq[j:j + merSize]
                    if kmer in kmer_dict:
                        temp_list = kmer_dict[kmer]
                        temp_list[0] += 1
                        if pos:
                            temp_list[1] -= (j + 1) / len(refseq)
                        if file_name not in temp_list: temp_list.append(file_name)
                    else:
                        temp_list = [1, -(j + 1) / len(refseq), file_name] if pos else [1, 0, file_name]
                        kmer_dict[sys.intern(kmer)] = temp_list  # sys,inter利用了地址信息，保证相同元素出现在同一地址
            if print_info:
                message1=" " * 100
                My_Log_Recorder(message=message1,keep_quiet=keep_quiet,path="",stage="",printout_flag=True,stored_flag=False,make_a_newline=False,use_cutting_line=False)
                message2="Memory:{0},G, Number of Seq:{1}".format(round(6 * sys.getsizeof(kmer_dict) / 1024 / 1024 / 1024, 4),gene_number)
                My_Log_Recorder(message=message2, keep_quiet=keep_quiet,path="", stage="", printout_flag=True, stored_flag=False,make_a_newline=False, use_cutting_line=False)


                # print(" " * 100, flush=True, end='\r')
                # print('Memory:', round(6 * sys.getsizeof(kmer_dict) / 1024 / 1024 / 1024, 4), 'G, Number of Seq',
                #       gene_number, end='\r', flush=True)  # 6倍率
        infile.close()
    return kmer_dict




def filter_reads(_dict, _kmer, stepSize, readseq, get_rc=False):
    tmplist = []
    for j in range(0, len(readseq) - _kmer + 1 - stepSize, stepSize):
        kmer = readseq[j:j + _kmer]
        if kmer in _dict:
            tmplist.extend(_dict.get(kmer)[2:])  # 剔除kmercount 和 位置信息

    if get_rc:  # data1的反向互补并不是data2 最正确的是过滤data1,过滤data2,取交集---------实际上根据data1过滤后的title定位data2也是可以的
        # temp_seq = Seq(readseq).reverse_complement()   #Bio.Seq(seq).reverse_complement()
        # temp_seq = Seq.reverse_complement(readseq)       #Bio.Seq.reverse_complement(seq) 快
        temp_seq = reverse_complement_limit(readseq)  # 最快
        for j in range(0, len(temp_seq) - _kmer + 1 - stepSize, stepSize):
            kmer = temp_seq[j:j + _kmer]
            if kmer in _dict:
                tmplist.extend(_dict.get(kmer)[2:])
    return set(tmplist)


def filter_paired_reads(_dict, _kmer, setpSize, file_1, file_2, out_dir, t_id, t_count, get_rc=False, data_size='all',keep_quiet=False):
    t1, t2, reads_count = time.time(), 0, 0
    infile_1 = gzip.open(file_1, 'r') if file_1[-3:].lower() == ".gz" else open(file_1, 'rb')
    infile_2 = gzip.open(file_2, 'r') if file_2[-3:].lower() == ".gz" else open(file_2, 'rb')
    temp_rec1 = [infile_1.readline(), infile_1.readline(), infile_1.readline(), infile_1.readline()]
    temp_rec2 = [infile_2.readline(), infile_2.readline(), infile_2.readline(), infile_2.readline()]
    data_size = str(data_size)
    data_size = 1000000000000 if data_size.lower() == 'all' else int(data_size)
    for _ in infile_1:
        reads_count += 1
        if reads_count % t_count == t_id:
            for file_name in filter_reads(_dict, _kmer, setpSize, temp_rec1[1].decode('utf-8'),
                                          get_rc):  # 如果过滤成功返回一个去重后的列表，该列表存放着满足过滤条件的基因名（文件名）
                with open(os.path.join(out_dir, file_name + "_" + str(t_id) + "_" + ".fasta"), "a+") as outfile:
                    outfile.writelines(['>', temp_rec1[0].decode('utf-8')[1:], temp_rec1[1].decode('utf-8')])
                    outfile.writelines(['>', temp_rec2[0].decode('utf-8')[1:], temp_rec2[1].decode('utf-8')])
        temp_rec1 = [_, infile_1.readline(), infile_1.readline(), infile_1.readline()]
        temp_rec2 = [infile_2.readline(), infile_2.readline(), infile_2.readline(), infile_2.readline()]
        if reads_count * t_count % 1000000 == 0:
            t2 = time.time()
            t1, t2 = t2, t2 - t1
            # print(" " * 50, flush=True, end='\r')
            # info = "handled{0:>4}m reads, {1:>4}s/m reads".format(reads_count * t_count // 1000000, round(t2, 2))
            # print(info, end='\r', flush=True)
            message1 = " " * 50
            My_Log_Recorder(message=message1, keep_quiet=keep_quiet, path="", stage="", printout_flag=True,
                            stored_flag=False, make_a_newline=False, use_cutting_line=False)
            message2 = "handled{0:>4}m reads, {1:>4}s/m reads".format(reads_count * t_count // 1000000, round(t2, 2))
            My_Log_Recorder(message=message2, keep_quiet=keep_quiet, path="", stage="", printout_flag=True,
                            stored_flag=False, make_a_newline=False, use_cutting_line=False)

        if reads_count * t_count >= data_size:
            break
    infile_1.close()
    infile_2.close()


def filter_single_reads(_dict, _kmer, setpSize, file_1, out_dir, t_id, t_count, get_rc=False, data_size='all',keep_quiet=False):
    t1, t2, reads_count = time.time(), 0, 0
    infile_1 = gzip.open(file_1, 'r') if file_1[-3:].lower() == ".gz" else open(file_1, 'rb')
    temp_rec1 = [infile_1.readline(), infile_1.readline(), infile_1.readline(), infile_1.readline()]
    data_size = str(data_size)
    data_size = 10000000000 if data_size.lower() == 'all' else int(data_size)
    for _ in infile_1:
        reads_count += 1
        if reads_count % t_count == t_id:
            for file_name in filter_reads(_dict, _kmer, setpSize, temp_rec1[1].decode('utf-8'),
                                          get_rc):  # 如果过滤成功返回一个去重后的列表，该列表存放着满足过滤条件的基因名（文件名）
                with open(os.path.join(out_dir, file_name + "_" + str(t_id) + "_" + ".fasta"), "a+") as outfile:
                    outfile.writelines(['>', temp_rec1[0].decode('utf-8')[1:], temp_rec1[1].decode('utf-8')])
        temp_rec1 = [_, infile_1.readline(), infile_1.readline(), infile_1.readline()]
        if reads_count * t_count % 1000000 == 0:
            t2 = time.time()
            t1, t2 = t2, t2 - t1
            # print(" " * 50, flush=True, end='\r')
            # info = "handled{0:>4}m reads, {1:>4}s/m reads".format(reads_count * t_count // 1000000, round(t2, 2))
            # print(info, end='\r', flush=True)
            #
            message1=" " * 50
            My_Log_Recorder(message=message1, keep_quiet=keep_quiet, path="", stage="", printout_flag=True,stored_flag=False, make_a_newline=False, use_cutting_line=False)
            message2="handled{0:>4}m reads, {1:>4}s/m reads".format(reads_count * t_count // 1000000, round(t2, 2))
            My_Log_Recorder(message=message2, keep_quiet=keep_quiet, path="", stage="", printout_flag=True,stored_flag=False, make_a_newline=False, use_cutting_line=False)

        if reads_count * t_count >= data_size:
            break

    infile_1.close()


def filter_paired_reads_single_thread(_dict, _kmer, setpSize, file_1, file_2, out_dir, get_rc=False, data_size='all',keep_quiet=False):
    t1, t2, reads_count = time.time(), 0, 0
    infile_1 = gzip.open(file_1, 'r') if file_1[-3:].lower() == ".gz" else open(file_1, 'rb')
    infile_2 = gzip.open(file_2, 'r') if file_2[-3:].lower() == ".gz" else open(file_2, 'rb')
    temp_rec1 = [infile_1.readline(), infile_1.readline(), infile_1.readline(), infile_1.readline()]
    temp_rec2 = [infile_2.readline(), infile_2.readline(), infile_2.readline(), infile_2.readline()]
    data_size = str(data_size)
    data_size = 1000000000000 if data_size.lower() == 'all' else int(data_size)
    for _ in infile_1:
        reads_count += 1
        for file_name in filter_reads(_dict, _kmer, setpSize, temp_rec1[1].decode('utf-8'),
                                      get_rc):  # 如果过滤成功返回一个去重后的列表，该列表存放着满足过滤条件的基因名（文件名）
            with open(os.path.join(out_dir, file_name + ".fasta"), "a+") as outfile:
                outfile.writelines(['>', temp_rec1[0].decode('utf-8')[1:], temp_rec1[1].decode('utf-8')])
                outfile.writelines(['>', temp_rec2[0].decode('utf-8')[1:], temp_rec2[1].decode('utf-8')])
        temp_rec1 = [_, infile_1.readline(), infile_1.readline(), infile_1.readline()]
        temp_rec2 = [infile_2.readline(), infile_2.readline(), infile_2.readline(), infile_2.readline()]
        if reads_count % 1000000 == 0:
            t2 = time.time()
            t1, t2 = t2, t2 - t1

            # print(" " * 50, flush=True, end='\r')
            # info = "handled{0:>4}m reads, {1:>4}s/m reads".format(reads_count // 1000000, round(t2, 2))
            # print(info, end='\r', flush=True)

            message1 = " " * 50
            My_Log_Recorder(message=message1, keep_quiet=keep_quiet, path="", stage="", printout_flag=True,
                            stored_flag=False, make_a_newline=False, use_cutting_line=False)
            message2 = "handled{0:>4}m reads, {1:>4}s/m reads".format(reads_count // 1000000, round(t2, 2))
            My_Log_Recorder(message=message2, keep_quiet=keep_quiet, path="", stage="", printout_flag=True,
                            stored_flag=False, make_a_newline=False, use_cutting_line=False)

        if reads_count >= data_size:
            break
    infile_1.close()
    infile_2.close()


def filter_single_reads_single_thread(_dict, _kmer, setpSize, file_1, out_dir, get_rc=False, data_size='all',keep_quiet=False):
    t1, t2, reads_count = time.time(), 0, 0
    infile_1 = gzip.open(file_1, 'r') if file_1[-3:].lower() == ".gz" else open(file_1, 'rb')
    temp_rec1 = [infile_1.readline(), infile_1.readline(), infile_1.readline(), infile_1.readline()]
    data_size = str(data_size)
    data_size = 10000000000 if data_size.lower() == 'all' else int(data_size)

    for _ in infile_1:
        reads_count += 1
        for file_name in filter_reads(_dict, _kmer, setpSize, temp_rec1[1].decode('utf-8'),
                                      get_rc):  # 如果过滤成功返回一个去重后的列表，该列表存放着满足过滤条件的基因名（文件名）
            with open(os.path.join(out_dir, file_name + ".fasta"), "a+") as outfile:
                outfile.writelines(['>', temp_rec1[0].decode('utf-8')[1:], temp_rec1[1].decode('utf-8')])
        temp_rec1 = [_, infile_1.readline(), infile_1.readline(), infile_1.readline()]
        if reads_count % 1000000 == 0:
            t2 = time.time()
            t1, t2 = t2, t2 - t1
            # print(" " * 50, flush=True, end='\r')
            # info = "handled{0:>4}m reads, {1:>4}s/m reads".format(reads_count // 1000000, round(t2, 2))
            # print(info, end='\r', flush=True)
            message1 = " " * 50
            My_Log_Recorder(message=message1, keep_quiet=keep_quiet, path="", stage="", printout_flag=True,
                            stored_flag=False, make_a_newline=False, use_cutting_line=False)
            message2 = "handled{0:>4}m reads, {1:>4}s/m reads".format(reads_count // 1000000, round(t2, 2))
            My_Log_Recorder(message=message2, keep_quiet=keep_quiet, path="", stage="", printout_flag=True,
                            stored_flag=False, make_a_newline=False, use_cutting_line=False)
        if reads_count >= data_size:
            break
    infile_1.close()


####################################################33
######################################################
'''
多线程写入合并文件 combine file
'''


##########################################################
#########################################################

def combine_file(path, path_list, gene_name):
    with open(os.path.join(path, gene_name + ".fasta"), "a", newline="\n") as f1:
        for i in path_list:
            with open(i, "r", newline="\n") as f2:
                lines = f2.readlines()
                f1.writelines(lines)

    for i in path_list:
        os.remove(i)
    return gene_name


def combine_file_parallel(path, path_dict, thread):
    if not path_dict:
        return 0
    executor = ProcessPoolExecutor(max_workers=thread)
    task_pool = []
    results = []
    for key, value in path_dict.items():
        task_pool.append(executor.submit(combine_file, path, value, key))

    for i in task_pool:
        results.append(i.result())
    executor.shutdown()


def get_filtered_path_for_combine(path):
    files = get_file_list(path)
    path_dict = defaultdict(list)
    for i in files:
        name = get_basename(i) + ".fasta"
        suffix = re.findall(r"_\d{1,}_.fasta", name)[0]
        name = name.split(suffix)[0]
        if name not in path_dict:
            path_dict[name] = [i]
        else:
            path_dict[name].append(i)
    return path_dict


######################################################################
#######################################################################
'''
re-filter
'''
######################################################################
#######################################################################
'''
1.1获取哪些序列需要重过滤 #第一次的评判标准依据文件大小
'''


def get_re_filter_path(filtered_out, filtered_limit_size):
    files = get_file_list(filtered_out)
    re_filter_files = []
    no_re_filter_files = []

    if files != []:
        for i in files:
            if get_file_physical_size(i) > filtered_limit_size:
                re_filter_files.append(i)
            else:
                no_re_filter_files.append(i)
    return re_filter_files, no_re_filter_files


'''
1.2.1获得ref 平均长度
返回 {filename:length}
'''


def get_ref_info(reference_list_path, ref_length_dict):
    length_all = 0
    ref_count = 0
    for file in reference_list_path:
        file_name = get_basename(file)
        for rec in SeqIO.parse(file, "fasta"):
            length = len(rec.seq)
            length_all += length
            ref_count += 1
        average_length = int(length_all / ref_count)
        ref_length_dict[file_name] = average_length
    return ref_length_dict


'''
1.2.2

获取filtered reads信息 平均长度==》总长度
'''

def get_reads_length(file):
    infile = gzip.open(file, 'rt') if file[-3:].lower() == ".gz" else open(file, 'r')
    temp_rec = [infile.readline(), infile.readline(), infile.readline(), infile.readline()]
    reads_length = []
    line_number = 0
    line_limit = 100
    for _ in infile:
        if line_number < line_limit:
            reads_length.append(len(temp_rec[1].strip()))
            line_number += 1
            temp_rec = [_, infile.readline(), infile.readline(),
                        infile.readline()]  # 每次读取的第一行 如果换成infile.readline(),相当于在此基础上又读了一行
        else:
            break
    infile.close()
    length = int(sum(reads_length) / len(reads_length))

    return length


'''
按照reads平均长度计算filtered_out 总bp数目，不影响结果但比详细计算bp数目快10倍
'''


def get_filtered_reads_bp_number(filtered_out_path_list, reads_length, filtered_reads_whole_bp_dict):
    '''
    :param filtered_out_path_list:  filtered_out大文件结果
    :param reads_length:
    :param filtered_reads_whole_bp_dict:
    :return:
    '''
    for file in filtered_out_path_list:
        file_name = get_basename(file)
        f = open(file, "rb")
        number = 0
        for _ in f:
            number += 1
        number = int(number / 2)
        whole_bp = reads_length * number
        filtered_reads_whole_bp_dict[file_name] = whole_bp
    return filtered_reads_whole_bp_dict


'''
重过滤
'''


def re_filter_reads(_dict, _kmer, stepSize, readseq, get_rc=False):
    '''
    :param _dict:     ref中一个基因的hashdict
    :param _kmer:     wordsize
    :param stepSize:  step
    :param readseq:   read
    :param get_rc:    是否需要反向互补，filter中已经反向互补了，并写进一个文件。ref已经反向互补了，此时不需要反向互补了
    :return:
    '''
    tmplist = []
    for j in range(0, len(readseq) - _kmer + 1 - stepSize, stepSize):
        kmer = readseq[j:j + _kmer]
        if kmer in _dict:
            tmplist.extend(_dict.get(kmer)[2:])  # 剔除kmercount 和 位置信息

    if get_rc:
        temp_seq = reverse_complement_limit(readseq)
        for j in range(0, len(temp_seq) - _kmer + 1 - stepSize, stepSize):
            kmer = temp_seq[j:j + _kmer]
            if kmer in _dict:
                tmplist.extend(_dict.get(kmer)[2:])

    return set(tmplist)


def do_re_filter_reads(parameter_information, used_dict, kmer, big_reads_path, re_filtered_reads_whole_bp_dict,
                       get_rc=False):
    '''
    :param parameter_information: 参数表
    :param used_dict: ref中一个基因的hashdict
    :param kmer: kmer
    :param big_reads_path: filtered_out超过500层或大于10M，big reads不发生变化
    :param re_filtered_reads_whole_bp_dict: 用于记录新生成的re-filter 过滤的reads bp总数
    :param get_rc: 是否需要反向互补
    :return:
    '''
    stepSize = parameter_information["step"]
    out_dir = parameter_information["out"]
    infile = open(big_reads_path, "rb")

    for _ in infile:
        line = _
        if not line.strip():  # 规避空行
            continue
        else:
            temp_rec = [_, infile.readline()]
        for file_name in re_filter_reads(used_dict, kmer, stepSize, temp_rec[1].decode('utf-8'), get_rc):
            re_filtered_reads_whole_bp_dict[file_name] += len(temp_rec[1].decode('utf-8'))
            with open(os.path.join(out_dir, file_name + ".fasta"), "a+", newline="\n") as f:  # newline="\n" 统一不同操作系统换行符
                f.writelines(['>', temp_rec[0].decode('utf-8')[1:], temp_rec[1].decode('utf-8')])
    infile.close()


def do_re_filter_loop(parameter_information, re_filtered_reads_whole_bp_dict, ref_length_dict, re_filter_dict,
                      thread_id, thread,keep_quiet):
    '''
    :param parameter_information:
    :param re_filtered_reads_whole_bp_dict:
    :param ref_length_dict: 记录参考序列的平均长度，用于计算richness
    :param re_filter_dict: 记录需要重过滤基因 大于500层或大于10M
    :param thread_id:  标记线程id
    :param thread:     线程总数
    :param keep_quiet:     是否保持沉默
    :return:
    '''

    current = 0
    for i in re_filter_dict:
        current += 1
        if current % thread == thread_id:
            big_reads_path = i["big_reads_path"]
            filtered_reads_path = i["filtered_reads_path"]
            filtered_reads_path_backup = i["filtered_reads_path_backup"]
            ref_path = i["ref_path"]
            gene_name = i["gene_name"]
            kmer = parameter_information["wordsize"]
            kmer_add = 2
            re_filter_number = 0

            filter_csv_path = parameter_information["filter_csv_path"]
            k1 = parameter_information["wordsize"]

            while True:
                old_kmer = kmer  # 记录上一次kmer
                kmer += kmer_add
                new_kmer = kmer  # 记录本次kmer值

                if kmer > 127:
                    sth = [gene_name, k1, old_kmer]
                    mylog(filter_csv_path, sth)
                    break

                # 备份数据 （每次开始优先备份数据，覆盖式的备份; 达到循环退出条件时候，删除备份）
                if re_filter_number == 0:
                    shutil.copy(big_reads_path, filtered_reads_path_backup)
                else:
                    shutil.copy(filtered_reads_path, filtered_reads_path_backup)
                re_filter_number += 1

                used_dict = gethashdict(ref_path, kmer, get_rc=True, pos=False,print_info=False,keep_quiet=keep_quiet)  # filter  re-filter 的策略都是 参考序列反向互补,reads不反向互补
                with open(filtered_reads_path, "w") as f:  # 清空序列
                    pass
                re_filtered_reads_whole_bp_dict[gene_name] = 0  # 重置
                do_re_filter_reads(parameter_information, used_dict, kmer, big_reads_path,
                                   re_filtered_reads_whole_bp_dict,
                                   False)
                richness = int(re_filtered_reads_whole_bp_dict[gene_name] / ref_length_dict[gene_name])
                size = get_file_physical_size(filtered_reads_path)

                # 相关参数
                sharp_cutoff_threshold = 512 * 1024  # 锐减阈值
                tolerate_threshold = 10 * 1024 * 1024  # 容忍阈值
                richness_level_1 = 512
                richness_level_2 = 1024
                richness_level_3 = 2048
                kmer_add_1 = 2
                kmer_add_2 = 4
                kmer_add_3 = 8

                # for test
                # sharp_cutoff_threshold = 200 * 1024  # 锐减阈值
                # tolerate_threshold = 250 * 1024  # 容忍阈值
                # richness_level_1 = 100
                # richness_level_2 = 1024
                # richness_level_3 = 2048
                # kmer_add_1 = 6
                # kmer_add_2 = 10
                # kmer_add_3 = 14

                # 打印输出
                if size >= sharp_cutoff_threshold:
                    # print(" " * 50, flush=True, end='\r')  # 之前的信息比较长，冲刷掉
                    # print("gene:{} kmer:{} richness:{} size:{}KB".format(gene_name, kmer, richness, int(size / 1024)),
                    #       flush=True, end='\r')

                    message1 = " " * 50
                    My_Log_Recorder(message=message1, keep_quiet=keep_quiet, path="", stage="", printout_flag=True,
                                    stored_flag=False, make_a_newline=False, use_cutting_line=False)
                    message2 ="gene:{} kmer:{} richness:{} size:{}KB".format(gene_name, kmer, richness, int(size / 1024))
                    My_Log_Recorder(message=message2, keep_quiet=keep_quiet, path="", stage="", printout_flag=True,
                                    stored_flag=False, make_a_newline=False, use_cutting_line=False)





                # kmer太大，结果锐减且突破最低临界值 ---保留上一次的结果
                if size < sharp_cutoff_threshold:
                    # 删除备份，保留上次结果
                    os.remove(filtered_reads_path)
                    shutil.copy(filtered_reads_path_backup, filtered_reads_path)
                    os.remove(filtered_reads_path_backup)
                    # log
                    sth = [gene_name, k1, kmer]
                    mylog(filter_csv_path, sth)
                    break
                elif size >= sharp_cutoff_threshold and size <= tolerate_threshold:
                    # 删除备份
                    os.remove(filtered_reads_path_backup)
                    # log
                    sth = [gene_name, k1, kmer]
                    mylog(filter_csv_path, sth)
                    break
                else:
                    pass

                # 确定Kmer增大幅度
                if richness <= richness_level_1 or size <= tolerate_threshold:
                    # 删除备份
                    os.remove(filtered_reads_path_backup)
                    # log
                    sth = [gene_name, k1, kmer]
                    mylog(filter_csv_path, sth)
                    break
                elif kmer > 127:
                    # 删除备份
                    os.remove(filtered_reads_path_backup)
                    # log
                    sth = [gene_name, k1, old_kmer]
                    mylog(filter_csv_path, sth)
                    break

                elif richness >= richness_level_1 and richness < richness_level_2:
                    kmer_add = kmer_add_1
                    os.remove(filtered_reads_path_backup)
                elif richness >= richness_level_2 and richness <= richness_level_3:
                    kmer_add = kmer_add_2
                    os.remove(filtered_reads_path_backup)
                else:
                    kmer_add = kmer_add_3
                    os.remove(filtered_reads_path_backup)
                del used_dict
                gc.collect()  # 释放内存


def my_filter_main(configuration_information):
    # 初始化操作
    data1 = configuration_information["data1"]
    data2 = configuration_information["data2"]
    single = configuration_information["single"]
    thread = configuration_information["thread_number"]
    wordsize = configuration_information["k1"]
    out_dir = configuration_information["out"]
    step = configuration_information["step_length"]
    reference = configuration_information["reference"]
    data_size = configuration_information["data_size"]
    quiet=configuration_information["quiet"]


    filter_csv_path = os.path.join(out_dir, "filter.csv")

    # 读写瓶颈 硬盘
    if thread >= 4:
        limit_thread = 4
    else:
        limit_thread = thread
    parameter_information = {"data1": data1, "data2": data2, "single": single,
                             "thread": thread, "wordsize": wordsize, "out": out_dir, "step": step,
                             "reference": reference, "data_size": data_size, "filter_csv_path": filter_csv_path
                             }

    # print('======================== Filter =========================')
    # cutting_line(" Filter ")
    My_Log_Recorder(message=" Filter ", keep_quiet=quiet, path="", stage="", printout_flag=False,
                    stored_flag=False, make_a_newline=False, use_cutting_line=True)

    hashdict = {}
    reads_length = 100  # 默认，后续会重新计算
    t1 = time.time()
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    hashdict = gethashdict(reference, wordsize, get_rc=True, pos=False, print_info=True,keep_quiet=quiet)  # 注意reference做了反向互补

    # print(" " * 50, flush=True, end='\r')  # 之前的信息比较长，冲刷掉
    # print("Hashdict has prepared")
    # print('Memory of hashdict:', round(6 * sys.getsizeof(hashdict) / 1024 / 1024 / 1024, 4), 'G')
    message1 = " " * 50
    My_Log_Recorder(message=message1, keep_quiet=quiet, path="", stage="", printout_flag=True,
                    stored_flag=False, make_a_newline=False, use_cutting_line=False)
    message2 = "Hashdict has prepared"
    My_Log_Recorder(message=message2, keep_quiet=quiet, path="", stage="", printout_flag=True,
                    stored_flag=False, make_a_newline=True, use_cutting_line=False)
    message3 = "Memory of hashdict:{}G".format(round(6 * sys.getsizeof(hashdict) / 1024 / 1024 / 1024, 4))
    My_Log_Recorder(message=message3, keep_quiet=quiet, path="", stage="", printout_flag=True,
                    stored_flag=False, make_a_newline=True, use_cutting_line=False)

    size = get_size(reference)  # 3.500
    ref_limit_size = 3.53 * 1024 * 1024  # 353*10*1000bp
    combine_flag = 1  #合并多线程输出的结果

    if sys.platform in ('darwin'):
        # multiprocessing.set_start_method('fork')
        multiprocessing.set_start_method('fork', force=True) #macos
    if data1 and data2:
        if sys.platform in ("win32") and size >= ref_limit_size or thread == 1:
            reads_length = get_reads_length(data1)
            filter_paired_reads_single_thread(hashdict, wordsize, step, data1, data2, out_dir, False, data_size,quiet)
            combine_flag = 0
        else:
            reads_length = get_reads_length(data1)
            thread_list = []
            for t in range(limit_thread):
                p = Process(target=filter_paired_reads,
                            args=(hashdict, wordsize, step, data1, data2, out_dir, t, limit_thread, False, data_size,quiet))
                thread_list.append(p)
            for t in thread_list:
                t.start()
            for t in thread_list:
                t.join()
    if single:
        if sys.platform in ("win32") and size >= ref_limit_size or thread == 1:
            reads_length = get_reads_length(single)
            filter_single_reads_single_thread(hashdict, wordsize, step, single, out_dir, get_rc=False,
                                              data_size=data_size,keep_quiet=quiet)
            combine_flag = 0
        else:
            reads_length = get_reads_length(single)
            thread_list = []
            for t in range(limit_thread):
                p = Process(target=filter_single_reads,
                            args=(hashdict, wordsize, step, single, out_dir, t, limit_thread, False, data_size,quiet))
                thread_list.append(p)
            for t in thread_list:
                t.start()
            for t in thread_list:
                t.join()

    ############################## combine file #############################
    try:
        if combine_flag:
            path_dict = get_filtered_path_for_combine(out_dir)
            combine_file_parallel(out_dir, path_dict, thread)
    except:
        pass

    t2 = time.time()
    # print(" " * 50, flush=True, end='\r')  # 之前的信息比较长，冲刷掉
    message=" " * 50
    My_Log_Recorder(message=message, keep_quiet=quiet, path="", stage="", printout_flag=True,
                    stored_flag=False, make_a_newline=False, use_cutting_line=False)


    used_temp_time = format((t2 - t1), ".2f")
    # print("Filter time used: {}s".format(used_temp_time))
    message="Filter time used: {}s".format(used_temp_time)
    My_Log_Recorder(message=message, keep_quiet=quiet, path="", stage="", printout_flag=True,
                    stored_flag=False, make_a_newline=True, use_cutting_line=False)


    ############################## re-filter #############################
    # cutting_line(" Re-Filter ")
    My_Log_Recorder(message=" Re-Filter ", keep_quiet=quiet, path="", stage="", printout_flag=False,
                    stored_flag=False, make_a_newline=False, use_cutting_line=True)

    # filtered_limit_size = 0.5 * 1024 * 1024 #for test
    filtered_limit_size = 10 * 1024 * 1024
    re_filter_files, no_re_filter_files = get_re_filter_path(out_dir, filtered_limit_size)
    header = ["gene", "k1", "re_k1"]
    mylog(filter_csv_path, header)
    if no_re_filter_files:
        for i in no_re_filter_files:
            gene = get_basename(i)
            sth = [gene, str(wordsize), "None"]
            mylog(filter_csv_path, sth)

    if re_filter_files == []:
        t3 = time.time()
        # print("Re-Filter time used: {}s".format(round((t3 - t2), 4)))
        used_temp_time = format((t3 - t2), ".2f")
        message="Re-Filter time used: {}s".format(used_temp_time)
        My_Log_Recorder(message=message, keep_quiet=quiet, path="", stage="", printout_flag=True,
                        stored_flag=False, make_a_newline=True, use_cutting_line=False)
        # print("Re-Filter time used: {}s".format(used_temp_time))
    else:
        re_filter_ref_files = []
        for i in re_filter_files:
            name = get_basename(i)
            re_filter_ref_path = os.path.join(reference, name + ".fasta")
            re_filter_ref_files.append(re_filter_ref_path)

        ref_length_dict = {}  # 参考序列的平均长度
        filtered_reads_whole_bp_dict = defaultdict(int)  ## 过滤reads总数
        re_filtered_reads_whole_bp_dict = defaultdict(int)  ## 重过滤reads总数
        ref_length_dict = get_ref_info(re_filter_ref_files, ref_length_dict)  # ref_length_dict  第一次过滤的参考信息
        filtered_reads_whole_bp_dict = get_filtered_reads_bp_number(re_filter_files, reads_length,
                                                                    filtered_reads_whole_bp_dict)

        # print(ref_length_dict)
        # print(filtered_reads_whole_bp_dict)

        if not os.path.isdir(os.path.join(out_dir, 'big_reads')):
            os.mkdir(os.path.join(out_dir, 'big_reads'))
        re_filter_dict = []
        for i in zip(ref_length_dict.items(), filtered_reads_whole_bp_dict.items()):
            temp = {}
            gene_name = i[0][0]
            filtered_reads_path = os.path.join(out_dir, gene_name + ".fasta")
            filtered_reads_path_backup = os.path.join(out_dir, gene_name + "_backup" + ".fasta")
            richness = int(i[1][1] / i[0][1])
            size = get_file_physical_size(filtered_reads_path)

            big_reads_path = os.path.join(out_dir, "big_reads", gene_name + ".fasta")
            if os.path.isfile(reference):
                ref_path = reference
            else:
                ref_path = os.path.join(reference, gene_name + ".fasta")
            shutil.move(filtered_reads_path, big_reads_path)
            temp["gene_name"] = gene_name
            temp["richness"] = richness
            temp["size"] = size
            temp["filtered_reads_path"] = filtered_reads_path
            temp["filtered_reads_path_backup"] = filtered_reads_path_backup
            temp["ref_path"] = ref_path
            temp["big_reads_path"] = big_reads_path
            re_filter_dict.append(temp)

        if len(re_filter_dict) > 0:
            thread_list = []
            for t in range(min(thread, len(re_filter_dict))):  # 减少进程开销
                p = Process(target=do_re_filter_loop,
                            args=(parameter_information, re_filtered_reads_whole_bp_dict, ref_length_dict, re_filter_dict,
                                t, min(thread, len(re_filter_dict)),quiet))
                thread_list.append(p)
            for t in thread_list:
                t.start()
            for t in thread_list:
                t.join()

        t3 = time.time()
        print(" " * 50, flush=True, end='\r')  # 之前的信息比较长，冲刷掉
        used_temp_time = format((t3 - t2), ".2f")
        message="Re-Filter time used: {}s".format(used_temp_time)
        My_Log_Recorder(message=message, keep_quiet=quiet, path="", stage="", printout_flag=True,
                        stored_flag=False, make_a_newline=True, use_cutting_line=False)
        # print("Re-Filter time used: {}s".format(used_temp_time))


if __name__ == '__main__':
    pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="filter",
                                   usage="%(prog)s <-1 -2| -s> <-r> <-o> [options]")
    pars.add_argument("-1", "--data1", dest="data1",
                      help="One end of the paired-end reads, support fastq/fastq.gz/fastq.bz2", metavar="")
    pars.add_argument("-2", "--data2", dest="data2",
                      help="Another end of the paired-end reads, support fastq/fastq.gz/fastq.bz2",
                      metavar="")
    pars.add_argument("-s", "--single", dest="single",
                      help="Single-read, support fastq/fastq.gz/fastq.bz2", metavar="")
    pars.add_argument("-o", "--out", dest="out", help="Specify the result folder",
                      metavar="", required=True)
    pars.add_argument("-k1", "--kmer1", dest="kmer1", help="length of a word-size [default=29]",
                      default=29, type=int, metavar="")
    pars.add_argument("-step_length", metavar="", dest="step_length", type=int,
                      help="the length of the sliding window on the reads, [default=4]", default=4)  # step length
    pars.add_argument("-t", "--thread", metavar="", dest="thread_number", help="Thread", type=int, default=4)
    pars.add_argument("-r", "--reference", metavar="", dest="reference", type=str, help="references", required=True)
    pars.add_argument("-d", "--data", metavar="", dest="data_size", help="data size", default='all')
    pars.add_argument("-quiet", dest="quiet", help="Do not write progress messages to stderr", default=False,action='store_true')

    args = pars.parse_args()

    # 初始化操作
    data1 = get_absolute(args.data1)
    data2 = get_absolute(args.data2)
    single = get_absolute(args.single)
    thread_number = args.thread_number
    k1 = args.kmer1
    out_dir = get_absolute(args.out)
    step_length = args.step_length
    reference = get_absolute(args.reference)
    data_size = args.data_size
    quiet=args.quiet


    filter_configuration_information = {"data1": data1, "data2": data2, "single": single,
                                        "thread_number": thread_number, "k1": k1,
                                        "out": out_dir,
                                        "step_length": step_length,
                                        "reference": reference, "data_size": data_size,
                                        "quiet":quiet
                                        }

    my_filter_main(filter_configuration_information)
