#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2022/3/31 9:32
# @Author  : xiepulin
# @File    : my_assemble.py
# @Software: PyCharm
import sys
import argparse
from Bio import SeqIO
import gzip
import os
import time
from collections import defaultdict
import platform
import gc
from concurrent.futures import ProcessPoolExecutor
import csv
import copy
import math

'''
获取绝对路径
'''


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


'''
获得目标文件夹下的fasta文件，仅第一层
'''


def get_fasta_file(path):
    path = get_absolute(path)
    files = os.listdir(path)
    suffix_list = ["fa", "fas", "fasta"]
    files_list = []
    for i in files:
        suffix = i.split(".")[-1].lower()
        file_path = os.path.join(path, i)
        if os.path.isfile(file_path) and suffix in suffix_list:
            files_list.append(file_path)
    return files_list


'''
获得参考序列高质量的平均长度
1.获得平均长度 2.剔除低于平均长度的序列 3.再取一次平均长度
'''


def get_ref_info(reference, ref_length_dict):
    files = get_file_list(reference)
    for file in files:
        length_list = []
        ref_count = 0
        file_name = get_basename(file)
        for rec in SeqIO.parse(file, "fasta"):
            length = len(rec.seq)
            length_list.append(length)
            ref_count += 1

        average_length = sum(length_list) / len(length_list)
        # 剔除过短的序列
        superior_length_list = [i for i in length_list if i >= average_length]
        if superior_length_list:
            superior_average_length = sum(superior_length_list) / len(superior_length_list)
            average_length = max(average_length, superior_average_length)
        ref_length_dict[file_name] = int(average_length)
    return ref_length_dict


# 获得序列长度 剔除N
def get_length(path):
    length_all = 0
    number = 0
    for rec in SeqIO.parse(path, "fasta"):
        seq = str(rec.seq).replace("N", "")
        length_all += len(seq)
        number += 1
    length = int(length_all / number)
    return length


# 获取reads长度
def get_reads_length(file, type, _number=100):
    reads_length = []
    number = 0
    number_limit = _number

    if type == "fastq":
        infile = gzip.open(file, 'rt') if file[-3:].lower() == ".gz" else open(file, 'r')
        for _ in infile:
            if number < number_limit:
                line = _
                if not line.strip():  # 规避空行
                    continue
                else:
                    temp_rec = [_.strip(), infile.readline().strip(), infile.readline().strip(),
                                infile.readline().strip()]
                    number += 1
                    reads_length.append(len(temp_rec[1]))
            else:
                break
        infile.close()
        length = int(sum(reads_length) / len(reads_length))

    elif type == "fasta":
        infile = gzip.open(file, 'rt') if file[-3:].lower() == ".gz" else open(file, 'r')
        for _ in infile:
            if number < number_limit:
                line = _
                if not line.strip():  # 规避空行
                    continue
                else:
                    temp_rec = [_.strip(), infile.readline().strip()]
                    number += 1
                    reads_length.append(len(temp_rec[1]))
            else:
                break
        infile.close()
        length = int(sum(reads_length) / len(reads_length))

    else:
        length = 150
    return length


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



# 反向互补
def reverse_complement_all(seq):
    return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]


# 简化反向互补
def reverse_complement_limit(seq):
    return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


# kmers生成器
def make_kmers(seq, k):
    for i in range(len(seq) - k + 1): yield seq[i:i + k]


# 补齐生成器
def forward(seq):
    for x in 'ACGT': yield seq[1:] + x


def reverse(seq):
    for x in 'TGCA': yield x + seq[:-1]


# 创建哈西字典, 将reference存入字典

'''
data_structure {"kmer1":[kmer1_count,pos1,gene_name1,gene_name2,gene_name3],
                "kmer2":[kmer2_count,pos2,gene_name1,gene_name2,gene_name3]
                }

{'CGTTTTGACTGTATCGC': [3, 0.0017985611510791366, 'matK'],
'GCGATACAGTCAAAACG': [3, -2.9712230215827335, 'matK']}
'''


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



'''
assemble的hashdict  记录kmercount频次，记录kmer在ref的平均位置信息
'''

'''
data_structure   [kmer: [kmercount,[location]]

'CCCTCAAGAAAAGGGCAACCAATGCCCACAA': [1, [0.927315357561547]] 
'TAACGAGGCAAAGGCAAGTGCAGCAGCACTT': [1, [0]]
'''


def make_assemble_hashdict(filtered_out, merSize, ref_dict, get_rc=False):
    files = get_file_list(filtered_out)
    if files == []:
        return 0
    kmer_dict = defaultdict(list)
    kmer_count = 0
    for file in files:
        infile = open(file, 'r', encoding='utf-8', errors='ignore')
        infile.readline()
        for line in infile:
            temp_str = []
            while line and line[0] != '>':
                temp_str.append(line)
                line = infile.readline()
            refseq_r = ''.join(filter(str.isalnum, ''.join(temp_str).upper()))  # 过滤器，返回迭代器对象
            if get_rc:
                for refseq in [refseq_r, reverse_complement_limit(refseq_r)]:
                    kmer_count += len(refseq) - merSize + 1
                    for j in range(0, len(refseq) - merSize + 1):
                        temp_list, temp_len, kmer = [], [0], refseq[j:j + merSize]
                        if kmer in kmer_dict:
                            temp_list = kmer_dict[kmer]
                            temp_list[0] += 1
                        else:
                            temp_list = [1, temp_len]
                            kmer_dict[kmer] = temp_list
                        if kmer in ref_dict:
                            temp_len[0] = 1 + ref_dict[kmer][1] / ref_dict[kmer][0] if ref_dict[kmer][1] < 0 else \
                                ref_dict[kmer][1] / ref_dict[kmer][0]

                        # ref_dict[kmer][1] / ref_dict[kmer][0]  累计位置信息/kmercount = 平均位置信息  之所以用1+，是因为考虑了反向互补的情况

            else:
                for refseq in [refseq_r]:
                    kmer_count += len(refseq) - merSize + 1
                    for j in range(0, len(refseq) - merSize + 1):
                        temp_list, temp_len, kmer = [], [0], refseq[j:j + merSize]
                        if kmer in kmer_dict:
                            temp_list = kmer_dict[kmer]
                            temp_list[0] += 1
                        else:
                            temp_list = [1, temp_len]
                            kmer_dict[kmer] = temp_list
                        if kmer in ref_dict:
                            temp_len[0] = 1 + ref_dict[kmer][1] / ref_dict[kmer][0] if ref_dict[kmer][1] < 0 else \
                                ref_dict[kmer][1] / ref_dict[kmer][0]

        infile.close()
    return kmer_dict, kmer_count


def get_dis(_pos, new_pos, weight=0.5):
    return 1 - abs(_pos - new_pos) if (_pos and new_pos) else weight


'''
权重赋分
当出现多歧，赋分差距越大，则越倾向于保守而不是延长序列。 a:10分   b：1分  b路径可能要延长很久才能超过a
'''


def get_weight(_count, _pos, average_pos, _dis=0.5):
    '''
    :param _count:
    :param _pos: 当前拼接到的位置
    :param average_pos:  kmer在ref上出现的平均位置,经make_assemble_hashdict处理后， 位置信息都在0,1之间 （考虑了反向互补）
    :param _weight:
    :return:
    '''
    if average_pos and _pos:
        dis = (1 - abs(_pos - average_pos))
    else:
        dis = _dis
    weight = _count ** dis
    return weight


def get_contig_forward(_dict, seed, iteration=1000):
    '''
    :param _dict:  filtered_out哈希字典    {["kmer1":[kmer1_count,[pos1]],["kmer2":[kmer2_count,[pos2]]  }       {'TATACCTGTCAACGGAT': [17, [0]]}
    :param seed:   seed
    :param iteration:  迭代次数 之前指最大递归次数
    :return:
    '''

    '''
    temp_list:  存储用于组装的kmer,每一次弹出最后的一个kmer,用于延申
    reads_list: 存储用于组装的kmer
    stack_list: 栈 将Kmer记录为节点，遍历树的方法，每一次到叶子节点就回溯到上一个分歧节点，直到栈为空（出现分歧的时候记录）
    best_kmc:  记录最大累计Kmercount
    cur_kmc: 记录当前累计kmercount
    best_pos 记录最佳位置 ，用于contig组装scaffold
    cur_pos   记录当前位置
    best_seq:记录最佳路径
    cur_seq:记录当前路径

    核心，kmer的组装是按照权重组装的，权重由reads的kmercout 和reads在ref上的位置决定

    _dict[i][0] ** get_dis(_pos, _dict[i][1][0], weight)  没有位置信息，将开根号处理（默认0.5）;    核心思想距离越近，kmercount越大的权重越大

    '''

    temp_list, reads_list, stack_list = [seed], [seed], []
    best_kmc, cur_kmc, best_pos, cur_pos, best_seq, cur_seq = [], [], [], [], [], []
    _pos, node_distance = 0, 0
    best_used_kmercount, cur_used_kmercount = [], []  # 增添用于assemble的kmer计数

    while True:
        # node = sorted([(i, _dict[i][1][0], _dict[i][0] ** get_dis(_pos, _dict[i][1][0], weight), _dict[i][0]) for i in forward(temp_list[-1]) if i in _dict], key=lambda _: _[2], reverse=True)
        node = sorted([(i, _dict[i][1][0], get_weight(_dict[i][0], _pos, _dict[i][1][0]), _dict[i][0]) for i in
                       forward(temp_list[-1]) if i in _dict], key=lambda _: _[2], reverse=True)
        # node  [ "kmer",cur_pos, weight，kmercount]  [('GTTTTGACTGTATCGCA', 0.001199040767386091, 24.995845787996004,146)]
        while node:
            if node[0][0] in temp_list:
                node.pop(0)
            else:
                break
        if not node:
            iteration -= 1
            if sum(cur_kmc) > sum(best_kmc):
                best_kmc, best_seq, best_pos = cur_kmc.copy(), cur_seq.copy(), cur_pos.copy()
                best_used_kmercount = cur_used_kmercount.copy()
            # [(temp_list.pop(), cur_kmc.pop(), cur_seq.pop(), cur_pos.pop()) for _ in range(node_distance)]
            [(temp_list.pop(), cur_kmc.pop(), cur_seq.pop(), cur_pos.pop(), cur_used_kmercount.pop()) for _ in
             range(node_distance)]  # 新增Kmercount计数
            if stack_list:
                node, node_distance, _pos = stack_list.pop()
            else:
                break
        if len(node) >= 2:
            stack_list.append((node[1:], node_distance, _pos))
            node_distance = 0
        if node[0][1] > 0:
            _pos = node[0][1]

        temp_list.append(node[0][0])
        reads_list.append(node[0][0])
        cur_pos.append(node[0][1])  # 要么有位置信息，要么未出现在ref上则为0  有位置信息的会被优先组装，随后在原有基础上向两边延申
        cur_kmc.append(node[0][2])
        cur_seq.append(node[0][0][-1])
        cur_used_kmercount.append(node[0][3])  # 新增kmercount计数

        node_distance += 1

        if not iteration: break

    return best_seq, reads_list, best_kmc, best_pos, best_used_kmercount


def get_contig(_dict, seed, iteration=1000, weight=0.5):
    '''
    :param _dict:  filtered_dict
    :param seed:  seed
    :param iteration:  迭代次数
    :param weight:   权重（kmer在ref上的位置+kmer频次）
    :return: contig contig在ref上的最小位置，contig在ref上的最大位置，kmer_list
    '''

    right, reads_list_1, k_1, pos_1, best_used_kmercount_1 = get_contig_forward(_dict, seed,
                                                                                iteration)  # best_seq, reads_list, best_kmc, best_pos,best_used_kmercount
    left, reads_list_2, k_2, pos_2, best_used_kmercount_2 = get_contig_forward(_dict, reverse_complement_limit(seed),
                                                                               iteration)

    # print("pos_1:{}".format(pos_1))   #[0.1,0.2,0.3...0.9, 0 , 0, 0, 0 .0 ,0 ,0 ]
    # print("pos_2:{}".format(pos_2))   #[0 , 0, 0, 0 .0 ,0 ,0 ]

    _pos = [x for x in pos_1 + pos_2 if x > 0]
    min_pos, max_pos = 0, 1
    if _pos: min_pos, max_pos = min(_pos), max(_pos)
    return reverse_complement_limit(''.join(left)) + seed + ''.join(
        right), min_pos, max_pos, reads_list_1 + reads_list_2, sum(best_used_kmercount_1) + sum(best_used_kmercount_2)


def get_scaffold(_dict, seed_list, limit, ref_length, iteration=1000, weight=0.5):
    mid_seed, min_dis = seed_list[0], 1
    # for seed in seed_list:
    #     if abs(seed[2] - 0.5) < min_dis:
    #         mid_seed, min_dis = seed, abs(seed[2] - 0.5)
    contig, min_pos, max_pos, read_list, best_used_kmercount = get_contig(_dict, mid_seed[0], iteration, weight)

    for i in read_list:
        if i in _dict[i]: del _dict[i]
    contigs = [(contig, min_pos, max_pos)]
    # print(seed_list) #[('TTTGGTGGAGTTGTTGGTG', 92, 79.03477756740409)] ["kmer",kmercount.pos_sum累计位置]

    # 左侧拼接 left
    new_seed = []
    for x in seed_list:
        if x[2] / x[1] < min_pos and x[1] > limit:
            new_seed = x
    while new_seed:
        # 将使用过的reads剔除
        contig, _min_pos, _max_pos, read_list, best_used_kmercount = get_contig(_dict, new_seed[0], iteration, weight)
        for i in read_list:
            if i in _dict[i]:
                del _dict[i]
        # print(new_seed[2] / new_seed[1],"1")
        # print(min_pos,max_pos,"pos")
        if new_seed[2] / new_seed[1] < _min_pos or new_seed[2] / new_seed[1] > _max_pos:
            break
        # print(_min_pos, _max_pos,contig,"left")
        if _max_pos < contigs[0][2]:  # 新contig完全在左侧  [(contig, min_pos, max_pos)]
            contigs.insert(0, (contig, _min_pos, _max_pos))
        else:
            break
        new_seed = []
        for x in seed_list:
            if x[2] / x[1] < _min_pos and x[1] > limit:
                new_seed = x

    # 右侧拼接
    new_seed = []
    for x in seed_list:
        if x[2] / x[1] > max_pos and x[1] > limit:
            new_seed = x
    while new_seed:
        contig, _min_pos, _max_pos, read_list, best_used_kmercount = get_contig(_dict, new_seed[0], iteration, weight)
        for i in read_list:
            if i in _dict[i]:
                del _dict[i]
        # print(new_seed[2] / new_seed[1],"1")
        # print(min_pos,max_pos,"pos")
        if new_seed[2] / new_seed[1] < _min_pos or new_seed[2] / new_seed[1] > _max_pos:
            break
        # print(_min_pos, _max_pos, contig,"right")
        if _min_pos > contigs[-1][1]:  # 完全在右侧 [(contig, min_pos, max_pos)]
            contigs.append((contig, _min_pos, _max_pos))
        else:
            break
        new_seed = []
        for x in seed_list:
            if x[2] / x[1] > _max_pos and x[1] > limit:
                new_seed = x

    scaffold = contigs[0][0]
    for x in range(1, len(contigs)):
        scaffold += max(2, int((contigs[x][1] - contigs[x - 1][2]) * ref_length)) * 'N' + contigs[x][
            0]  # contig 从左往右组装， 下一个最小contig最小位置减去上一个contig最大位置 能估算出gap空洞长度
    return scaffold


'''
利用交集确定种子，利用动态limit过滤低质量reads
'''


def get_seed(ref_dict, filtered_dict, ref_length, limit_count):
    # 利用交集确定种子
    for i in ref_dict:
        if i not in filtered_dict:
            ref_dict[i][0] = 0
    # 动态limit
    if limit_count < 0:
        limit_count = dynamic_limit_v2(filtered_dict, ref_length, 2)
    # 剔除低质量reads: 频次太低且未在ref中出现，注意频次很高未在ref中的reads是保留了的
    if limit_count > 0:
        _filter = [x for x in filtered_dict if filtered_dict[x][0] <= limit_count]
        for x in _filter:
            if filtered_dict[x][1][0] == 0:
                del filtered_dict[x]
        _filter = []
    seed_list = [(x, ref_dict[x][0]) for x in ref_dict if ref_dict[x][0] > 0]
    list.sort(seed_list, key=lambda x: x[1], reverse=True)
    return seed_list


def dynamic_limit_v1(_dict, smoother_level=4, smoother_size=64):
    # 计算kmercount 频次分布图,默认记录kmer最大频次为256
    count_list = [0] * smoother_size
    for x in _dict:
        if _dict[x][0] <= smoother_size:
            count_list[_dict[x][0] - 1] += 1

    # print(count_list, "dynamic_limit_v1")
    # 平滑器
    for i in range(smoother_level):
        for x in range(1, smoother_size - smoother_level + i):
            if count_list[x + smoother_level - 1 - i] < count_list[x - 1] and count_list[x] < count_list[
                x + smoother_level - i]:
                return min(x + 1, 64)
    return 2


# F0=F1+F2+length*rate rate取值为2的时候，比较宽松.
def dynamic_limit_v2(_dict, ref_length, assemble_kmer, reads_length, ref_rate=15, list_size=64):
    '''
    :param _dict:   filtered dict
    :param ref_length:   参考序列
    :param N:       filtered_kmer_count
    :param ref_rate:   参考序列倍率
    :param list_size:  记录kmer频次分布
    :return: (limit,deepth)
    理论支撑：
    F0-f1-f2= rate * ref_length*2  2代表2倍，因为计算了反向互补  rate 误差，取决于经验值（逻辑:基因组种类数约等于基因组长度）
    该方法对深度计算结果欠佳
    '''
    # 计算kmercount 频次分布图,默认记录kmer最大频次为256
    count_list = [0] * list_size
    for x in _dict:
        if _dict[x][0] <= list_size:
            count_list[
                _dict[x][0] - 1] += 1  # [46516, 5624, 1334, 376, 630, 108, 72, 44, 98, 94, 2,2,2,2,2,2,0,0,0,0,2]
    # 计算器
    KmerSum, F0, sum_f, sum_k, sum_k2 = 0, len(_dict), 0, 0, 0  # 总kmer数,kmer种类数,当前频次kmer数,累计kmer数
    for x in range(list_size):
        sum_f += count_list[x]
        sum_k += count_list[x]
        sum_k2 += count_list[x] * (x + 1)
        if (F0 - sum_k) / 2 < ref_length * ref_rate:  # ref_rate越大，limit越小，越是保守
            return max(2, x)  # 返回Limit,deepth

    return 2  # 返回Limit


# 自动识别浅层或者深层
def dynamic_limit_v3(_dict, ref_length, assemble_kmer, reads_length, ref_rate=15, list_size=65536):
    count_list = [0] * list_size
    for x in _dict:
        if _dict[x][0] <= list_size:
            count_list[
                _dict[x][0] - 1] += 1  # [46516, 5624, 1334, 376, 630, 108, 72, 44, 98, 94, 2,2,2,2,2,2,0,0,0,0,2]

    # print(count_list, "dynamic_limit_v3_count_list")

    level_1, level_2, level_max = 32, 64, list_size
    KmerSum_temp, KmerSum_level_1, KmerSum_level_2, KmerSum = 0, 0, 0, 0
    F0, sum_f, sum_k, sum_k2 = len(_dict), 0, 0, 0
    skimming_threshold_1 = 0.4  # 浅层数据kmercount会集中在64之内
    skimming_threshold_2 = 0.6
    tran_threshold_1 = 0.6  # 深层数据kmercount值会比较大，常常在200以上大范围集中
    tran_threshold_2 = 0.5
    kmer_threshold = 0.5 * reads_length
    for x in range(list_size):
        KmerSum += count_list[x] * (x + 1)
        KmerSum_temp += count_list[x] * (x + 1)
        if x == level_1 - 1:
            KmerSum_level_1 = KmerSum_temp
        if x == level_2 - 1:
            KmerSum_level_2 = KmerSum_temp

    # print(KmerSum_level_1/KmerSum,"k_l1")
    # print(KmerSum_level_2/KmerSum, "k_l2")
    # print(KmerSum,"k_all")
    # 当kmer过大，没有对应的kmercount
    if KmerSum == 0:
        limit = 2
        deepth = 0
        return limit, deepth

    '''
    1.kmer较小
    浅层kmercount64占比很大; 深层 kmercount64占比很小
    2.kmer较大
    浅层kmercount64占比更大; 深层 kmercount64占比增大，但仍不不是主要的
    '''
    if assemble_kmer <= kmer_threshold:
        # 浅层
        if KmerSum_level_1 / KmerSum >= skimming_threshold_1 or KmerSum_level_2 / KmerSum >= skimming_threshold_2:
            flag = 0
        # 深层
        elif (KmerSum - KmerSum_level_1) / KmerSum >= tran_threshold_1:
            flag = 1
        else:
            flag = 0
    else:
        # 浅层
        if KmerSum_level_1 / KmerSum >= skimming_threshold_1:
            flag = 0
        # 深层
        elif (KmerSum - KmerSum_level_1) / KmerSum >= tran_threshold_2:
            flag = 1
        else:
            flag = 0

    # 深层 limit用平滑曲线评估limit，由于数据质量高，用KmerSum计算深度
    if flag == 1:
        limit = dynamic_limit_v1(_dict, 4, 64)
        for x in range(limit + 1):
            sum_f += count_list[x]
            sum_k += count_list[x]
            sum_k2 += count_list[x] * (x + 1)
        deepth = round((KmerSum - sum_k2) * reads_length / (reads_length - assemble_kmer + 1) / ref_length / 2,
                       2)  # 适合深层

        # print("used_kmer:{}".format(KmerSum-sum_k2))


    # 浅层 用kmer种类的方法评估Limit.由于数据质量低，用F0计算深度
    else:
        limit_test = dynamic_limit_v1(_dict, 4, 64)  # 此处Limit用于计算深度，不作为真实的limit
        for x in range(limit_test + 1):
            sum_f += count_list[x]
            sum_k += count_list[x]
            sum_k2 += count_list[x] * (x + 1)
        deepth = round((KmerSum - sum_k2) * reads_length / (reads_length - assemble_kmer + 1) / ref_length / 2,
                       2)  # 适合深层
        # print("used_kmer:{}".format(KmerSum - sum_k2))
        limit = dynamic_limit_v2(_dict, ref_length, assemble_kmer, reads_length, ref_rate=15, list_size=64)
        # print("limit_dy2:{}".format(limit))

    return limit, deepth


def static_limit_v1(_dict, ref_length, assemble_kmer, limit, reads_length=150, ref_rate=15, list_size=65536):
    limit_static = copy.deepcopy(limit)
    count_list = [0] * list_size
    for x in _dict:
        if _dict[x][0] <= list_size:
            count_list[
                _dict[x][0] - 1] += 1  # [46516, 5624, 1334, 376, 630, 108, 72, 44, 98, 94, 2,2,2,2,2,2,0,0,0,0,2]

    # print(count_list, "dynamic_limit_v3_count_list")

    level_1, level_2, level_max = 32, 64, list_size
    KmerSum_temp, KmerSum_level_1, KmerSum_level_2, KmerSum = 0, 0, 0, 0
    F0, sum_f, sum_k, sum_k2 = len(_dict), 0, 0, 0
    skimming_threshold_1 = 0.4  # 浅层数据kmercount会集中在64之内
    skimming_threshold_2 = 0.6
    tran_threshold_1 = 0.6  # 深层数据kmercount值会比较大，常常在200以上大范围集中
    tran_threshold_2 = 0.5
    kmer_threshold = 0.5 * reads_length
    for x in range(list_size):
        KmerSum += count_list[x] * (x + 1)
        KmerSum_temp += count_list[x] * (x + 1)
        if x == level_1 - 1:
            KmerSum_level_1 = KmerSum_temp
        if x == level_2 - 1:
            KmerSum_level_2 = KmerSum_temp

    # print(KmerSum_level_1/KmerSum,"k_l1")
    # print(KmerSum_level_2/KmerSum, "k_l2")
    # print(KmerSum,"k_all")
    # 当kmer过大，没有对应的kmercount
    if KmerSum == 0:
        limit = 2
        deepth = 0
        return limit, deepth

    '''
    1.kmer较小
    浅层kmercount64占比很大; 深层 kmercount64占比很小
    2.kmer较大
    浅层kmercount64占比更大; 深层 kmercount64占比增大，但仍不不是主要的
    '''
    if assemble_kmer <= kmer_threshold:
        # 浅层
        if KmerSum_level_1 / KmerSum >= skimming_threshold_1 or KmerSum_level_2 / KmerSum >= skimming_threshold_2:
            flag = 0
        # 深层
        elif (KmerSum - KmerSum_level_1) / KmerSum >= tran_threshold_1:
            flag = 1
        else:
            flag = 0
    else:
        # 浅层
        if KmerSum_level_1 / KmerSum >= skimming_threshold_1:
            flag = 0
        # 深层
        elif (KmerSum - KmerSum_level_1) / KmerSum >= tran_threshold_2:
            flag = 1
        else:
            flag = 0

    # 深层 limit用平滑曲线评估limit，由于数据质量高，用KmerSum计算深度
    if flag == 1:
        limit = dynamic_limit_v1(_dict, 4, 64)
        for x in range(limit + 1):
            sum_f += count_list[x]
            sum_k += count_list[x]
            sum_k2 += count_list[x] * (x + 1)
        deepth = round((KmerSum - sum_k2) * reads_length / (reads_length - assemble_kmer + 1) / ref_length / 2,
                       2)  # 适合深层

    # 浅层 用kmer种类的方法评估Limit.由于数据质量低，用F0计算深度
    else:
        limit_test = dynamic_limit_v1(_dict, 4, 64)  # 此处Limit用于计算深度，不作为真实的limit
        for x in range(limit_test + 1):
            sum_f += count_list[x]
            sum_k += count_list[x]
            sum_k2 += count_list[x] * (x + 1)
        deepth = round((KmerSum - sum_k2) * reads_length / (reads_length - assemble_kmer + 1) / ref_length / 2,
                       2)  # 适合深层
        # limit = dynamic_limit_v2(_dict, ref_length, assemble_kmer, reads_length, ref_rate=15, list_size=64)
        # print("limit_dy2:{}".format(limit))
    return limit_static, deepth


def dynamic_weight(_dict):  # 匹配上的reads越多， 值越小。
    temp_count = 0
    for i in _dict:
        if _dict[i][1][0] > 0:  # _dict[i][1][0] 位置信息 大于0说明该条reads匹配到了ref
            temp_count += 1
    return 1 - temp_count / len(
        _dict)  # temp_count/len(_dict)代表相似度  1-相似度=差异度  weight=kmercount**weight 所以，差异度越大，ref不可信，则需要越依赖kmercount的大小


def do_assemble(ref_path, filtered_path, k2, limit_count, limit_min_length, limit_max_length,
                ref_length_dict, assembled_out_path, reads_length, unique_identity, change_seed=32,
                scaffold_or_not=True,keep_quiet=False):
    '''
    :param ref_path:
    :param filtered_path:
    :param k2:
    :param limit_count:
    :param limit_min_length:
    :param limit_max_length:
    :param ref_length_dict:
    :param assembled_out_path:
    :param reads_length:
    :param change_seed:
    :param scaffold_or_not:
    :param keep_quiet:


    ref_dict: [key1:[kmercount1,location1,gene1,gene2], key2:[kmercount2,location2,gene1,gene2],]
    filtered_dict : [key1:[kmercount1,location1,gene1,gene2], key2:[kmercount2,location2,gene1,gene2],]   eg: 'CCCTCAAGAAAAGGGCAACCAATGCCCACAA': [1, 0.9273153575615475, 'ETS']
    种子选择 (1)在ref 和reads 里面同时出现 （2）考虑ref频次与reads频次  score= x1 * log(x2+1) x1代表ref kmercount, x2代表 reads kmercount
    :return:
    '''

    gene_name = get_basename(ref_path)
    ref_dict = gethashdict(ref_path, k2, get_rc=True, pos=True, print_info=False,keep_quiet=keep_quiet)
    filtered_dict, filtered_kmer_count = make_assemble_hashdict(filtered_path, k2, ref_dict, get_rc=True)

    filtered_dict_back = copy.deepcopy(filtered_dict)  # 备份
    # print(filtered_dict)
    # print(ref_dict)

    ref_length = ref_length_dict[gene_name]
    contig_path = os.path.join(assembled_out_path, "contig", gene_name + ".fasta")
    short_contig_path = os.path.join(assembled_out_path, "short_contig", gene_name + ".fasta")
    scaffold_path = os.path.join(assembled_out_path, "scaffold", gene_name + ".fasta")
    ref_length_min = ref_length * limit_min_length
    ref_length_max = ref_length * limit_max_length

    write_contig = False
    best_seed = ["None"]  # 防止报错，只要能做出结果就会被覆盖掉

    # 利用交集确定种子
    for i in ref_dict:
        if i not in filtered_dict:
            ref_dict[i][0] = 0

    if limit_count == "auto":
        limit_count, deepth = dynamic_limit_v3(filtered_dict, ref_length, k2, reads_length, ref_rate=15,
                                               list_size=1024)  # ref_rate的选择非常关键
        # print("limit_count:{} deepth:{}".format(limit_count,deepth))
        while (limit_count):
            filtered_dict = copy.deepcopy(filtered_dict_back)  # 重新赋值
            _filter = [x for x in filtered_dict if filtered_dict[x][0] <= limit_count]
            for x in _filter:
                if filtered_dict[x][1][0] == 0:
                    del filtered_dict[x]
            _filter = []

            # seed_list = [(x, ref_dict[x][0], ref_dict[x][1]) for x in ref_dict if ref_dict[x][0] > 0 and ref_dict[x][1] > 0]  # 新增位置信息,仅考虑ref kmercount
            seed_list = [(x, ref_dict[x][0] * math.log10((filtered_dict[x][0] + 1)), ref_dict[x][1]) for x in ref_dict
                         if ref_dict[x][0] > 0 and ref_dict[x][1] > 0]  # 新增位置信息,同时考虑reads,ref kmercount

            if seed_list == []:
                return 0  # 找不到seed，一般是k2值太大
            list.sort(seed_list, key=lambda x: x[1], reverse=True)

            # print("limit_count:{0},deepth:{1}".format(limit_count,deepth))
            # 利用相似度计算权重 filtered reads 与ref 的大致距离
            cur_weight = dynamic_weight(filtered_dict) if filtered_dict else 0.5  # 0.52
            '''
            best_seed,best_contig,best_length 设定第一次默认最佳
            '''
            best_seed = seed_list[0]  # ('TTTGGTGGAGTTGTTGG', 95, 81.74894261950934)

            # 速度与允许迭代次数高度正相关，python最大递归1000次，我们改写为循环后，可突破上限

            best_contig, min_pos, max_pos, read_list, best_used_kmercount = get_contig(filtered_dict, best_seed[0],
                                                                                       iteration=1000,
                                                                                       weight=cur_weight)  # for test

            best_length = len(best_contig)

            # 对于过短的处理,更换seed

            if best_length < ref_length_min:

                # 更换seed
                change_time = 0
                last_seed = seed_list[0]
                for seed in seed_list[1:]:
                    if last_seed[0][1:] == seed[0][:-1] or last_seed[0][2:] == seed[0][:-2]:  # 移动1~2bp
                        last_seed = seed
                        continue
                    change_time += 1
                    last_seed = seed
                    new_contig, min_pos, max_pos, read_list, new_best_used_kmercount = get_contig(filtered_dict,
                                                                                                  last_seed[0],
                                                                                                  iteration=1000,
                                                                                                  weight=cur_weight)

                    # new_contig = get_contig(filtered_dict, best_seed[0], iteration=1000, weight=cur_weight)[0]
                    new_length = len(new_contig)
                    # print(" " * 100, flush=True, end='\r')  # 之前的信息比较长，冲刷掉
                    # print("{0}: use another seed: {1} ... {2:<2}".format(gene_name, last_seed[0][0:8], change_time),
                    #       flush=True,
                    #       end='\r')

                    message1 = " " * 100
                    My_Log_Recorder(message=message1, keep_quiet=keep_quiet, path="", stage="", printout_flag=True,
                                    stored_flag=False, make_a_newline=False, use_cutting_line=False)
                    message2 = "{0}: use another seed: {1} ... {2:<2}".format(gene_name, last_seed[0][0:8], change_time)
                    My_Log_Recorder(message=message2, keep_quiet=keep_quiet, path="", stage="", printout_flag=True,
                                    stored_flag=False, make_a_newline=False, use_cutting_line=False)

                    if new_length >= ref_length_min:
                        best_seed = last_seed
                        best_contig = new_contig
                        best_length = new_length
                        write_contig = True
                        break
                    elif change_time >= change_seed:  # 超过最大更换种子数目
                        # 种子为未达标情况下，最佳seed
                        break
                    else:  # 未达标情况下的最好结果
                        if new_length > best_length:
                            best_contig = new_contig
                            best_length = new_length
                            best_seed = last_seed

            # 对于过长的处理，不更换seed，改为增加limit

            elif best_length > ref_length_max:
                if limit_count >= 32:
                    limit_count_max = limit_count + 32
                else:
                    limit_count_max = 33

                # 加大limit
                for temp_limit in range(limit_count + 1, limit_count_max):
                    # 备份上一次的结果
                    old_seed = best_seed
                    old_contig = best_contig
                    old_length = len(old_contig)
                    # 本次的结果，后续会覆盖刷新
                    new_seed = best_seed
                    new_contig = best_contig
                    new_length = len(new_contig)

                    # print(new_length)

                    # 对filtered_dict再次精简
                    _filter = [x for x in filtered_dict if filtered_dict[x][0] <= temp_limit]
                    for x in _filter:
                        if filtered_dict[x][1][0] == 0:
                            del filtered_dict[x]
                    _filter = []
                    # 此处更新权重影响很小，不更新
                    # cur_weight = assemble.dynamic_weight(filtered_dict)
                    contig, min_pos, max_pos, read_list, new_best_used_kmercount = get_contig(filtered_dict,
                                                                                              best_seed[0],
                                                                                              iteration=1000,
                                                                                              weight=cur_weight)
                    new_seed = best_seed
                    new_contig = contig
                    new_length = len(contig)

                    # 最大长度小于或等于最大阈值的时候，退出循环  考虑突然变短的增大limit,突然过短的情况(<limit_min_length) 该情况会保留上一次的结果；
                    # 最大长度大于最大阈值的时候，持续增大limit
                    # if new_length < 100:  #for test
                    if new_length <= ref_length_max:
                        if new_length >= ref_length_min:
                            best_contig = new_contig
                            best_length = new_length
                            best_seed = new_seed
                            break
                        else:
                            best_contig = old_contig
                            best_length = old_length
                            best_seed = old_seed
                            break
                    else:
                        best_contig = new_contig
                        best_length = new_length
                        best_seed = new_seed

                    # print("temp_limit:{}".format(temp_limit))

                write_contig = True
            else:
                write_contig = True

            if write_contig:
                # print(" " * 100, flush=True, end='\r')  # 之前的信息比较长，冲刷掉
                # print("{0} best_seed: {1}... best_length: {2} ref_length: {3}".format(gene_name, best_seed[0][0:8],
                #                                                                       best_length,
                #                                                                       ref_length))

                message1 = " " * 100
                My_Log_Recorder(message=message1, keep_quiet=keep_quiet, path="", stage="", printout_flag=True,
                                stored_flag=False, make_a_newline=False, use_cutting_line=False)
                message2 = "{0} best_seed: {1}... best_length: {2} ref_length: {3}".format(gene_name, best_seed[0][0:8],best_length,ref_length)
                My_Log_Recorder(message=message2, keep_quiet=keep_quiet, path="", stage="", printout_flag=True,
                                stored_flag=False, make_a_newline=True, use_cutting_line=False)

                with open(contig_path, "w") as f:
                    f.write('>' + unique_identity+ "_" + gene_name + '_contig_k' + str(k2) + '_' + str(best_length) + '\n')
                    f.write(best_contig + '\n')

                break  # 跳出循环
            else:
                # print(" " * 100, flush=True, end='\r')  # 之前的信息比较长，冲刷掉
                # print("{0} limit_count: {1} best_seed: {2}... best_length: {3} ref_length: {4}".format(gene_name,
                #                                                                                        limit_count,
                #                                                                                        best_seed[0][
                #                                                                                        0:8],
                #                                                                                        best_length,
                #                                                                                        ref_length))

                message1 = " " * 100
                My_Log_Recorder(message=message1, keep_quiet=keep_quiet, path="", stage="", printout_flag=True,
                                stored_flag=False, make_a_newline=False, use_cutting_line=False)
                message2 = "{0} limit_count: {1} best_seed: {2}... best_length: {3} ref_length: {4}".format(gene_name,limit_count,best_seed[0][0:8],
                                                                                                       best_length,
                                                                                                       ref_length)
                My_Log_Recorder(message=message2, keep_quiet=keep_quiet, path="", stage="", printout_flag=True,
                                stored_flag=False, make_a_newline=True, use_cutting_line=False)




                with open(short_contig_path, "w") as f:
                    f.write('>' + unique_identity+ "_" + gene_name + '_short_contig_k' + str(k2) + '_' + str(best_length) + '\n')
                    f.write(best_contig + '\n')

                del filtered_dict
                gc.collect()  # 释放内存

                # 动态limit，优化空间还很大
                if limit_count > 5:
                    limit_count = 5
                elif limit_count >= 5:
                    limit_count -= 1
                elif limit_count > 2:
                    limit_count -= 1
                else:
                    break  # limit_count=2 跳出循环

    # 静态Limit
    else:
        limit_count = int(limit_count)
        limit_count, deepth = static_limit_v1(filtered_dict, ref_length, k2, limit_count, reads_length,
                                              ref_rate=15, list_size=1024)

        _filter = [x for x in filtered_dict if filtered_dict[x][0] <= limit_count]
        for x in _filter:
            if filtered_dict[x][1][0] == 0:
                del filtered_dict[x]
        _filter = []
        # seed_list = [(x, ref_dict[x][0], ref_dict[x][1]) for x in ref_dict if ref_dict[x][0] > 0 and ref_dict[x][1] > 0]  #  新增位置信息,仅考虑ref kmercount
        seed_list = [(x, ref_dict[x][0] * math.log10((filtered_dict[x][0] + 1)), ref_dict[x][1]) for x in ref_dict if
                     ref_dict[x][0] > 0 and ref_dict[x][1] > 0]  # 新增位置信息,同时考虑reads,ref kmercount

        if seed_list == []:
            return 0

        list.sort(seed_list, key=lambda x: x[1], reverse=True)

        cur_weight = dynamic_weight(filtered_dict) if filtered_dict else 0.5  # 0.52
        best_seed = seed_list[0]  # ('TTTGGTGGAGTTGTTGG', 95, 81.74894261950934)
        best_contig, min_pos, max_pos, read_list, best_used_kmercount = get_contig(filtered_dict, best_seed[0],
                                                                                   iteration=1000,
                                                                                   weight=cur_weight)  # for test

        # real_deepth = best_used_kmercount * reads_length / (reads_length - assemble_kmer + 1) / ref_length
        # print("best_used_kmercount:{0}  real_length:{1}".format(best_used_kmercount,real_deepth))

        best_length = len(best_contig)

        # 对于过短的处理,更换seed
        if best_length < ref_length_min:
            # 更换seed
            change_time = 0
            last_seed = seed_list[0]
            for seed in seed_list[1:]:
                if last_seed[0][1:] == seed[0][:-1] or last_seed[0][2:] == seed[0][:-2]:  # 移动1~2bp
                    last_seed = seed
                    continue
                change_time += 1
                last_seed = seed
                # print("last_seed:{}".format(last_seed))
                new_contig, min_pos, max_pos, read_list, new_best_used_kmercount = get_contig(filtered_dict,
                                                                                              last_seed[0],
                                                                                              iteration=1000,
                                                                                              weight=cur_weight)

                # new_contig = get_contig(filtered_dict, best_seed[0], iteration=1000, weight=cur_weight)[0]
                new_length = len(new_contig)

                # print(" " * 100, flush=True, end='\r')  # 之前的信息比较长，冲刷掉
                # print("{0}: use another seed: {1} ... {2:<2}".format(gene_name, last_seed[0][0:8], change_time),
                #       flush=True,
                #       end='\r')

                message1 = " " * 100
                My_Log_Recorder(message=message1, keep_quiet=keep_quiet, path="", stage="", printout_flag=True,
                                stored_flag=False, make_a_newline=False, use_cutting_line=False)
                message2 = "{0}: use another seed: {1} ... {2:<2}".format(gene_name, last_seed[0][0:8], change_time)
                My_Log_Recorder(message=message2, keep_quiet=keep_quiet, path="", stage="", printout_flag=True,
                                stored_flag=False, make_a_newline=False, use_cutting_line=False)

                if new_length >= ref_length_min:
                    best_seed = last_seed
                    best_contig = new_contig
                    best_length = new_length
                    write_contig = True
                    break
                elif change_time >= change_seed:  # 超过最大更换种子数目
                    # 种子为未达标情况下，最佳seed
                    break
                else:  # 未达标情况下的最好结果
                    if new_length > best_length:
                        best_contig = new_contig
                        best_length = new_length
                        best_seed = last_seed


        # 对于过长的处理，不更换seed，改为增加limit
        elif best_length > ref_length_max:
            if limit_count >= 32:
                limit_count_max = limit_count + 32
            else:
                limit_count_max = 33
            # 加大limit
            for temp_limit in range(limit_count + 1, limit_count_max):
                # 备份上一次的结果
                old_seed = best_seed
                old_contig = best_contig
                old_length = len(old_contig)
                # 本次的结果，后续会覆盖刷新
                new_seed = best_seed
                new_contig = best_contig
                new_length = len(new_contig)

                # 对filtered_dict再次精简
                _filter = [x for x in filtered_dict if filtered_dict[x][0] <= temp_limit]
                for x in _filter:
                    if filtered_dict[x][1][0] == 0:
                        del filtered_dict[x]
                _filter = []
                # 此处更新权重影响很小，不更新
                # cur_weight = assemble.dynamic_weight(filtered_dict)
                contig, min_pos, max_pos, read_list, new_best_used_kmercount = get_contig(filtered_dict, best_seed[0],
                                                                                          iteration=1000,
                                                                                          weight=cur_weight)

                new_seed = best_seed
                new_contig = contig
                new_length = len(contig)

                # 最大长度小于或等于最大阈值的时候，退出循环  考虑突然变短的增大limit,突然过短的情况(<limit_min_length) 该情况会保留上一次的结果；
                # 最大长度大于最大阈值的时候，持续增大limit
                # if new_length < 100:  #for test
                if new_length <= ref_length_max:
                    if new_length >= ref_length_min:
                        best_contig = new_contig
                        best_length = new_length
                        best_seed = new_seed
                        break
                    else:
                        best_contig = old_contig
                        best_length = old_length
                        best_seed = old_seed
                        break
                else:
                    best_contig = new_contig
                    best_length = new_length
                    best_seed = new_seed
            write_contig = True
        else:
            write_contig = True

        if write_contig:
            # print(" " * 100, flush=True, end='\r')  # 之前的信息比较长，冲刷掉
            # print("{0} best_seed: {1}... best_length: {2} ref_length: {3}".format(gene_name, best_seed[0][0:8],
            #                                                                       best_length,
            #                                                                       ref_length))

            message1 = " " * 100
            My_Log_Recorder(message=message1, keep_quiet=keep_quiet, path="", stage="", printout_flag=True,
                            stored_flag=False, make_a_newline=False, use_cutting_line=False)

            message2 = "{0} best_seed: {1}... best_length: {2} ref_length: {3}".format(gene_name, best_seed[0][0:8],
                                                                                  best_length,
                                                                                  ref_length)
            My_Log_Recorder(message=message2, keep_quiet=keep_quiet, path="", stage="", printout_flag=True,
                            stored_flag=False, make_a_newline=True, use_cutting_line=False)

            with open(contig_path, "w") as f:
                f.write('>'  + unique_identity+ "_" + gene_name + '_contig_k' + str(k2) + '_' + str(best_length) + '\n')
                f.write(best_contig + '\n')
        else:
            # print(" " * 100, flush=True, end='\r')  # 之前的信息比较长，冲刷掉
            # print("{0} best_seed: {1}... best_length: {2} ref_length: {3}".format(gene_name, best_seed[0][0:8],
            #                                                                       best_length,
            #                                                                       ref_length))

            message1 = " " * 100
            My_Log_Recorder(message=message1, keep_quiet=keep_quiet, path="", stage="", printout_flag=True,
                            stored_flag=False, make_a_newline=False, use_cutting_line=False)

            message2 ="{0} best_seed: {1}... best_length: {2} ref_length: {3}".format(gene_name, best_seed[0][0:8],
                                                                                  best_length,
                                                                                  ref_length)
            My_Log_Recorder(message=message2, keep_quiet=keep_quiet, path="", stage="", printout_flag=True,
                            stored_flag=False, make_a_newline=True, use_cutting_line=False)


            with open(short_contig_path, "w") as f:
                f.write('>' + unique_identity + "_" + gene_name + '_short_contig_k' + str(k2) + '_' + str(best_length) + '\n')
                f.write(best_contig + '\n')

    if scaffold_or_not:
        filtered_dict = copy.deepcopy(filtered_dict_back)
        _filter = [x for x in filtered_dict if filtered_dict[x][0] <= limit_count]
        for x in _filter:
            if filtered_dict[x][1][0] == 0:
                del filtered_dict[x]
        _filter = []

        seed_list = [(x, ref_dict[x][0], ref_dict[x][1]) for x in ref_dict if
                     ref_dict[x][0] > 0 and ref_dict[x][1] > 0]  # 新增位置信息
        if seed_list == []:
            return 0
        list.sort(seed_list, key=lambda x: x[1], reverse=True)

        # 利用相似度计算权重 filtered reads 与ref 的大致距离
        cur_weight = dynamic_weight(filtered_dict) if filtered_dict else 0.5  # 0.52
        scaffold = get_scaffold(filtered_dict, seed_list, limit_count, ref_length, weight=cur_weight)
        with open(scaffold_path, 'w') as out:
            scaffold_length = len(scaffold.replace("N", ""))
            if scaffold_length >= ref_length_min:
                out.write('>' + unique_identity+ "_" +gene_name + '_scaffold_k' + str(k2) + '_' + str(scaffold_length) + '\n')
                out.write(scaffold + '\n')
                message="{}: try scaffold. scaffold length: {}".format(gene_name, scaffold_length)
                My_Log_Recorder(message=message, keep_quiet=keep_quiet, path="", stage="", printout_flag=True,
                                stored_flag=False, make_a_newline=True, use_cutting_line=False)

                # print("{}: try scaffold. scaffold length: {}".format(gene_name, scaffold_length))
            else:
                message = "{}: try scaffold. failed".format(gene_name)
                My_Log_Recorder(message=message, keep_quiet=keep_quiet, path="", stage="", printout_flag=True,
                                stored_flag=False, make_a_newline=True, use_cutting_line=False)
                # print("{}: try scaffold. failed".format(gene_name))

    if is_exist(contig_path) or is_exist(scaffold_path):
        if is_exist(short_contig_path):
            os.remove(short_contig_path)

    # 记录assemble csv文件
    assemble_csv_path = os.path.join(assembled_out_path, "assemble.csv")
    #     header = ["gene", "richness", "limit","seed","k2","ref_length", "short_contig_length", "contig_length", "scaffold_length"]
    _gene = gene_name
    _richness = round(deepth, 2)
    _limit = limit_count
    _seed = best_seed[0]  # ('TTTGGTGGAGTTGTTGG', 95, 81.74894261950934)
    _k2 = k2
    _ref_length = ref_length

    if is_exist(short_contig_path):
        _short_contig_length = get_length(short_contig_path)
    else:
        _short_contig_length = "None"
    if is_exist(contig_path):
        _contig_length = get_length(contig_path)
    else:
        _contig_length = "None"
    if is_exist(scaffold_path):
        _scaffold_length = get_length(scaffold_path)
    else:
        _scaffold_length = "None"
    sth = [_gene, _richness, _limit, _seed, _k2, ref_length, _short_contig_length, _contig_length, _scaffold_length]
    mylog(assemble_csv_path, sth)


def do_assemble_parallel(ref_all_path, filtered_all_path, k2, limit_count, limit_min_length,
                         limit_max_length, ref_length_dict, assembled_out_path, reads_length, thread_number,
                         unique_identity, change_seed=32,
                         scaffold_or_not=True,keep_quiet=False):
    files = get_fasta_file(filtered_all_path)
    if files == []:
        return 0
    task_pool = []
    results = []
    if not os.path.isdir(assembled_out_path):
        os.mkdir(assembled_out_path)
    short_contig_path = os.path.join(assembled_out_path, "short_contig")
    contig_path = os.path.join(assembled_out_path, "contig")
    scaffold_path = os.path.join(assembled_out_path, "scaffold")

    if not os.path.isdir(short_contig_path):
        os.mkdir(short_contig_path)
    if not os.path.isdir(contig_path):
        os.mkdir(contig_path)
    if not os.path.isdir(scaffold_path) and scaffold_or_not:
        os.mkdir(scaffold_path)

    assemble_csv_path = os.path.join(assembled_out_path, "assemble.csv")
    header = ["gene", "richness", "limit", "seed", "k2", "ref_length", "short_contig_length", "contig_length",
              "scaffold_length"]
    mylog(assemble_csv_path, header)

    executor = ProcessPoolExecutor(max_workers=thread_number)
    for i in files:
        name = get_basename(i)
        ref_path = os.path.join(ref_all_path, name + ".fasta")
        task_pool.append(
            executor.submit(do_assemble, ref_path, i, k2, limit_count, limit_min_length, limit_max_length,
                            ref_length_dict, assembled_out_path, reads_length, unique_identity, change_seed,
                            scaffold_or_not,keep_quiet))
    for i in task_pool:
        results.append(i.result())
    executor.shutdown()


def my_assemble_main(configuration_information):
    thread_number = configuration_information["thread_number"]
    k2 = configuration_information["k2"]
    assembled_out_path = configuration_information["assembled_out_path"]
    reference = configuration_information["reference"]
    filtered_path = configuration_information["input"]
    limit_count = configuration_information["limit_count"]
    limit_min_length = configuration_information["limit_min_length"]
    limit_max_length = configuration_information["limit_max_length"]
    change_seed = configuration_information["change_seed"]
    scaffold_or_not = configuration_information["scaffold_or_not"]
    out_dir = configuration_information["out_dir"]
    quiet=configuration_information["quiet"]

    reads_length = 150  # 默认值
    if out_dir:
        unique_identity = os.path.split(out_dir)[-1]
    else:
        unique_identity = "GM"

        # print('======================== Assemble =========================')
    # cutting_line(" Assemble ")
    My_Log_Recorder(message=" Assemble ", keep_quiet=quiet, path="", stage="", printout_flag=True,
                    stored_flag=False, make_a_newline=False, use_cutting_line=True)

    t1 = time.time()
    ref_length_dict = {}
    ref_length_dict = get_ref_info(reference, ref_length_dict)
    files = get_fasta_file(filtered_path)

    if files != []:
        reads_length = get_reads_length(files[0], type="fasta", _number=100)

    do_assemble_parallel(reference, filtered_path, k2, limit_count, limit_min_length, limit_max_length, ref_length_dict,
                         assembled_out_path, reads_length, thread_number, unique_identity, change_seed, scaffold_or_not,quiet)

    t2 = time.time()
    used_temp_time = format((t2 - t1), ".2f")
    message="Assemble time used: {}s".format(used_temp_time)
    My_Log_Recorder(message=message, keep_quiet=quiet, path="", stage="", printout_flag=True,
                    stored_flag=False, make_a_newline=True, use_cutting_line=False)

    # print("Assemble time used: {}s".format(used_temp_time))


if __name__ == '__main__':
    # -i -r 中的文件名必须是xx.fasta 同名
    pars = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="assemble",
                                   usage="%(prog)s <-i> <-r> <-o> [options]")

    pars.add_argument("-o", "--out", dest="out", help="Specify the result folder <dir>",
                      metavar="", required=True)

    pars.add_argument("-k2", "--kmer2", dest="kmer2", help="length of kmer[default=31]",
                      default=31, type=int, metavar="")
    pars.add_argument("-t", "--thread", metavar="", dest="thread_number", help="Thread", type=int, default=1)
    pars.add_argument("-r", "--reference", metavar="", dest="reference", type=str, help="references  <dir>",
                      required=True)
    pars.add_argument("-i", "--input", metavar="", dest="input", type=str, help="filtered reads  <dir>", required=True)

    pars.add_argument('-limit_count', metavar='', dest='limit_count', help='''limit of kmer count,[default=auto]''',
                      required=False,
                      default='auto')
    pars.add_argument('-limit_min_length', metavar='', dest='limit_min_length', type=float,
                      help='''limit of contig length''',
                      required=False, default=0.5)
    pars.add_argument('-limit_max_length', metavar='', dest='limit_max_length', type=float,
                      help='''limit of contig length''',
                      required=False, default=2)
    pars.add_argument("-change_seed", metavar="", dest="change_seed", type=int,
                      help='''times of changing seed [default=32]''', required=False,
                      default=32)
    pars.add_argument('-scaffold', metavar="", dest="scaffold", type=str, help='''make scaffold''', default=False)
    pars.add_argument("-quiet", dest="quiet", help="Do not write progress messages to stderr", default=False,action='store_true')

    args = pars.parse_args()

    # 初始化操作
    thread = args.thread_number
    k2 = args.kmer2
    assembled_out_path = args.out
    reference = args.reference
    filtered_path = args.input
    limit_count = args.limit_count
    limit_min_length = args.limit_min_length
    limit_max_length = args.limit_max_length
    change_seed = args.change_seed
    scaffold_or_not = args.scaffold
    quiet=args.quiet
    reads_length = 150  # 默认值

    True_set = ["yes", "y", "Yes", "true", "True", "1", 1]
    if scaffold_or_not in True_set:
        scaffold_or_not = True
    else:
        scaffold_or_not = False

    assemble_configuration_information = {
        "thread_number": thread, "k2": k2,
        "assembled_out_path": assembled_out_path,
        "reference": reference,
        "input": filtered_path,
        "limit_count": limit_count,
        "limit_min_length": limit_min_length, "limit_max_length": limit_max_length,
        "change_seed": change_seed,
        "scaffold_or_not": scaffold_or_not,
        "out_dir": "",
        "quiet":quiet

    }
    my_assemble_main(assemble_configuration_information)
