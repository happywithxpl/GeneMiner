#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2021/10/20 16:47
# @Author  : xiepulin
# @File    : basic.py
# @Software: PyCharm


import sys
import subprocess
import os
from Bio import SeqIO
from Bio import pairwise2
from Bio.SeqRecord import  SeqRecord
from Bio.Seq import  Seq
import platform
import shutil
import re
import time
from collections import  defaultdict

import global_var as gv


################################################################
##################################################################
'''
流程部分
'''
###################################################################
###################################################################

'''
检测终止信号，优雅退出
'''
def signal_handler(signal, frame):
    print("GeneMiner has been terminated")
    gv.set_value("my_gui_flag", 0)
    sys.exit(0)
'''
判断文件是否为fasta格式
'''
def is_fasta(filename):
    try:
        with open(filename, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            flag=any(fasta)  # any:如果都不是fasta 会提示false。但凡有一个是fasta,会提示True。
            # 但每次都是对一个文件进行判定，所以功能得以实现
    except:
        flag=False
    return flag


'''
判断文件是否为Genbank格式
'''
def is_gb(filename):
    try:
        with open(filename, "r") as handle:
            fasta = SeqIO.parse(handle, "gb")
            flag= any(fasta)
    except:
        flag=False
    return flag


'''
判断文件是否是txt文本文件，布尔值
'''
def is_txt_file(file):
    try:
        flag=1
        suffix=os.path.splitext(file)[-1]
        if suffix!=".txt":
            flag=0
    except:
        flag=0
    return flag






'''
为了防止数据质量太差，导致某一步骤做不出来结果，无法接续运行
为了提高程序程序鲁棒性，判断文件是否存在，如果不存在或者文件大小为0，认定为无效返回0  有效返回1 
'''

def is_exist(file):
    if os.path.isfile(file):
        if os.path.getsize(file) > 0:
            flag = 1
        else:
            flag = 0

    elif os.path.isdir(file):
        files = get_files(file)
        if files==[]:
            flag=0
        else:
            flag = 1
            for i in files:
                if os.path.getsize(i) > 0:
                    continue
                else:
                    flag = 0
                    break
    else:
        flag=0

    return flag

'''
文件存在 ：文件存在且大于0
文件夹存在：文件夹存在且不为空
'''
def is_exist_simple(file):
    if os.path.isfile(file):
        if os.path.getsize(file) > 0:
            flag = 1
        else:
            flag = 0

    elif os.path.isdir(file):
        files = get_files(file)
        if files==[]:
            flag=0
        else:
            flag = 1
    else:
        flag=0
    return flag




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
新建文件夹
'''
def dir_make(out_dir):
    if os.path.exists(out_dir) and len(os.listdir(out_dir)) > 0:
        print("\nerror! {}：the folder has already existed,and there are files under the folder.".format(out_dir))
        sys.exit()
    #针对于用户先创建文件夹，再使用
    elif os.path.exists(out_dir) and len(os.listdir(out_dir))==0:
        pass
    else:
        os.makedirs(out_dir)



#分割线 文字内容不得大于分割线
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




def get_basename(file):
    if is_exist(file):
        basename=os.path.basename(file)
        if ".fasta" in basename:
            basename = basename.split(".fasta")[0]
        elif ".fas" in basename:
            basename = basename.split(".fas")[0]
        elif ".fa" in basename:
            basename = basename.split(".fa")[0]
        else:
            basename=basename
        return basename


'''
判断输入参考基因组是文件还是文件夹，0代表文件夹，1代表文件
'''
def file_or_directory(ref):
    if os.path.isdir(ref):
        files = os.listdir(ref)
        if len(files) == 0:
            print("the input  reference genome folder is empty")
            sys.exit()
        else:
            flag = 0
            # print("this is a dir")
    elif os.path.isfile(ref):
        size = os.path.getsize(ref)
        if size == 0:
            print("The input reference genome file is empty")
            sys.exit()
        else:
            flag = 1
            # print("this is a file")
    else:
        sys.exit()
    return flag

# '''
# 获得目标陆路径 所有非空文件
# '''
# def get_file_list(path):
#     file_list=[]
#     if os.path.isdir(path):
#         files=get_files(path)
#         for i in files:
#             if os.path.getsize(i)==0:
#                 print("{} is empty".format(i))
#                 sys.exit()
#             else:
#                 file_list.append(i)
#     elif os.path.isfile(path):
#         size = os.path.getsize(path)
#         if size == 0:
#             print("{}：is empty".format(path))
#             sys.exit()
#         else:
#             file_list.append(path)
#     else:
#         sys.exit()
#     return file_list

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



def check_true(sth):
    true_set=["yes","y","YES","Yes","true","t","TRUE","True","1",True]
    false_set=["no","n","NO","No","false","f","FALSE","False","0",False,None,"None"]
    if sth in true_set:
        sth=True
    elif sth in false_set:
        sth=False
    else:
        sth=0  #实际上0 就是False 这样写有点问题
    return  sth




##############################################################
'''
bootstrap
'''
###########################################################


'''
双序列全局比对，采用动态规划算法。返回双序列全局比对得分，可以用来判断 序列方向和参考序列方向 是否一致
'''
def pairwise(seq1, seq2):
    seq1=str(seq1)
    seq2=str(seq2)

    seq1 = str.upper(seq1)
    seq2 = str.upper(seq2)
    alignment = pairwise2.align.globalms(seq1, seq2, 1, -1, -3, -2)       #同blast准则 75%保守度
    # alignment = pairwise2.align.localms(seq1, seq2, 2, -1, -2, -0.5)
    score = alignment[0][2]
    return score


def cut_align(seq, ref):
    seq=str(seq) #防止输入的是Seq格式（biopython）
    ref=str(ref)
    seq = str.upper(seq)
    ref = str.upper(ref)
    alignment = pairwise2.align.localms(seq, ref, 1, -1, -3, -2)
    seq_align = alignment[0][
        0]  # 第一条双序列比对结果  PMSK_align ,此时和ref_align相同长度  Alignment(seqA='-TCGAAAAAAAAAAAAAAA', seqB='ATCGATCG-----------', score=8.0, start=1, end=5)
    start_site = alignment[0][3]
    end_site = alignment[0][4]
    cut_align_sequence = seq_align[start_site:end_site]

    cut_align_sequence=cut_align_sequence.replace("-","")
    ref_length = len(ref)
    if len(cut_align_sequence) < 0.8 * ref_length:
        trimmed_contigs = ""
    elif len(cut_align_sequence) >= 0.8 * ref_length and len(cut_align_sequence) <= 1.2 * ref_length:  # 剪切的序列既不能太短了，又不能在剪切的时候引入太多空位
        trimmed_contigs = cut_align_sequence
        # trimmed_contigs = trimmed_contigs.replace("-", "")
    else:
        trimmed_contigs = ""
    return trimmed_contigs






###############################################
'''
计算一致度和覆盖度     平均时间 0.42290181159973145s  bootstrap
'''
#################################333


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
def get_seq(fasta_file,max_seq_number=100,seq_count_limit=False):
    infile = open(fasta_file, 'r', encoding='utf-8', errors='ignore')
    seq, name = "", ""
    my_list = []
    seq_number = 0
    while True:
        line = infile.readline()
        line = line.strip()
        line=line.replace("N","")  #特定对于scaffold
        if (line.startswith('>') or not line) and name:  # 保证最后一条序列能能在保存后退出
            temp = {}
            temp = {name: seq}
            my_list.append(temp)
            seq_number += 1
        if line.startswith('>'):
            name = line[1:]
            seq = ''
        else:
            seq += line
        if not line:
            break
        if seq_count_limit and max_seq_number:
            if seq_number >= max_seq_number:
                break
    infile.close()
    return my_list


def get_gap_information(alignment_seq):
    gap_number_list = re.findall(r"-", alignment_seq)  # 统计空位数量
    continuous_gap_list = re.findall(r"-{2,}", alignment_seq)  # 统计连续空位
    gap_number = len(gap_number_list)
    gap_extend = 0  # 将连续gap记为一个gap,被压缩gap数量=gap长度 - 1
    if continuous_gap_list == []:
        gap_extend = gap_extend
    else:
        for i in continuous_gap_list:
            gap_extend += len(i) - 1

    gap_open = gap_number - gap_extend

    continuous_gap_inter_list = re.findall(r"[ACGT]-{2,}[ACGT]", alignment_seq)  # 记录中间 连续gap数量
    single_gap_inter_list = re.findall(r"[ACGT]-[ACGT]", alignment_seq)  # 记录中间 单一gap数量
    continuous_gap_inter_list_new = []
    single_gap_inter_list_new = []
    if continuous_gap_inter_list != []:
        for i in continuous_gap_inter_list:
            temp = re.sub(r"[ACGT]", "", i)
            continuous_gap_inter_list_new.append(temp)
    if single_gap_inter_list != []:
        for i in single_gap_inter_list:
            temp = re.sub(r"[ACGT]", "", i)
            single_gap_inter_list_new.append(temp)
    gap_extend_inter = 0
    if continuous_gap_inter_list_new == []:
        gap_extend_inter = 0
    else:
        for i in continuous_gap_inter_list_new:
            gap_extend_inter += len(i) - 1

    continuous_gap_inter = "".join(continuous_gap_inter_list_new)
    single_gap_inter = "".join((single_gap_inter_list_new))
    continuous_gap_inter_number = len(continuous_gap_inter)
    single_gap_inter_number = len(single_gap_inter)

    start_end_gap_number = gap_number - continuous_gap_inter_number - single_gap_inter_number  # 首尾gap = 总gap - 中间gap(单一+连续)

    # print("con:{}".format(continuous_gap_inter_list))
    # print("con_new:{}".format(continuous_gap_inter_list_new))
    # print("sin:{}".format(single_gap_inter_list))
    # print("sin_new:{}".format(single_gap_inter_list_new))

    return [gap_number, gap_open, gap_extend, gap_extend_inter, start_end_gap_number]


'''
不使用库，gap压缩，free end
'''
def get_identity_v1(seq1, seq2):  # seq1 seq2(ref)
    #打分矩阵 Global alignment with free end gaps   5/-4/-12/-3 (65%)   1/-1/-3/-2 (75%)
    match_socre = 5
    mismatich_score = -4.0
    gap_open_socre = -12
    gap_extend_score = -3
    alignments = pairwise2.align.globalms(seq1, seq2, match_socre, mismatich_score, gap_open_socre,
                                          gap_extend_score,penalize_end_gaps=False,one_alignment_only=True)  # 全局比对，相同的残基就给1分，不同和gap不扣分

    alignment_seq1 = alignments[0][0]
    alignment_seq2 = alignments[0][1]
    length = len(alignment_seq1)
    gap_information1 = get_gap_information(alignment_seq1)
    gap_information2 = get_gap_information(alignment_seq2)
    gap_number, gap_open, gap_extend, gap_extend_inter, start_end_gap_number = gap_information2 if gap_information1[0] <= \
                                                                                 gap_information2[
                                                                                     0] else gap_information1
    match_number = 0
    mismatch_number = 0

    for i in zip(alignment_seq1,alignment_seq2):
        if "-" in i:
            continue
            # insert+=1 #不科学，是两条序列累计的插入，而非一条的插入
        elif i[0] == i[1] and i[0] != "-" and i[1] != "-":
            match_number += 1
        elif i[0] != i[1] and i[0] != "-" and i[1] != "-":
            mismatch_number += 1
        else:
            pass
    identity = round( (match_number / (length - gap_extend_inter-start_end_gap_number)),4)
    return identity



def get_identity(seq1, seq2):  # seq1 seq2(ref)
    #打分矩阵 Global alignment with free end gaps   5/-4/-12/-3 (65%)   1/-1/-3/-2 (75%)
    match_socre = 5
    mismatich_score = -4.0
    gap_open_socre = -12
    gap_extend_score = -3
    alignments = pairwise2.align.globalms(seq1, seq2, match_socre, mismatich_score, gap_open_socre,
                                          gap_extend_score)  # 全局比对，相同的残基就给1分，不同和gap不扣分


    alignment_seq1 = alignments[0][0]
    alignment_seq2 = alignments[0][1]
    score = alignments[0][2]  # 比对得分
    length = len(alignment_seq1)
    gap_information1 = get_gap_information(alignment_seq1)
    gap_information2 = get_gap_information(alignment_seq2)
    gap_number, gap_open, gap_extend, gap_extend_inter, start_end_gap_number = gap_information2 if gap_information1[0] <= \
                                                                                 gap_information2[
                                                                                     0] else gap_information1
    # print( gap_number, gap_open, gap_extend, gap_extend_inter, start_end_gap_number)
    # print("alignment1:{}".format(alignment_seq1))
    # print("alignment2:{}".format(alignment_seq2))

    # x = Symbol("x")  # match_number
    # y = Symbol("y")  # mismatch_number
    # expression = [x + y + gap_number - length,
    #               x * match_socre + y * mismatich_score + gap_open * gap_open_socre + gap_extend * gap_extend_score - score]
    # result = solve(expression, [x, y])
    # match_number = int(result[x])
    match_number=  ( (score -gap_open * gap_open_socre - gap_extend * gap_extend_score) - (mismatich_score * (length-gap_number))  )/(match_socre-mismatich_score)
    mismatch_number = ( match_socre * (length-gap_number) + gap_open * gap_open_socre + gap_extend * gap_extend_score - score ) / (match_socre-mismatich_score)
    match_number= int(match_number)
    mismatch_number=int(mismatch_number)
    # print(match_number,mismatch_number)
    identity = round((match_number / (length - gap_extend_inter-start_end_gap_number)), 4)
    return identity



'''
gap compress and free end
get mutation model
'''
def get_identity_and_mutate_model(seq1,seq2):
    # 打分矩阵 Global alignment with free end gaps   5/-4/-12/-3 (65%)   1/-1/-3/-2 (75%)
    match_socre = 5
    mismatich_score = -4.0
    gap_open_socre = -12
    gap_extend_score = -3
    alignments = pairwise2.align.globalms(seq1, seq2, match_socre, mismatich_score, gap_open_socre,
                                          gap_extend_score, penalize_end_gaps=False,
                                          one_alignment_only=True)  # 全局比对，相同的残基就给1分，不同和gap不扣分

    alignment_seq1 = alignments[0][0]
    alignment_seq2 = alignments[0][1]
    length = len(alignment_seq1)
    gap_information1 = get_gap_information(alignment_seq1)
    gap_information2 = get_gap_information(alignment_seq2)
    gap_number, gap_open, gap_extend, gap_extend_inter, start_end_gap_number = gap_information2 if gap_information1[
                                                                                                       0] <= \
                                                                                                   gap_information2[
                                                                                                       0] else gap_information1

    match_number = 0
    mismatch_number = 0
    # mutate_set=["AA","AT","AC","AG","TA","TT,","TC","TG","CA","CT","CC","CG","GA","GT","GC","GG"]
    mutate_number_dict=defaultdict(int)
    for i in zip(alignment_seq1, alignment_seq2):
        my_set=i[0]+i[1]
        if "-" in i:
            continue
            # insert+=1 #不科学，是两条序列累计的插入，而非一条的插入
        elif i[0] == i[1] and i[0] != "-" and i[1] != "-":
            match_number += 1
            mutate_number_dict[my_set]+=1

        elif i[0] != i[1] and i[0] != "-" and i[1] != "-":
            mismatch_number += 1
            mutate_number_dict[my_set] += 1
        else:
            pass

    mutate_a_number =mutate_number_dict["AA"]+mutate_number_dict["AT"]+mutate_number_dict["AC"]+mutate_number_dict["AG"]
    mutate_t_number = mutate_number_dict["TA"] + mutate_number_dict["TT"] + mutate_number_dict["TC"] + mutate_number_dict["TG"]
    mutate_c_number = mutate_number_dict["CA"] + mutate_number_dict["CT"] + mutate_number_dict["CC"] + mutate_number_dict["CG"]
    mutate_g_number = mutate_number_dict["GA"] + mutate_number_dict["GT"] + mutate_number_dict["GC"] + mutate_number_dict["GG"]

    ##需要谨慎考虑ATCG不全的情况
    if mutate_a_number !=0:
        mutate_aa,mutate_at,mutate_ac,mutate_ag = mutate_number_dict["AA"]/mutate_a_number, (mutate_number_dict["AA"]+mutate_number_dict["AT"])/mutate_a_number,(mutate_number_dict["AA"]+mutate_number_dict["AT"]+mutate_number_dict["AC"])/mutate_a_number,1
    else:
        mutate_aa,mutate_at,mutate_ac,mutate_ag = 0,0,0,1
    if mutate_t_number !=0:
        mutate_ta,mutate_tt,mutate_tc,mutate_tg=mutate_number_dict["TA"]/mutate_t_number,(mutate_number_dict["TA"] + mutate_number_dict["TT"])/mutate_t_number,( mutate_number_dict["TA"] + mutate_number_dict["TT"] + mutate_number_dict["TC"])/mutate_t_number,1
    else:
        mutate_ta, mutate_tt, mutate_tc, mutate_tg=0,0,0,1
    if mutate_c_number!=0:
        mutate_ca, mutate_ct, mutate_cc, mutate_cg=mutate_number_dict["CA"]/mutate_c_number,(mutate_number_dict["CA"] + mutate_number_dict["CT"])/mutate_c_number,(mutate_number_dict["CA"] + mutate_number_dict["CT"] + mutate_number_dict["CC"])/mutate_c_number,1
    else:
        mutate_ca, mutate_ct, mutate_cc, mutate_cg=0,0,0,1
    if mutate_g_number!=0:
        mutate_ga, mutate_gt, mutate_gc, mutate_gg= mutate_number_dict["GA"]/mutate_g_number,(mutate_number_dict["GA"] + mutate_number_dict["GT"])/mutate_g_number,(mutate_number_dict["GA"] + mutate_number_dict["GT"] + mutate_number_dict["GC"])/mutate_g_number,1
    else:
        mutate_ga, mutate_gt, mutate_gc, mutate_gg=0,0,0,1


    my_mutate_model={"A":[ mutate_aa,mutate_at,mutate_ac,mutate_ag],
                     "T":[ mutate_ta,mutate_tt,mutate_tc,mutate_tg],
                     "C":[mutate_ca, mutate_ct, mutate_cc, mutate_cg],
                     "G":[ mutate_ga, mutate_gt, mutate_gc, mutate_gg]


    }
    identity = match_number /(length - gap_extend_inter - start_end_gap_number)

    return (identity,my_mutate_model)


'''
得到系统类型
'''
def get_platform():
    current_os = platform.system().lower()   #linux windows darwin
    return  current_os


'''
得到绝对路径
'''
def get_absolute(path):
    if path==None:
        return None
    else:
        if os.path.isabs(path):
            return path
        else:
            path = os.path.abspath(path)
            return path


'''
获得目标文件夹下的fasta文件，仅第一层
'''
def get_fasta_file(path):
    path=get_absolute(path)
    files=os.listdir(path)
    suffix_list=["fa","fas","fasta"]
    files_list=[]
    for i in files:
        suffix=i.split(".")[-1].lower()
        file_path= os.path.join(path,i)
        if os.path.isfile(file_path) and suffix in suffix_list :
            files_list.append(file_path)
    return  files_list


def run_command(cmd):
    subprocess.call(cmd,shell=True)

# def run_command(cmd):
#     os.system(cmd)

# def run_command(cmd):
#     subprocess.Popen(cmd,shell=True)











##################################
'''
GUI换行输出
'''
def ML_str_re(in_str=""):
    in_str = in_str.replace("\r", "\r\n")
    global del_count
    in_str_list = in_str.split("\n")
    out_str = ""
    for i in in_str_list:
        if i.endswith('\r') == False:
            out_str += i + "\n"
    return out_str
#################################3







#记录日志信息
import csv
def mylog_v1(file_path,sth):
    with open(file_path,"a",newline='') as f:  #不同操作系统换行符不一样 linux:\n  windows:\r\n  mac \r newline参数可以统一换行符  否则中间会出现默认空白行
        writer=csv.writer(f)
        writer.writerow(sth)

def mylog(file_path,sth,row=True):
    if row==True:
        with open(file_path,"a",newline='') as f:  #不同操作系统换行符不一样 linux:\n  windows:\r\n  mac \r newline参数可以统一换行符  否则中间会出现默认空白行
            writer=csv.writer(f)
            writer.writerow(sth)
    else:
        with open(file_path, "a",newline='') as f:  # 不同操作系统换行符不一样 linux:\n  windows:\r\n  mac \r newline参数可以统一换行符  否则中间会出现默认空白行
            writer = csv.writer(f)
            writer.writerows(sth)


























