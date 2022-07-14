#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2021/10/20 16:53
# @Author  : xiepulin
# @File    : verify_parameters.py
# @Software: PyCharm

import multiprocessing
from basic import *
from global_var import get_value ,set_value, get_init
###############################################

'''
检测终止信号，优雅退出
'''
def signal_handler(signal, frame):
    print("GeneMiner has been terminated")
    set_value("my_gui_flag", 0)
    sys.exit(0)

'''
检测系统环境部分
(1)python版本
'''
##############################################
def check_python_version():
    this_python = sys.version_info[:2]
    min_version = (3, 6)
    if this_python < min_version:
        message_parts = [
            "GeneMiner does not work on Python {}.{}.".format(*this_python),            #*可以用来接受元组
            "The minimum supported Python version is {}.{}.".format(*min_version),
            "Please download the eligible version from https://www.python.org/.".format(*this_python)]
        print("ERROR: " + " ".join(message_parts))
        set_value("my_gui_flag", 0)
        sys.exit()

##############################################################





##############################################################
'''
第一部分 检测各项参数是否正确,正确后打印参数
1.检测线程数量默认为所有线程数量的1/4 默认不大于4
2限定kmer的大小 k1过滤  k2拼接
3 限定步长
4 限定limit count
5 限定limit length
6 限定change_seed
7 检测软边界 默认为75 范围为0~200之间
8 检测长度
9 参考基因组格式是否对应  
10 检测数据量大小 -d
11 检测输入
12.检测自展参数
13.检测输出文件夹 out

'''
#############################################################


'''
1.检测线程数量默认为所有线程数量的1/4 默认不大于4
'''
def check_threads_number(thread):
    thread_number = 1  # 初始化，至少一个
    thread_number_all = multiprocessing.cpu_count()
    if thread == "auto":
        thread_number = int(thread_number_all * 0.25)
        if thread_number == 0:
            thread_number = thread_number + 1
        elif thread_number>=4:
            thread_number=4
        else:
            thread_number=thread_number
    else:
        thread=str(thread)
        if thread.isdigit():         #判断是否为纯数字
            thread = int(thread)
            if thread > thread_number_all:
                print("Number of threads exceed the  maximum, please check the -t parameter")
                set_value("my_gui_flag", 0)
                sys.exit()
            elif thread <= 0:
                print("Number of threads shuold exceed  0, please check the -t parameter")
            else:
                thread_number = thread

        else:  #防止输入乱七八糟的东西
            print("The number of threads must be an integer greater than 0, please check the -t parameter")
            set_value("my_gui_flag", 0)
            sys.exit()

    return thread_number

'''
2限定kmer的大小 k1过滤  k2拼接
2.1 filter k1 [17,127]
'''
def check_k1(k1):
    max = 127
    min = 17
    if k1 < min or k1 > max:
        print("K1 ranges from 17 to 127, please check the -k1 parameter")
        set_value("my_gui_flag", 0)
        sys.exit()
    else:
        k2 = k1
    return k2

'''
2.2 assemble k2 [17,127]
'''
def check_k2(k2):
    max=127
    min=17
    if k2<min or k2>max:
        print("K2 ranges from 17 to 127, please check the -k2 parameter")
        set_value("my_gui_flag", 0)
        sys.exit()
    else:
        k2=k2
    return k2


'''
3 限定步长
'''
def check_step_length(step_length):
    max=10
    min=1
    if step_length<min or step_length>max:
        print("Step_length ranges from 1 to 10, please check the -step_length parameter")
        set_value("my_gui_flag", 0)
        sys.exit()
    else:
        step_length=step_length
    return step_length

'''
4 限定limit count
'''
def check_limit_count(limit_count):
    if limit_count.isdigit(): #数字型
        limit_count=int(limit_count)
        min = 0
        if limit_count < min:
            print("Limit_count cannot be less than 0, please check the -limit_count parameter")
            set_value("my_gui_flag", 0)
            sys.exit()

    elif isinstance(limit_count,str): #字符型
        limit_count=str(limit_count)
        set = ["auto", "AUTO", "Auto"]
        if limit_count not in set:
            print("Unsupported value encountered, please check the -limit_count parameter")
            set_value("my_gui_flag", 0)
            sys.exit()
        else:
            limit_count="auto"
    else:
        print("Unsupported value encountered, please check the -limit_count parameter")
        set_value("my_gui_flag", 0)
        sys.exit()

    return limit_count



'''
5 限定limit length
'''
def check_limit_length(limit_min_length,limit_max_length):
    min = 0
    if limit_min_length <= min:
        print("Limit_min_length is greater than 0, please check the -limit_min_length parameter")
        set_value("my_gui_flag", 0)
        sys.exit()
    else:
        limit_min_length=limit_min_length

    if limit_max_length <= min:
        print("Limit_max_length is greater than 0, please check the -limit_max_length parameter")
        set_value("my_gui_flag", 0)
        sys.exit()
    else:
        limit_max_length=limit_max_length
    return limit_min_length,limit_max_length


'''
6 限定change_seed
'''
def check_change_seed(change_seed):
    min=0
    if change_seed<min :
        print("Change_seed cannot be less than 0, please check the -change_seed parameter")
        set_value("my_gui_flag", 0)
        sys.exit()
    else:
        change_seed=change_seed
    return change_seed


def check_scaffold(scaffold):
    scaffold_or_not=check_true(scaffold)
    if  scaffold_or_not == True or scaffold_or_not == False:
        scaffold_or_not=scaffold_or_not
    else:
        print("Unsupported value encountered, please check the -scaffold parameter")
        set_value("my_gui_flag", 0)
        sys.exit()
    return scaffold_or_not



'''
7 检测软边界 默认为75 范围为0~200之间
'''
def check_soft_boundary(soft_boundary):
    soft_boundary = int(soft_boundary)
    max = 200
    min = 0
    if soft_boundary > max or soft_boundary < min:
        print(
            "The length of the soft boundary ranges from 0 to 200, and the recommended length is 0.5 * reads_length, please check the -b parameter")
        set_value("my_gui_flag", 0)
        sys.exit()
    else:
        soft_boundary = soft_boundary
    return soft_boundary

'''
8 检测长度
'''
def check_max_min_length(max_length,min_length):
    if max_length<= min_length:
        print("The maximum gene length should be greater than the minimum gene length, please check the -max parameter ")
        set_value("my_gui_flag", 0)
        sys.exit()

    if max_length<=0:
        print("The maximum gene length should be greater than zero, please check the -max parameter")
        set_value("my_gui_flag", 0)
        sys.exit()

    if min_length <= 0:
        print("The minimum gene length should not be less than zero, please check the -min parameter")
        set_value("my_gui_flag", 0)
        sys.exit()




'''
9 参考基因组格式是否对应  
(1)二选一
（2）路径是否真实存在
（3）后缀格式
'''
def check_reference(target_reference_fa, target_reference_gb):
    #二选一
    temp=[target_reference_fa,target_reference_gb]
    ref_list = [i for i in temp if i!=None]

    if ref_list==[]:
        print("Please input reference,check the -rtfa,-rtgb ")
        set_value("my_gui_flag", 0)
        sys.exit()
    if len(ref_list)==2:
        print("Please choose only one parameter from -rtfa and -rtgb")
        set_value("my_gui_flag", 0)
        sys.exit()




    #路径是否真实存在
    for i in ref_list:
        if not is_exist_simple(i):  #没用is_exist 用的is_exit_simple
            print("{}: xx the file does not exist".format(i))
            set_value("my_gui_flag", 0)
            sys.exit()
    # 格式是否正确   0代表文件夹 1代表文件
    if target_reference_fa in ref_list:
        Nonconforming_file = []  # 记录不合格的文件
        files=get_file_list(target_reference_fa)
        for  file in  files:
            if not is_fasta(file):
                Nonconforming_file.append(file)
        if Nonconforming_file!=[]:
            print("references should be in fasta format,please check -rtfa parameter")
            Nonconforming_file = [os.path.basename(file) for file in Nonconforming_file]
            if len(Nonconforming_file) == 1:
                Nonconforming_file = "".join(Nonconforming_file)
                print("{} is not in fasta format".format(Nonconforming_file))
            else:
                Nonconforming_file = ",".join(Nonconforming_file)
                print("{} are not in fasta format".format(Nonconforming_file))
            set_value("my_gui_flag", 0)
            sys.exit()
    elif target_reference_gb in ref_list:
        Nonconforming_file = []  # 记录不合格的文件
        files = get_file_list(target_reference_gb)
        for file in files:
            if not is_gb(file):
                Nonconforming_file.append(file)
        if Nonconforming_file !=[]:
            print("references should be in GenBank-format, please check -rtgb parameter")
            Nonconforming_file = [os.path.basename(file) for file in Nonconforming_file]
            if len(Nonconforming_file) == 1:
                Nonconforming_file = "".join(Nonconforming_file)
                print("{} is not in GenBank-format".format(Nonconforming_file))
            else:
                Nonconforming_file = ",".join(Nonconforming_file)
                print("{} are not in GenBank-format".format(Nonconforming_file))
            set_value("my_gui_flag", 0)
            sys.exit()
    else:
        pass



'''
10 检测数据量大小 -d
'''
def check_datasize(data_size):
    data_size=str(data_size) # 先把数字或者"all" 都转为字符串处理   防止'int' object has no attribute 'isdigit'

    min_data_size = 1000000
    if data_size.lower()=="all":
        ultimate_data_size=data_size.lower()
        return ultimate_data_size

    elif data_size.isdigit():
        data_size=int(data_size)
        if data_size < min_data_size:
            print(
                "Please check -n parameter. for better results, input data should not be less than {0} lines.".format(
                    min_data_size))
            set_value("my_gui_flag", 0)
            sys.exit()
        else:
            ultimate_data_size = data_size - data_size % 100000  # 10w起步，保证是4的倍数
        return ultimate_data_size
    else:
        print("Please check -n parameter, it must be an integer greater than {} or 'all' ".format(min_data_size))
        set_value("my_gui_flag", 0)
        sys.exit()


'''
11 检测输入
 -1 -2   -single 是否真实存在，是否正确选择
'''
def check_input(data1, data2, single):
    #是否存在
    if (data1 == None or data2 == None) and single == None:
        print("Please choose only one parameter from <-1 -2> and <-s>")
        set_value("my_gui_flag", 0)
        sys.exit()

    if data1 and data2  and single:
        print("Please choose only one parameter from <-1 -2> and <-s>")
        set_value("my_gui_flag", 0)
        sys.exit()

    if (data1 == None and data2 != None) or (data1 != None and data2 == None):
        print("You must choose both the -1 and -2 ")
        set_value("my_gui_flag", 0)
        sys.exit()

    #后缀是否正确
    temp=[data1,data2,single]
    input_data=[i for i in temp if i!=None]

    for i in input_data:
        if not is_exist(i):
            print("{}:the file does not exist".format(i))
            set_value("my_gui_flag", 0)
            sys.exit()
        if ".fastq.gz" in i or ".fq.gz" in i or ".fq" in i  or ".fastq" in i:
            pass
        else:
            print("{}:the file need .fq/.fastq/.fq.gz/.fastq.gz as the suffix".format(i))
            set_value("my_gui_flag", 0)
            sys.exit()







'''
12.检测自展参数 -bn
'''
def check_bootstrap_parameter(bootstrap_number):
    if bootstrap_number==None or bootstrap_number=="None":
        flag = "No"
        bootstrap_number = "None"
        bootstrap_information = [flag, bootstrap_number]
    else:
        if  bootstrap_number <= 0:
            print("The number must be greater than 0, please check the -bn parameter")
            set_value("my_gui_flag", 0)
            sys.exit()
        else:
            flag = "Yes"
            bootstrap_number = bootstrap_number
            bootstrap_information = [flag, bootstrap_number]

    return bootstrap_information



'''
13.检测输出文件夹 out
'''
def check_out_dir(out):
    if os.path.isdir(out):
        if len(os.listdir(out)) == 0:  # 文件夹存在但里面为空是能够使用的
            out_dir_name = out
        else:
            print("{} already exists and there are files under the folder, please check the -o parameter".format(out))
            set_value("my_gui_flag", 0)
            sys.exit(0)
    else:
        out_dir_name = out

    return out_dir_name



'''
参数校验合格之后，打印参数表
在extract_raw_data之前
'''
def print_parameter_information(parameter_information_dict):
    temp={}
    for key,value in parameter_information_dict.items():
        if parameter_information_dict[key] != None:
            temp[key] = value
    message_all = []
    for key, value in temp.items():
        message = "{0:<22}:{1:<}".format(key, value) #左对齐 22位宽  target_reference最长 16位宽
        message_all.append(message)

    symbol="-"
    print(symbol*22,flush=True)
    header="GeneMiner: a software for extracting phylogenetic markers from next generation sequencing data\n" \
           "Version: 1.0\n" \
           "Copyright (C) 2022 Pulin Xie\n" \
           "Please contact <xiepulin@stu.edu.scu.cn>, if you have any bugs or questions"

    print(header,flush=True)
    for i in message_all:
        # i=i.replace("_"," ")
        i=i.capitalize()#首字母大写
        print(i, flush=True)
    print(symbol * 22, flush=True)






