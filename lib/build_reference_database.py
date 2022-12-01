#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2021/10/20 16:48
# @Author  : xiepulin
# @File    : build_reference_database.py
# @Software: PyCharm

import re
import os
from Bio import  SeqIO
from Bio.Seq import  Seq
import  time
from collections import  defaultdict
from  Bio.SeqRecord import SeqRecord
from concurrent.futures import ProcessPoolExecutor
from lib.basic import get_file_list,get_basename
# from basic import get_file_list,get_basename

##########################################################
##########################################################
'''
第三部分
提取出参考基因的fasta序列，用于构建参考基因库。主要分为两类：一类从gb格式中根据基因名提取，一类从fasta格式中提取
'''
###########################################################
###########################################################



class Extract_reference():
    def __init__(self,configuration_information):

        self.configuration_information=configuration_information  #包含各级文件名字信息
        self.out_dir=configuration_information["out"]  #最大一级的输出文件夹
        self.rtfa=configuration_information["rtfa"]                   #参考基因组
        self.rtgb = configuration_information["rtgb"]  # 参考基因组
        # self.soft_boundary=configuration_information["soft_boundary"] #软边界
        self.soft_boundary=0  #2022-10-31 在my_assemble中新增更广义的软边界，即处理组装结果而不是参考序列
        self.gene_max_length=configuration_information["max_length"]  #基因最大长度
        self.gene_min_length=configuration_information["min_length"]   #基因最小长度
        self.reference_database=configuration_information["reference_database"]  #"reference_database"
        self.thread_number=configuration_information["thread_number"]


    def add_soft_boundary(self,start,end,start_all,end_all):
        soft_boundary=self.soft_boundary
        gene_min_length=self.gene_min_length
        gene_max_length=self.gene_max_length
        if (end - start < gene_min_length) and (end - start > gene_max_length):
            return (start,end)
        soft_start = start - soft_boundary
        soft_end = end + soft_boundary
        if soft_start <= start_all:  #  Out of Left Boundary
            start = start
        else:
            start = soft_start
        if soft_end >= end_all:  # Out of Right Boundary
            end = end
        else:
            end = soft_end
        return  (start,end)


    def write_fasta_file(self,record,path):
        SeqIO.write(record,path,"fasta")

    def extract_reference_from_gb_parallel(self):
        out_dir = self.out_dir
        ref = self.rtgb
        soft_boundary = self.soft_boundary
        gene_max_length = self.gene_max_length
        gene_min_length = self.gene_min_length
        reference_database = self.reference_database
        thread_number=self.thread_number
        reference_database_path = os.path.join(out_dir, reference_database)
        if not os.path.isdir(reference_database_path):
            os.makedirs(reference_database_path)
        files = get_file_list(ref)

        task_pool=[]
        results=[]
        executor = ProcessPoolExecutor(max_workers=thread_number)
        for file in files:
            task_pool.append(executor.submit(self.extract_reference_from_gb,file))
        for i in task_pool:
            results.append(i.result())
        if not results:
            return 0

        '''
        Unified format
        '''
        my_records=defaultdict(list)   # [ All_Records1,All_Records2 ] -->  [ {"gene1":[SeqRecord1,Seqrecord2],"gene2":[SeqRecord1,Seqrecord2]}, {"gene2":[SeqRecord1,Seqrecord2],"gene3":[SeqRecord1,Seqrecord2]} ]
        for i in results:
            for key,value in i.items():
                # print(key) #psbA
                # print(value) #[SeqRecord(),SeqRecord(),SeqRecord().....]
                if key not in my_records:
                    my_records[key]=value
                else:
                    my_records[key].extend(value)  # value is list ,so use "extend" not "append"

        '''
        Processing Alias
        '''
        uniform_name_dict = defaultdict(dict)  #defaultdict(<class 'dict'>, {'acc': {'ACC': 2, 'acc': 3, 'aCc': 1} ...)
        Alias2Real_dict = defaultdict()  # defaultdict(None, {'acc': 'acc', 't': 't', 'ss': 'SS'})
        real_name_list = []  # ['acc', 't', 'SS']

        for key , value in my_records.items():
            number=len(value)
            new_name = str(key).lower()
            if new_name not in uniform_name_dict:
                uniform_name_dict[new_name] = {key: number}
            else:
                uniform_name_dict[new_name].update({key: number})
        for key, value in uniform_name_dict.items():
            real_name = max(value, key=lambda gene_name: value[gene_name])
            Alias2Real_dict[key] = real_name
            real_name_list.append(real_name)

        # print(uniform_name_dict)
        # print(Alias2Real_dict)
        '''
        Unified format again
        '''
        my_records_ultimate=defaultdict(list)
        for key,value in my_records.items():
            if key in real_name_list:
                my_records_ultimate[key].extend(value)
            else:
                key_ultimate= Alias2Real_dict[str(key).lower()]
                my_records_ultimate[key_ultimate].extend(value)


        ####write
        task_pool2=[]
        results2=[]
        executor2=ProcessPoolExecutor(max_workers=thread_number)
        for key,value in my_records_ultimate.items():
            path = os.path.join(reference_database_path, key + ".fasta")
            task_pool2.append(executor2.submit(self.write_fasta_file,value,path))
        for i in task_pool2:
            results2.append(i.result())





    def extract_reference_from_gb(self,file):
        gene_min_length=self.gene_min_length
        gene_max_length=self.gene_max_length

        All_Records=defaultdict(list)  #all rec

        for rec in SeqIO.parse(file, "gb"):
            Records=defaultdict(list)# rec      {"ycf1":[{"gene_name":xx,"gene_sequence":xx},{"gene_name":xx,"gene_sequence":xx}]， "matk":[{"gene_name":xx,"gene_sequence":xx}]   }
            repeated_gene = []
            multi_fragment_complex_gene = []  # (rps12)
            appeared_gene=[]

            crossed_origin_gene = ["psbA", "trnH-GUG"]
            sequence = rec.seq
            temp = [i.strand for i in rec.features if i.type == "source"]
            strand_all = temp[0] if temp != [] else 1  # genbank默认为正义链
            start_all = 1
            end_all = len(rec.seq)
            organism = rec.annotations["organism"].replace(" ", "_")  #  species name  eg. Ligusticum_chuanxiong
            id = rec.id  #id = accession + version eg. NC_057131.1
            identifier = organism + "_" + id

            for feature in rec.features:
                # OrderedDict([('gene', ['matK']), ('locus_tag', ['KQ413_pgp084']), ('db_xref', ['GeneID:65316243'])])
                # OrderedDict([('locus_tag', ['A4330_gr002']), ('db_xref', ['GeneID:27214299'])])鲁棒性
                #feature,type == gene
                if feature.type == "gene" and "gene" in feature.qualifiers.keys() and feature.qualifiers["gene"][0] in crossed_origin_gene:
                    gene_information={}
                    seq=feature.location.extract(sequence)
                    gene_name = feature.qualifiers["gene"][0].replace(" ", "_")
                    gene_information["gene_name"] = feature.qualifiers["gene"][0]
                    gene_information["gene_sequence"] = seq
                    gene_information["identifier"] = identifier
                    gene_information["organism"] = organism
                    gene_information["id"] = id
                    gene_information["length"]=len(str(seq))
                    gene_information["start"]=feature.location.start    #ExactPosition(84940)
                    gene_information["end"] = feature.location.end


                    # crossed_origin_gene  must be recorded
                    if gene_name not in appeared_gene:
                        appeared_gene.append(gene_name)
                    else:
                        pass

                    if gene_name not in Records:
                        Records[gene_name]=[gene_information]      #[gene_information] --> [{}]
                    else:
                        Records[gene_name].append(gene_information)  #[gene_information1,gene_information2 ] --> [{},{}]


                elif feature.type == "gene" and "gene" in feature.qualifiers.keys():
                    gene_information = {}
                    gene_name = feature.qualifiers["gene"][0].replace(" ", "_")
                    location=feature.location
                    if "join" in str(location):
                        multi_fragment_complex_gene.append(gene_name)
                        continue

                    strand=feature.strand
                    start=feature.location.start
                    end=feature.location.end
                    start,end=self.add_soft_boundary(int(start),int(end),start_all,end_all)  #ExactPosition(84940) --> int
                    if strand==strand_all:
                        seq=sequence[start:end]
                        gene_information["gene_sequence"] = seq
                    else:
                        seq=sequence[start:end].reverse_complement()
                        gene_information["gene_sequence"] = seq
                    gene_information["gene_name"] = feature.qualifiers["gene"][0]
                    gene_information["identifier"] = identifier
                    gene_information["organism"] = organism
                    gene_information["id"] = id
                    gene_information["length"] = len(str(seq))
                    gene_information["start"] = feature.location.start
                    gene_information["end"] = feature.location.end


                    appeared_gene.append(gene_name)
                    if gene_name not in Records:
                        Records[gene_name] = [gene_information]
                    else:
                        Records[gene_name].append(gene_information)
                else:
                    pass

            if not Records:
                continue

            for key,value in Records.items():
                number=0
                #Inverted Repeat Zone
                if len(Records[key]) >=2:
                    for i in range(len(value)):
                        number+=1
                        Records[key][i]["identifier"] = Records[key][i]["identifier"]+"_Repeat_"+str(number)

            # {"ycf1":[{"gene_name":xx,"gene_sequence":xx},{"gene_name":xx,"gene_sequence":xx}]， "matk":[{"gene_name":xx,"gene_sequence":xx}]   }


            for key, value in Records.items():
                for i in range(len(value)):
                    start=int(value[i]["start"])
                    end=int(value[i]["end"])
                    if (end - start  >=  gene_min_length) and (end - start  <= gene_max_length):  #Length Limit
                        temp=SeqRecord(id=value[i]["identifier"],seq=value[i]["gene_sequence"], description="")
                        if key not in All_Records:
                            All_Records[key]=[temp]
                        else:
                            All_Records[key].append(temp)  #temp is "SeqRecord" ,not list. so use "append" not "extend"
        # print(All_Records)   # { "gene1": [SeqRecord1,SeqRecord2] ,"gene2": [SeqRecord1,SeqRecord2]  }

        # print(len(All_Records["atpB"]))

        return All_Records

    def extract_reference_from_fasta(self):
        ref = self.rtfa
        out_dir = self.out_dir
        reference_database = self.reference_database
        reference_databese_path = os.path.join(out_dir, reference_database)
        thread_number=self.thread_number

        if not os.path.isdir(reference_databese_path):
            os.makedirs(reference_databese_path)

        files=get_file_list(ref)
        task_pool = []
        results = []

        executor = ProcessPoolExecutor(max_workers=thread_number)
        for file in files:
            file_name = get_basename(file) + ".fasta"
            path = os.path.join(reference_databese_path, file_name)
            task_pool.append(executor.submit(self.get_pure_fasta_format_sequence,file,path))
        for i in task_pool:
            results.append(i.result())


    # 剔除ACGTU之外的序列，如？-等等;剔除空行 使用Biopython
    def get_pure_fasta_format_sequence(self, file, output):
        my_records = []
        base=["A","T","C","G","U"]
        for rec in SeqIO.parse(file,"fasta"):
            name=rec.name
            description=rec.description
            seq=str(rec.seq).upper()
            seq=[i for i in seq if i in base]
            seq="".join(seq)
            if not seq:
                continue
            # record=[">"+name+" "+description+"\n",seq+"\n"]
            # my_records.extend(record)
            record=SeqRecord(id=name,seq=Seq(seq),description=description)
            my_records.append(record)

        if my_records:
            SeqIO.write(my_records,output,"fasta")
            # with open(output,"w") as f:
            #     f.writelines(my_records)


def my_bulid_reference_database_pipeline(configuration_information):
    rtfa=configuration_information["rtfa"]
    rtgb=configuration_information["rtgb"]

    if rtfa:
        my_get_reference_from_fasta=Extract_reference(configuration_information)
        my_get_reference_from_fasta.extract_reference_from_fasta()
    elif rtgb:
        my_get_reference_from_gb=Extract_reference(configuration_information)
        my_get_reference_from_gb.extract_reference_from_gb_parallel()
    else:
        pass


# if __name__ == '__main__':
#     t1 = time.time()
#     out = r"E:\Computer\python\GeneMiner\eeeeeeeee10 重构bootstrap\example\shallow"
#     rtgb=r"E:\Computer\python\GeneMiner\eeeeeeeee10 重构bootstrap\example\shallow_ref.gb"
#     # rtgb=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeeee10 重构bootstrap\example\demo.gb"
#     # rtgb=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeeee10 重构bootstrap\example\demo2.gb"
#     # rtgb=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeeee10 重构bootstrap\example\demo2"
#
#
#     rtfa=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeeee10 重构bootstrap\example\ITS_ref.fasta"
#     # rtfa = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeeee10 重构bootstrap\example\demo1"
#     rtfa=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeeee10 重构bootstrap\example\ref_Brassicaceae"
#
#
#     soft_boundary = 20
#     max_length = 5000
#     min_length = 0
#     reference_database = "reference_database"
#     thread_number=4
#
#     configuration_information = {"out": out, "rtfa": rtfa, "rtgb": rtgb, "soft_boundary": soft_boundary,
#                                  "max_length": max_length, "min_length": min_length,
#                                  "reference_database": reference_database,
#                                  "thread_number":thread_number}
#
#     # my_bulid_reference_database_pipeline(configuration_information)
#
#     target1 = Extract_reference(configuration_information)
#     target1.extract_reference_from_gb_parallel()
#
#     # target1 = Extract_reference(configuration_information)
#     # target1.extract_reference_from_fasta()
#
#
#
#     t2 = time.time()
#     print(t2 - t1)
































