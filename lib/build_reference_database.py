#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2021/10/20 16:48
# @Author  : xiepulin
# @File    : build_reference_database.py
# @Software: PyCharm

from basic import *
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
        self.out_dir=configuration_information["out_dir"]  #最大一级的输出文件夹
        self.rtfa=configuration_information["rtfa"]                   #参考基因组
        self.rtgb = configuration_information["rtgb"]  # 参考基因组
        self.soft_boundary=configuration_information["soft_boundary"] #软边界
        self.gene_max_length=configuration_information["max_length"]  #基因最大长度
        self.gene_min_length=configuration_information["min_length"]   #基因最小长度
        self.reference_database=configuration_information["reference_database"]  #"reference_database"

    '''
    extract_location_from_gb() identify_location() collect_fasta_information_from_gb_v2() 三个函数服务于extract_reference_from_gb（）
    '''

    '''
    提取出具体的基因位置
    非复合基因  [297:2890](-)  复合基因 join{[69127:69541](-), [97170:97602](-), [96808:98834](-)}
    将biopython格式的基因表示格式转为一个列表 如: ['35348', '35576', '33277', '34495', '51520', '51542', '221825', '222220', '220475', '220625']
    '''

    def extract_location_from_gb(self, feature, gene_location):
        if "join" not in gene_location:
            # 软边界
            start = int(feature.location.start)
            end = int(feature.location.end)

            location = [start, end]
            return location
        else:
            location = re.findall(r"\d+:\d+",
                                  gene_location)  # ['69427:69541', '97370:97602', '96808:96834']
            location = ",".join(location).replace(":", ",")  # 69427,69541,97370,97602,96808,96834
            location = location.split(",")  # ['69427', '69541', '97370', '97602', '96808', '96834']
            return location
    '''
    增添软边界 增添长度限制
    example: boundary=0 则['35348', '35576', '33277', '34495']----------[35348, 35576, 33277, 34495]  
    '''
    def identify_location(self, gene_location, soft_boundary, start_all, end_all, gene_max_length,
                          gene_min_length):
        # 过滤起始位点 和终止位点。
        # (1)加上和减去软边界，不得超过0 或者148512 （起点和终点） (2)长度在要求范围中间
        new_location = []
        for i in range(0, len(gene_location), 2):
            start = int(gene_location[i])
            end = int(gene_location[i + 1])
            soft_start = start - soft_boundary
            soft_end = end + soft_boundary
            if soft_start <= start_all:  # 不能超出最左的边界 0
                start = start
            else:
                start = soft_start
            if soft_end >= end_all:  # 不能超出右边界 全长148512
                end = end
            else:
                end = soft_end
            if (end - start >= gene_min_length) and (end - start <= gene_max_length):
                new_location.append(start)
                new_location.append(end)
        return new_location

    '''
    提取出基因的所有信息  
    [{'gene_name': 'orf56b', 'gene_sequence': Seq('ATGCCTCTAGATGATTTCGACGTTATATTGGGCATTGATTTCCTATTGGCAGCT...TAA'), 'identifier': 'Daucus_carota_subsp._sativus_NC_017856', 'organism': 'Daucus_carota_subsp._sativus', 'accession': 'NC_017856'}]
    对于有多段的复合基因 gene_1 gene_2 gene_3 作为区分标识
    '''
    def collect_reference_information_from_gb_v2(self,gene_feature_location_key,gene_feature_location_value,soft_boundary,start_all,end_all, gene_max_length,gene_min_length,rec, strand_all, identifier, organism,accession):
        '''
        :param gene_feature_location_key:  基因名  gene_feature_location的键
        :param gene_feature_location_value:  feature与location存储信息  gene_feature_location的值
        :param start_all:
        :param end_all:
        :param gene_max_length:
        :param gene_min_length:
        :param rec:
        :param strand_all:
        :param identifier:
        :param organism:
        :param accession:
        :return:
        '''
        #数据结构  {gene":[ [feature1,location1],[feature2,location2] ]  }
        gene_information_all = []
        '''
        将位置信息重新存为列表
        '''
        new_gene_feature_location=defaultdict(list)
        for i in gene_feature_location_value:
            gene_location = str(i[1])  # 原始的location有特定的数据结构，需要强制转换为str   #join{[69427:69541](-), [96808:97602](-)}
            feature = i[0]
            gene_location = self.extract_location_from_gb(feature, gene_location)  #['69427', '69541', '96808', '97602'] 单纯的将biopython数据结构类型的位置信息转为列表
            gene_location = self.identify_location(gene_location, soft_boundary, start_all, end_all, gene_max_length,gene_min_length) #[69427, 69541, 96808, 97602]  #添加软边界，添加长度限制

            if gene_location==[]:
                continue
            if new_gene_feature_location[gene_feature_location_key] ==[]:
                new_gene_feature_location[gene_feature_location_key] = [[feature, gene_location]]
            else:
                new_gene_feature_location[gene_feature_location_key].append([feature, gene_location])

        if new_gene_feature_location==[]:  #如果都不符合要求，返回
            return gene_information_all

        '''
        记录基因信息，包括location,seq,description,index等等 特别是对多基因多片段的处理
        '''
        gene_number=len(new_gene_feature_location[gene_feature_location_key])  #基因出现的次数 一般100多个基因有10多个会出现两次
        sequence = rec.seq
        index=0  #区分基因的多个片段，甚至多个基因的多个片段
        if gene_number == 0:
            pass
        elif gene_number==1:  #[[feature1,location1]]
            gene_information = {} #临时
            feature = new_gene_feature_location[gene_feature_location_key][0][0]       #[[feature1,location1]] --feature
            gene_location = new_gene_feature_location[gene_feature_location_key][0][1] #[[feature1,location1]] --location
            strand = feature.strand
            gene_information["gene_name"] = feature.qualifiers["gene"][0]
            gene_information["gene_sequence"] = sequence[gene_location[0]:gene_location[1]] if strand_all == strand else sequence[gene_location[0]:gene_location[-1]].reverse_complement()
            gene_information["identifier"] = identifier
            gene_information["organism"] = organism
            gene_information["accession"] = accession
            gene_information_all = [gene_information]

        else:  #[[feature1,location1], [feature2,location2]]
            gene_feature_location=new_gene_feature_location[gene_feature_location_key]
            for i in range(gene_number):
                feature = gene_feature_location[i][0]
                gene_location = gene_feature_location[i][1]
                strand = feature.strand
                for i in range(0, len(gene_location), 2):
                    index += 1
                    gene_information = {}
                    gene_information["gene_name"] = feature.qualifiers["gene"][0] + "_" + str(index)  # 用于区分
                    gene_information["gene_sequence"] = sequence[gene_location[i]:gene_location[
                        i + 1]] if strand_all == strand else sequence[
                                                             gene_location[i]:gene_location[i + 1]].reverse_complement()
                    gene_information["identifier"] = identifier
                    gene_information["organism"] = organism
                    gene_information["accession"] = accession
                    gene_information_all.append(gene_information)


        # print(gene_information_all)
        return gene_information_all




    def extract_reference_from_gb(self):
        out_dir = self.out_dir
        ref = self.rtgb
        soft_boundary = self.soft_boundary
        gene_max_length = self.gene_max_length
        gene_min_length = self.gene_min_length
        reference_database = self.reference_database

        reference_database_path = os.path.join(out_dir, reference_database)
        dir_make(reference_database_path)
        files = get_file_list(ref)

        Genes_information = []
        # omit_genes = ["rps12"]  # 分段太多，省略
        omit_genes = ["rps12"]  # 分段太多，省略
        special_genes = ["psbA", "trnH-GUG"]  # 两个可能跨越原点的基因：trnH-GUG 和 psbA
        for file in files:
            for rec in SeqIO.parse(file, "gb"):
                sequence = rec.seq
                temp = [i.strand for i in rec.features if i.type == "source"]
                strand_all = temp[0] if temp != [] else 1  # genbank默认为正义链

                start_all = 1
                end_all = len(rec.seq)
                organism = rec.annotations["organism"].replace(" ", "_")  # 物种名  Ligusticum_chuanxiong
                accession = rec.name.replace(" ", "_")  # 样本名/genbank id
                identifier = organism + "_" + accession
                genes_feature_location = defaultdict(list)
                gene_name_list = []

                for feature in rec.features:
                    # OrderedDict([('gene', ['matK']), ('locus_tag', ['KQ413_pgp084']), ('db_xref', ['GeneID:65316243'])])
                    # OrderedDict([('locus_tag', ['A4330_gr002']), ('db_xref', ['GeneID:27214299'])])鲁棒性
                    #feature,type == gene 可能是不存在的
                    if feature.type == "gene" and "gene" in feature.qualifiers.keys() and feature.qualifiers["gene"][0] in special_genes:
                        gene_information = {}
                        gene_name = feature.qualifiers["gene"][0]
                        if gene_name not in gene_name_list:
                            gene_name_list.append(gene_name)

                        gene_information["gene_name"] = feature.qualifiers["gene"][0]
                        gene_information["gene_sequence"] = feature.location.extract(sequence)
                        gene_information["identifier"] = identifier
                        gene_information["organism"] = organism
                        gene_information["accession"] = accession
                        Genes_information.append(gene_information)
                    elif feature.type == "gene" and "gene" in feature.qualifiers.keys() and feature.qualifiers["gene"][0] in omit_genes:
                        continue

                    elif  feature.type == "gene" and "gene" in feature.qualifiers.keys():
                        gene_name = feature.qualifiers["gene"][0]
                        if gene_name not in gene_name_list:
                            location = feature.location
                            genes_feature_location[gene_name] = [[feature, location]]
                            gene_name_list.append(gene_name)
                        else:
                            location = feature.location
                            genes_feature_location[gene_name].append([feature, location])
                    else:
                        pass

                for i in genes_feature_location:   #genes_feature_location[i] 传得值
                    gene_information_all=self.collect_reference_information_from_gb_v2(i,genes_feature_location[i],soft_boundary,start_all,end_all, gene_max_length,gene_min_length,rec, strand_all, identifier, organism,accession)
                    if gene_information_all != []:
                        for i in gene_information_all:
                            Genes_information.append(i)


        # print(Genes_information)     #[{'gene_name': 'matK', 'gene_sequence': Seq('ATGGAGGAATTCCAAAGATATTTAAAGCTAAATATATCTCAACAACACTACTTC...TAA'), 'identifier': 'Ligusticum_chuanxiong_NC_038088', 'organism': 'Ligusticum_chuanxiong', 'accession': 'NC_038088'}]
        '''
        写fasta格式的文件
        '''
        # gene_name_list 只记录基因名，没有考虑join情况下 (gene_1 gene_2的情况)
        new_gene_name_list = []
        for i in Genes_information:
            gene_name = i["gene_name"]
            new_gene_name_list.append(gene_name)
        # print(new_gene_name_list)
        for i in new_gene_name_list:
            my_records = []
            for j in Genes_information:
                if j["gene_name"] == i and len(j["gene_sequence"]) >= gene_min_length and len(
                        j["gene_sequence"]) <= gene_max_length:
                    my_record = SeqRecord(seq=j["gene_sequence"], id=j["accession"],
                                          # 文件名为基因名  >后面的名字为GenBank号码  description为物种名
                                          description=j["organism"])
                    my_records.append(my_record)
                    Genes_information.remove(j)
            # print(my_records)
            # print(len(my_records))
            if my_records != []:
                SeqIO.write(my_records, os.path.join(reference_database_path, i + ".fasta"), "fasta")


        #剔除杂质序列
        files=get_file_list(reference_database_path)
        for file in files:
            file_name = get_basename(file) + ".fasta"
            path = os.path.join(reference_database_path, file_name)
            self.get_pure_fasta_format_sequence(file, path)





    def extract_reference_from_fasta(self):
        ref = self.rtfa
        out_dir = self.out_dir
        reference_database = self.reference_database
        reference_databese_path = os.path.join(out_dir, reference_database)

        dir_make(reference_databese_path)
        files=get_file_list(ref)
        for file in files:
            file_name=get_basename(file)+".fasta"
            path = os.path.join(reference_databese_path, file_name)
            self.get_pure_fasta_format_sequence(file, path)


    # 剔除ACGTU之外的序列，如？-等等
    def get_pure_fasta_format_sequence(self, file, output):
        my_records = []
        base=["A","T","C","G"]
        infile = open(file, 'r', encoding='utf-8', errors='ignore')
        for line in infile:
            if line[0]==">":
                my_records.append(line)
            else:
                # line=line.upper()
                # line=[i for i in line if i in base]
                # line="".join(line)+"\n"
                # my_records.append(line)
                #新增剔除空行
                line=line.strip()
                if line:
                    line = line.upper()
                    line = [i for i in line if i in base]
                    line = "".join(line) + "\n"
                    my_records.append(line)

        with open(output,"w") as f:
            f.writelines(my_records)


def my_bulid_reference_database_pipeline(configuration_information):
    rtfa=configuration_information["rtfa"]
    rtgb=configuration_information["rtgb"]

    if rtfa:
        my_get_reference_from_fasta=Extract_reference(configuration_information)
        my_get_reference_from_fasta.extract_reference_from_fasta()
    elif rtgb:
        my_get_reference_from_gb=Extract_reference(configuration_information)
        my_get_reference_from_gb.extract_reference_from_gb()
    else:
        pass


if __name__ == '__main__':
    t1 = time.time()
    out_dir = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeee9 重构filter\example\AT_cp"
    # rtfa = r"example\gene_353"
    # rtfa=r"example\demo.fa"
    rtgb=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeee9 重构filter\example\AT_ref"
    # rtgb = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeee9 重构filter\example\ref_gb"

    #剔除空行
    rtfa=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeee9 重构filter\example\ETS.fasta"

    soft_boundary = 0
    max_length = 5000
    min_length = 300
    reference_database = "reference_database"

    configuration_information = {"out_dir": out_dir, "rtfa": rtfa, "rtgb": rtgb, "soft_boundary": soft_boundary,
                                 "max_length": max_length, "min_length": min_length,
                                 "reference_database": reference_database}

    target1 = Extract_reference(configuration_information)
    # target1.extract_reference_from_fasta()
    target1.extract_reference_from_gb()

    t2 = time.time()
    print(t2 - t1)
































