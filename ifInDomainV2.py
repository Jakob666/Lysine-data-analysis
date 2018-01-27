#-*- coding:utf-8 -*-
'''
Description:
    只考虑所有显著蛋白在 motif区域内的情况，且将每个蛋白的7种修饰数据合并，
    定位到 domain和非domain区域，并使用 bootstrap检验是否显著差异。
==================================
@author:hbs
@date:2018-1-26
@version:1.0
'''
import pandas as pd
from bootstrapTest import Bootstrap_test
from ifInDomain import IfInDomain
import os
import multiprocessing as mp
import warnings


class ifInDomainV2(IfInDomain):

    def __init__(self, cancer, sequence_data, domain_data):
        self.cancer = cancer
        self.sequences = sequence_data
        self.domain = domain_data



    def load_motif_data(self, mod_type):
        # motif数据的导入，从依据modification分类后的文件中导入 protein_id(uniprot isoform),Canonical,cancer_type,
        # mutated position,from,to,sample_num（目前还没计数）,direct mutation
        in_motif_file = "/data1/hbs/all_cancer_analysis/all_cancer_sig_prot_motif_modify/"+self.cancer+"/"+mod_type+"_in_motif.txt"
        try:
            in_motif_data = pd.read_csv(in_motif_file, sep="\t",
                                    usecols=["protein_id", "Canonical", "cancer_type", "mutated position",
                                             "sample_num"], encoding="utf-8")
        except pd.io.common.EmptyDataError:
            in_motif_data = pd.DataFrame(columns=["protein_id", "Canonical", "cancer_type", "mutated position",
                                             "sample_num"])
        self.in_motif = pd.concat([self.in_motif, in_motif_data], axis=0)
        return None

    def merge_motif_data(self):
        self.in_motif = pd.DataFrame(columns=["protein_id", "Canonical", "cancer_type",
                                              "mutated position", "sample_num"])

        mod_list = ["ace", "gly", "mal", "met", "suc", "sum", "ubi"]
        map(self.load_motif_data, mod_list)
        self.in_motif.drop_duplicates(inplace=True)
        self.in_motif.set_index("Canonical", drop=False, inplace=True)
        return None

    def merge_data(self):
        data = pd.merge(left=self.in_motif, right=self.domain, left_index=True, right_index=True)
        data.drop_duplicates(inplace=True)

        # interproscan预测某个蛋白序列很可能没有domain，没有domain的蛋白是不会出现在interproscan文件中
        # 此时就要找出这些蛋白，定义为 protein_without_domain，并最终归总在 self.out_domain中
        no_domain_protein = self.in_motif.drop(set(data.index), axis=0)
        self.protein_without_domain = pd.DataFrame(
            columns=["idx", "Uniprot Accession", "protein length", "domain start", "domain end", "non_domain length",
                     "position", "sample_num", "sample_count", "in_domain"])
        no_domain_protein.apply(self.no_domain_protein, axis=1)

        self.protein_without_domain.set_index("idx", drop=False, inplace=True)
        self.protein_without_domain = pd.merge(left=self.protein_without_domain, right=self.sequences, left_index=True,
                                               right_index=True)
        # 获取到没有domain的蛋白的序列长度
        self.protein_without_domain = self.protein_without_domain.apply(self.get_seq_length, axis=1)
        self.protein_without_domain.drop_duplicates(inplace=True)
        return data


if __name__ == "__main__":
    warnings.filterwarnings("ignore")

    domain = pd.DataFrame(columns=["protein_id", "Uniprot Accession", "protein length", "start", "end"])

    def load_domain_data(mod_type):
        global domain
        mod_dic = {"ace": "Acetylation", "gly": "Glycation", "mal": "Malonylation", "met": "Methylation",
                   "suc": "Succinylation", "sum": "Sumoylation", "ubi": "Ubiquitination"}
        domain_file = "/data1/hbs/all_cancer_analysis/all_domain_classify_by_modify/"+ mod_dic[mod_type] +"_simplified.txt"
        domain_data = pd.read_csv(domain_file, sep="\t", encoding="utf-8")

        # 计算domain长度，此时dataframe中会有两个长度，一个是整个蛋白序列的长度，另一个是单独的domain区域的长度
        domain_data["domain length"] = domain_data["end"] - domain_data["start"] + 1
        domain = pd.concat([domain, domain_data], axis=0)
        return None

    def merge_domain_data():
        global domain
        i = IfInDomain("")
        mod_list = ["ace", "gly", "mal", "met", "suc", "sum", "ubi"]

        map(load_domain_data, mod_list)
        domain.drop_duplicates(inplace=True)
        domain.set_index("Uniprot Accession", drop=False, inplace=True)
        # 将存在overlap的domain合并
        domain = i.merge_overlap_domain(domain)
        return None

    def execute_class_ifindomainV2(cancer,domain_data,sequence_data,in_domain_file,out_domain_file,in_domain_calc,out_domain_calc):
        '''
        :param files: 是一个列表，里面存放相关的蛋白质信息文件和输出的文件地址
        '''
        i = ifInDomainV2(cancer, sequence_data, domain_data)
        i.merge_motif_data()
        data = i.merge_data()
        i.form_new_data(data)
        i.to_doc(in_domain_file, in_domain_calc,out_domain_file,out_domain_calc)
        return None

    pool = mp.Pool(7)

    cancer_kind = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP',
                   'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM',
                   'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']
    # cancer_kind = ["PCPG"]

    sequence_file = "/data1/hbs/total_fastr/Total.elm"
    sequence_data = pd.read_csv(sequence_file, sep="\t", usecols=[0, 1, 4], encoding="utf-8")
    sequence_data.set_index("Uniprot Accession", inplace=True)

    merge_domain_data()
    # print(domain)
    # exit()

    result_dir = "/data1/hbs/all_cancer_analysis/all_cancer_modification_in_domain_or_not/seven_modify_test_result/"

    for c in cancer_kind:
        cancer_dir = result_dir + c
        # os.system("mkdir " + cancer_dir)
        in_domain_file = cancer_dir + "/" + "in_domain.txt"
        out_domain_file = cancer_dir + "/" + "out_domain.txt"
        in_domain_calc = cancer_dir + "/" + "in_domain_calc.txt"
        out_domain_calc = cancer_dir + "/" + "out_domain_calc.txt"
        # execute_class_ifindomainV2(c,domain,sequence_data,in_domain_file,out_domain_file,in_domain_calc,out_domain_calc)

        pool.apply_async(execute_class_ifindomainV2, (c,domain,sequence_data,in_domain_file,out_domain_file,in_domain_calc,out_domain_calc))
    
    pool.close()
    pool.join()