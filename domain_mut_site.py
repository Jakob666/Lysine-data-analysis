#-*- coding:utf-8 -*-
'''
Description:
    将蛋白质突变位点匹配到相应的基因组坐标上（annovar获得），
    文件存放位置 /data1/data/zengyanru/LysineTCGA/find_mutsite_wholegenome
    将⑧中得到的位点与之前的domain文件相结合，得到 每种癌症、每种修饰的显著性蛋白
    在domain区域的突变位点是哪些，为的是进一步使用phastCons100软件进行保守性分析
===================================================================
@author: hbs
@date: 2018-1-18（由于原始文件有问题，暂时停滞）
@version: 1.0
===================================================================
@update: 由于原始文件有所变动，之前的代码全部推掉
@version:2.0
'''
import pandas as pd
import warnings
import multiprocessing as mp
import os


class MutSiteOnGenome:

    def __init__(self):
        pass

    @staticmethod
    def load_genome_data(genome_file, mod_type):
        '''
        将annovar生成注释的genome文件读入内存中，文件中比较的列都比较重要，不做选读。
        文件中大多数的行是这种形式存在的
        ACC	O95104	in	ENST00000286835	O95104-1	no_ptm_in	no_mutsite	no_mutsite	no_mutsite	nonsense	nonsense
        ACC	P40937	in	ENST00000454402	P40937-1	Sumoylation	7	*METSALKQQEQPA	E10K	12	118016855
        ACC	Q92556	in	ENST00000310758	Q92556-1	no_ptm_in	no_mutsite	no_mutsite	no_mutsite	nonsense	nonsense
        需要去除nonsense的行，这些属于无用信息
        :param genome_file: 基因组文件的路径
        :param mod_type: 翻译后修饰的类型，如果是all，则是对7种修饰不进行区分
        :return: genome data in dataframe format
        '''
        #将数据文件读入内存，形成dataframe
        genome_data = pd.read_csv(genome_file, sep="\t", encoding="utf-8")
        #过滤无效信息
        genome_data = genome_data[genome_data["Start"] != "nonsense"]
        if mod_type != "all":
            genome_data = genome_data[genome_data["ptm_type"] == mod_type]
        del genome_data["ptm_type"]
        genome_data.drop_duplicates(inplace=True)
        genome_data.fillna("emp", inplace=True)
        return genome_data

    @staticmethod
    def load_domain_data(domain_file):
        '''
        将存放in domain突变位点的数据读入dataframe，具体选取的列见下方注释
        :param domain_file: domain文件的路径
        :return: domain文件相关的dataframe
        '''
        #读取之前进行分类的 in domain的突变的数据，读取的三列分别为 Uniprot ID（精确到亚形）、position和sample num
        domain_data = pd.read_csv(domain_file, sep="\t", usecols=[0, 5, 6], encoding="utf-8")
        domain_data.drop_duplicates(subset=["position", "sample_num"], inplace=True)
        return domain_data

    def match_mutation_on_genome(self, domain_data, genome_data):
        '''
        将得到的 in domain突变的数据和 genome突变的数据整合在一起，确定蛋白上
        突变位点在基因组上的位置。
        :param domain_data: in domain突变位点的数据，dataframe形式
        :param genome_data: 基因组位点的数据，dataframe形式
        :return:
        '''
        #创建索引的过程
        domain_data["position"] = pd.Series(domain_data["position"].map(str))
        domain_data["idx"] = domain_data["Uniprot Accession"] + ":" + domain_data["position"]
        genome_data["idx"] = ""
        genome_data = genome_data.apply(self.__create_genome_index, axis=1)

        domain_data.set_index("idx", inplace=True)
        genome_data.set_index("idx", inplace=True)

        data = pd.merge(left=domain_data, right=genome_data, left_index=True, right_index=True)
        return data

    def __create_genome_index(self, i):
        isoform = i["Isoform"]
        position = i["mutation_site"][1:-1]
        if isoform != "emp":
            idx = isoform + ":" + position
        idx = i["protein"] + ":" + position
        i["idx"] = idx
        return i

    def to_doc(self, data, output_file):
        '''
        将匹配到基因组位点的in domain 突变数据写入文件
        :param data:需要写入数组的的数据
        :param output_file:
        :return: None
        '''
        data.to_csv(output_file, sep="\t", mode="w", encoding="utf-8")
        return None

    def main(self, cancer, genome_file, output_file, mode_type):
        '''将之前定义的方法组合形成一整套的流程'''
        if mode_type == "all":
            mod_list = ["ace", "gly", "mal", "met", "suc", "sum", "ubi"]
        else:
            mod_list = [mode_type]
        domain_dir = "/data1/hbs/all_cancer_analysis/all_cancer_mutation_site_conservation/all_cancer_in_domain_mut_site/"
        domain_data = pd.DataFrame(columns=["Uniprot Accession", "position", "sample_num"])
        for m in mod_list:
            domain_file = domain_dir + cancer + m + "_in_domain.txt"
            data = MutSiteOnGenome.load_domain_data(domain_file)
            domain_data = pd.concat([domain_data, data], axis=0)
        domain_data.drop_duplicates(inplace=True)

        genome_data = MutSiteOnGenome.load_genome_data(genome_file, mode_type)

        matched_data = self.match_mutation_on_genome(domain_data, genome_data)
        self.to_doc(matched_data, output_file)


if __name__ == "__main__":
    warnings.filterwarnings("ignore")

    def get_genome_location(cancer, genome_file, output_file, mode_type="all"):
        m = MutSiteOnGenome()
        m.main(cancer, genome_file, output_file, mode_type)
        return None

    pool = mp.Pool(8)
    cancer_kind = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP',
                   'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM',
                   'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']
    result_dir = "/data1/hbs/all_cancer_analysis/all_cancer_mutation_site_conservation/all_cancer_mut_site_genome_location/"

    for c in cancer_kind:
        cancer_dir = result_dir + c
        os.system("mkdir " + cancer_dir)
        genome_file = "/data1/data/zengyanru/LysineTCGA/find_mutsite_wholegenome/%s_protein_mutsite_wholegenome.txt" % c
        output_file = cancer_dir + "/" + c + "_in_domain_mut_genome_location.txt"
        pool.apply_async(get_genome_location, (c, genome_file, output_file))

    pool.close()
    pool.join()