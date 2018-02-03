#-*- coding:utf-8 -*-
'''
Description:
    从曾严儒生成的文件中匹配所有 in domain 突变位点的保守性打分，此次打分使用的是 phyloP100way和 phastCons100，
    目前已有的文件内容为，
Cancer	protein	mut_in_out	ENST	Isoform	        ptm_type	ptm_site	ptm_motif	mutation_site	Chr	Start	phyloP100way_vertebrate	phyloP100way_vertebrate_rankscore	phastCons100way_vertebrate	phastCons100way_vertebrate_rankscore
KICH	P10809	    in	ENST00000345042	P10809-1	Succinylation	269	EIANAHRKPLVII	H267L	        2	197493393	7.962	                0.874	                            1.0	                        0.715
KICH	P00491	    in	ENST00000361505	nan	        Malonylation	265	EEVLAAGKQAAQKLE	G264S	        14	20476521	1.743	                0.377	                            0.996	                    0.391
KICH	P20936	    in	ENST00000274376	P20936-1	Malonylation	834	SCELSPSKLEKNEDV	E839K	        5	87379762	9.841	                0.983	                            1.0	                        0.715
KICH	P10809	    in	ENST00000345042	P10809-1	Malonylation	269	EIANAHRKPLVIIAE	H267L	        2	197493393	7.962	                0.874	                            1.0	                        0.715
    其中能用到的是protein、Isoform、mutation_site、phyloP100way_vertebrate和phastCons100way_vertebrate，
    用之前生成的存放于/data1/hbs/all_cancer_analysis/all_cancer_mutation_site_conservation/all_cancer_in_domain_mut_site中
所记录的数据所记录的点与曾严儒生成的文件中的点进行匹配，得到如下四列信息，
    protein（精确到亚形）  mutation_site   phyloP100way_vertebrate     phastCons100way_vertebrate
生成文件后用这里的数据进行最终的保守性打分的 cdf图绘制。
=============================================
@author:hbs
@date:2018-1-30
@version:1.0
'''
import pandas as pd
import os
import multiprocessing as mp
import warnings


class GetConserveScore:
    def __init__(self):
        self.mod_list = {"ace":"Acetylation", "gly":"Glycation", "mal":"Malonylation", "met":"Methylation",
                         "suc":"Succinylation", "sum":"Sumoylation", "ubi":"Ubiquitination"}
        self.in_domain_mut_site = pd.DataFrame(columns=["Uniprot Accession", "position"])

    def load_mut_conserve_score_data(self, mut_conserve_score_file):
        try:
            mut_score_data = pd.read_csv(mut_conserve_score_file, sep="\t", usecols=[1, 4, 8, 11, 13], encoding="utf-8")

        except:
            mut_score_data = pd.DataFrame(["protein", "Isoform", "phyloP100way_vertebrate", "mutation_site", "phastCons100way_vertebrate"])
        mut_score_data.fillna("emp", inplace=True)
        mut_score_data = mut_score_data.apply(self.preprocess,axis=1)
        self.mut_score_data = mut_score_data.drop_duplicates(subset=["Isoform", "mutation_site"])
        return None

    def load_in_domain_mut_file(self, in_domain_mut_file):
        try:
            in_domain_mut_data = pd.read_csv(in_domain_mut_file, sep="\t", usecols=[0, 5], encoding="utf-8")

        except:
            in_domain_mut_data = pd.DataFrame(columns=["Uniprot Accession", "position"])
        self.in_domain_mut_site = pd.concat([self.in_domain_mut_site, in_domain_mut_data], axis=0)
        self.in_domain_mut_site.drop_duplicates(inplace=True)
        return None

    def preprocess(self, i):
        '''score文件中的位点的形式是 XdddY的形式，其中X、Y是氨基酸，ddd代表的是位点数值型字符串，需要将这个位点提取出来'''
        if i["Isoform"] == "emp":
            i["Isoform"] = i["protein"]
        pos = i["mutation_site"][1:-1]
        i["mutation_site"] = float(pos)
        return i

    def match_with_conserve_score(self):
        self.mut_conserve_data = pd.merge(left=self.in_domain_mut_site, right=self.mut_score_data,left_on=["Uniprot Accession", "position"],
                        right_on=["Isoform", "mutation_site"])
        return None

    def to_doc(self, output_file):
        self.mut_conserve_data.to_csv(output_file, sep="\t", index=False, columns=["Uniprot Accession", "mutation_site",
                                    "phyloP100way_vertebrate","phastCons100way_vertebrate"], encoding="utf-8")
        return None


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    def get_phastCons_score(cancer_type, mut_conserve_score_file, in_domain_dir, output_file):
        c = GetConserveScore()
        c.load_mut_conserve_score_data(mut_conserve_score_file)
        mod = ["ace", "gly", "mal", "met", "suc", "sum", "ubi"]
        for m in mod:
            in_domain_mut_file = in_domain_dir + cancer_type + "/" + m + "_in_domain.txt"
            c.load_in_domain_mut_file(in_domain_mut_file)
        c.match_with_conserve_score()
        c.to_doc(output_file)

    pool = mp.Pool(8)
    cancer_kind = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP',
                   'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM',
                   'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']

    in_domain_dir = "/data1/hbs/all_cancer_analysis/all_cancer_mutation_site_conservation/all_cancer_in_domain_mut_site/"
    conserve_dir = "/data/zengyanru/LysineTCGA/find_mutsite_wholegenome/find_extra_score/"
    result_dir = "/data/hupus/mut_site_conservation/"
    for c in cancer_kind:
        conserve_file = conserve_dir + c + "_wholegenome_extra_score.txt"
        output_file = result_dir + c + "_mut_site_cobservation.txt"
        pool.apply_async(get_phastCons_score, (c, conserve_file, in_domain_dir, output_file))

    pool.close()
    pool.join()
