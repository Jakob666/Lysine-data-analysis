#-*- coding:utf-8 -*-
'''
Description:
    将 domain_mut_site.py文件生成的在domain区域的突变文件有效信息提取
    生成保守性打分需要的 PhastCons输入的 bed格式文件，该格式文件主要的内容
    是，
        chr   start   end（这三列必须）  name   score   strand  ……
    bed文件作为输入文件进行保守性打分。
'''
import pandas as pd
import os
import multiprocessing as mp
import warnings


class BedFormer:
    def __init__(self):
        pass

    def load_data(self, prot_mut_genome_file):
        '''
        :param prot_mut_genome_file: domain_mut_site.py文件生成的蛋白突变位点相应的基因组坐标的文件
        :return:
        '''
        #读入突变位点的基因组坐标文件，读入的内容 Uniprot Accession、chr、start
        site_data = pd.read_csv(prot_mut_genome_file, sep="\t", usecols=[0, 11, 12])
        return site_data

    def form_bed_file(self, site_data):
        site_data["Chr"] = pd.Series(site_data["Chr"].map(str))
        site_data["Chr"] = "chr" + site_data["Chr"]
        site_data["Chr_start"] = site_data["Start"] - 7
        site_data["Chr_end"] = site_data["Start"] + 7
        return site_data

    def to_doc(self, site_data, bed_file):
        site_data.to_csv(bed_file, sep="\t", index=False, columns=["Chr", "Chr_start", "Chr_end"], header=False, mode="w", encoding="utf-8")
        return None


if __name__ == "__main__":
    warnings.filterwarnings("ignore")

    def form_bed_file(prot_mut_genome_file, bed_file):
        b = BedFormer()
        data = b.load_data(prot_mut_genome_file)
        site_data = b.form_bed_file(data)
        b.to_doc(site_data, bed_file)

    pool = mp.Pool(8)
    prot_mut_dir = "/data1/hbs/all_cancer_analysis/all_cancer_mutation_site_conservation/all_cancer_mut_site_genome_location/"
    cancer_kind = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP',
                   'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM',
                   'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']
    target_dir = "/data1/hbs/all_cancer_analysis/all_cancer_mutation_site_conservation/all_cancer_bed_file/"
    for c in cancer_kind:
        cancer_dir = target_dir + c
        os.system("mkdir " + cancer_dir)
        prot_mut_genome_file = prot_mut_dir + c + "/" + c + "_in_domain_mut_genome_location.txt"
        bed_file = cancer_dir + "/" + c + "_mut_region.bed"
        pool.apply_async(form_bed_file, (prot_mut_genome_file, bed_file))
    pool.close()
    pool.join()