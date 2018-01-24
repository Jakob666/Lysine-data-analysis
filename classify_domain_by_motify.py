#-*- coding:utf-8 -*-
'''
Description:
    将interproscan处理fasta序列得到的domain，按照七种修饰进行分类。
    其中一些有用的列进行保留，没有用的就不提取出来。
===============================================
@author:hbs
@date:2018-1-19
@version:1.0
'''
import pandas as pd
from multiprocessing import Pool


def domain_classify(mod, file):
    output_file = file.split(".")[0] + "_simplified.txt"
    data = pd.read_csv(file, sep="\t", header=None, usecols=[0, 2, 6, 7], encoding="utf-8")
    data.columns = ["title", "protein length", "start", "end"]
    #protein_id这里是PLMD ID
    data["protein_id"] = ""
    #Uniprot Accession没有精确到是什么亚形
    data["Uniprot Accession"] = ""
    data = data.apply(extract_modification, axis=1)
    data.drop("title", inplace=True)
    data.to_csv(output_file, sep="\t",index=False, columns=["protein_id","Uinprot Accession","protein length","start","end"], mode="w", encoding="utf-8")

def extract_modification(i):
    '''从类似于 PLMD-3268@P04908@Succinylation@Homo 的id中提取PLMD id、Uniprot ID和修饰类型信息'''
    info = i["protein_id"].split("@")
    i["protein_id"] = info[0]
    i["Uniprot Accession"] = info[1]


if __name__ == "__main__":
    pool = Pool(4)
    mod_list = ["Glycation", "Malonylation", "Methylation", "Succinylation", "Sumoylation"]
    for mod in mod_list:
        f = "/data1/hbs/all_cancer_analysis/all_domain_classify_by_modify/" + mod + ".tsv"
        pool.apply_async(domain_classify, (mod, f))

    pool.close()
    pool.join()
    print("finish")