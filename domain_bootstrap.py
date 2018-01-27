# -*- coding:utf-8 -*-
'''
Description:
    统计sig的domain区和非domain区的突变数目，之后进行bootstrap检验，
    检验为 单尾检验，α=0.05，对于突变的数目需要进行标准化处理，公式如下：
                in_domain_mut_num = 统计的数目/蛋白domain区长度  *1000
                out_domain_mut_num = 统计的数目/蛋白非domain区长度  *1000
                其中 非domain长度 = sequence length - domain区长度
    检验目的：检验显著蛋白的domain和非domain的突变数目是否存在显著性差异。
    前期工作：目前已经将不同癌症的蛋白质修饰精确分类为如下几类，
    insignificant/significant protein
        |-----------> in domain ------->in motif mutation----->central mutation
        |                       |                       |----->surround mutation
        |                       |------>out motif mutation
        |
        |-----------> out domain------->in motif mutation----->central mutatu=ion
                               |                        |----->surround mutation
                               |------->out motif mutation
    可进行如下检验：
    1.某一癌症显著性蛋白的某种修饰 in domain和 out domain是否显著差异。（此时要整合in和out的三种文件）
    2.某一癌症显著性蛋白的某种修饰 in domain in motif和 out domain in motif是否显著差异。（需要整合in motif的两种文件）
    3.某一癌症显著性蛋白的某种修饰 in domain in motif和 in domain out motif是否显著差异。
    4.某一癌症显著性蛋白的某种修饰 in domain central mutation和 out domain central mutation是否显著差异。
======================================================
@author: hbs
@date: 2018-1-13
@version 1.0
======================================================
@update: 2018-1-23
    将所有33种癌症分别进行不同修饰的bootstrap检验
@version 1.1
'''
import warnings
import pandas as pd
from bootstrapTest import Bootstrap_test
import os
from multiprocessing import Pool


class DomainBootstrapTest:
    def __init__(self, in_domain_file, out_domain_file, mode_type):
        self.in_domain_file = in_domain_file
        self.out_domain_file = out_domain_file
        self.mod = mode_type

    def load_data(self, f):
        data = pd.read_csv(f, sep="\t", usecols=[0, 3], encoding="utf-8")
        return data

    def bootstrap(self, in_domain_data, out_domain_data):
        b = Bootstrap_test(time=10000, side="one-side")
        x_data = list(in_domain_data["mut_per_k_aa"])
        y_data = list(out_domain_data["mut_per_k_aa"])
        res, p_val, p_ref, x_mean, y_mean, x_std, y_std = b.main(x_data, y_data)
        return res, p_val, p_ref, x_mean, y_mean, x_std, y_std

    def judge(self,res, p_val, p_ref, x_mean, y_mean, x_std, y_std):
        content = "%s in %s mutations counts between domain and non-domain area on significant proteins' sequences,"%(res, self.mod)
        if p_val > p_ref:
            content += "p value %s > %s\n" % (p_val, p_ref)
        else:
            content += "p value %s < %s\n" % (p_val, p_ref)
        content += "average mutations per k aa in domain area: %s\naverage mutations per k aa in non-domain area: %s\n" % (x_mean, y_mean)
        content += "standard variance in domain area: %s\nstandard variance in non-domain area: %s\n" % (x_std, y_std)
        if res == "Significant difference" and x_mean > y_mean:
            content += "the occurance of mutations in domain area is significant higher than which in non-domain area.\n\n"
        elif res == "Significant difference" and x_mean < y_mean:
            content += "the occurance of mutations in domain area is significant lower than which in non-domain area.\n\n"
        return content

    def main(self):
        in_domain_data = self.load_data(self.in_domain_file)
        out_domain_data = self.load_data(self.out_domain_file)
        res, p_val, p_ref, x_mean, y_mean, x_std, y_std = self.bootstrap(in_domain_data, out_domain_data)
        content = self.judge(res, p_val, p_ref, x_mean, y_mean, x_std, y_std)
        return content


if __name__ == "__main__":
    warnings.filterwarnings("ignore")

    def domain_test(in_domain_file, out_domain_file, cancer, mod_type, output_doc):
        d = DomainBootstrapTest(in_domain_file, out_domain_file, mod_type)
        content = d.main()
        if content.find("higher") != -1:
            print(cancer + ":" + mod_type)
        with open(output_doc, "w") as f:
            f.write(content)
            f.close()
        return None

    pool = Pool(7)
    # 33种癌症
    cancer_kind =['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH','KIRC', 'KIRP',
                  'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM',
                  'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']
    # cancer_kind = ["ACC"]
    # 7种修饰的全称
    mod_list = {"ace": "Acetylation", "gly": "Glycation", "mal": "Malonylation", "met": "Methylation",
                "suc": "Succinylation", "sum": "Sumoylation", "ubi": "Ubiquitination"}
    # 7种修饰的简写
    mod_abbr_list = ["ace", "gly", "mal", "met", "suc", "sum", "ubi"]
    # mod_abbr_list = ["ace"]

    calc_dir = "/data1/hbs/all_cancer_analysis/all_cancer_modification_in_domain_or_not/calc_result/"
    test_result_dir = "/data1/hbs/all_cancer_analysis/all_cancer_modification_in_domain_or_not/test_result/"
    for c in cancer_kind:
        cancer_dir = test_result_dir + c
        os.system("mkdir " + cancer_dir)
        for mod in mod_abbr_list:
            in_domain_file = calc_dir + c + "/" + mod + "_in_domain_calc.txt"
            out_domain_file = calc_dir + c + "/" + mod + "_out_domain_calc.txt"
            output_doc = test_result_dir + c + "/" + mod_list[mod] + "_test_result.txt"
            pool.apply_async(domain_test, (in_domain_file, out_domain_file, c, mod, output_doc))
    pool.close()
    pool.join()