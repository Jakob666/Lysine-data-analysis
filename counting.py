# -*- coding:utf-8 -*-
'''
    Description:
    该模块用于将修饰数据进行整理，模型预测出来的显著和非显著性蛋白的各种位置的修饰样本总数，
    当同一个样本因为motif和序列不同归类在 sig=yes和sig≠yes中时
            motif优先级 sig=yes -> sig≠ yes
    当相同的样本因为选取的蛋白序列不同分到 central和surround 两类时，优先保留 central的
            修饰位置优先级 central -> surround -> out
=======================
@author: hbs
@date: 2018-1-3
@version: 1.3
'''
import os
import pandas as pd
import multiprocessing as mp
import warnings


class Counting:
    def __init__(self, sig_motif_dir, insig_motif_dir, cancer_type):
        '''
        :param sig_motif_dir: 存放显著蛋白motif结果的文件目录
        :param insig_motif_dir: 存放不显著蛋白结果的文件目录
        :param cancer_type: 是本次统计的癌症类型
        '''
        self.target_dir = sig_motif_dir
        self.target_dir2 = insig_motif_dir
        self.cancer = cancer_type
        #创建两个文件夹分别存储同一种癌症、不同类型的修饰的计数文件
        os.system("mkdir /data1/hbs/all_cancer_analysis/all_cancer_modification_calculate/sig_prot/"+ self.cancer)
        os.system("mkdir /data1/hbs/all_cancer_analysis/all_cancer_modification_calculate/insig_prot/"+ self.cancer)

    def load_data(self, f):
        '''用于将文件中的数据进行读取， motif数据中读取的数据'''
        #读取到的列为protein_id（Uniprot详细到是否是isoform）、Canonical、cancer_type、left_flank、right_flank、mutated position、from、to、sample_num（此时只是sample号码，并没有计数）、direct mutation
        data = pd.read_csv(f, sep="\t", usecols=[0, 1, 2, 5, 6, 7, 8, 9, 10, 11])
        data.drop_duplicates(inplace=True)  #去除完全相同的重复项
        data = data.applymap(str)
        return data

    def preprocess_with_motif(self, sig_in, sig_out, insig_in, insig_out):
        '''当同一个修饰因为motif和序列不同归类在 sig=yes和sig≠yes中时
            motif优先级 sig=yes -> sig≠ yes
                       in motif -> out motif
        :param sig_in: 是显著蛋白in motif的数据，传入是dataframe形式
        :param sig_out: 是显著蛋白out motif的数据，传入的是dataframe形式
        :param insig_in: 是非显著蛋白in motif的数据，传入的是dataframe形式
        :param insig_out: 是非显著蛋白out motif的数据，传入的是dataframe形式
        '''
        #传入的时候各个数据的dataframe都是没有索引的，使用common_index方法设置索引（设置的索引足以区分不同的蛋白）
        sig_in, sig_out, insig_in, insig_out = map(Counting.common_index, [sig_in, sig_out, insig_in, insig_out])
        #获取经过处理后的索引集合
        sig_in_tags = set(sig_in.index)
        sig_out_tags = set(sig_out.index)
        insig_in_tags = set(insig_in.index)
        insig_out_tags = set(insig_out.index)

        #从优先级最低的 insig_out里面开始提取要清除的蛋白
        insig_out_need_drop = (sig_in_tags | sig_out_tags | insig_in_tags) & insig_out_tags
        # print(len(insig_out_need_drop))
        #从优先级次低的 insig_in中提取要清除的蛋白
        insig_in_need_drop = (sig_in_tags | sig_out_tags) & insig_in_tags
        # print(len(insig_in_need_drop))
        #从四个蛋白分级中优先级第二的 sig_out中提取要清除的蛋白
        sig_out_need_drop = sig_in_tags & sig_out_tags
        # print(len(sig_out_need_drop))

        sig_out.drop(sig_out_need_drop, axis=0, inplace=True)
        insig_in.drop(insig_in_need_drop, axis=0, inplace=True)
        insig_out.drop(insig_out_need_drop, axis=0, inplace=True)
        return sig_in, sig_out, insig_in, insig_out

    def preprocess_with_mut_location(self, sig_in, sig_out, insig_in, insig_out):
        '''此时数据已经在sig和insig区域的重复已经去除，进一步要做的是去除 in motif和out motif之间的重复
                修饰位置优先级 central -> surround
        :param sig_in: 是显著蛋白in motif的数据，传入是dataframe形式
        :param sig_out: 是显著蛋白out motif的数据，传入的是dataframe形式
        :param insig_in: 是非显著蛋白in motif的数据，传入的是dataframe形式
        :param insig_out: 是非显著蛋白out motif的数据，传入的是dataframe形式
        '''
        #提取出in motif显著蛋白central和surround的修饰
        sig_in_central_mut = sig_in[sig_in["direct mutation"] == "yes"]
        sig_in_surround_mut = sig_in[sig_in["direct mutation"] == "no"]
        #显著蛋白out motif而言不必进行进一步分类
        sig_out_surround_mut = sig_out

        #非显著蛋白的 in motif修饰分成 central和 surround两类
        insig_in_central_mut = insig_in[insig_in["direct mutation"] == "yes"]
        insig_in_surround_mut = insig_in[insig_in["direct mutation"] == "no"]
        #非显著蛋白out motif不用进一步分类
        insig_out_surround_mut = insig_out

        #显著蛋白 in motif修饰中 central和surround重复时除去 surround中的记录
        sig_in_surround_need_drop = set(sig_in_central_mut.index) & set(sig_in_surround_mut.index)
        #非显著蛋白in motif修饰中 central和surround重复时除去 surround中的记录
        insig_in_surround_need_drop = set(insig_in_central_mut.index) & set(insig_in_surround_mut.index)

        sig_in_surround_mut.drop(sig_in_surround_need_drop, axis=0, inplace=True)
        insig_in_surround_mut.drop(insig_in_surround_need_drop, axis=0, inplace=True)

        final_classify = {"sig_in_central_mut":sig_in_central_mut, "sig_in_surround_mut":sig_in_surround_mut,
                          "sig_out_mut":sig_out_surround_mut, "insig_in_central_mut":insig_in_central_mut,
                          "insig_in_surround_mut":insig_in_surround_mut, "insig_out_mut":insig_out_surround_mut}

        return final_classify

    @staticmethod
    def common_index(df):
        '''传入的是不同数据的dataframe，将protein_id（Uniprot ID）、样本号、突变位置结合，设置为索引'''
        df["tag"] = df["protein_id"] + ":" + df["sample_num"] + ":" + df["mutated position"]
        df.set_index("tag", inplace=True)
        return df

    def counting(self, classify, mod_type):
        '''将上面两部去重复产生的数据都加上用于计数的一列 sample count，传入到write_to_doc中写入文件'''
        sig_in_central = classify["sig_in_central_mut"]
        sig_in_central["sample_count"] = 1
        sig_in_surround = classify["sig_in_surround_mut"]
        sig_in_surround["sample_count"] = 1
        sig_out = classify["sig_out_mut"]
        sig_out["sample_count"] = 1
        insig_in_central = classify["insig_in_central_mut"]
        insig_in_central["sample_count"] = 1
        insig_in_surround = classify["insig_in_surround_mut"]
        insig_in_surround["sample_count"] = 1
        insig_out = classify["insig_out_mut"]
        insig_out["sample_count"] = 1

        self.write_in_to_doc(sig_in_central, mod_type)
        self.write_in_to_doc(sig_in_surround, mod_type, position="_in_motif_surround")
        self.write_in_to_doc(sig_out, mod_type,position="_out_motif")
        self.write_in_to_doc(insig_in_central, mod_type, sig="no")
        self.write_in_to_doc(insig_in_surround, mod_type, sig="no", position="_in_motif_surround")
        self.write_in_to_doc(insig_out, mod_type,sig="no", position="_out_motif")

    def write_in_to_doc(self, data, mod_type, sig="yes", position="_in_motif_central"):
        '''

        :param data: 经过counting方法处理得到的dataframe，有一列 sample count用于计数
        :param mod_type: 具体的某种修饰
        :param sig: 显著性是显著还是非显著
        :param position: 修饰的三种位置 in_motif_central、in_motif_surround和out_motif
        :return:
        '''

        data.reset_index(drop=True)
        #只有UniprotID、sample count两列的计数文件
        res = data.groupby(["protein_id"], as_index=False)["sample_count"].sum()
        #较为详细信息UniprotId、样本号、从什么到什么突变的计数
        detailed_res = data.groupby(["protein_id", "sample_num", "from", "to"], as_index=False)["sample_count"].sum()
        if sig == "yes":
            path = "/data1/hbs/all_cancer_analysis/all_cancer_modification_calculate/sig_prot/" + self.cancer +"/"
        else:
            path = "/data1/hbs/all_cancer_analysis/all_cancer_modification_calculate/insig_prot/" + self.cancer + "/"
        detailed_res.to_csv(path + mod_type + "/" + position + "DetailMut.txt", sep="\t", index=False, mode="w", encoding="utf-8")
        res.to_csv(path + mod_type + "/" + position + "Mut.txt", sep="\t", index=False, mode="w", encoding="utf-8")


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    mod_list = ["ace", "gly", "mal", "met", "suc", "sum", "ubi"]
    # cancer_types = ['ACC','BLCA','BRCA','CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH','KIRC','KIRP',
    #                'LAML','LGG','LIHC','LUAD','LUSC','MESO','OV','PAAD','PCPG','PRAD','READ','SARC','SKCM',
    #                'STAD','TGCT','THCA','THYM','UCEC','UCS','UVM']
    cancer_types = ["PAAD"]

    def main(mod, cancer):
        sig_motif_dir = "/data1/hbs/all_cancer_analysis/all_cancer_sig_prot_motif_modify/"
        insig_motif_dir = "/data1/hbs/all_cancer_analysis/all_cancer_insig_prot_motif_modify/"
        c = Counting(sig_motif_dir, insig_motif_dir, cancer)

        in_motif = mod + "_in_motif.txt"
        out_motif = mod + "_out_motif.txt"
        sig_in_motif =  sig_motif_dir + cancer + "/" + in_motif
        sig_out_motif = sig_motif_dir + cancer + "/" + out_motif
        insig_in_motif = insig_motif_dir + cancer + "/" + in_motif
        insig_out_motif = insig_motif_dir + cancer + "/" + out_motif

        sig_in_motif_data = c.load_data(sig_in_motif)
        sig_out_motif_data = c.load_data(sig_out_motif)
        insig_in_motif_data = c.load_data(insig_in_motif)
        insig_out_motif_data = c.load_data(insig_out_motif)

        sig_in, sig_out, insig_in, insig_out = c.preprocess_with_motif(sig_in_motif_data, sig_out_motif_data, insig_in_motif_data,insig_out_motif_data)
        classify = c.preprocess_with_mut_location(sig_in, sig_out, insig_in, insig_out)
        c.counting(classify, mod)

    pool = mp.Pool(8)
    for cancer in cancer_types:
        for m in mod_list:
            pool.apply_async(main, (m, cancer))

    pool.close()
    # c = Counting("/data1/hbs/modify_in_motif(sig=yes)", "/data1/hbs/modify_in_motif(signotyes)")
    # c.main()