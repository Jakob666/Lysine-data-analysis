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
import pandas as pd
import warnings


class Counting:
    def __init__(self, sig_motif_dir, insig_motif_dir):
        self.target_dir = sig_motif_dir
        self.target_dir2 = insig_motif_dir

    def load_data(self, f):
        data = pd.read_csv(f, sep="\t", usecols=[0, 1, 2, 5, 6, 7, 8, 9])
        data.drop_duplicates(inplace=True)  #去除完全相同的重复项
        data = data.applymap(str)
        return data

    def preprocess_with_motif(self, sig_in, sig_out, insig_in, insig_out):
        '''当同一个修饰因为motif和序列不同归类在 sig=yes和sig≠yes中时
            motif优先级 sig=yes -> sig≠ yes
        '''
        sig_in, sig_out, insig_in, insig_out = map(Counting.common_index, [sig_in, sig_out, insig_in, insig_out])

        sig_in_tags = set(sig_in.index)
        sig_out_tags = set(sig_out.index)
        insig_in_tags = set(insig_in.index)
        insig_out_tags = set(insig_out.index)

        insig_out_need_drop = (sig_in_tags | sig_out_tags | insig_in_tags) & insig_out_tags
        # print(len(insig_out_need_drop))
        insig_in_need_drop = (sig_in_tags | sig_out_tags) & insig_in_tags
        # print(len(insig_in_need_drop))
        sig_out_need_drop = sig_in_tags & sig_out_tags
        # print(len(sig_out_need_drop))

        sig_out.drop(sig_out_need_drop, axis=0, inplace=True)
        insig_in.drop(insig_in_need_drop, axis=0, inplace=True)
        insig_out.drop(insig_out_need_drop, axis=0, inplace=True)
        return sig_in, sig_out, insig_in, insig_out

    def preprocess_with_mut_location(self, sig_in, sig_out, insig_in, insig_out):
        '''当相同的样本因为选取的蛋白序列不同分到 central和surround 两类时，优先保留 central的'''
        sig_in_central_mut = sig_in[sig_in["direct mutation"] == "yes"]
        sig_in_surround_mut = sig_in[sig_in["direct mutation"] == "no"]
        sig_out_surround_mut = sig_out
        insig_in_central_mut = insig_in[insig_in["direct mutation"] == "yes"]
        insig_in_surround_mut = insig_in[insig_in["direct mutation"] == "no"]
        insig_out_surround_mut = insig_out

        sig_in_surround_need_drop = set(sig_in_central_mut.index) & set(sig_in_surround_mut.index)
        insig_in_surround_need_drop = set(insig_in_central_mut.index) & set(insig_in_surround_mut.index)

        sig_in_surround_mut.drop(sig_in_surround_need_drop, axis=0, inplace=True)
        insig_in_surround_mut.drop(insig_in_surround_need_drop, axis=0, inplace=True)

        final_classify = {"sig_in_central_mut":sig_in_central_mut, "sig_in_surround_mut":sig_in_surround_mut,
                          "sig_out_mut":sig_out_surround_mut, "insig_in_central_mut":insig_in_central_mut,
                          "insig_in_surround_mut":insig_in_surround_mut, "insig_out_mut":insig_out_surround_mut}

        return final_classify

    @staticmethod
    def common_index(df):
        df["tag"] = df["protein_id"] + ":" + df["sample_num"] + ":" + df["mutated position"]
        df.set_index("tag", inplace=True)
        return df

    def counting(self, classify, mod_type):
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
        data.reset_index(drop=True)
        res = data.groupby(["protein_id"])["sample_count"].sum()
        detailed_res = data.groupby(["protein_id", "sample_num", "from", "to"])["sample_count"].sum()
        if sig == "yes":
            path = "/data1/hbs/centralOrNotSigIsyes/"
        else:
            path = "/data1/hbs/centralOrNotSigNotyes/"
        detailed_res.to_csv(path + mod_type + position + "DetailMut.txt", sep="\t", mode="w", encoding="utf-8")
        res.to_csv(path + mod_type + position + "Mut.txt", sep="\t", mode="w", encoding="utf-8")

    def main(self):
        mod_list = ["ace", "gly", "mal", "met", "suc", "sumo", "ubi"]
        for mod in mod_list:
            in_motif = mod + "_in_motif.txt"
            out_motif = mod + "_out_motif.txt"
            sig_in_motif = "/data1/hbs/modify_in_motif(sig=yes)/" + in_motif
            sig_out_motif = "/data1/hbs/modify_in_motif(sig=yes)/" + out_motif
            insig_in_motif = "/data1/hbs/modify_in_motif(signotyes)/" + in_motif
            insig_out_motif = "/data1/hbs/modify_in_motif(signotyes)/" + out_motif

            sig_in_motif_data = self.load_data(sig_in_motif)
            sig_out_motif_data = self.load_data(sig_out_motif)
            insig_in_motif_data = self.load_data(insig_in_motif)
            insig_out_motif_data = self.load_data(insig_out_motif)

            sig_in, sig_out, insig_in, insig_out = self.preprocess_with_motif(sig_in_motif_data, sig_out_motif_data, insig_in_motif_data,insig_out_motif_data)
            classify = self.preprocess_with_mut_location(sig_in, sig_out, insig_in, insig_out)
            self.counting(classify, mod)


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    # c = Counting("/data1/hbs/modify_in_motif(sig=yes)", "/data1/hbs/modify_in_motif(signotyes)")
    # c.main()