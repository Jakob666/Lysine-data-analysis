#-*- coding:utf-8 -*-
'''
Description:
    根据文件 all_sig_pro1.txt中的数据，挑选出在某个癌症中出现两次及以上的蛋白，该文件的形式为
-----------------------------------------------------------------------------------------------
Q0VDD8-1	4	[['LUAD_sumo_all_tes_result', '0.05'], ['LUAD_ubi_all_tes_result', '0.04885897435897436'], ['SKCM_sumo_all_tes_result', '0.0'], ['STAD_ubi_all_tes_result', '0.042928571428571434']]	DNAH14	F	F
Q9NPG3-1	3	[['CHOL_suc_all_tes_result', '0.0'], ['ESCA_sumo_all_tes_result', '0.0'], ['KIRC_suc_all_tes_result', '0.0']]	UBN1	F	F
P38606-1	1	[['BLCA_gly_all_tes_result', '0.0']]	ATP6V1A	F	F
------------------------------------------------------------------------------------------------
比如文件的第一行，Q0VDD8-1 这个蛋白在LUAD这个癌症中出现两次，将这种蛋白视为高频蛋白；而像Q9NPG3-1和P38606-1即使出现在癌症中，但是同一个癌症中没有出现多次，则视为低频蛋白。
选出这些高频蛋白的突变点在 in motif和 out motif中出现次数（当同一个位点因为选取的序列不同，同时出现在in motif 和 out motif中时，保留 in motif的记录）
==============================
@author:hbs
@date:2018-2-1
@version:1.1
'''
import pandas as pd
import numpy as np
import re
import os
import warnings


class HighFrequencyProt:
    def __init__(self, target_file):
        self.target_file = target_file

    def load_data(self):
        data = pd.read_csv(self.target_file, sep="\t", usecols=[0, 1, 2], encoding="utf-8")
        data.columns = ["Uniprot Accession", "related_mut", "record"]
        data = data[data["related_mut"] > 1]
        data["high_frequency"] = ""
        return data

    def extract_related_cancer(self, i):
        pattern = re.compile("[A-Z]+", re.S)
        record = i["record"]
        res = re.findall(pattern, record)
        if len(res) != len(set(res)):
            i["high_frequency"] = "yes"
        else:
            i["high_frequency"] = "no"
        return i

    def screen_high_frequency_prot(self, data):
        data = data.apply(self.extract_related_cancer, axis=1)
        high_freq_prot = data[data["high_frequency"] == "yes"]
        return high_freq_prot

    def main(self):
        data = self.load_data()
        high_freq_prot = self.screen_high_frequency_prot(data)
        high_freq_prot = high_freq_prot["Uniprot Accession"]
        return high_freq_prot


class GetMotifAndNonMotifLength:
    def __init__(self, high_freq_prot):
        self.high_freq_prot = high_freq_prot
        self.motif_length = pd.DataFrame(columns=["protein_id", "motif length", "non_motif length"])

    def load_data(self, f):
        try:
            data = pd.read_csv(f, sep="\t", usecols=[0, 4, 5])
        except:
            data = pd.DataFrame(columns=["protein_id", "motif length", "non_motif length"])
        return data

    def main(self):
        target_dir = "/data1/hbs/all_cancer_analysis/all_cancer_motif_significant_prot_test/in_motif_detail/"
        files = os.listdir(target_dir)
        for f in files:
            f = target_dir + f
            d = self.load_data(f)
            self.motif_length = pd.concat([self.motif_length, d], axis=0)
        self.motif_length.drop_duplicates(inplace=True)
        high_freq_prot = list(set(self.motif_length.index)&set(self.high_freq_prot))
        self.motif_length = self.motif_length.ix[high_freq_prot]


    class GetMutSiteNumber:
        def __init__(self, high_freq_prot):

            self.cancer_kind = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC',
                       'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ',
                       'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']
            self.mod = ["ace", "gly", "mal", "met", "suc", "sum", "ubi"]
            self.target_dir = "/data1/hbs/all_cancer_analysis/all_cancer_sig_prot_motif_modify/"
            self.high_freq_prot = high_freq_prot

        def load_data(self):
            in_motif_data = pd.DataFrame(columns=["protein_id", "mutated position"])
            out_motif_data = pd.DataFrame(columns=["protein_id", "mutated position"])
            for c in self.cancer_kind:
                for m in self.mod:
                    in_motif_file = self.target_dir + c + "/" + m + "_in_motif.txt"
                    out_motif_file = self.target_dir + c + "/" + m + "_out_motif.txt"
                    try:
                        in_mot = pd.read_csv(in_motif_file, sep="\t", usecols=[0, 7], encoding="utf-8")
                    except:
                        in_mot = pd.DataFrame(columns=["protein_id", "mutated position"])

                    try:
                        out_mot = pd.read_csv(out_motif_file, sep="\t", usecols=[0, 7], encoding="utf-8")
                    except:
                        out_mot = pd.DataFrame(columns=["protein_id", "mutated position"])
                    in_motif_data = pd.concat([in_motif_data, in_mot], axis=0)
                    out_motif_data = pd.concat([out_motif_data, out_mot], axis=0)
            in_motif_data["mutated position"] = pd.Series(in_motif_data["mutated position"].values.astype("String"))
            out_motif_data["mutated position"] = pd.Series(out_motif_data["mutated position"].values.astype("String"))
            in_motif_data["idx"] = in_motif_data["protein_id"] + in_motif_data["mutated position"]
            in_motif_data.set_index("idx", inplace=True)
            out_motif_data["idx"] = out_motif_data["protein_id"] + out_motif_data["mutated position"]
            out_motif_data.set_index("idx", inplace=True)
            in_motif_data.drop_duplicates(inplace=True)
            out_motif_data.drop_duplicates(inplace=True)

            common_prot = set(in_motif_data.index) & set(out_motif_data.index)
            out_motif_data.drop(common_prot, axis=0, inplace=True)

            return in_motif_data, out_motif_data

        def counting(self, in_motif_data, out_motif_data):
            in_motif_res = pd.value_counts(in_motif_data["protein_id"])
            out_motif_res = pd.value_counts(out_motif_data["protein_id"])
            high_freq_prot = list(set(in_motif_res.index) & set(self.high_freq_prot))
            in_motif_res = in_motif_res.ix[high_freq_prot]

            high_freq_prot = list(set(out_motif_res.index) & set(self.high_freq_prot))
            out_motif_res = out_motif_res.ix[high_freq_prot]
            return in_motif_res, out_motif_res


if __name__ == "__main__":
    warnings.filterwarnings("ignore")






    # h = HighFrequencyProt("C:/Users/hbs/Desktop/all_sig_pro1.txt")
    # high_freq_prot = h.main()
    # print(high_freq_prot)