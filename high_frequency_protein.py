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
            data = pd.read_csv(f, sep="\t", usecols=["protein_id", "motif length", "non_motif length"])
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
        self.motif_length.set_index("protein_id", inplace=True)
        # high_freq_prot = list(set(self.motif_length["protein_id"].values)&set(self.high_freq_prot))
        motif_length = self.motif_length.ix[self.high_freq_prot]
        self.protein_length = motif_length["motif length"] + motif_length["non_motif length"]
        # self.motif_length.dropna(inplace=True)
        self.protein_length.drop_duplicates(inplace=True)

class GetMutSiteNumber:
    def __init__(self, high_freq_prot):
        self.cancer_kind = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC',
                   'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ',
                   'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']
        self.cancer = {"O60506-1":"UCEC", "Q9UPN3-1":"UCEC","F5H2G6":"BRCA","Q8IVF2-1":"LAML","Q8NCE2-1":"GBM","P49753-1":"LUSC","Q9NRY4":"TGCT",
                        "Q7L2E3-1":"COAD","P07602-1":"COAD","Q9Y4P8-1":"UCEC","P62258-1":"STAD","Q8NBS9-1":"LIHC","Q99757":"UCEC", "P31946-1":"LUSC",
                        "Q14CX7-1":"CESC", "P23246-1":"THYM"}
        
        self.mod = ["ace", "gly", "mal", "met", "suc", "sum", "ubi"]
        self.target_dir = "/data1/hbs/all_cancer_analysis/all_cancer_sig_prot_motif_modify/"
        self.high_freq_prot = list(high_freq_prot)

    def load_data(self):
        in_motif_data = pd.DataFrame(columns=["protein_id", "cancer_type", "K position", "left flank", "right flank", "mutated position"])
        out_motif_data = pd.DataFrame(columns=["protein_id", "cancer_type", "K position", "left flank", "right flank", "mutated position"])
        for c in self.cancer_kind:
            for i in self.high_freq_prot:
                # c = self.cancer[i]
                for m in self.mod:
                    in_motif_file = self.target_dir + c + "/" + m + "_in_motif.txt"
                    out_motif_file = self.target_dir + c + "/" + m + "_out_motif.txt"
                    try:
                        in_mot = pd.read_csv(in_motif_file, sep="\t", usecols=[0, 3,5,6, 7], encoding="utf-8")
                        in_mot = in_mot[in_mot["protein_id"]==i]
                        in_mot["cancer_type"] = c
                    except:
                        in_mot = pd.DataFrame(columns=["protein_id", "K position", "left flank", "right flank", "mutated position","cancer_type"])

                    try:
                        out_mot = pd.read_csv(out_motif_file, sep="\t", usecols=[0, 3,5,6,7], encoding="utf-8")
                        out_mot = out_mot[out_mot["protein_id"]==i]
                        out_mot["cancer_type"] = c
                    except:
                        out_mot = pd.DataFrame(columns=["protein_id", "K position", "left flank", "right flank", "mutated position", "cancer_type"])
                    in_motif_data = pd.concat([in_motif_data, in_mot], axis=0)
                    out_motif_data = pd.concat([out_motif_data, out_mot], axis=0)

        in_motif_data["start"] = in_motif_data["K position"] - in_motif_data["left flank"]
        in_motif_data["end"] = in_motif_data["K position"] + in_motif_data["right flank"]
        in_motif_data.set_index("protein_id", drop=False, inplace=True)
        in_motif_data = in_motif_data.ix[list(self.high_freq_prot)]
        in_motif_data = in_motif_data[["protein_id", "cancer_type", "mutated position","start","end"]]
        in_motif_data.drop_duplicates(subset=["protein_id","mutated position"],inplace=True)
        in_motif_data["mutated position"] = pd.Series(in_motif_data["mutated position"].map(str))
        in_motif_data["idx"] = in_motif_data["protein_id"] + in_motif_data["mutated position"]
        in_motif_data.set_index("idx", inplace=True)

        out_motif_data["start"] = out_motif_data["K position"] - out_motif_data["left flank"]
        out_motif_data["end"] = out_motif_data["K position"] + out_motif_data["right flank"]
        out_motif_data.set_index("protein_id", drop=False, inplace=True)
        out_motif_data = out_motif_data.ix[list(self.high_freq_prot)]
        out_motif_data = out_motif_data[["protein_id", "cancer_type", "mutated position","start","end"]]
        out_motif_data.drop_duplicates(subset=["protein_id","mutated position"],inplace=True)
        out_motif_data["mutated position"] = pd.Series(out_motif_data["mutated position"].map(str))
        out_motif_data["idx"] = out_motif_data["protein_id"] + out_motif_data["mutated position"]
        out_motif_data.set_index("idx", inplace=True)

        common_prot = set(in_motif_data.index) & set(out_motif_data.index)
        out_motif_data.drop(common_prot, axis=0, inplace=True)

        in_motif_data.set_index("protein_id", drop=False, inplace=True)
        # print(in_motif_data)
        out_motif_data.set_index("protein_id", drop=False, inplace=True)
        # print("================================")
        # in_motif_data = self.merge_overlap_domain(in_motif_data)
        # print(in_motif_data)
        # in_motif_data["mutated position"] = pd.Series(in_motif_data["mutated position"].map(str))
        # out_motif_data["mutated position"] = pd.Series(out_motif_data["mutated position"].map(str))
        # in_motif_data["idx"] = in_motif_data["protein_id"] + in_motif_data["mutated position"]
        # in_motif_data.reindex(list(in_motif_data["idx"]))
        # out_motif_data["idx"] = out_motif_data["protein_id"] + out_motif_data["mutated position"]
        # out_motif_data.reindex(list(out_motif_data["idx"]))
        
        # in_motif_data.drop_duplicates(inplace=True)
        # out_motif_data.drop_duplicates(inplace=True)

        # common_prot = set(in_motif_data.index) & set(out_motif_data.index)
        # 

        return in_motif_data, out_motif_data

    def merge_overlap_domain(self, domain_data):
        '''将重叠的domain区域合并'''
        new_domain_data = pd.DataFrame(columns=["protein_id", "cancer_type", "start", "end"])
        idxs = set(domain_data.index)

        for idx in idxs:
            try:
                intervals = domain_data.ix[idx]
            except:
                continue
            # mut_site = intervals["mutated position"]
            # index = i["protein_id"] + str(i["mutated position"])
            if "Series" in str(type(intervals)):  #Series说明该蛋白对应的domain只有一个
                d = pd.DataFrame([[idx, intervals["cancer_type"], intervals["start"], intervals["end"]]], columns=["protein_id", "cancer_type", "start", "end"])
                new_domain_data = pd.concat([new_domain_data, d], axis=0)

            elif "DataFrame" in str(type(intervals)):   #蛋白上有多个domain
                intervals.sort_values("start", inplace=True)  #起始位点从小到大排序
                domain_area = list(zip(list(intervals["start"]), list(intervals["end"])))
                cancer_type = list(intervals["cancer_type"])[0]
                need_drop = []
                for i in range(len(domain_area)-1):
                    interval1 = domain_area[i]
                    interval2 = domain_area[i+1]
                    overlapped = self.if_overlapped(interval1, interval2)

                    if overlapped:
                        domain_area[i+1] = [min([interval1[0], interval2[0]]), max([interval1[1], interval2[1]])]
                        need_drop.append(interval1)
                for j in need_drop:
                    domain_area.remove(j)
                domain_area = list(map(list, domain_area))
                d = pd.DataFrame(domain_area, columns=["start", "end"])
                d["protein_id"] = idx
                d["cancer_type"] = cancer_type
                new_domain_data = pd.concat([new_domain_data, d], axis=0)
        new_domain_data.set_index("protein_id", inplace=True)
        # print(new_domain_data)
        return new_domain_data

    def if_overlapped(self, interval1, interval2):
        '''判断同一个蛋白上两个domain区间是否重叠
        :param interval1: 传入第一个区间
        :param interval2: 传入第二个区间
        '''
        overlapped = False
        if max([interval1[0], interval2[0]]) <= min([interval1[1], interval2[1]]):  # 相交的情况1
            overlapped = True

        elif not (interval1[0] > interval2[1] or interval1[1] < interval2[0]):  # 相交的情况2
            overlapped = True
        return overlapped

    def counting(self, in_motif_data, out_motif_data):
        in_motif = pd.DataFrame(columns=["protein_id", "cancer_type", "K position", "left flank", "right flank", "mutated position"])
        out_motif = pd.DataFrame(columns=["protein_id", "cancer_type", "K position", "left flank", "right flank", "mutated position"])
        for i in self.high_freq_prot:
            in_mot = in_motif_data[(in_motif_data["protein_id"]==i) & (in_motif_data["cancer_type"] == self.cancer[i])]
            in_motif = pd.concat([in_motif, in_mot], axis=0)
            out_mot = out_motif_data[(out_motif_data["protein_id"]==i) & (out_motif_data["cancer_type"] == self.cancer[i])]
            out_motif = pd.concat([out_motif, out_mot], axis=0)
        in_motif_res = pd.value_counts(in_motif["protein_id"])
        out_motif_res = pd.value_counts(out_motif["protein_id"])
        high_freq_prot = list(set(in_motif_res.index) & set(self.high_freq_prot))
        in_motif_res = in_motif_res.ix[high_freq_prot]

        # high_freq_prot = list(set(out_motif_res.index) & set(self.high_freq_prot))
        out_motif_res = out_motif_res.ix[high_freq_prot]
        in_motif_res = pd.DataFrame(in_motif_res)
        out_motif_res = pd.DataFrame(out_motif_res)
        return in_motif_res, out_motif_res


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    seq_length = pd.read_csv("/data1/hbs/total_fastr/Total.elm", sep="\t", usecols=[1, 4], encoding="utf-8")
    seq_length.drop_duplicates(inplace=True)
    seq_length["protein_length"] = 0
    def get_len(i):
        i["protein_length"] = len(i["Sequence"])
        return i
    seq_length = seq_length.apply(get_len, axis=1)
    seq_length = seq_length[["Uniprot Accession", "protein_length"]]

    h = HighFrequencyProt("/data1/hbs/all_cancer_analysis/significant_prot_for_IBS/all_sig_pro1.txt")
    high_freq_prot = h.main()

    gm = GetMutSiteNumber(high_freq_prot)
    in_motif_data, out_motif_data = gm.load_data()

    in_motif_merged = gm.merge_overlap_domain(in_motif_data)
    in_motif_merged["motif length"] = in_motif_merged["end"] - in_motif_merged["start"] + 1
    in_motif_merged = pd.DataFrame(in_motif_merged.groupby([in_motif_merged.index])["motif length"].sum())

    in_motif_common = set(in_motif_data.index)&set(in_motif_merged.index)

    in_motif_data = in_motif_data.ix[in_motif_common]
    in_motif_res,out_motif_res = gm.counting(in_motif_data, out_motif_data)

    df = pd.merge(right=in_motif_merged, left=in_motif_res,left_index=True,right_index=True)
    df.columns=["in motif mut", "motif length"]
    df = pd.merge(left=df, right=out_motif_res, left_index=True, right_index=True)
    df.columns= ["in motif mut", "motif length", "out motif mut"]
    df.fillna(0,inplace=True)
    df["canonical"] = ""
    df["protein_id"] = ""
    def get_canonical(i):
        i["protein_id"] = i.name
        i["canonical"] = i.name.split("-")[0]
        return i

    df = df.apply(get_canonical, axis=1)
    df = pd.merge(left=df, right=seq_length,left_on="canonical",right_on="Uniprot Accession")
    df["non_motif length"] = df["protein_length"] - df["motif length"]
    df["cancer_type"]=""
    def get_cancer(i):
        i["cancer_type"] = gm.cancer[i["protein_id"]]
        return i

    df = df.apply(get_cancer, axis=1)
    df = df[["protein_id","cancer_type", "protein_length","in motif mut", "motif length", "out motif mut", "non_motif length"]]
    df.to_csv("/data1/hbs/all_cancer_analysis/significant_prot_for_IBS/high_freq_protV2.txt", sep="\t", index=False, mode="w", encoding="utf-8")
    # out_motif_merged = gm.merge_overlap_domain(out_motif_data)
    # out_motif_merged["motif length"] = out_motif_merged["end"] - out_motif_merged["start"] + 1
    # print(in_motif_data.ix["Q8IVF2-1"])
    # print(in_motif_merged)
