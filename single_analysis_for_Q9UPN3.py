#-*- coding:utf-8 -*-
'''
Description:
    单独的对Q9UPN3-1这个蛋白进行分析，统计该蛋白的突变点在 in motif和 out motif中出现次数，以及in motif和 out motif
    区域的长度。
==================================
@author:hbs
@date:2018-2-5
@version:1.1
'''
import os
import pandas as pd
from high_frequency_protein import GetMutSiteNumber


class AnalysisQ9UPN3(GetMutSiteNumber):

    def load_data(self):
        high_freq_prot = []
        in_motif_data = pd.DataFrame(
            columns=["protein_id", "cancer_type", "K position", "left flank", "right flank", "mutated position"])
        out_motif_data = pd.DataFrame(
            columns=["protein_id", "cancer_type", "K position", "left flank", "right flank", "mutated position"])
        for i in range(33):
            pro = "Q9UPN3-1"
            prot = "Q9UPN3-1" + str(i)
            high_freq_prot.append(prot)
            c = self.cancer_kind[i]
            for m in self.mod:
                in_motif_file = self.target_dir + c + "/" + m + "_in_motif.txt"
                out_motif_file = self.target_dir + c + "/" + m + "_out_motif.txt"
                try:
                    in_mot = pd.read_csv(in_motif_file, sep="\t", usecols=[0, 3, 5, 6, 7], encoding="utf-8")
                    in_mot = in_mot[in_mot["protein_id"] == pro]
                    in_mot["cancer_type"] = c
                    in_mot["protein_id"] = prot
                except:
                    in_mot = pd.DataFrame(
                        columns=["protein_id", "K position", "left flank", "right flank", "mutated position",
                                 "cancer_type"])

                try:
                    out_mot = pd.read_csv(out_motif_file, sep="\t", usecols=[0, 3, 5, 6, 7], encoding="utf-8")
                    out_mot = out_mot[out_mot["protein_id"] == pro]
                    out_mot["cancer_type"] = c
                    out_mot["protein_id"] = prot
                except:
                    out_mot = pd.DataFrame(
                        columns=["protein_id", "K position", "left flank", "right flank", "mutated position",
                                 "cancer_type"])
                in_motif_data = pd.concat([in_motif_data, in_mot], axis=0)
                out_motif_data = pd.concat([out_motif_data, out_mot], axis=0)
        in_motif_data["start"] = in_motif_data["K position"] - in_motif_data["left flank"]
        in_motif_data["end"] = in_motif_data["K position"] + in_motif_data["right flank"]
        in_motif_data.set_index("protein_id", drop=False, inplace=True)
        in_motif_data = in_motif_data.ix[list(high_freq_prot)]
        in_motif_data = in_motif_data[["protein_id", "cancer_type", "mutated position", "start", "end"]]
        in_motif_data.drop_duplicates(subset=["protein_id", "mutated position"], inplace=True)
        in_motif_data["mutated position"] = pd.Series(in_motif_data["mutated position"].map(str))
        in_motif_data["idx"] = in_motif_data["protein_id"] + in_motif_data["mutated position"]
        in_motif_data.set_index("idx", inplace=True)

        out_motif_data["start"] = out_motif_data["K position"] - out_motif_data["left flank"]
        out_motif_data["end"] = out_motif_data["K position"] + out_motif_data["right flank"]
        out_motif_data.set_index("protein_id", drop=False, inplace=True)
        out_motif_data = out_motif_data.ix[list(high_freq_prot)]
        out_motif_data = out_motif_data[["protein_id", "cancer_type", "mutated position", "start", "end"]]
        out_motif_data.drop_duplicates(subset=["protein_id", "mutated position"], inplace=True)
        out_motif_data["mutated position"] = pd.Series(out_motif_data["mutated position"].map(str))
        out_motif_data["idx"] = out_motif_data["protein_id"] + out_motif_data["mutated position"]
        out_motif_data.set_index("idx", inplace=True)

        common_prot = set(in_motif_data.index) & set(out_motif_data.index)
        out_motif_data.drop(common_prot, axis=0, inplace=True)

        in_motif_data.set_index("protein_id", drop=False, inplace=True)
        # print(in_motif_data)
        out_motif_data.set_index("protein_id", drop=False, inplace=True)
        return in_motif_data, out_motif_data


if __name__ == '__main__':
    seq_length = pd.read_csv("/data1/hbs/total_fastr/Total.elm", sep="\t", usecols=[1, 4], encoding="utf-8")
    seq_length.drop_duplicates(inplace=True)
    seq_length["protein_length"] = 0

    def get_len(i):
        i["protein_length"] = len(i["Sequence"])
        return i

    seq_length = seq_length.apply(get_len, axis=1)
    seq_length = seq_length[["Uniprot Accession", "protein_length"]]

    high_freq_prot = []
    a = AnalysisQ9UPN3(high_freq_prot)
    in_motif_data, out_motif_data = a.load_data()

    in_motif_merged = a.merge_overlap_domain(in_motif_data)
    in_motif_merged["motif length"] = in_motif_merged["end"] - in_motif_merged["start"] + 1
    in_motif_merged = pd.DataFrame(in_motif_merged.groupby([in_motif_merged.index])["motif length"].sum())

    in_motif_common = set(in_motif_data.index)&set(in_motif_merged.index)

    in_motif_data = in_motif_data.ix[in_motif_common]
    in_motif_res,out_motif_res = a.counting(in_motif_data, out_motif_data)

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
        i["cancer_type"] = a.cancer[i["protein_id"]]
        return i

    df = df.apply(get_cancer, axis=1)
    df = df[["protein_id","cancer_type", "protein_length","in motif mut", "motif length", "out motif mut", "non_motif length"]]
    df.to_csv("/data1/hbs/all_cancer_analysis/significant_prot_for_IBS/Q9UPN3.txt", sep="\t", index=False, mode="w", encoding="utf-8")