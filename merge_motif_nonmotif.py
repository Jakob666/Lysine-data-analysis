#-*- coding:utf-8 -*-
'''
Description:
    用于将蛋白的motif区域与非motif区域分别进行整合，存入文件，
    用于之后的保守型打分
=========================
@author: hbs
@date: 2018-1-16
@version: 1.0
'''
import warnings
import pandas as pd
import os
import re


class Merger:
    def __init__(self, icgc_file, tcga_file, cosmic_file, seq_file, motif_file, non_motif_file):
        '''
        :param icgc_file: 某种癌症的icgc文件
        :param tcga_file: 某种癌症的tcga文件
        :param cosmic_file: 某种癌症的cosmic文件
        :param seq_file: 存储所有蛋白序列的文件
        :param motif_file: 用于存储motif区域序列的文件
        :param non_motif_file: 用于存储非motif区域序列的文件
        '''
        self.icgc = icgc_file
        self.tcga = tcga_file
        self.cos = cosmic_file
        self.mod_list = ["Acetylation", "Glycation", "Malonylation", "Methylation",
                    "Succinylation", "Sumoylation", "Ubiquitination"]
        self.seq_file = seq_file
        self.pattern = re.compile(".*?->(.*)", re.S)
        self.motif = pd.DataFrame(columns=["protein_id", "isoform", "protein_length","central", "start", "end", "Sequence"])
        self.non_motif = pd.DataFrame(columns=["protein_id", "isoform", "start", "end"])
        self.motif_file = motif_file
        self.non_motif_file = non_motif_file

    def load_motif_data(self, file):
        try:
            data = pd.read_csv(file, sep="\t", usecols=[1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18], encoding="utf-8")
            data.columns = ["protein_id", "Ensembl_Transcript_Isoform", "protein_length", "Acetylation_in", "Glycation_in", "Malonylation_in",
                        "Methylation_in", "Succinylation_in", "Sumoylation_in", "Ubiquitination_in", "Acetylation_out",
                        "Glycation_out", "Malonylation_out", "Methylation_out", "Succinylation_out", "Sumoylation_out",
                        "Ubiquitination_out"]
            data.fillna("", inplace=True)
            data = data.apply(self.extract_isoform, axis=1)
        except:
            data = pd.DataFrame(columns=["protein_id", "Ensembl_Transcript_Isoform", "protein_length", "Acetylation_in",
                        "Glycation_in", "Malonylation_in","Methylation_in", "Succinylation_in", "Sumoylation_in",
                        "Ubiquitination_in", "Acetylation_out", "Glycation_out", "Malonylation_out", "Methylation_out",
                        "Succinylation_out", "Sumoylation_out", "Ubiquitination_out"])
        return data

    def extract_isoform(self, i):
        '''读取到的数据中蛋白质很可能是亚形中存在修饰而经典类中不存在修饰，
            如果不是亚形被修饰则注释为空，如果是亚形的修饰注释的形式为：
                ENST00000357254 -> P17027-1
            所要提取的是后半部分，即 P17027-1。
        '''
        try:
            protein_id = i["protein_id"]
            isoform = i["Ensembl_Transcript_Isoform"]
            if isoform != "":
                res = re.findall(self.pattern, isoform)[0]
                res = res.strip()
                i["Ensembl_Transcript_Isoform"] = res
            else:
                i["Ensembl_Transcript_Isoform"] = protein_id
            return i
        except:
            print(i)

    def load_sequence_data(self):
        '''将序列的信息从文件中提取出来'''
        data = pd.read_csv(self.seq_file, sep="\t", usecols=[1, 4], encoding="utf-8")
        data.set_index("Uniprot Accession", inplace=True)
        data.drop_duplicates(inplace=True)
        return data

    def merge_data(self):
        '''该方法用于把motif的信息和序列信息进行整合'''
        #将三个motif文件进行读取，读取后暂时先不进行整合
        icgc = self.load_motif_data(self.icgc)
        tcga = self.load_motif_data(self.tcga)
        cosmic = self.load_motif_data(self.cos)
        #合并后去重
        total_data = pd.concat([icgc, tcga, cosmic], axis=0)
        total_data.drop_duplicates(inplace=True)

        seq = self.load_sequence_data()
        #生成 self.motif
        total_data.apply(self.extract_motif_seq, axis=1)

        self.motif.drop_duplicates(inplace=True)
        self.motif.set_index("protein_id", inplace=True)
        self.motif.astype({"start":"int", "end":"int"})
        self.motif = self.merge_motif_seq()
        self.motif = pd.merge(left=self.motif, right=seq, left_index=True, right_index=True)
        self.motif = self.motif.apply(self.pick_motif_seq, axis=1)

        #根据蛋白长度、motif区域得到non_motif区域的范围，结果为 self.non_motif
        self.merge_non_motif_seq(self.motif)
        self.non_motif.set_index("protein_id", drop=False, inplace=True)
        self.non_motif = pd.merge(left=self.non_motif, right=seq, left_index=True, right_index=True)
        self.non_motif = self.non_motif.apply(self.pick_motif_seq, axis=1)

    def merge_motif_seq(self):
        '''检查各个蛋白motif的区域是否有overlap，有overlap则将两个motif合并
            避免最终将所有motif合并的时候合并出错。'''
        new_interval = []
        idxs = set(self.motif.index)
        for idx in idxs:
            data = self.motif.ix[idx]
            if "DataFrame" in str(type(data)):
                data.sort_values("start", inplace=True)
                need_drop = []
                intervals = list(zip(data["start"], data["end"]))
                for i in range(len(intervals)-1):
                    interval1 = intervals[i]
                    interval2 = intervals[i+1]
                    overlapped = self.if_overlapped(interval1, interval2)
                    if overlapped:
                        intervals[i+1] = [min([interval1[0], interval2[0]]), max([interval1[1], interval2[1]])]
                        need_drop.append(intervals[i])
                for j in need_drop:
                    intervals.remove(j)
                annotate = [idx, list(data["isoform"])[0], list(data["protein_length"])[0]]
            elif "Series" in str(type(data)):
                intervals = [[data["start"], data["end"]]]
                annotate = [idx, data["isoform"], data["protein_length"]]
            intervals = list(map(list, intervals))
            intervals = [annotate + v for v in intervals]
            new_interval += intervals
        motif_seq = pd.DataFrame(new_interval, columns=["protein_id", "isoform", "protein_length", "start", "end"])
        motif_seq.set_index("protein_id", drop=False, inplace=True)
        return motif_seq


    def extract_motif_seq(self, i):
        '''将各个序列中的motif的起始位点终止位点和对应的序列抽提出来。'''
        motifs = []
        location = ["_in", "_out"]
        for mod in self.mod_list:
            for l in location:
                motif = mod + l
                info = i[motif]
                if info != "":
                    seqs = info.strip().split("|")
                    for seq in seqs:
                        sites = seq.split(",")
                        central = sites[0] #中心位点
                        aa_seq = sites[4]  #序列
                        left_flank = sites[2] #左侧aa数目
                        right_flank = sites[3]#右侧aa数目
                        start = int(central) - int(left_flank) -1 #再减一是将其变为索引值
                        end = int(central) + int(right_flank)   #end不减一是便于之后选取序列
                        useful = [central, start, end, aa_seq]
                        motifs.append(useful)
        df = pd.DataFrame(motifs,columns=["central", "start", "end", "Sequence"])
        df["protein_id"] = i["protein_id"]
        df["isoform"] = i["Ensembl_Transcript_Isoform"]
        df["protein_length"] = i["protein_length"]
        self.motif = pd.concat([self.motif, df], axis=0)
        return i

    def create_index(self, i):
        i["protein_id"] = i["isoform"].split("-")[0]
        return i

    def if_overlapped(self, interval1, interval2):
        overlapped = False
        try:
            if int(max((interval1[0], interval2[0]))) <= int(min((interval1[1]), int(interval2[1]))):
                overlapped = True

            elif not (int(interval1[0]) > int(interval2[1]) or int(interval1[1]) < int(interval2[0])):
                overlapped = True
            return overlapped
        except:
            print("i1",interval1)
            print("i2",interval2)

    def merge_non_motif_seq(self, motif_seq):
        '''根据merge_motif_seq方法得到的motif_seq处理得到non_motif_seq，
           需要注意的是 motif_seq中start是索引值，end是索引值+1，之前之所以
           这样是便于截取序列'''
        non_motif_seq = pd.DataFrame(columns=["protein_id", "isoform", "start", "end"])
        idxs = set(motif_seq.index)

        for idx in idxs:
            non_motif_area = []
            data = motif_seq.ix[idx]
            if "DataFrame" in str(type(data)):
                protein_length = list(data["protein_length"])[0]
                protein_id = list(data["protein_id"])[0]
                isoform = list(data["protein_id"])[0]
                motif_area = list(zip(data["start"], data["end"]))
                for i in range(len(motif_area)):
                    if i == 0 and motif_area[i][0] > 0:
                        non_motif_area.append([protein_id, isoform, 0, motif_area[i][0]])
                    elif i == len(motif_area)-1 and motif_area[i][1] < protein_length-1:
                        non_motif_area.append([protein_id, isoform, motif_area[i][1], protein_length])
                    else:
                        non_motif_area.append([protein_id, isoform, motif_area[i-1][1], motif_area[i][0]])

            elif "Series" in str(type(data)):
                protein_length = data["protein_length"]
                motif_area = [data["start"], data["end"]]
                if motif_area[0] > 0:
                    non_motif_area.append([data["protein_id"], data["isoform"], 0, motif_area[0]])
                if motif_area[1] != protein_length-1:
                    non_motif_area.append([data["protein_id"], data["isoform"], motif_area[1], protein_length])
            non_motif = pd.DataFrame(non_motif_area, columns=["protein_id", "isoform", "start", "end"])
            self.non_motif = pd.concat([self.non_motif, non_motif], axis=0)

    def pick_motif_seq(self, i):
        i["Sequence"] = i["Sequence"][i["start"]:i["end"]]
        return i

    def connect_all_sequence(self):
        '''将self.motif和 self.non_motif的序列进行如下操作，
            同一个蛋白的domain序列连在一起， 同一蛋白的非domain序列连在一起
        '''
        self.motif = self.motif.groupby(["protein_id", "isoform"], as_index=False)["Sequence"].sum()
        self.non_motif = self.non_motif.groupby(["protein_id", "isoform"], as_index=False)["Sequence"].sum()

    def to_doc(self):
        self.motif.to_csv(self.motif_file, sep="\t", index=False, mode="w", encoding="utf-8")
        self.non_motif.to_csv(self.non_motif_file, sep="\t", index=False, mode="w", encoding="utf-8")


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    tcga = "/data1/data/zengyanru/LysineTCGA/find_mutation_in_motif/all_info_V3/UCEC_all_info_v4.txt"
    icgc = "/data1/data/zengyanru/LysineTCGA/ICGC/cancer_annovar_correct/info_v3/UCEC-US_annotated__all_info_v2.txt"
    cosmic = ""
    m = Merger(icgc, tcga, cosmic, seq, "/data1/hbs/conservation_analysis/total_motif.txt", "/data1/hbs/conservation_analysis/total_non_motif.txt")
    m.merge_data()
    m.connect_all_sequence()
    m.to_doc()

    # icgc = "C:/Users/hbs/Desktop/ucec/UCEC_all_info_v4.txt"
    # icgc = "C:/Users/hbs/Desktop/motif_for_test.txt"
    # tcga = "C:/Users/hbs/Desktop/ucec/UCEC-US_annotated__all_info_v2.txt"
    # cosmic = ""
    # seq = "C:/Users/hbs/Desktop/Total.elm"

    # m = Merger(icgc, tcga, cosmic, seq, "m.txt", "nm.txt")
    # m.merge_data()
    # m.connect_all_sequence()
    # m.to_doc()