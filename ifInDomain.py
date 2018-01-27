#-*- coding:utf-8 -*-
'''
Description:
    目前所有的motif的样本突变已被分为central、surround和out三类
    需要进一步将他们分为 in domain和out domain两类。domain的文件
    是interproscan处理过后生成在 /data1/hbs/all_cancer_analysis/all_domain_classify_by_modify
    目录中按照不同的修饰将所有蛋白的domain写入到不同文件中。
=======================================================
@author:hbs
@date:2018-1-20
@version:1.2
'''
import pandas as pd
import warnings
import os
import multiprocessing as mp


class IfInDomain:

    def __init__(self, sequence_data):
        self.sequences = sequence_data

    def load_data(self, domain_file, in_motif_file, out_motif_file):
        '''将记录某一种修饰的domain文件中interproscan产生的数据读取出来，读取protein_id(PLMD),uniprot id、modification,seq length,start,end'''
        domain_data = pd.read_csv(domain_file, sep="\t", encoding="utf-8")
        domain_data.set_index("Uniprot Accession", drop=False, inplace=True)
        # 将存在overlap的domain合并
        domain_data = self.merge_overlap_domain(domain_data)
        #计算domain长度，此时dataframe中会有两个长度，一个是整个蛋白序列的长度，另一个是单独的domain区域的长度
        domain_data["domain length"] = domain_data["end"] - domain_data["start"] + 1

        #motif数据的导入，从依据modification分类后的文件中导入 protein_id(uniprot isoform),Canonical,cancer_type,mutated position,from,to,sample_num（目前还没计数）,direct mutation
        in_motif_data = pd.read_csv(in_motif_file, sep="\t", usecols=["protein_id","Canonical","cancer_type",	"mutated position","sample_num"], encoding="utf-8")
        in_motif_data.set_index("Canonical", drop=False, inplace=True)

        out_motif_data = pd.read_csv(out_motif_file, sep="\t", usecols=["protein_id", "Canonical", "cancer_type", "mutated position","sample_num"], encoding="utf-8")
        out_motif_data.set_index("Canonical", drop=False, inplace=True)

        motif_data = pd.concat([in_motif_data, out_motif_data], axis=0)
        motif_data.drop_duplicates(inplace=True)

        #将domain信息与motif的信息进行合并（取交集），二者索引相同
        data = pd.merge(left=motif_data, right=domain_data, left_index=True, right_index=True)
        data.drop_duplicates(inplace=True)

        #获取到序列信息，从文件中读取的内容包括PLMD ID、Uniprot ID和sequence序列

        #interproscan预测某个蛋白序列很可能没有domain，没有domain的蛋白是不会出现在interproscan文件中
        #此时就要找出这些蛋白，定义为 protein_without_domain，并最终归总在 self.out_domain中
        no_domain_protein = motif_data.drop(set(data.index), axis=0)
        self.protein_without_domain = pd.DataFrame(columns=["idx", "Uniprot Accession", "protein length", "domain start", "domain end", "non_domain length", "position", "sample_num", "sample_count", "in_domain"])
        no_domain_protein.apply(self.no_domain_protein, axis=1)

        self.protein_without_domain.set_index("idx", drop=False, inplace=True)
        self.protein_without_domain = pd.merge(left=self.protein_without_domain, right=self.sequences, left_index=True, right_index=True)
        #获取到没有domain的蛋白的序列长度
        self.protein_without_domain = self.protein_without_domain.apply(self.get_seq_length, axis=1)
        self.protein_without_domain.drop_duplicates(inplace=True)
        # self.protein_without_domain["protein_id"] = ""
        # self.protein_without_domain["position"] = self.protein_without_domain["position"].map(str)
        # self.protein_without_domain["idx"] = self.protein_without_domain["idx"] + ":" + self.protein_without_domain["position"]
        # self.protein_without_domain.set_index("idx", inplace=True)
        return data

    # def create_index(self, i):
    #     '''对domain_data的protein id进行处理，目前的protein ID是详细到亚形isoform的，需要处理
    #         处理成Caonical形式，因为domain中是Canonical作为索引'''
    #     info = i["protein_id"].split("-")
    #     i["Canonical"] = info[0]

    def get_seq_length(self, i):
        '''self.protein_without_domain中的蛋白一开始没有从domain文件中获取到蛋白质长度，所以需要对照序列文件获取蛋白长度
        ，同时从sequence文件中获取到PLMD ID号'''

        seq_len = len(i["Sequence"])    #获取序列长度，之后的bootstrap检验会用到
        i["protein length"] = seq_len
        return i

    def merge_overlap_domain(self, domain_data):
        '''将重叠的domain区域合并'''
        new_domain_data = pd.DataFrame(columns=["protein_id", "protein length", "start", "end"])
        idxs = set(domain_data.index)
        for idx in idxs:
            intervals = domain_data.ix[idx]
            protein_length = intervals["protein length"]
            if "Series" in str(type(intervals)):  #Series说明该蛋白对应的domain只有一个
                d = pd.DataFrame([[idx, protein_length, intervals["start"], intervals["end"]]], columns=["protein_id", "protein length", "start","end"])
                new_domain_data = pd.concat([new_domain_data, d], axis=0)

            elif "DataFrame" in str(type(intervals)):   #蛋白上有多个domain
                intervals.sort_values("start", inplace=True)  #起始位点从小到大排序
                domain_area = list(zip(list(intervals["start"]), list(intervals["end"])))
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
                d["protein length"] = list(protein_length)[0]
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

    def form_new_data(self, data):
        '''根据载入的数据生成在domain和不在domain的位点的dataframe，
            同一个位点因为选取的domain不同 in_domain的判断同时出现yes和no时，
            优先保留yes的情况，删除out_domain中的数据'''
        self.in_domain = pd.DataFrame(columns=["Uniprot Accession", "protein length", "domain start", "domain end", "domain length", "position", "sample_num","sample_count", "in_domain"])
        self.out_domain = pd.DataFrame(columns=["Uniprot Accession", "protein length", "domain start", "domain end", "non_domain length", "position", "sample_num","sample_count", "in_domain"])
        data["in_domain"] = ""
        #判断是否是在domain区域中
        data.apply(self.if_in_domain, axis=1)

        self.in_domain.drop_duplicates(subset=["Uniprot Accession","sample_num", "position"], inplace=True)
        self.out_domain.drop_duplicates(subset=["Uniprot Accession","sample_num", "position"], inplace=True)

        self.in_domain["position"] = self.in_domain["position"].map(str)
        #首先统计同一个区域中突变样本的数目
        self.in_domain = self.in_domain.groupby(["Uniprot Accession", "protein length","domain start","domain length","domain end", "in_domain"],as_index=False)["sample_count"].sum()
        #将所有的domain区域加和
        self.in_domain = self.in_domain.groupby(["Uniprot Accession", "protein length", "in_domain"], as_index=False)["domain length", "sample_count"].sum()
        self.in_domain.set_index("Uniprot Accession", drop=False, inplace=True)

        #将不在domain中的蛋白信息和蛋白上没有domain的数据结合在一起
        self.out_domain = self.out_domain.groupby(["Uniprot Accession","protein length","in_domain","non_domain length"],as_index=False)["sample_count"].sum()
        self.protein_without_domain = self.protein_without_domain.groupby(["Uniprot Accession","protein length","in_domain","non_domain length"], as_index=False)["sample_count"].sum()
        self.out_domain = pd.concat([self.out_domain, self.protein_without_domain], axis=0)
        self.out_domain.set_index("Uniprot Accession", drop=False, inplace=True)

        #一个位点由于选取的序列不同被同时分在了 in domain和 out domain中时，此时舍去out domain中的条目，避免之后统计过程中的重复计数
        self.common_prot = set(self.in_domain.index) & set(self.out_domain.index)
        self.out_domain = self.out_domain.apply(self.calc_non_domain, axis=1)
        return None

    def if_in_domain(self, i):
        '''判断突变位点是否在domain区域中'''
        domain_range = [i["start"], i["end"]]
        domain_length = i["end"] - i["start"] + 1
        pos = int(i["mutated position"])
        # for num in range(len(protein_mod_position) - 1):
        #     pos = int(protein_mod_position[num])
        if pos >= domain_range[0] and pos <= domain_range[1]:
            data = pd.DataFrame([[i["protein_id"], i["protein length"], i["start"], i["end"], domain_length, pos, i["sample_num"], 1, "yes"]])
            data.columns=["Uniprot Accession", "protein length", "domain start", "domain end", "domain length", "position", "sample_num", "sample_count", "in_domain"]
            self.in_domain = pd.concat([self.in_domain, data], axis=0)
        elif pos < domain_range[0] or pos > domain_range[1]:
            data = pd.DataFrame([[i["protein_id"], i["protein length"], 0, 0, 0, pos, i["sample_num"], 1, "no"]])
            data.columns=["Uniprot Accession", "protein length", "domain start", "domain end", "non_domain length", "position", "sample_num","sample_count", "in_domain"]
            self.out_domain = pd.concat([self.out_domain, data], axis=0)
        return None

    def no_domain_protein(self, i):
        '''对protein上没有domain的蛋白质处理，默认该蛋白上所有修饰位点都是out domain的。
           默认domain起始和终止位点都是0，长度为0。加入到 self.protein_without_domain中'''
        pos = i["mutated position"]
        # idx = i.index + ":" + i["protein_id"]
        data = pd.DataFrame([[i["Canonical"], i["protein_id"], 0, 0, 0, 0, pos, i["sample_num"], 1, "no"]])
        data.columns = ["idx", "Uniprot Accession", "protein length", "domain start", "domain end", "non_domain length", "position", "sample_num", "sample_count", "in_domain"]
        self.protein_without_domain = pd.concat([self.protein_without_domain, data], axis=0)

    def calc_non_domain(self, i):
        idx = i["Uniprot Accession"]
        if idx in self.common_prot:
            non_domain_length = i["protein length"] - self.in_domain.ix[idx]["domain length"]
            i["non_domain length"] = non_domain_length
        else:
            i["non_domain length"] = i["protein length"]
        return i

    def to_doc(self, in_domain_file, in_domain_calc, out_domain_file, out_domain_calc):
        '''
        :param in_domain_file: in domain的突变的详细信息写入的文件
        :param in_domain_calc: in domain的统计计数信息写入的文件，最终用于bootstrap检验
        :param out_domain_file: out domain的突变的详细信息写入的文件
        :param out_domain_calc: out domain的统计计数信息写入的文件，最终用于bootstrap检验
        :return:
        '''
        self.in_domain.to_csv(in_domain_file, sep="\t", index=False, mode="w", encoding="utf-8")
        in_domain_res = self.in_domain.groupby(["Uniprot Accession", "domain length"], as_index=False)["sample_count"].sum()
        in_domain_res["mut_per_k_aa"] = in_domain_res["sample_count"]/in_domain_res["domain length"] *1000
        in_domain_res.to_csv(in_domain_calc, sep="\t", index=False, mode="w", encoding="utf-8")
        # print(in_domain_res)
        self.out_domain.to_csv(out_domain_file, sep="\t", index=False, mode="w", encoding="utf-8")
        out_domain_res = self.out_domain.groupby(["Uniprot Accession", "non_domain length"],as_index=False)["sample_count"].sum()
        out_domain_res["mut_per_k_aa"] = out_domain_res["sample_count"] / out_domain_res["non_domain length"] *1000
        out_domain_res.to_csv(out_domain_calc, sep="\t", index=False, mode="w", encoding="utf-8")
        # print(out_domain_res)   
        return None      


if __name__ == "__main__":
    warnings.filterwarnings("ignore")

############################################################################################################
    def execute_class_ifindomain(domain_file,in_motif_file,out_motif_file,sequence_data,in_domain_file,out_domain_file,in_domain_calc,out_domain_calc):
        '''
        :param files: 是一个列表，里面存放相关的蛋白质信息文件和输出的文件地址
        '''
        i = IfInDomain(sequence_data)
        data = i.load_data(domain_file, in_motif_file, out_motif_file)
        i.form_new_data(data)
        i.to_doc(in_domain_file,in_domain_calc, out_domain_file,out_domain_calc)
        return None
    #33种癌症
    cancer_kind = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP',
                   'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM',
                   'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']

    #7种修饰的全称
    mod_list = {"ace":"Acetylation", "gly":"Glycation", "mal":"Malonylation", "met":"Methylation",
                         "suc":"Succinylation", "sum":"Sumoylation", "ubi":"Ubiquitination"}
    #7种修饰的简写
    mod_abbr_list = ["ace", "gly", "mal", "met", "suc", "sum", "ubi"]

    #interproscan处理后按照不同修饰类型存放的蛋白的domain的文件目录
    interproscan_dir = "/data1/hbs/all_cancer_analysis/all_domain_classify_by_modify/"

    #motif的计数文件
    total_protein_data = "/data1/hbs/all_cancer_analysis/all_cancer_sig_prot_motif_modify/"

    #三种不同位置的突变
    mod_in_motif = ["_in_motif.txt", "_out_motif.txt"]
    file_list = []

    sequence_file = "/data1/hbs/total_fastr/Total.elm"
    sequence_data = pd.read_csv(sequence_file, sep="\t", usecols=[0, 1, 4], encoding="utf-8")
    sequence_data.set_index("Uniprot Accession", inplace=True)

    result_dir = "/data1/hbs/all_cancer_analysis/all_cancer_modification_in_domain_or_not/calc_result/"

    pool = mp.Pool(8)

    for c in cancer_kind:
        #创建每一种癌症的文件目录
        os.system("mkdir "+result_dir+c)
        cancer_dir = total_protein_data + c
        for m in mod_abbr_list:
            #每种癌症都分为七种修饰
            #经过interproscan处理后的各类修饰蛋白的domain文件
            domain_file = interproscan_dir + mod_list[m] +"_simplified.txt"
            #某类癌症在某一位置的某种修饰的motif文件
            in_motif_file = total_protein_data + c+ "/" + m + mod_in_motif[0]
            out_motif_file = total_protein_data + c+ "/" + m + mod_in_motif[1]

            #详尽的in domain和 out domain的记录文件
            in_domain_file = result_dir+ c +"/" + m+ "_in_domain.txt" 
            out_domain_file = result_dir+ c +"/" + m+ "_out_domain.txt" 
            #统计某种癌症、某一修饰是否在domain中的样本计数文件
            in_domain_calc = result_dir+ c +"/" + m+"_in_domain_calc.txt"
            out_domain_calc = result_dir+ c +"/" + m+"_out_domain_calc.txt"
            # execute_class_ifindomain(domain_file, in_motif_file, out_motif_file, sequence_data, in_domain_file, out_domain_file, in_domain_calc,
            # out_domain_calc)
            pool.apply_async(execute_class_ifindomain,(domain_file,in_motif_file,out_motif_file,sequence_data,in_domain_file,out_domain_file,in_domain_calc,out_domain_calc))
    pool.close()
    pool.join()

#################################################################################################
#上面这一大串不能封装在类中的原因：
    #pool方法都使用了queue.Queue将task传递给工作进程。multiprocessing必须将数据序列化以在进程间传递。
    #方法只有在模块的顶层时才能被序列化，跟类绑定的方法不能被序列化，就会出现上面的异常。   
#意思就是，在类中的方法无法使线程和进程形成Queue。

    # pool = mp.Pool(8)
    # # pool.map(execute_class_ifindomain,file_list)  这种和下面的apply_async方式都可以
    # for i in file_list:
    #     print(i)
    #     # pool.apply_async(execute_class_ifindomain, (i,))    #传入参数的形式是tuple
    # print("executing...")
    # pool.close()    #其他子进程无法加入池中
    # pool.join()     #等待所有子线程结束
    # print("finish!")

#测试代码，不用理会
    # "C:/Users/hbs/Desktop/domain/", "C:/Users/hbs/Desktop/motif/"
    # sequence_data = pd.read_csv("C:/Users/hbs/Desktop/Total.elm", sep="\t", usecols=[0, 1, 4], encoding="utf-8")
    # sequence_data.set_index("Uniprot Accession", inplace=True)
    # i = IfInDomain(sequence_data)
    # data = i.load_data("C:/Users/hbs/Desktop/in_out_domain/Glycation_simplified.txt", "C:/Users/hbs/Desktop/in_out_domain/gly_in_motif.txt", "C:/Users/hbs/Desktop/in_out_domain/gly_out_motif.txt")
    #
    # i.form_new_data(data)
    # i.to_doc("in.txt","in_calc.txt", "out.txt", "out_calc.txt")
