#-*- coding:utf-8 -*-
import pandas as pd
import warnings
import os
import multiprocessing as mp


class IfInDomain:

    def load_data(self, domain_file, motif_file):

        domain_data = pd.read_csv(domain_file, sep="\t", header=None, usecols=[0, 6, 7], encoding="utf-8")
        domain_data.columns = ["protein_id", "domain_start", "domain_end"]
        domain_data["protein_id"] = domain_data.apply(self.create_index, axis=1)
        domain_data.set_index("protein_id", inplace=True)
        domain_data = self.merge_overlap_domain(domain_data)
        domain_data["domain_length"] = domain_data["domain_end"] - domain_data["domain_start"] + 1

        motif_data = pd.read_csv(motif_file, sep="\t", usecols=[0, 1, 2, 3, 5], encoding="utf-8")
        motif_data["sequence"] = motif_data.apply(self.get_seq_length, axis=1)
        motif_data["idx"] = motif_data["protein_id"] + ":" + motif_data["Uniprot Accession"]
        motif_data.set_index("idx", inplace=True)

        data = pd.merge(left=motif_data, right=domain_data, left_index=True, right_index=True)
        data.drop_duplicates(inplace=True)
        
        #interproscan预测某个蛋白序列很可能没有domain，没有domain的蛋白是不会出现在interproscan文件中
        #此时就要找出这些蛋白，定义为 protein_without_domain，并最终归总在 self.out_domain中
        no_domain_protein = motif_data.drop(set(data.index), axis=0)
        self.protein_without_domain = pd.DataFrame(columns=["idx", "protein_id", "Uniprot Accession", "sequence_length", "domain_start", "domain_end", "domain_length", "position", "in_domain"])
        no_domain_protein.apply(self.no_domain_protein, axis=1)
        return data

    def create_index(self, i):
        '''对domain_data的protein id进行处理'''
        info = i["protein_id"].split("@")
        idx = info[0] + ":" + info[1]
        return idx

    def get_seq_length(self, i):
        seq_len = len(i["sequence"])    #获取序列长度，之后的bootstrap检验会用到
        return seq_len

    def merge_overlap_domain(self, domain_data):
        '''将重叠的domain区域合并'''
        new_domain_data = pd.DataFrame(columns=["protein_id", "domain_start", "domain_end"])
        idxs = set(domain_data.index)
        for idx in idxs:
            intervals = domain_data.ix[idx]
            if "Series" in str(type(intervals)):  #Series说明该蛋白对应的domain只有一个
                d = pd.DataFrame([[idx, intervals["domain_start"], intervals["domain_end"]]], columns=["protein_id", "domain_start","domain_end"])
                new_domain_data = pd.concat([new_domain_data, d], axis=0)

            elif "DataFrame" in str(type(intervals)):   #蛋白上有多个domain
                intervals.sort_values("domain_start", inplace=True)  #起始位点从小到大排序
                domain_area = list(zip(list(intervals["domain_start"]), list(intervals["domain_end"])))
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
                d = pd.DataFrame(domain_area, columns=["domain_start", "domain_end"])
                d["protein_id"] = idx
                new_domain_data = pd.concat([new_domain_data, d], axis=0)
        new_domain_data.set_index("protein_id", inplace=True)
        # print(new_domain_data)
        return new_domain_data

    def if_overlapped(self, interval1, interval2):
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
        self.in_domain = pd.DataFrame(columns=["idx", "protein_id", "Uniprot Accession", "sequence_length", "domain_start", "domain_end", "domain_length", "position", "in_domain"])
        self.out_domain = pd.DataFrame(columns=["idx", "protein_id", "Uniprot Accession", "sequence_length", "domain_start", "domain_end", "domain_length", "position", "in_domain"])
        data["in_domain"] = ""
        data.apply(self.if_in_domain, axis=1)
        self.in_domain.set_index("idx", inplace=True)
        self.out_domain.set_index("idx", inplace=True)
        need_drop = set(self.in_domain.index) & set(self.out_domain.index)
        self.out_domain.drop(need_drop, axis=0, inplace=True)

        #将不在domain中的蛋白信息和蛋白上没有domain的数据结合在一起
        self.out_domain = pd.concat([self.out_domain, self.protein_without_domain], axis=0) 
        return None

    def if_in_domain(self, i):
        '''判断突变位点是否在domain区域中'''
        domain_range = [i["domain_start"], i["domain_end"]]
        protein_mod_position = i["position"].split(";")
        idx = i["protein_id"] + ":" + i["Uniprot Accession"]
        for num in range(len(protein_mod_position) - 1):
            pos = int(protein_mod_position[num])
            if pos >= domain_range[0] and pos <= domain_range[1]:
                data = pd.DataFrame([[idx, i["protein_id"], i["Uniprot Accession"], i["sequence"], i["domain_start"], i["domain_end"], i["domain_length"], pos, "yes"]], columns=["idx", "protein_id", "Uniprot Accession", "sequence_length", "domain_start", "domain_end", "domain_length", "position", "in_domain"])
                self.in_domain = pd.concat([self.in_domain, data], axis=0)
            else:
                data = pd.DataFrame([[idx, i["protein_id"], i["Uniprot Accession"], i["sequence"], i["domain_start"], i["domain_end"], i["domain_length"], pos, "no"]], columns=["idx", "protein_id", "Uniprot Accession", "sequence_length", "domain_start", "domain_end", "domain_length", "position", "in_domain"])
                self.out_domain = pd.concat([self.out_domain, data], axis=0)
        return None

    def no_domain_protein(self, i):
        '''对protein上没有domain的蛋白质处理，默认该蛋白上所有修饰位点都是out domain的。
           默认domain起始和终止位点都是0，长度为0。加入到 self.protein_without_domain中'''
        protein_mod_position = i["position"].split(";")
        idx = i["protein_id"] + ":" + i["Uniprot Accession"]
        for num in range(len(protein_mod_position) - 1):
            pos = int(protein_mod_position[num])
            data = pd.DataFrame([[idx, i["protein_id"], i["Uniprot Accession"], i["sequence"], 0, 0, 0, pos, "no"]], columns=["idx", "protein_id", "Uniprot Accession", "sequence_length", "domain_start", "domain_end", "domain_length", "position", "in_domain"])
            self.protein_without_domain = pd.concat([self.protein_without_domain, data], axis=0)

    def to_doc(self, in_domain_file, out_domain_file):
        self.in_domain.to_csv(in_domain_file, sep="\t", index=False, mode="w", encoding="utf-8")
        self.out_domain.to_csv(out_domain_file, sep="\t", index=False, mode="w", encoding="utf-8")
        return None      


if __name__ == "__main__":
    warnings.filterwarnings("ignore")

############################################################################################################
    mod_list = ["Acetylation", "Glycation", "Malonylation", "Methylation",
                         "Succinylation", "Sumoylation", "Ubiquitination"]
    interproscan_dir = "/data1/hbs/processedByInterproscan/"
    total_protein_data = "/data1/hbs/mutation_separate(sig!=yes)/"
    protein_significance = ["_sig", "_insig"]
    file_list = []

    domain_files = os.listdir(interproscan_dir)
    sig_insig_files = os.listdir(total_protein_data)
    for i in range(len(mod_list)):
        mod = mod_list[i]
        for j in protein_significance:
            domain_file = interproscan_dir + mod + j + ".fasta.tsv"
            motif_file = total_protein_data + mod + j + ".txt"
            in_domain_file = os.path.dirname(os.path.realpath(__file__)) + "/" + mod + j + "_in_domain.txt"
            out_domain_file = os.path.dirname(os.path.realpath(__file__)) + "/" + mod + j + "_out_domain.txt"
            file_list.append([domain_file, motif_file, in_domain_file, out_domain_file])

    def execute_class_ifindomain(files):
        '''
        :param files: 是一个列表，里面存放相关的蛋白质信息文件和输出的文件地址
        '''
        domain_file = files[0]
        motif_file = files[1] 
        in_domain_file = files[2] 
        out_domain_file = files[3]
        
        i = IfInDomain()
        data = i.load_data(domain_file, motif_file)
        print(data.head())
        i.form_new_data(data)
        i.to_doc(in_domain_file, out_domain_file)
        return None 
#################################################################################################
#上面这一大串不能封装在类中的原因：
    #pool方法都使用了queue.Queue将task传递给工作进程。multiprocessing必须将数据序列化以在进程间传递。
    #方法只有在模块的顶层时才能被序列化，跟类绑定的方法不能被序列化，就会出现上面的异常。   
#意思就是，在类中的方法无法使线程和进程形成Queue。

    pool = mp.Pool(8)
    # pool.map(execute_class_ifindomain,file_list)  这种和下面的apply_async方式都可以
    for i in file_list:
        pool.apply_async(execute_class_ifindomain, (i,))    #传入参数的形式是tuple
    print("executing……")
    pool.close()    #其他子进程无法加入池中
    pool.join()     #等待所有子线程结束
    print("finish!")

#测试代码，不用理会
    # i = IfInDomain("C:/Users/hbs/Desktop/domain/", "C:/Users/hbs/Desktop/motif/")
    # data = i.load_data("C:/Users/hbs/Desktop/domain/Succinylation_sig.fasta.tsv", "C:/Users/hbs/Desktop/motif/Succinylation_sig.txt")
    # i.form_new_data(data)