#-*- coding:utf-8 -*-
import os
import re
import sys
import getopt
import pandas as pd
import warnings

class Mut_Prot:
    def __init__(self, mut_prot_file, if_sig):
        '''
        :param mut_prot: 突变蛋白相关文件目录
        :param if_sig: yes、no-yes和no三种
        '''
        self.mp = mut_prot_file
        self.if_sig = if_sig

    def load_data(self):
        # data = pd.DataFrame(columns=["mutated_prot", "if_sig"])
        data = pd.read_csv(self.mp, sep="\t",header=None,usecols=[0,1])
        data.columns = ["mutated_prot", "if_sig"]
        data.drop_duplicates(inplace=True)
        return data

    def preprocess(self, data):
        data = data[data["if_sig"] == self.if_sig]
        return data

class Site_Annotated:
    def __init__(self, icgc, tcga, cosmic, mod_type,kind):
        '''
        :param icgc: icgc位点注释文件目录
        :param tcga: tcga位点注释文件目录
        :param cosmic: cosmic位点注释文件目录
        :param mod_type: 突变位点的修饰类型
        '''
        self.icgc = icgc
        self.tcga = tcga
        self.cos = cosmic
        mod_dict = {"ace":["Acetylation_in","Acetylation_out"],"gly":["Glycation_in","Glycation_out"],
                    "mal":["Malonylation_in","Malonylation_out"],"met":["Methylation_in","Methylation_out"],
                    "suc":["Succinylation_in","Succinylation_out"], "sumo":["Sumoylation_in","Sumoylation_out"],
                    "ubi":["Ubiquitination_in","Ubiquitination_out"]}
        self.modify = mod_dict[mod_type]
        self.cancer_type = kind

    def load_data(self, path):
        files = os.listdir(path)
        cols = ["cancer_type","protein_id","Ensembl_Transcript","Ensembl_Transcript_Isoform","protein_length",
                "Acetylation_in", "Glycation_in", "Malonylation_in", "Methylation_in",
                "Succinylation_in", "Sumoylation_in", "Ubiquitination_in", "Acetylation_out",
                "Glycation_out", "Malonylation_out", "Methylation_out", "Succinylation_out",
                "Sumoylation_out",	"Ubiquitination_out"]
        use_cols = [0, 1, 3]
        columns = ["cancer_type", "Canonical","isoform","in","out"]
        for mod in self.modify:
            mod_index = cols.index(mod)
            use_cols.append(mod_index)
        data = pd.DataFrame(columns=columns)
        for file in files:
            file = path + '/' + file
            try:
                d = pd.read_csv(file, sep='\t',header=None,usecols=use_cols)
                d.columns = columns
                data = pd.concat([data, d], axis=0)
            except:
                continue
        return data

    def annotated_data(self):
        paths = [self.icgc, self.tcga, self.cos]
        data = pd.DataFrame(columns=["cancer_type","Canonical", "isoform","in","out"])
        for path in paths:
            d = self.load_data(path)
            data = pd.concat([data, d], axis=0)
        # data.set_index("protein_id",inplace=True)
        data["isoform"].fillna('',inplace=True)
        data["isoform"] = data["Canonical"] + data["isoform"]
        data["isoform"] = data["isoform"].apply(self.standardize_prot_name)
        data = data[data["cancer_type"]==self.cancer_type]
        return data

    def standardize_prot_name(self, i):
        '''有的基因名称类似于 ENST00000262510 -> Q86WI3-1 预处理后的结果为 Q86WI3-1'''
        pattern = re.compile('.*?->.*?([A-Z,a-z,0-9]+-\d{1})', re.S)
        res = re.findall(pattern, i)
        if res!= []:
            return res[0]
        else:
            return i

class Mutation_Sites:
    def __init__(self, mut_prot, if_sig,icgc, tcga, cosmic, mod_type,kind):
        self.mp = mut_prot
        self.if_sig = if_sig
        self.icgc = icgc
        self.tcga = tcga
        self.cos = cosmic
        self.mod = mod_type
        self.cancer_type = kind

    def mutated_proteins(self):
        m = Mut_Prot(self.mp, self.if_sig)
        data = m.load_data()
        data = m.preprocess(data)
        proteins = set(data.mutated_prot)
        return proteins

    def annotated_sites(self):
        a = Site_Annotated(self.icgc, self.tcga, self.cos, self.mod,self.cancer_type)
        data = a.annotated_data()
        data.set_index("isoform",drop=False,inplace=True)
        return data

    def mut_sites(self):
        proteins = self.mutated_proteins()
        annotation = self.annotated_sites()
        annotation.fillna('',inplace=True)
        annotation[["in","out"]] = annotation[["in","out"]].applymap(self.analysis)   #每个蛋白注释信息解读返回list
        common_protein = proteins & set(annotation.isoform)
        if len(common_protein) != 0:
            # data_in = pd.DataFrame(["protein_id","cancer","isoform","seq","pos","from","to","sample_num","direct_mutation"])    #in domain的修饰
            # data_out = pd.DataFrame(["protein_id","cancer","isoform","seq","pos","from","to","sample_num","direct_mutation"])   #out domain的修饰
            entry_in = []
            entry_out = []
            for pro in common_protein:
                try:
                    data = annotation.ix[pro]
                    d_in = data[["cancer_type","Canonical","in"]]
                    if "Series" in str(type(d_in)):
                        if d_in["in"]!="":
                            l = d_in["in"].split(';')
                            for i in l:
                                i = i.split('\t')
                                i.insert(0,pro)
                                i.insert(1,d_in["Canonical"])
                                i.insert(2,d_in["cancer_type"])
                                if i[3]==i[5]:
                                    i.append ("yes")
                                else:
                                    i.append("no")
                                entry_in.append(i)
                    elif "DataFrame" in str(type(d_in)):
                        info = list(d_in["in"])
                        info3 = list(d_in["Canonical"])
                        info2 = list(d_in["cancer_type"])
                        for i in range(len(info)):
                            if info[i] != '':
                                ls = info[i].split(';')
                                for l in ls:
                                    l = l.strip().split('\t')
                                    l.insert(0,pro)
                                    l.insert(1,info3[i])
                                    l.insert(2,info2[i])
                                    if l[3] == l[5]:
                                        l.append("yes")
                                    else:
                                        l.append("no")
                                    entry_in.append(l)
                except KeyError:
                    continue

            for pro in proteins:
                try:
                    data = annotation.ix[pro]
                    d_in = data[["cancer_type","Canonical","out"]]
                    if "Series" in str(type(d_in)):
                        if d_in["out"]!="":
                            l = d_in["out"].split(';')
                            for i in l:
                                i = i.split('\t')
                                i.insert(0,pro)
                                i.insert(1,d_in["Canonical"])
                                i.insert(2,d_in["cancer_type"])
                                if i[3]==i[5]:
                                    i.append ("yes")
                                else:
                                    i.append("no")
                                entry_out.append(i)
                    elif "DataFrame" in str(type(d_in)):
                        info = list(d_in["out"])
                        info3 = list(d_in["Canonical"])
                        info2 = list(d_in["cancer_type"])
                        for i in range(len(info)):
                            if info[i] != '':
                                ls = info[i].split(';')
                                for l in ls:
                                    l = l.strip().split('\t')
                                    l.insert(0,pro)
                                    l.insert(1,info3[i])
                                    l.insert(2,info2[i])
                                    if l[3] == l[5]:
                                        l.append("yes")
                                    else:
                                        l.append("no")
                                    entry_out.append(l)
                except KeyError:
                    continue
            columns=["protein_id","Canonical","cancer_type","K position","seq","mutated position","from","to","sample_num","direct mutation"]
            in_domain_data = pd.DataFrame(entry_in)
            in_domain_data.columns = columns
            out_domain_data = pd.DataFrame(entry_out)
            out_domain_data.columns = columns
            in_domain_file = os.path.dirname(os.path.realpath(__file__)) + '/' +self.mod + "_in_domain.txt"
            out_domain_file = os.path.dirname(os.path.realpath(__file__)) + '/'+self.mod + "_out_domain.txt"
            in_domain_data.to_csv(in_domain_file,sep='\t',index=False,mode='w',encoding="utf-8")
            out_domain_data.to_csv(out_domain_file,sep='\t',index=False,mode='w',encoding="utf-8")
        else:
            in_domain_file = self.mod + "_in_domain.txt"
            out_domain_file = self.mod + "_out_domain.txt"
            with open(in_domain_file,'w') as f:
                f.write("")
            f.close()
            with open(out_domain_file,'w') as f:
                f.write("")
            f.close()

    def analysis(self,i):
        '''将关键修饰信息从文件中提取'''
        data = []
        if i != '':
            info = i.strip().split("|")     #分割得到这个蛋白中有几个突变的序列
            pattern = re.compile('(\d+),\d+,\d,\d,(.*?),(.*)',re.S)
            for item in info:
                res = re.findall(pattern,item)
                pos_and_sample = res[0][2]
                pos_and_sample = pos_and_sample.split(';')
                for j in pos_and_sample:
                    pos = j.split('-')[0]
                    sample_num = str(len(j.split(',')))
                    try:
                        mut = pos[1:-1]
                        f = pos[0]
                        t = pos[-1]
                        result = [res[0][0],res[0][1],mut,f,t,sample_num]
                        result = "\t".join(result)
                        data.append(result)
                    except:
                        print(j)
            data = ';'.join(data)
        else:
            data = ''
        return data


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    try:
        opts,args = getopt.getopt(sys.argv[1::],"m:s:i:t:c:k:",["mod=","sig=","icgc=","tcga=","cosmic=","kind="])
    except:
        exit(0);
    mod = 'ace'
    sig = 'yes'
    icgc = '/data1/data/zengyanru/LysineTCGA/ICGC/cancer_annovar_correct/info_v3'
    tcga = '/data1/data/zengyanru/LysineTCGA/find_mutation_in_motif/all_info_V3'
    cosmic = '/data1/data/zengyanru/LysineTCGA/COSMIC/correct_protein'
    mutation_file = '/data1/hbs/UCECresult/UCEC'+mod+'_all.txt2.txt'
    kind = 'UCEC'

    for opt,arg in opts:
        if opt in ('-m','--mod'):
            mod = arg
        elif opt in ('-s','--sig'):
            sig = args
        elif opt in ('-i','--icgc'):
            icgc = arg
        elif opt in ('-t','--tcga'):
            tcga = arg
        elif opt in ('-c','--cosmic'):
            cosmic = arg
        elif opt in ('-k','--kind'):
            kind = arg
        else:
            print("invalid input")
            exit(1)

    m = Mutation_Sites(mutation_file,sig,icgc,tcga,cosmic,mod,kind)
    p = m.mut_sites()