#-*- coding:utf-8 -*-
import os
import re
import sys
import getopt
import pandas as pd
import warnings
import multiprocessing as mp

class Mut_Prot:
    def __init__(self, mut_prot_file, if_sig):
        '''
        :param mut_prot: 突变蛋白相关文件目录
        :param if_sig: yes、no-yes和no三种，目前文件中使用p值来定义显著与非显著，p<0.05被认为是显著
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
        data = data[data["if_sig"] < self.if_sig]
        return data

class Site_Annotated:
    def __init__(self, icgc, tcga, cosmic, mod_type,kind):
        '''
        每种注释文件的区别仅仅是 数据来源不同（有的来自icgc，有的是tcga和cosmic），但是列都是一样的，
        cancer_type	protein_id	Ensembl_Transcript	Ensembl_Transcript_Isoform	protein_length  七种modify_in  七种modify_out
        由于文件都不算大，所以每次读取的都是所有癌症的数据，根据传入的 kind参数从中挑选本次所要提取的癌症类型。
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
        '''读取某一个文件目录下（icgc、tcga和cosmic的三个目录之一）全部文件的数据，'''
        files = os.listdir(path)
        cols = ["cancer_type","protein_id","Ensembl_Transcript","Ensembl_Transcript_Isoform","protein_length",
                "Acetylation_in", "Glycation_in", "Malonylation_in", "Methylation_in",
                "Succinylation_in", "Sumoylation_in", "Ubiquitination_in", "Acetylation_out",
                "Glycation_out", "Malonylation_out", "Methylation_out", "Succinylation_out",
                "Sumoylation_out",	"Ubiquitination_out"]
        use_cols = [0, 1, 3]
        columns = ["cancer_type", "Canonical","isoform","in","out"]

        for mod in self.modify: #self.modify是 [(mode_type)_in, (mode_type)_out] ，其中mode_type是传入的修饰类型
            mod_index = cols.index(mod)
            use_cols.append(mod_index)  #指定修饰在cols中对应的列，也是在文件中对应的列加入到 use_cols中，方便之后读取
        data = pd.DataFrame(columns=columns)
        for file in files:
            file = path + '/' + file
            #这里的try是为了防止读到不对的文件使得程序无法进行
            try:
                d = pd.read_csv(file, sep='\t',header=None,usecols=use_cols)
                d.columns = columns
                data = pd.concat([data, d], axis=0)
            except:
                continue
        return data

    def annotated_data(self):
        paths = [self.icgc, self.tcga, self.cos]    #读取icgc、tcga和cosmic的所有文件，提取某一修饰这5种信息
        data = pd.DataFrame(columns=["cancer_type","Canonical", "isoform","in","out"])
        for path in paths:
            d = self.load_data(path)
            data = pd.concat([data, d], axis=0)
        # data.set_index("protein_id",inplace=True)
        data["isoform"].fillna('',inplace=True)
        data["isoform"] = data["Canonical"] + data["isoform"]
        #如果有亚形，下面一行代码能够获取到亚形的确切名称
        data["isoform"] = data["isoform"].apply(self.standardize_prot_name)
        #在全部癌症中选取本次要研究的癌症
        data = data[data["cancer_type"]==self.cancer_type]
        # data.insert()
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
        '''data是Mut_prot类产生的 sig=self.if_sig的dataframe，
            用这个dataframe得到所有的相应的显著性的蛋白名称，用集合形式返回。'''
        m = Mut_Prot(self.mp, self.if_sig)
        data = m.load_data()
        data = m.preprocess(data)
        proteins = set(data.mutated_prot)
        return proteins

    def annotated_sites(self):
        '''创建一个 Site_annotated对象，读取到某一癌症，某种修饰的 icgc、tcga和cosmic的值
            返回如下信息：cancer_type,Canonical, isoform,in,out的dataframe'''
        a = Site_Annotated(self.icgc, self.tcga, self.cos, self.mod,self.cancer_type)
        data = a.annotated_data()
        data.set_index("isoform",drop=False,inplace=True)
        return data

    def mut_sites(self):
        #所有显著或者非显著的蛋白（集合形式）
        proteins = self.mutated_proteins()
        #某种癌症，某一个修饰的注释信息（dataframe形式）
        annotation = self.annotated_sites()
        annotation.fillna('',inplace=True)
        #每个蛋白注释信息解读返回list
        annotation[["in","out"]] = annotation[["in","out"]].applymap(self.analysis)
        #取出所有蛋白与某一癌症的某种修饰的蛋白的交集
        common_protein = proteins & set(annotation.isoform)
        if len(common_protein) != 0:
            # data_in = pd.DataFrame(["protein_id","cancer","isoform","seq","pos","from","to","sample_num","direct_mutation"])    #in domain的修饰
            # data_out = pd.DataFrame(["protein_id","cancer","isoform","seq","pos","from","to","sample_num","direct_mutation"])   #out domain的修饰
            entry_in = []
            entry_out = []
            for pro in common_protein:
                try:
                    data = annotation.ix[pro]
                    # in motif修饰的蛋白
                    d_in = data[["cancer_type","Canonical","in"]]
                    #series类型，即icgc、tcga和cosmic总共对该蛋白质有一条记录
                    if "Series" in str(type(d_in)):
                        if d_in["in"]!="":
                            l = d_in["in"].split(';')
                            for i in l:
                                i = i.split('\t')
                                i.insert(0,pro)
                                i.insert(1,d_in["Canonical"])
                                i.insert(2,d_in["cancer_type"])
                                #判断是否是中心位点修饰
                                if i[3]==i[8]:
                                    i.append ("yes")
                                else:
                                    i.append("no")
                                entry_in.append(i)
                    #dataframe类型，即存在多条记录
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
                                    if l[3] == l[8]:
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
                                if i[3]==i[8]:
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
                                    if l[3] == l[8]:
                                        l.append("yes")
                                    else:
                                        l.append("no")
                                    entry_out.append(l)
                except KeyError:
                    continue
            columns=["protein_id","Canonical","cancer_type","K position","position in motif","left flank","right flank","seq","mutated position","from","to","sample_num","direct mutation"]
            in_domain_data = pd.DataFrame(entry_in)
            in_domain_data.columns = columns
            out_domain_data = pd.DataFrame(entry_out)
            out_domain_data.columns = columns
            #每种癌症创建自己的目录并写入每种癌症每个修饰的文件
            os.system("mkdir "+ os.path.dirname(os.path.realpath(__file__)) + '/' + self.cancer_type + "/")
            in_domain_file = os.path.dirname(os.path.realpath(__file__)) + '/' + self.cancer_type + "/" + self.mod + "_in_motif.txt"
            out_domain_file = os.path.dirname(os.path.realpath(__file__)) + '/' + self.cancer_type + "/" +self.mod + "_out_motif.txt"
            in_domain_data.to_csv(in_domain_file,sep='\t',index=False,mode='w',encoding="utf-8")
            out_domain_data.to_csv(out_domain_file,sep='\t',index=False,mode='w',encoding="utf-8")
        else:
            in_domain_file = os.path.dirname(os.path.realpath(__file__)) + '/' + self.cancer_type + "/" +self.mod + "_in_motif.txt"
            out_domain_file = os.path.dirname(os.path.realpath(__file__)) + '/' + self.cancer_type + "/" +self.mod + "_out_motif.txt"
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
            #文件中用“|”分隔开不同的序列突变，split后得到这个蛋白中有几个突变的序列
            info = i.strip().split("|")
            #pattern中所提取的信息是中心位点位置、修饰位置和左、右的侧链长度、序列和突变的详细信息（第几位从什么突变到什么）
            pattern = re.compile('(\d+),(\d+),(\d),(\d),(.*?),(.*)',re.S)
            for item in info:
                res = re.findall(pattern,item)
                pos_and_sample = res[0][-1]
                #不同的详细信息是由“；”分隔开的，split后可以逐条提取
                pos_and_sample = pos_and_sample.split(';')
                for j in pos_and_sample:
                    #突变位点（蛋白序列上第几个aa）
                    pos = j.split('-')[0]
                    #所有的样本id都存在于一个中括号中[]并通过“，”逗号分隔，split获取到某个突变的所有样本
                    p = re.compile("\[(.*?)\]", re.S)
                    samples = re.findall(p, j)[0]
                    samples = samples[1:-1].split(',')
                    for s in samples:
                        try:
                            mut = pos[1:-1]  #是突变在aa序列的位置
                            f = pos[0]       #从什么aa突变而来
                            t = pos[-1]      #突变为什么aa
                            #result中是 中心位点位置、修饰位置和左、右的侧链长度、序列、mut，f,t,样本号
                            result = [res[0][0],res[0][1],res[0][2],res[0][3],mut,f,t,s]
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
    #小于0.05的被认为是显著蛋白；反之，是非显著蛋白。
    sig = 0.05
    icgc = '/data1/data/zengyanru/LysineTCGA/ICGC/cancer_annovar_correct/info_v3'
    tcga = '/data1/data/zengyanru/LysineTCGA/find_mutation_in_motif/all_info_V3'
    cosmic = '/data1/data/zengyanru/LysineTCGA/COSMIC/correct_protein'
    # mutation_file = '/data1/hbs/all_cancer_analysis/all_cancer_modification/UCEC'+mod+'_all.txt2.txt'
    cancer_kind = ['ACC','BLCA','BRCA','CESC','CHOL','COAD','DLBC','ESCA','GBM','HNSC','KICH','KIRC','KIRP',
                   'LAML','LGG','LIHC','LUAD','LUSC','MESO','OV','PAAD','PCPG','PRAD','READ','SARC','SKCM',
                   'STAD','TGCT','THCA','THYM','UCEC','UCS','UVM']

    for opt,arg in opts:
        # if opt in ('-m','--mod'):
        #     mod = arg
        if opt in ('-s','--sig'):
            sig = args
        elif opt in ('-i','--icgc'):
            icgc = arg
        elif opt in ('-t','--tcga'):
            tcga = arg
        elif opt in ('-c','--cosmic'):
            cosmic = arg
        # elif opt in ('-k','--kind'):
        #     kind = arg
        else:
            print("invalid input")
            exit(1)

    pool = mp.Pool(14)

    def motif_mod_doc(mutation_file, sig, icgc, tcga, cosmic, mod, kind):
        m = Mutation_Sites(mutation_file, sig, icgc, tcga, cosmic, mod, kind)
        p = m.mut_sites()

    files = os.listdir("/data1/hbs/all_cancer_analysis/all_cancer_modification")
    for kind in cancer_kind:
        pattern = re.compile(kind + "_", re.S)
        target_files = []
        for f in files:
            if re.match(pattern, f)!=None:
                target_files.append(f)  #得到某一种癌症的全部文件
        print(files)
        # for tf in target_files:
        #     mod = tf.split("_")[1][0:3].lowercase
        #     tf = "/data1/hbs/all_cancer_analysis/all_cancer_modification" + tf
        #     pool.apply_async(motif_mod_doc, (tf, sig, icgc, tcga, cosmic,mod, kind))