#-*- coding:utf-8 -*-
import os
import sys
import pandas as pd
import re
import warnings


class TotalFastr:
    def __init__(self, total_file):
        self.total_file = total_file

    def load_data(self):
        data = pd.read_csv(self.total_file, sep='\t', header=None)
        columns = ["protein_id","Uniprot Accession", "position", "mod_type", "sequence","Species","PMIDs"]
        data.columns = columns
        return data


class SignificantProt:
    def __init__(self, sig_prot_file):
        self.target_file = sig_prot_file

    def load_data(self):
        data = pd.read_csv(self.target_file, sep='\t', header=None, usecols=[0, 1])
        data.columns = ["protein", "Canonical"]
        return data


class InsignificantProt:
    def __init__(self, target_dir, total_file):
        self.mod_dict = {"ace": "Acetylation", "gly": "Glycation", "mal": "Malonylation", "met": "Methylation",
                         "suc": "Succinylation", "sumo": "Sumoylation", "ubi": "Ubiquitination"}
        self.target_dir = target_dir
        self.total_file = total_file

    def select_prot(self):
        files = os.listdir(self.target_dir)
        pattern = re.compile("([a-z]+)_.*",re.S)
        t = TotalFastr(self.total_file)
        total = t.load_data()
        for f in files:
            mod_type = re.findall(pattern, f)
            if mod_type != []:
                mod_type = mod_type[0]
                mod_type = self.mod_dict[mod_type]
                f = self.target_dir + '/' + f
                sig_prot = SignificantProt(f)
                sig_prot = sig_prot.load_data()  #sig=yes的蛋白id的集合
                relevant_total = total[total["mod_type"] == mod_type]
                relevant_total.set_index("Uniprot Accession", drop=False, inplace=True)
                all_proteins = set(relevant_total["Uniprot Accession"])
                common_proteins = set(sig_prot["protein"]) & all_proteins
                for p in common_proteins:
                    need_drop = relevant_total[relevant_total["Uniprot Accession"]==p]
                    relevant_total.drop(need_drop,axis=0,inplace=True)
                    need_drop.to_csv(os.path.dirname(os.path.realpath(__file__))+"/"+mod_type+"_sig.txt", sep='\t',header = None,index=False,mode='a',encoding="utf-8")
                relevant_total.to_csv(os.path.dirname(os.path.realpath(__file__))+"/"+mod_type+"_insig.txt", sep='\t',index=False,mode='w',encoding="utf-8")

if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    # t = TotalFastr("/data1/data/chenli/lys_analysis/second_analysis/Total.elm")
    # data = t.load_data()
    # print(data)
    # s = SignificantProt("/data1/hbs/mutation_separate(sig=yes)/ace_in_domain.txt")
    # data = s.load_data()
    # print(data)
    i = InsignificantProt('/data1/hbs/mutation_separate(sig=yes)', '/data1/data/chenli/lys_analysis/second_analysis/Total.elm')
    i.select_prot()