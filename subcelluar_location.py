#-*- coding:utf-8 -*-
'''
Description:
    使用权威文献中的亚细胞数据对模型预测的显著性蛋白进行亚细胞定位，
    数据集中一个蛋白可能对应多个亚细胞位点，按照数据集中的 main location
    作为最终定位的结果。
========================================
@author:hbs
@date:2018-2-5
@version:1.1
'''
import pandas as pd
import warnings


class SubCellarLoc:
    def __init__(self, subcellar_dataset, sig_prot_dir, output_file):
        self.subcellar = subcellar_dataset
        self.sig_prot_dir = sig_prot_dir
        self.cancer_kind = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH',
                            'KIRC','KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD',
                            'READ','SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']
        self.mod = ["ace", "gly", "mal", "met", "suc", "sum", "ubi"]
        self.output = output_file

    def load_subcellar_data(self):
        self.subcell = pd.read_excel(self.subcellar, sheetname="Protein location results", usecols=[2, 34])
        return None

    def load_sig_prot_data(self):
        self.sig_prot = pd.DataFrame(columns=["protein_id", "Canonical"])
        for c in self.cancer_kind:
            d = self.sig_prot_dir + c
            for m in self.mod:
                in_f = d + "/" + m + "_in_motif.txt"
                out_f = d + "/" + m + "_out_motif.txt"
                try:
                    in_data = pd.read_csv(in_f, sep="\t", usecols=[0, 1], encoding="utf-8")
                    out_data = pd.read_csv(out_f, sep="\t", usecols=[0, 1], encoding="utf-8")
                except:
                    continue
                self.sig_prot = pd.concat([self.sig_prot, in_data, out_data], axis=0)
                self.sig_prot.drop_duplicates(inplace=True)
        return None

    def match_subcellar_location(self):
        self.protein_location = pd.merge(left=self.sig_prot, right=self.subcellar, left_on="Canonical", right_on="Uniprot")
        self.location_classification = pd.DataFrame(columns=["protein_id", "subcell_structure"])
        self.protein_location.apply(self.multi_location_separate, axis=1)
        return None

    def multi_location_separate(self, i):
        main_location = i["IF main protein location"]
        location = main_location.split(";")
        for l in location:
            if l in ["Nucleoli", "Fibrillar center", "Rim of uncleoli"]:
                l = "Nucleoli"
            elif l in ["Golgi", "apparatus"]:
                l = "Glogi"
            elif l in ["Vesicles", "Lipid droplets"]:
                l = "Vesicles"
            elif l in ["Plasma membrane", "Cell junction"]:
                l = "Plasma menbrane"
            elif l in ["Cytosol", "Cytoplasmic bodies", "Aggresome", "Rods and rings"]:
                l = "Cytosol"
            elif l in ["Intermediate", "filaments"]:
                l = "Intermediate"
            elif l in ["Microtubules", "Microtuble end", "Cytokinetic bridge", "Midbody", "Midbody ring", "Mitotic spindle"]:
                l = "Microtubules"
            elif l in ["Centrosome", "MTOC"]:
                l = "Centrosome"
            elif l in ["Actin filaments", "Focal adhesions"]:
                l = "Actin filaments"
            else:
                continue
            d = pd.DataFrame(data=[i["protein_id"], l], columns=["protein_id", "subcell_structure"])
            self.location_classification = pd.concat([self.location_classification, d], axis=0)
        return None

    def to_doc(self):
        self.location_classification.drop_duplicates(inplace=True)
        self.location_classification.to_csv(self.output, sep="\t", index=False, mode="w", encoding="utf-8")
        return None

    def main(self):
        self.load_sig_prot_data()
        self.load_subcellar_data()
        self.match_subcellar_location()
        self.to_doc()


if __name__ == '__main__':
    warnings.filterwarnings("ignore")
    s = SubCellarLoc("/data1/hbs/all_cancer_analysis/sig_prot_subcellar_location/location.xlsx",
                     "/data1/hbs/all_cancer_analysis/all_cancer_sig_prot_motif_modify"
                     , "/data1/hbs/all_cancer_analysis/sig_prot_subcellar_location/sig_prot_location.txt")
    s.main()