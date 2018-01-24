#-*- coding:utf-8 -*-
import pandas as pd

data = pd.read_csv("C:/Users/hbs/Desktop/in_domain_in_motif_surroundDetailMut.txt", sep="\t")
data.set_index("Uniprot Accession", drop=False, inplace=True)
l = data.groupby(["Uniprot Accession","domain length"], as_index=False)["sample_count"].sum()
print(l)