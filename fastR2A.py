#-*- coding:utf-8 -*-
import pandas as pd
import os
import re

class FastR2A:
    def __init__(self, fastr_dir):
        self.fastr_dir = fastr_dir

    def load_data(self, f):
        data = pd.read_csv(f, sep='\t', encoding="utf-8")
        return data

    def formFastr(self):
        files = os.listdir(self.fastr_dir)
        for f in files:
            outputFile = f.split(".")[0]
            f = self.fastr_dir + '/' + f
            data = self.load_data(f)
            data["fasta"] = ">" + data["protein_id"] + "@" + data["Uniprot Accession"] + "@" + data["mod_type"] + "@" + data["Species"] +"@" + data["position"] + "\n" \
                            + data["sequence"]
            data = data["fasta"]
            outputFile = "/data1/hbs/total_fasta/" + outputFile + ".fasta"
            data.to_csv(outputFile, sep="\t", index=False, mode="w", encoding="utf-8")
            content = open(outputFile,"r").read()
            pattern = re.compile('\"', re.S)
            content = re.sub(pattern, "", content)
            f = open(outputFile, "w")
            f.write(content)
            f.close()


if __name__ == "__main__":
    f = FastR2A("/data1/hbs/mutation_separate(sig!=yes)")
    f.formFastr()
