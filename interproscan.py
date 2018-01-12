#-*- coding:utf-8 -*-
'''将生成的fasta文件通过interproscan得到蛋白质的domain'''
import os

def run_interproscan(fasta_dir):
    files = os.listdir(fasta_dir)
    for file in files:
        interproscanFile = fasta_dir + "/" + file
        print(interproscanFile)
        outputFile = "/data1/hbs/processedByInterproscan/" + file
        os.system('/bin/sh  /home/targetge/public_software/my_interproscan/interproscan-5.25-64.0/interproscan.sh -appl Pfam,SMART -iprlookup -goterms -i ' + interproscanFile + ' -f tsv -dp -o ' + outputFile + '.tsv')
    print("finish")

if __name__ == "__main__":
    run_interproscan("/data1/hbs/total_fasta")