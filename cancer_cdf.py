#-*- coding:utf-8 -*-
'''
Description:
    分别使用get_phastCons_score.py文件得到的 phyloP100way_vertebrate和phastCons100way_vertebrate两个值
    绘制33种癌症的cdf图。
    两种分值哪个效果好最终使用哪个。
=====================================================
@author:hbs
@date:2018-1-30
@version:1.1
'''
import pandas as pd
import numpy as np
from draw_cdf import CDFDrawer
import multiprocessing as mp
from matplotlib import pyplot as plt
import warnings


class CancerCdfGraphs(CDFDrawer):
    def __init__(self, cancer_type, fig_name):
        self.cancer = cancer_type
        self.fig_name = fig_name

    def load_data(self, mut_conservation_file):
        self.conservation_score = pd.read_csv(mut_conservation_file, sep="\t", usecols=[2, 3], encoding="utf-8")
        try:
            self.conservation_score = self.conservation_score[self.conservation_score["phyloP100way_vertebrate"] != "None"]
            self.conservation_score = self.conservation_score[self.conservation_score["phastCons100way_vertebrate"] != "None"]
            self.conservation_score = self.conservation_score[self.conservation_score["phyloP100way_vertebrate"] != "emp"]
            self.conservation_score = self.conservation_score[self.conservation_score["phastCons100way_vertebrate"] != "emp"]

        except:
            pass
        return None

    def draw_with_phylop_score(self):
        c = CDFDrawer("continuous")
        values = self.conservation_score["phyloP100way_vertebrate"].values.astype("float")
        title = self.cancer + " cdf graph"
        c.main(values, x_label="phylop conservation score", y_label="probability", color_and_mark="r-", graph_title=title, fig_name=self.fig_name)

    def draw_with_phastCons_score(self):
        c = CDFDrawer("continuous")
        values = self.conservation_score["phastCons100way_vertebrate"].values.astype("float")
        title = self.cancer + " cdf graph"
        c.main(values, x_label="phastCons conservation score", y_label="probability", color_and_mark="b-", graph_title=title, fig_name=self.fig_name)

    def draw_four_lines(self, file_list, x_range="auto", y_range=[0,1.1], use_score="phastcons"):
        x_axis_plots = []
        y_axis_plots = []
        for f in file_list:
            self.load_data(f)
            if use_score=="phastcons":
                values = self.conservation_score["phastCons100way_vertebrate"].values.astype("float")
            elif use_score == "phylop":
                values = self.conservation_score["phyloP100way_vertebrate"].values.astype("float")
            else:
                raise ValueError("invalid input of use_score parament")
            values, test_num, distribution_seq = self.contiunous_data(values, distribution_seq=None, distribution_func=None)
            x_plots, y_plots = self.create_plots(values, distribution_seq)
            x_axis_plots.append(x_plots)
            y_axis_plots.append(y_plots)

        if x_range == "auto":
            x_range = [0,1.1]

        fig = plt.figure(figsize=(4, 4))
        plt.xlim(x_range[0], x_range[1])
        plt.ylim(y_range[0], y_range[1])
        plt.plot([x_range[0], x_range[1]], [1, 1], color="black", linestyle="--", alpha=0.4)
        plt.plot([1,1],[y_range[0], y_range[1]], color="black", linestyle="--", alpha=0.4)

        colors = ["red", "blue", "darkgreen", "orange"]
        linestyles = ["-", "-.", ":", "--"]
        labels = ["lysine related in domain mutation", "lysine related out domain mutation", "lysine irrelated in domain mutation", "lysine irrelated out domain mutation"]

        for i in range(len(x_axis_plots)):
            x_plots = x_axis_plots[i]
            y_plots = y_axis_plots[i]
            plt.plot(x_plots, y_plots, color=colors[i], linestyle=linestyles[i], label=labels[i])
        plt.xlabel("conservation_score")
        plt.ylabel("probability")
        plt.legend(fontsize=5, loc="upper left")

        plt.title(self.cancer+" CDF graph")

        fig.savefig(self.fig_name)

    def combine_all_cancer_data(self, y_range=[0,1.05]):
        related_in_domain_data = pd.DataFrame(columns=["phyloP100way_vertebrate", "phastCons100way_vertebrate"])
        related_out_domain_data = pd.DataFrame(columns=["phyloP100way_vertebrate", "phastCons100way_vertebrate"])
        irrelated_in_domain_data = pd.DataFrame(columns=["phyloP100way_vertebrate", "phastCons100way_vertebrate"])
        irrelated_out_domain_data = pd.DataFrame(columns=["phyloP100way_vertebrate", "phastCons100way_vertebrate"])

        cancer_kind = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC',
                       'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ',
                       'SARC','SKCM','STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']
        data_source = "C:/Users/hbs/Desktop/conservation_score/"
        for c in cancer_kind:
            related_in_domain_file = data_source + "related_in_domain/" + c + "_mut_site_conservation.txt"
            related_out_domain_file = data_source + "related_out_domain/" + c + "_mut_site_conservation.txt"
            irrelated_in_domain_file = data_source + "irrelated_in_domain/" + c + "_mut_site_conservation.txt"
            irrelated_out_domain_file = data_source + "irrelated_out_domain/" + c + "_mut_site_conservation.txt"
            d1 = self.load_all_data(related_in_domain_file)
            d2 = self.load_all_data(related_out_domain_file)
            d3 = self.load_all_data(irrelated_in_domain_file)
            d4 = self.load_all_data(irrelated_out_domain_file)
            related_in_domain_data = pd.concat([related_in_domain_data, d1], axis=0)
            related_out_domain_data = pd.concat([related_out_domain_data, d2], axis=0)
            irrelated_in_domain_data = pd.concat([irrelated_in_domain_data, d3], axis=0)
            irrelated_out_domain_data = pd.concat([irrelated_out_domain_data, d4], axis=0)
        #生成索引
        related_in_domain_data["idx"] = ""
        related_in_domain_data = related_in_domain_data.apply(self.create_idx, axis=1)
        related_out_domain_data["idx"] = ""
        related_out_domain_data = related_out_domain_data.apply(self.create_idx, axis=1)
        irrelated_in_domain_data["idx"] = ""
        irrelated_in_domain_data = irrelated_in_domain_data.apply(self.create_idx, axis=1)
        irrelated_out_domain_data["idx"] = ""
        irrelated_out_domain_data = irrelated_out_domain_data.apply(self.create_idx, axis=1)
        #下方的四行去除组内重复
        related_in_domain_data.drop_duplicates(inplace=True)
        related_out_domain_data.drop_duplicates(inplace=True)
        irrelated_in_domain_data.drop_duplicates(inplace=True)
        irrelated_out_domain_data.drop_duplicates(inplace=True)
        #下方的代码用于去除组间重复
        irrelated_out_need_drop = set(related_in_domain_data.index)|set(related_out_domain_data.index)|set(irrelated_in_domain_data.index) & set(irrelated_out_domain_data.index)
        irrelated_in_need_drop = set(related_in_domain_data.index)|set(related_out_domain_data.index)& set(irrelated_in_domain_data.index)
        related_out_need_drop = set(related_in_domain_data.index)&set(related_out_domain_data.index)

        irrelated_out_domain_data.drop(irrelated_out_need_drop, axis=0, inplace=True)
        irrelated_in_domain_data.drop(irrelated_in_need_drop, axis=0, inplace=True)
        related_out_domain_data.drop(related_out_need_drop, axis=0, inplace=True)

        data_list = [related_in_domain_data, irrelated_in_domain_data, related_out_domain_data, irrelated_out_domain_data]
        x_axis_plots = []
        y_axis_plots = []
        for i in data_list:
            values = i["phastCons100way_vertebrate"].values.astype("float")
            values, test_num, distribution_seq = self.contiunous_data(values, distribution_seq=None,
                                                                      distribution_func=None)
            x_plots, y_plots = self.create_plots(values, distribution_seq)
            x_axis_plots.append(x_plots)
            y_axis_plots.append(y_plots)

        x_range = [0, 1.05]
        fig = plt.figure(figsize=(8, 8), frameon=False)
        plt.xlim(x_range[0], x_range[1])
        plt.ylim(y_range[0], y_range[1])
        plt.plot([x_range[0], x_range[1]], [1, 1], color="black", linestyle="--", alpha=0.4)
        plt.plot([1, 1], [y_range[0], y_range[1]], color="black", linestyle="--", alpha=0.4)

        colors = ["red", "darkgreen", "blue", "orange"]
        linestyles = ["-", ":", "-.", "--"]
        labels = ["lysine related in domain mutation", "lysine irrelated in domain mutation",
                  "lysine related out domain mutation", "lysine irrelated out domain mutation"]

        for i in range(len(x_axis_plots)):
            x_plots = x_axis_plots[i]
            y_plots = y_axis_plots[i]
            plt.plot(x_plots, y_plots, color=colors[i], linestyle=linestyles[i], label=labels[i])
        plt.xlabel("conservation_score")
        plt.ylabel("probability")
        plt.xticks([0,0.25,0.5,0.75,1.0])
        plt.yticks([0, 0.25, 0.5, 0.75,1.0])
        plt.legend(fontsize=12, loc="upper left")

        plt.title("All cancer CDF graph")

        fig.savefig("C:/Users/hbs/Desktop/conservation_cdf/lysine_related_irrelated/all_cancer.pdf")

    def load_all_data(self, f):
        data = pd.read_csv(f, sep="\t", encoding="utf-8")
        try:
            data = data[data["phyloP100way_vertebrate"] != "None"]
            data = data[data["phastCons100way_vertebrate"] != "None"]
            data = data[data["phyloP100way_vertebrate"] != "emp"]
            data = data[data["phastCons100way_vertebrate"] != "emp"]

        except:
            pass
        return data

    def create_idx(self, i):
        idx=i["Uniprot Accession"] + str(i["mutation_site"])
        i["idx"] = idx
        return i

    def only_in_out_motif(self, y_range=[0,1.05]):
        in_motif_data = pd.DataFrame(columns=["phyloP100way_vertebrate", "phastCons100way_vertebrate"])
        out_motif_data = pd.DataFrame(columns=["phyloP100way_vertebrate", "phastCons100way_vertebrate"])

        cancer_kind = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC',
                       'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ',
                       'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']
        data_source = "C:/Users/hbs/Desktop/conservation_score/"
        for c in cancer_kind:
            related_in_domain_file = data_source + "related_in_domain/" + c + "_mut_site_conservation.txt"
            related_out_domain_file = data_source + "related_out_domain/" + c + "_mut_site_conservation.txt"
            irrelated_in_domain_file = data_source + "irrelated_in_domain/" + c + "_mut_site_conservation.txt"
            irrelated_out_domain_file = data_source + "irrelated_out_domain/" + c + "_mut_site_conservation.txt"
            d1 = self.load_all_data(related_in_domain_file)
            d2 = self.load_all_data(related_out_domain_file)
            d3 = self.load_all_data(irrelated_in_domain_file)
            d4 = self.load_all_data(irrelated_out_domain_file)
            in_motif_data = pd.concat([in_motif_data, d1], axis=0)
            in_motif_data = pd.concat([in_motif_data, d2], axis=0)
            out_motif_data = pd.concat([out_motif_data, d3], axis=0)
            out_motif_data = pd.concat([out_motif_data, d4], axis=0)
        # 生成索引
        in_motif_data["idx"] = ""
        in_motif_data = in_motif_data.apply(self.create_idx, axis=1)
        out_motif_data["idx"] = ""
        out_motif_data = out_motif_data.apply(self.create_idx, axis=1)

        in_motif_data.drop_duplicates(inplace=True)
        out_motif_data.drop_duplicates(inplace=True)
        out_need_drop = set(in_motif_data.index)&set(out_motif_data.index)
        out_motif_data.drop(out_need_drop, axis=0, inplace=True)

        data_list = [in_motif_data, out_motif_data]
        x_axis_plots = []
        y_axis_plots = []
        for i in data_list:
            values = i["phastCons100way_vertebrate"].values.astype("float")
            values, test_num, distribution_seq = self.contiunous_data(values, distribution_seq=None,
                                                                      distribution_func=None)
            x_plots, y_plots = self.create_plots(values, distribution_seq)
            x_axis_plots.append(x_plots)
            y_axis_plots.append(y_plots)

        x_range = [0, 1.05]
        fig = plt.figure(figsize=(8, 8), frameon=False)
        plt.xlim(x_range[0], x_range[1])
        plt.ylim(y_range[0], y_range[1])
        plt.plot([x_range[0], x_range[1]], [1, 1], color="black", linestyle="--", alpha=0.4)
        plt.plot([1, 1], [y_range[0], y_range[1]], color="black", linestyle="--", alpha=0.4)

        colors = ["red", "blue", "orange"]
        linestyles = ["-", "--", ":", "-."]
        labels = ["lysine related mutation", "lysine irrelated mutation"]

        for i in range(len(x_axis_plots)):
            x_plots = x_axis_plots[i]
            y_plots = y_axis_plots[i]
            plt.plot(x_plots, y_plots, color=colors[i], linestyle=linestyles[i], label=labels[i])
        plt.xlabel("conservation_score")
        plt.ylabel("probability")
        plt.xticks([0, 0.25, 0.5, 0.75, 1.0])
        plt.yticks([0, 0.25, 0.5, 0.75, 1.0])
        plt.legend(fontsize=12, loc="upper left")

        plt.title("All cancer CDF graph")

        fig.savefig("C:/Users/hbs/Desktop/conservation_cdf/lysine_related_irrelated/all_cancer_in_out_motif.pdf")

    def only_in_out_domain(self, y_range=[0,1.05]):
        in_domain_data = pd.DataFrame(columns=["phyloP100way_vertebrate", "phastCons100way_vertebrate"])
        out_domain_data = pd.DataFrame(columns=["phyloP100way_vertebrate", "phastCons100way_vertebrate"])

        cancer_kind = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC',
                       'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ',
                       'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']
        data_source = "C:/Users/hbs/Desktop/conservation_score/"
        for c in cancer_kind:
            related_in_domain_file = data_source + "related_in_domain/" + c + "_mut_site_conservation.txt"
            related_out_domain_file = data_source + "related_out_domain/" + c + "_mut_site_conservation.txt"
            irrelated_in_domain_file = data_source + "irrelated_in_domain/" + c + "_mut_site_conservation.txt"
            irrelated_out_domain_file = data_source + "irrelated_out_domain/" + c + "_mut_site_conservation.txt"
            d1 = self.load_all_data(related_in_domain_file)
            d2 = self.load_all_data(related_out_domain_file)
            d3 = self.load_all_data(irrelated_in_domain_file)
            d4 = self.load_all_data(irrelated_out_domain_file)
            in_domain_data = pd.concat([in_domain_data, d1], axis=0)
            out_domain_data = pd.concat([out_domain_data, d2], axis=0)
            in_domain_data = pd.concat([in_domain_data, d3], axis=0)
            out_domain_data = pd.concat([out_domain_data, d4], axis=0)
        # 生成索引
        in_domain_data["idx"] = ""
        in_domain_data = in_domain_data.apply(self.create_idx, axis=1)
        out_domain_data["idx"] = ""
        out_domain_data = out_domain_data.apply(self.create_idx, axis=1)

        in_domain_data.drop_duplicates(inplace=True)
        out_domain_data.drop_duplicates(inplace=True)
        out_need_drop = set(in_domain_data.index)&set(out_domain_data.index)
        out_domain_data.drop(out_need_drop, axis=0, inplace=True)

        data_list = [in_domain_data, out_domain_data]
        x_axis_plots = []
        y_axis_plots = []
        for i in data_list:
            values = i["phastCons100way_vertebrate"].values.astype("float")
            values, test_num, distribution_seq = self.contiunous_data(values, distribution_seq=None,
                                                                      distribution_func=None)
            x_plots, y_plots = self.create_plots(values, distribution_seq)
            x_axis_plots.append(x_plots)
            y_axis_plots.append(y_plots)

        x_range = [0, 1.05]
        fig = plt.figure(figsize=(8, 8), frameon=False)
        plt.xlim(x_range[0], x_range[1])
        plt.ylim(y_range[0], y_range[1])
        plt.plot([x_range[0], x_range[1]], [1, 1], color="black", linestyle="--", alpha=0.4)
        plt.plot([1, 1], [y_range[0], y_range[1]], color="black", linestyle="--", alpha=0.4)

        colors = ["red", "blue", "orange"]
        linestyles = ["-",  "-.","--", ":"]
        labels = ["lysine in domain mutation", "lysine out domain mutation"]

        for i in range(len(x_axis_plots)):
            x_plots = x_axis_plots[i]
            y_plots = y_axis_plots[i]
            plt.plot(x_plots, y_plots, color=colors[i], linestyle=linestyles[i], label=labels[i])
        plt.xlabel("conservation_score")
        plt.ylabel("probability")
        plt.xticks([0, 0.25, 0.5, 0.75, 1.0])
        plt.yticks([0, 0.25, 0.5, 0.75, 1.0])
        plt.legend(fontsize=12, loc="upper left")

        plt.title("All cancer CDF graph")

        fig.savefig("C:/Users/hbs/Desktop/conservation_cdf/lysine_related_irrelated/all_cancer_in_out_domain.pdf")


if __name__ == "__main__":
    warnings.filterwarnings("ignore")


    def get_phylop_cdf_graph(f, cancer, fig_name):
        fig_name = "C:/Users/hbs/Desktop/conservation_cdf/phylop/" + fig_name
        c = CancerCdfGraphs(cancer, fig_name)
        c.load_data(f)
        c.draw_with_phylop_score()

    def get_phastcons_cdf_graph(f, cancer, fig_name):
        fig_name = "C:/Users/hbs/Desktop/conservation_cdf/phastcons/" + fig_name
        c = CancerCdfGraphs(cancer, fig_name)
        c.load_data(f)
        c.draw_with_phastCons_score()

    def drawing(f, cancer, fig_name):
        c = CancerCdfGraphs(cancer, fig_name)
        c.draw_four_lines(f)

    # cancer_kind = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP',
    #                'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM',
    #                'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']
    # c = "UVM"
    # data_source_dir = "C:/Users/hbs/Desktop/conservation_score/"
    # sig_in_domain = data_source_dir + "related_in_domain/" + c + "_mut_site_conservation.txt"
    # sig_out_domain = data_source_dir + "related_out_domain/" + c + "_mut_site_conservation.txt"
    # insig_in_domain = data_source_dir + "irrelated_in_domain/" + c + "_mut_site_conservation.txt"
    # insig_out_domain = data_source_dir + "irrelated_out_domain/" + c + "_mut_site_conservation.txt"
    # f = [sig_in_domain,sig_out_domain,insig_in_domain,insig_out_domain]
    # fig_name = "C:/Users/hbs/Desktop/conservation_cdf/lysine_related_irrelated/" + c + "_conservation_cdf.pdf"
    # drawing(f, c, fig_name)
    c = CancerCdfGraphs("","")
    # c.combine_all_cancer_data()
    # c.only_in_out_motif()
    c.only_in_out_domain()

    # 测试代码不必理会
    # f = ["C:/Users/hbs/Desktop/conservation_score/sig_prot_in_domain/BRCA_mut_site_cobservation.txt",
    #      "C:/Users/hbs/Desktop/conservation_score/sig_prot_out_domain/BRCA_mut_site_conservation.txt",
    #      "C:/Users/hbs/Desktop/conservation_score/insig_prot_in_domain/BRCA_mut_site_conservation.txt",
    #      "C:/Users/hbs/Desktop/conservation_score/insig_prot_out_domain/BRCA_mut_site_conservation.txt"]
    # c = CancerCdfGraphs("BRCA", "")
    # c.draw_four_lines(f)
    # c.load_data(f)
    # c.draw_with_phylop_score()
    # c.draw_with_phastCons_score()