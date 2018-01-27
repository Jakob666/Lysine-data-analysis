#-*- coding:utf-8 -*-
import pandas as pd
from if_overlapped import if_overlapped as io
# from boostrapTest import Boostrap_test
import warnings


class Motif_bootstrap_test:

    def load_data(self, file):
        data = pd.read_csv(file, sep="\t", usecols=[0, 1, 2, 3, 4, 5, 8], encoding="utf-8")
        data.set_index("Canonical", inplace=True)
        return data

    def sequence_length(self, sequence_file):
        sequences = pd.read_csv(sequence_file, sep="\t", usecols=[1, 4], encoding="utf-8")
        sequences.drop_duplicates(inplace=True)
        sequences["prot_length"] = 0
        sequences = sequences.apply(self.calc_length, axis=1)
        sequences.drop("Sequence", axis=1, inplace=True)
        sequences.set_index("Uniprot Accession", inplace=True)
        return sequences

    def calc_length(self, i):
        i["prot_length"] = len(i["Sequence"])
        return i

    def calc_motif_start_end(self, i):
        start = i["K position"] - i["left flank"]
        end = i["K position"] + i["right flank"]
        if start < 0:
            start = 0
        if end > i["prot_length"]:
            end = i["prot_length"]
        i["start"] = start
        i["end"] = end
        return i

    def merge_overlap(self, df):
        new_dataframe = pd.DataFrame(columns=["protein_id","Canonical","prot_length","start","end","mutated position"])
        idxs = set(df.index)
        for idx in idxs:
            intervals = df.ix[idx]
            protein_length = intervals["prot_length"]
            if "Series" in str(type(intervals)):  # Series说明该蛋白对应的domain只有一个
                d = pd.DataFrame([[intervals["protein_id"],idx,protein_length,intervals["start"],intervals["end"],intervals["mutated position"]]],
                                 columns=["protein_id","Canonical","prot_length","start","end","mutated position"])
                new_dataframe = pd.concat([new_dataframe, d], axis=0)

            elif "DataFrame" in str(type(intervals)):  # 蛋白上有多个domain
                intervals.sort_values("start", inplace=True)  # 起始位点从小到大排序
                motif_area = list(zip(list(intervals["start"]), list(intervals["end"])))
                need_drop = []
                for i in range(len(motif_area) - 1):
                    interval1 = motif_area[i]
                    interval2 = motif_area[i + 1]
                    overlapped = io(interval1, interval2)

                    if overlapped:
                        motif_area[i + 1] = [min([interval1[0], interval2[0]]), max([interval1[1], interval2[1]])]
                        need_drop.append(interval1)
                for j in need_drop:
                    motif_area.remove(j)
                    motif_area = list(map(list, motif_area))
                d = pd.DataFrame(motif_area, columns=["start", "end"])
                d["protein_id"] = list(intervals["protein_id"])[0]
                d["Canonical"] = idx
                d["prot_length"] = list(protein_length)[0]
                d["mutated position"] = list(intervals["mutated position"])[0]
                new_dataframe = pd.concat([new_dataframe, d], axis=0)
        new_dataframe.set_index("protein_id", drop=False, inplace=True)
        # print(new_domain_data)
        return new_dataframe

    def calc_motif_length(self, df):
        df["motif length"] = df["end"] - df["start"] + 1
        return df

    def preprocess(self, motif_data, length_data):
        motif_data = pd.merge(left=motif_data, right=length_data, left_index=True, right_index=True)
        motif_data["start"] = 0
        motif_data["end"] = 0
        motif_data.apply(self.calc_motif_start_end, axis=1)
        motif_data = self.merge_overlap(motif_data)
        motif_data = self.calc_motif_length(motif_data)

    def calc_non_motif_length(self, i):
        idx = i["protein_id"]
        try:
            non_motif_length = self.in_motif.ix[idx]["non_motif length"]
        except:
            # non_motif_length = i["prot_length"]
            return i
        i["non_motif length"] = non_motif_length
        return i

    # def bootstrap(self, l1, l2):
    #     bs = Boostrap_test(side="one-side")
    #     res, p_val, p_ref, x_mean, y_mean, x_std, y_std = bs.main(l1, l2)
    #
    #     return res, p_val, p_ref, x_mean, y_mean, x_std, y_std

    def main(self, mod_type, cancer_type, sig):
        # path = "/data1/hbs/all_cancer_analysis/all_cancer_modification_calculate/" + sig + "/" + cancer_type + "/" + mod_type
        # in_motif_central = "_in_motif_centralDetailMut.txt"
        # in_motif_surround = "_in_motif_surroundDetailMut.txt"
        # out_motif = "_out_motifDetailMut.txt"
        #
        # in_central_file = path + "/" + in_motif_central
        # in_surround_file = path + "/" + in_motif_surround
        # out_motif_file = path + "/" + out_motif
        in_central_file = "C:/Users/hbs/Desktop/in_out_motif _test/_in_motif_centralDetailMut.txt"
        central_calc_file = "C:/Users/hbs/Desktop/in_out_motif _test/_in_motif_centralMut.txt"
        in_surround_file = "C:/Users/hbs/Desktop/in_out_motif _test/_in_motif_surroundDetailMut.txt"
        surround_calc_file = "C:/Users/hbs/Desktop/in_out_motif _test/_in_motif_surroundMut.txt"
        out_motif_file = "C:/Users/hbs/Desktop/in_out_motif _test/_out_motifDetailMut.txt"
        out_calc_file = "C:/Users/hbs/Desktop/in_out_motif _test/_out_motifMut.txt"

        central_data = self.load_data(in_central_file)
        central_calc = pd.read_csv(central_calc_file, sep="\t",encoding="utf-8")

        surround_data = self.load_data(in_surround_file)
        surround_calc = pd.read_csv(surround_calc_file, sep="\t", encoding="utf-8")

        out_data = self.load_data(out_motif_file)
        out_calc = pd.read_csv(out_calc_file, sep="\t", encoding="utf-8")

        in_calc = pd.concat([central_calc, surround_calc],axis=0)
        in_calc = in_calc.groupby(["protein_id"],as_index=False)["sample_count"].sum()

        # prot_length = self.sequence_length("/data1/hbs/total_fastr/Total.elm")
        prot_length = self.sequence_length("C:/Users/hbs/Desktop/Total.elm")
        central_data = pd.merge(left=central_data, right=prot_length, left_index=True, right_index=True)
        central_data["start"] = 0
        central_data["end"] = 0
        central_data = central_data.apply(self.calc_motif_start_end, axis=1)
        in_central = self.merge_overlap(central_data)

        surround_data = pd.merge(left=surround_data, right=prot_length, left_index=True, right_index=True)
        surround_data["start"] = 0
        surround_data["end"] = 0
        surround_data = surround_data.apply(self.calc_motif_start_end, axis=1)
        in_surround = self.merge_overlap(surround_data)

        out_data = pd.merge(left=out_data, right=prot_length, left_index=True, right_index=True)

        in_motif = pd.concat([central_data, surround_data], axis=0)
        in_motif = self.merge_overlap(in_motif)
        in_motif["motif length"] = in_motif["end"] - in_motif["start"] + 1
        in_motif = in_motif.groupby(["protein_id", "Canonical", "prot_length"], as_index=False)["motif length"].sum()

        in_motif.set_index("Canonical", drop=False, inplace=True)
        in_motif["non_motif length"] = in_motif["prot_length"] - in_motif["motif length"]

        self.in_motif = in_motif.copy()
        in_motif = pd.merge(left=in_motif,right=in_calc,left_on="protein_id",right_on="protein_id")
        in_motif["mut_per_k_aa"] = in_motif["sample_count"] / in_motif["non_motif length"] * 1000
        #
        out_data["non_motif length"] = out_data["prot_length"]
        out_data = out_data.apply(self.calc_non_motif_length, axis=1)
        out_data = out_data[["protein_id","non_motif length"]]

        out_data.drop_duplicates(inplace=True)
        out_data = pd.merge(left=out_data, right=out_calc, left_on="protein_id", right_on="protein_id")
        out_data["mut_per_k_aa"] = out_data["sample_count"] / out_data["non_motif length"] * 1000
        #
        # data1 = list(in_motif["mut_per_k_aa"])
        # data2 = list(out_data["mut_per_k_aa"])
        #
        # res, p_val, p_ref, x_mean, y_mean, x_std, y_std = self.bootstrap(data1, data2)
        #
        # content = "%s in %s mutations counts between motif and non-motif area on significant proteins' sequences, " % (
        # res, mod)
        # if p_val > p_ref:
        #     content += "p value %s > %s\n" % (p_val, p_ref)
        # else:
        #     content += "p value %s < %s\n" % (p_val, p_ref)
        # content += "average mutations in motif area: %s\naverage mutations in non-motif area: %s\n" % (x_mean, y_mean)
        # content += "standard variance in motif area: %s\nstandard variance in non-motif area: %s\n" % (x_std, y_std)
        # if res == "Significant difference" and x_mean > y_mean:
        #     content += "the occurance of mutations in motif area is significant higher than which in non-motif area.\n\n"
        # elif res == "Significant difference" and x_mean < y_mean:
        #     content += "the occurance of mutations in motif area is significant lower than which in non-motif area.\n\n"
        # print(content)


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    in_central_file = "C:/Users/hbs/Desktop/in_out_motif _test/_in_motif_centralDetailMut.txt"
    in_surround_file = "C:/Users/hbs/Desktop/in_out_motif _test/_in_motif_surroundDetailMut.txt"
    out_motif_file = "C:/Users/hbs/Desktop/in_out_motif _test/_out_motifDetailMut.txt"
    d = Motif_bootstrap_test()
    d.main("Ace", "BRCA","")

    # ms = Motif_bootstrap_test()
    # central_data = ms.load_data(in_central_file)
    # surround_data = ms.load_data(in_surround_file)
    # out_data = ms.load_data(out_motif_file)
    # prot_length = ms.sequence_length("C:/Users/hbs/Desktop/Total.elm")
    #
    # central_data = pd.merge(left=central_data, right=prot_length, left_index=True, right_index=True)
    # central_data["start"] = 0
    # central_data["end"] = 0
    # central_data = central_data.apply(ms.calc_motif_start_end, axis=1)
    # in_central = ms.merge_overlap(central_data)
    #
    # surround_data = pd.merge(left=surround_data, right=prot_length, left_index=True, right_index=True)
    # surround_data["start"] = 0
    # surround_data["end"] = 0
    # surround_data = surround_data.apply(ms.calc_motif_start_end, axis=1)
    # in_surround = ms.merge_overlap(surround_data)
    #
    # out_data = pd.merge(left=out_data, right=prot_length, left_index=True, right_index=True)
    # # out_data["start"] = 0
    # # out_data["end"] = 0
    # # out_data = out_data.apply(ms.calc_motif_start_end, axis=1)
    # # out = ms.merge_overlap(out_data)
    #
    # in_motif = pd.concat([central_data,surround_data], axis=0)
    # in_motif = ms.merge_overlap(in_motif)
    # in_motif["motif length"] = in_motif["end"] - in_motif["start"] + 1
    # in_motif = in_motif.groupby(["protein_id","Canonical","prot_length"], as_index=False)[["sample_count","motif length"]].sum()
    # in_motif.set_index("Canonical", drop=False, inplace=True)
    # in_motif["non_motif length"] = in_motif["prot_length"] - in_motif["motif length"]
    # ms.in_motif = in_motif.copy()
    # in_motif["mut_per_k_aa"] = in_motif["sample_count"]/in_motif["non_motif length"] *1000
    #
    # out_data["non_motif length"] = out_data["prot_length"]
    # out_data = out_data.apply(ms.calc_non_motif_length, axis=1)
    # out_data = out_data.groupby(["protein_id","non_motif length"], as_index=False)["sample_count"].sum()
    # out_data["mut_per_k_aa"] = out_data["sample_count"]/out_data["non_motif length"] *1000
    # print(out_data)
