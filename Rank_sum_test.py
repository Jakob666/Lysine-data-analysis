# -*- coding:utf-8 -*-
import pandas as pd
import numpy as np
import re
# from scipy import stats
import warnings


class Rank_sum_test:
    def __init__(self, datasets):
        '''形式很像是 r × c 的列联表，column是不同总体，row是总体的不同分类
           首先传入的是每个总体相应类别的文件，通过dict的形式传入'''
        self.datasets = datasets
        self.population_list = list(self.datasets.keys()) #不同的总体
        self.classifications = list(self.datasets[self.population_list[0]].keys())  #每个总体基本分类相同

    def load_data(self):
        for p in self.population_list:
            for c in self.classifications:
                f = self.datasets[p][c]
                data = pd.read_csv(f, sep="\t", header=None)
                data.columns = ["protein_id", "sample_count"]
                counting = np.array(list(data["sample_count"])).sum()   #计算每类的样本数目
                self.datasets[p][c] = counting
        return None

    def form_contingency_table(self):
        self.contingency_table = pd.DataFrame(self.datasets)
        self.contingency_table["row_sum"] = self.contingency_table.apply(sum, axis=1)   #每一级别的样本数目
        row_sum = list(self.contingency_table["row_sum"])
        ave_rank = pd.DataFrame(self.ave_rank(row_sum))
        self.contingency_table = pd.merge(left=self.contingency_table,right=ave_rank, left_index=True, right_index=True)
        self.contingency_table, rank_sum_cols = self.rank_sum(self.contingency_table)
        self.contingency_table.loc["col_sum"] = self.contingency_table.apply(sum, axis=0)   #每列数据求和

        total_samples = self.contingency_table.ix["col_sum"]["row_sum"] #总样本数目

        score = self.calc_score(rank_sum_cols, total_samples)

        p_val = self.get_p_val(score)
        if p_val < 0.05:
            res = "Significant difference"
        else:
            res = "Insignificant difference"
        return res, p_val

    def ave_rank(self, l):
        '''求每个分类的平均秩次，传入的 l 的元素是每个分类所包含的样本总数（此时不考虑总体）'''
        result = dict()
        current_index = list(self.contingency_table.index)
        for i in range(len(l)):
            if i == 0:
                ave_rank = (1.0 + l[0]) / len(self.population_list)
                result[current_index[0]] = ave_rank
            else:
                ave_rank = (l[i-1] + 1 + l[i]) / len(self.population_list)
                result[current_index[i]] = ave_rank
        result = {"ave_rank":result}
        return result

    def rank_sum(self, df):
        indexs = list(df.index)
        rank_sum_cols = []
        rank_sum = dict()
        for p in self.population_list:
            col_name = p + "_rank_sum"
            rank_sum_cols.append(col_name)
            rank_sum[col_name] = dict()
            for i in indexs:
                rank_sum[col_name][i] = df.ix[i]["ave_rank"] * df.ix[i][p]
        rank_sum = pd.DataFrame(rank_sum)
        df = pd.merge(left=self.contingency_table, right=rank_sum, left_index=True, right_index=True)
        return df, rank_sum_cols

    def calc_score(self, rank_sum_cols, total_samples):
        if len(self.population_list) < 2:
            return "error"

        elif len(self.population_list) == 2:
            samples_list = list(self.contingency_table.ix["col_sum"][self.population_list]) #两个总体的样本总数
            T_list = list(self.contingency_table.ix["col_sum"][rank_sum_cols])  #两组样本的秩和
            T = min(T_list)     #选最小的秩和作为统计量
            idx = T_list.index(T)

            numerator = np.abs(T - float(samples_list[idx]) * (total_samples+1)/2) - 0.5
            denominator = np.sqrt(samples_list[0] * samples_list[1] * float(total_samples + 1)/12)
            z_score = numerator/denominator

            c_val = 0
            for c in self.classifications:
                tj = self.contingency_table.ix[c]["row_sum"]
                print(tj)
                c_val += tj ** 3 - tj
            c_val = 1.0 - float(c_val) / (total_samples ** 3 - total_samples)

            return z_score/np.sqrt(c_val)

        elif len(self.population_list) > 2: #当样本组大于等于两个的时候
            sum = 0
            for i in range(len(self.population_list)):
                Ti = self.contingency_table.ix["col_sum"][self.population_list[i]]
                ni = self.contingency_table.ix["col_sum"][rank_sum_cols[i]]
                sum += Ti**2/ni
            h_score = 12.0/total_samples*(total_samples+1) * sum - 3*(total_samples + 1)

            c_val = 0
            for c in self.classifications:
                tj = self.contingency_table.ix[c]["row_sum"]
                c_val += tj**3 - tj
            c_val = 1.0 - float(c_val)/(total_samples**3 - total_samples)

            return h_score/c_val

    def get_p_val(self, score):
        if len(self.population_list) == 2:
            z_standard = 1.96   #α=0.05，df=2时
            if score > z_standard:
                p_val = 0.0
            else:
                p_val = 0.06
        # elif len(self.population_list) > 2:
        #     df = len(self.population_list) - 1
        #     p_val = 1.0 - stats.chi2.cdf(score, df)

        return p_val


class Testing:
    def __init__(self, target_dir, population_list,classification):
        self.mod_list = ["ace", "gly", "mal", "met", "suc", "sumo", "ubi"]
        self.target_dir = target_dir
        self.population_list = population_list
        self.classification = classification

    def form_datasets(self, mod):
        datasets = dict()
        for p in self.population_list:
            datasets[p] = dict()
            for c in self.classification:
                datasets[p][c] = self.target_dir[p] + mod + c + "Mut.txt"
        return datasets

    def main(self):
        content = ""
        for mod in self.mod_list:
            datasets = self.form_datasets(mod)
            r = Rank_sum_test(datasets)
            r.load_data()
            res, p_val = r.form_contingency_table()
            content += "%s in %s modification counts p value %s\n"%(res, mod, p_val)
        f = open("/data1/hbs/UCEC_significant_difference/sig_insig_in_motif_comparison(rs).txt", "w")
        f.write(content)
        f.close()

if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    # datasets = {"sig":{"central":"C:/Users/hbs/Desktop/sig/gly_in_motif_centralMut.txt",
    #                   "surround":"C:/Users/hbs/Desktop/sig/gly_in_motif_surroundMut.txt",
    #                   "out":"C:/Users/hbs/Desktop/sig/gly_out_motifMut.txt"},
    #            "insig":{"central":"C:/Users/hbs/Desktop/insig/gly_in_motif_centralMut.txt",
    #                     "surround":"C:/Users/hbs/Desktop/insig/gly_in_motif_surroundMut.txt",
    #                     "out":"C:/Users/hbs/Desktop/insig/gly_out_motifMut.txt"}}
    # r = Rank_sum_test(datasets)
    # r.load_data()
    # res, p_val = r.form_contingency_table()
    # print(res, p_val)
    target_dir = {"sig":"/data1/hbs/centralOrNotSigIsyes",
                  "insig":"/data1/hbs/centralOrNotSigNotyes"}
    population_list = ["sig", "insig"]
    classification = ["_in_motif_central", "_in_motif_surround", "_out_motif"]
    t = Testing(target_dir, population_list, classification)
    t.main()