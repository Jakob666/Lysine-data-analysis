#-*- coding:utf-8 -*-
'''
Description:
    参考网上大牛的统计学方法，相关链接 https://web.as.uky.edu/statistics/users/pbreheny/621/F12/notes/10-11.pdf
    实现 bootstrap方法，对模型预测出的显著性和非显著性蛋白在特定位置修饰的频率等进行差异分析。
==========================
@author: hbs
@date: 2017-1-6
@version: 1.0
==========================
@update: 对最终的显著性检验的方式进行优化
@version: 1.1
'''
import pandas as pd
import numpy as np
# from matplotlib import pyplot as plt
import os
import warnings
# from t_test import T_test


class Boostrap_test:

    def __init__(self, Confidence=0.95, time=1000):
        '''
        :param Confidence: 置信区间
        :param time: boostrap抽样的次数
        '''
        self.conf = Confidence
        self.times = time
        self.t_val_list = []

    def preprocess(self, x_data, y_data):
        '''求出x、y样本的合并均值z，将x、y样本的值进行相应的替换：
            xi'= xi - x_mean + z
            yi'= yi - y_mean + z
        '''
        conbimed = x_data + y_data
        z = np.mean(conbimed)
        x_mean = np.mean(x_data)
        x_std = np.std(x_data)
        y_mean = np.mean(y_data)
        y_std = np.std(y_data)
        #计算 observed数据的t值
        t_obs = Boostrap_test.calc_t_val(x_mean, y_mean, len(x_data), len(y_data), x_std, y_std)

        x_data = map(lambda i: i - x_mean + z, x_data)
        y_data = map(lambda i: i - y_mean + z, y_data)

        return x_data, y_data, t_obs, x_mean, y_mean

    def boostrap(self, x_data, y_data, t_obs):
        '''
        模拟boostrap抽样，每次抽样 size个样本
        :return:
        '''
        x_size = len(x_data)
        y_size = len(y_data)
        for t in range(self.times):
            x_samples = np.random.choice(x_data, x_size)  #choice方法可以模拟又放回的抽取，size的值可以大于len(data)
            y_samples = np.random.choice(y_data, y_size)

            x_samples_mean = np.mean(x_samples)
            x_samples_std = np.std(x_samples)
            y_samples_mean = np.mean(y_samples)
            y_samples_std = np.std(y_samples)

            t_val = Boostrap_test.calc_t_val(x_samples_mean, y_samples_mean, x_size, y_size, x_samples_std, y_samples_std)
            self.t_val_list.append(t_val)
        p_val = self.calc_p_val(t_obs, self.t_val_list)

        p_ref = 1.0 - self.conf
        if p_val < p_ref:
            res = "Significant difference"
        else:
            res = "Insignificant difference"

        return res, p_val, p_ref

    @staticmethod
    def calc_t_val(x_mean, y_mean, x_size, y_size, x_std, y_std):
        '''计算 bootstrap抽样样本的 t值，此处两组样本的t值使用的是不具备方差齐性时的 t'检验公式'''
        numerator = np.abs(x_mean - y_mean)
        denominator = np.sqrt(np.square(x_std)/x_size + np.square(y_std)/y_size)
        t_val = numerator/denominator
        return t_val

    def calc_p_val(self, t_obs, t_boot_list):
        '''通过bootstrap抽样所得的数据求得的t值，与observed数据测得的t值比较，求出p值'''
        count = 0.0
        for t in range(self.times):
            if t_obs < t_boot_list[t]:
                count += 1.0
        p_val = count/self.times
        return p_val

    def main(self, x_data, y_data):
        x_data, y_data, t_obs, x_mean, y_mean = self.preprocess(x_data, y_data)
        res, p_val, p_ref = self.boostrap(x_data, y_data, t_obs)
        return res, p_val, p_ref, x_mean, y_mean


class Significant_difference_test:
    def __init__(self):
        self.mod_list = ["ace", "gly", "mal", "met", "suc", "sumo", "ubi"]
        self.sig_dir = "/data1/hbs/centralOrNotSigIsyes/"
        self.insig_dir = "/data1/hbs/centralOrNotSigNotyes/"
        self.b = Boostrap_test()

    def compare_in_motif(self):
        f = open("/data1/hbs/UCEC_significant_difference/sig_insig_in_motif_comparison(bo).txt", "w")
        content = ""
        for mod in self.mod_list:
            sig_central_mod_file = self.sig_dir + mod + "_in_motif_centralMut.txt"
            sig_central_mod_data = Significant_difference_test.load_data(sig_central_mod_file)
            sig_surround_mod_file = self.sig_dir + mod + "_in_motif_surroundMut.txt"
            sig_surround_mod_data = Significant_difference_test.load_data(sig_surround_mod_file)
            sig_in_motif = pd.concat([sig_central_mod_data, sig_surround_mod_data], axis=0)
            sig_in_motif = pd.DataFrame(sig_in_motif.groupby(["protein_id"])["sample_count"].sum())
            sig_data = sorted(list(sig_in_motif["sample_count"]))
            # sig_data_size = len(sig_data)

            insig_central_mod_file = self.insig_dir + mod + "_in_motif_centralMut.txt"
            insig_central_mod_data = Significant_difference_test.load_data(insig_central_mod_file)
            insig_surround_mod_file = self.insig_dir + mod + "_in_motif_surroundMut.txt"
            insig_surround_mod_data = Significant_difference_test.load_data(insig_surround_mod_file)
            insig_in_motif = pd.concat([insig_central_mod_data, insig_surround_mod_data], axis=0)
            insig_in_motif = pd.DataFrame(insig_in_motif.groupby(["protein_id"])["sample_count"].sum())
            insig_data = sorted(list(insig_in_motif["sample_count"]))
            # insig_data_size = len(insig_data)

            res, p_val, p_ref, sig_mean, insig_mean = self.b.main(sig_data, insig_data)
            if res == "Significant difference":
                content += "%s in %s modification counts in motif area, p value %s < %s\n"%(res, mod, p_val, p_ref)
            else:
                content += "%s in %s modification counts in motif area, p value %s > %s\n"%(res, mod, p_val, p_ref)
            content += "average modification in significant proteins: %s\n"\
                       "average modification in insignificant proteins: %s\n\n"%(sig_mean, insig_mean)
        f.write(content)
        f.close()

    def compare_out_motif(self):
        f = open("/data1/hbs/UCEC_significant_difference/sig_insig_out_motif_comparison(bo).txt", "w")
        content = ""
        for mod in self.mod_list:
            sig_out_motif_file = self.sig_dir + mod + "_out_motifMut.txt"
            insig_out_motif_file = self.insig_dir + mod + "_out_motifMut.txt"
            sig_out_motif_data = Significant_difference_test.load_data(sig_out_motif_file)
            insig_out_motif_data = Significant_difference_test.load_data(insig_out_motif_file)
            sig_data = sorted(list(sig_out_motif_data["sample_count"]))
            # sig_data_size = len(sig_data)
            insig_data = sorted(list(insig_out_motif_data["sample_count"]))
            # insig_data_size = len(insig_data)

            res, p_val, p_ref, sig_mean, insig_mean = self.b.main(sig_data, insig_data)
            if res == "Significant difference":
                content += "%s in %s modification counts out of motif area, p value %s < %s\n"%(res, mod, p_val, p_ref)
            else:
                content += "%s in %s modification counts out of motif area, p value %s > %s\n"%(res, mod, p_val, p_ref)
            content += "average modification in significant proteins: %s\n" \
                       "average modification in insignificant proteins: %s\n\n"%(sig_mean, insig_mean)
        f.write(content)
        f.close()

    def compare_in_motif_central(self):
        f = open("/data1/hbs/UCEC_significant_difference/sig_insig_central_motif_comparison(bo).txt", "w")
        content = ""
        for mod in self.mod_list:
            sig_central_motif_file = self.sig_dir + mod + "_in_motif_centralMut.txt"
            insig_central_motif_file = self.insig_dir + mod + "_in_motif_centralMut.txt"
            sig_central_motif_data = Significant_difference_test.load_data(sig_central_motif_file)
            insig_central_motif_data = Significant_difference_test.load_data(insig_central_motif_file)
            sig_data = sorted(list(sig_central_motif_data["sample_count"]))
            # sig_data_size = len(sig_data)
            insig_data = sorted(list(insig_central_motif_data["sample_count"]))
            # insig_data_size = len(insig_data)

            res, p_val, p_ref, sig_mean, insig_mean = self.b.main(sig_data, insig_data)
            if res == "Significant difference":
                content += "%s in %s modification counts in central of motif area, p value %s < %s\n"%(res, mod, p_val, p_ref)
            else:
                content += "%s in %s modification counts in central of motif area, p value %s > %s\n"%(res, mod, p_val, p_ref)
            content += "average modification in significant proteins: %s\n" \
                       "average modification in insignificant proteins: %s\n\n"%(sig_mean, insig_mean)
        f.write(content)
        f.close()

    def compare_in_motif_surround(self):
        f = open("/data1/hbs/UCEC_significant_difference/sig_insig_surround_motif_comparison(bo).txt", "w")
        content = ""
        for mod in self.mod_list:
            sig_surround_motif_file = self.sig_dir + mod + "_in_motif_surroundMut.txt"
            insig_surround_motif_file = self.insig_dir + mod + "_in_motif_surroundMut.txt"
            sig_surround_motif_data = Significant_difference_test.load_data(sig_surround_motif_file)
            insig_surround_motif_data = Significant_difference_test.load_data(insig_surround_motif_file)
            sig_data = sorted(list(sig_surround_motif_data["sample_count"]))
            # sig_data_size = len(sig_data)
            insig_data = sorted(list(insig_surround_motif_data["sample_count"]))
            # insig_data_size = len(insig_data)

            res, p_val, p_ref, sig_mean, insig_mean = self.b.main(sig_data, insig_data)
            if res == "Significant difference":
                content += "%s in %s modification counts in motif flanking sequence area, p value %s < %s\n"%(res, mod, p_val, p_ref)
            else:
                content += "%s in %s modification counts in motif flanking sequence area, p value %s > %s\n"%(res, mod, p_val, p_ref)
            content += "average modification in significant proteins: %s\n" \
                       "average modification in insignificant proteins: %s\n\n"%(sig_mean, insig_mean)
        f.write(content)
        f.close()

    @staticmethod
    def load_data(f):
        data = pd.read_csv(f, sep="\t", header=None)
        data.columns = ["protein_id", "sample_count"]
        return data

if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    # c = Counting("/data1/hbs/modify_in_motif(sig=yes)", "/data1/hbs/modify_in_motif(signotyes)")
    # c.main()
    # b = Boostrap("C:/Users/hbs/Desktop/gly_in_motif_centralMut.txt", 0.95, 0.2, time=1000)
    # data = b.load_data("C:/Users/hbs/Desktop/gly_in_motif_centralMut.txt")
    # interval, q1, q3 = b.boostrap(data)
    # print(interval)
    # print(q1)
    # print(q3)
    s = Significant_difference_test()
    s.compare_in_motif()
    s.compare_out_motif()
    s.compare_in_motif_central()
    s.compare_in_motif_surround()