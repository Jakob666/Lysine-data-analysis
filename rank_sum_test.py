#-*- coding:utf-8 -*-
'''
Description:
    使用SciPy的stats库中kruskal方法和ranksuns方法对
    in_motif,out_motif 的位点保守性分数进行差异显著性检验
    in_motif_in_domain, in_motif_out_domain, out_motif_in_domain和out_motif_out_domain的位点保守性进行差异显著性检验
    因为样本不服从正态分布，且方差不具齐性所以选用非参数的秩和检验，
    两组样本选用Wilcoxon检验方法，多组样本选用Kruskal-Wallis检验
@Author：hbs
@Date：2018-5-11
@version：1.0
'''
import numpy as np
from scipy import stats
import pandas as pd
import os
import warnings


class RankSumTest:
    def related_sites(self):
        '''
        从目录中提取赖氨酸修饰相关的in_motif的位点
        :return: 已排序的in_domain, out_domain和全部related位点的一维数组
        '''
        related_in = "D:\\deserve\\conservation_score\\related_in_domain"
        related_out = "D:\\deserve\\conservation_score\\related_out_domain"
        related_in_data = np.concatenate([self.load_data(related_in+"/"+f) for f in os.listdir(related_in)])
        related_out_data = np.concatenate([self.load_data(related_out+"/"+f) for f in os.listdir(related_out)])
        related_data = np.concatenate([related_in_data, related_out_data])
        return np.sort(related_in_data), np.sort(related_out_data), np.sort(related_data)


    def irrelated_sites(self):
        '''
        从目录中提取赖氨酸修饰相关的out_motif的位点
        :return: 已排序的in_domain, out_domain和全部irrelated位点的一维数组
        '''
        irrelated_in = "D:\\deserve\\conservation_score\\irrelated_in_domain"
        irrelated_out = "D:\\deserve\\conservation_score\\irrelated_out_domain"
        irrelated_in_data = np.concatenate([self.load_data(irrelated_in + "/" + f) for f in os.listdir(irrelated_in)])
        irrelated_out_data = np.concatenate([self.load_data(irrelated_out + "/" + f) for f in os.listdir(irrelated_out)])
        irrelated_data = np.concatenate([irrelated_in_data, irrelated_out_data])
        return np.sort(irrelated_in_data), np.sort(irrelated_out_data), np.sort(irrelated_data)


    def load_data(self, file):
        '''
        将文件中的PhastCons数据读入Series
        :param file: 文件的绝对路径
        :return: 一维数组
        '''
        data = pd.read_csv(file, sep="\t", usecols=[3])
        try:
            emp = data[data["phastCons100way_vertebrate"] == "emp"]
            non = data[data["phastCons100way_vertebrate"] == "None"]
            data.drop(emp.index, inplace=True)
            data.drop(non.index, inplace=True)
        except TypeError as e:
            pass
        except ValueError:
            print(file)
            print(data.head())
        # 拉直成一维数组并返回
        data = data.values
        return np.ravel(data.astype("float"))


    def rank_sum_test(self):
        related_in, relate_out, related = self.related_sites()
        irrelated_in, irrelate_out, irrelated = self.irrelated_sites()
        print(stats.ranksums(related, irrelated))
        print(stats.kruskal(related_in, relate_out, irrelated_in, irrelate_out))


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    r = RankSumTest()
    r.rank_sum_test()
