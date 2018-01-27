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
==========================
@update: 纠正了其中一处语法错误
@version: 1.2
==========================
@update: 加了更详尽的注释，因为实际项目中出现易混淆点
@version: 1.3
==========================
@update: 对bootstrap的均值、方差的计算做了较大程度上的变化
@version: 2.0
==========================
@update: 纠正一处关键的拼写错误
@version: 2.1.1
==========================
@update: 使用矢量计算代替了之前的一处线性计算，提高效率
@version: 2.1.2
==========================
@update: 增加了对初始输入值的判断，传入空值的情况直接报错并退出。
@version: 2.1.3
'''
import numpy as np
# from matplotlib import pyplot as plt
import os
import warnings
# from t_test import T_test


class Bootstrap_test:

    def __init__(self, Confidence=0.95, time=1000, side="two-side"):
        '''
        :param Confidence: 置信区间
        :param time: boostrap抽样的次数
        :param side: 默认是two-side的双尾检验，如果是单尾检验则为"one-side"
        '''
        self.conf = Confidence
        self.times = time
        self.side = side
        self.t_val_list = []
        #经过self.time次抽样随机抽取到的样本，也就是模拟得到的总体
        self.x_total_samples = np.array([])
        self.y_total_samples = np.array([])

    def preprocess(self, x_data, y_data):
        '''
        求取x、y两组观测值的 t'值
        :param x_data: 参数类型是列表，由样本组一的样本值组成
        :param y_data: 参数类型是列表，由样本组二的样本值组成
        '''
        if len(x_data) == 0 | len(y_data) == 0:
            print("invalid input, at least one empty set in two datasets.")
            exit()
        #未经处理的x样本组数据均值
        x_mean = np.mean(x_data)
        x_std = np.std(x_data)
        #未经处理的y样本组数据均值
        y_mean = np.mean(y_data)
        y_std = np.std(y_data)
        #计算 observed数据的t值
        if x_mean < y_mean:
            t_obs = Bootstrap_test.calc_t_val(x_mean, y_mean, len(x_data), len(y_data), x_std, y_std, self.side)
        elif x_mean >= y_mean:
            t_obs = Bootstrap_test.calc_t_val(y_mean, x_mean, len(y_data), len(x_data), y_std, x_std, self.side)
        # print(t_obs)
        return t_obs

    def bootstrap(self, x_data, y_data, t_obs):
        '''
        模拟boostrap抽样，每次抽样 size个样本（有放回抽样）
        求出x、y样本抽样后形成的新的样本的合并均值z，将x、y新样本的值进行相应的替换：
            xi'= xi - x_mean + z
            yi'= yi - y_mean + z
        :param x_data: 参数类型是列表，由样本组一的样本值组成
        :param y_data: 参数类型是列表，由样本组二的样本值组成
        :param t_obs: 是一个数值，是preprocess方法中计算出来的两个样本组原始值的 t检验值
        :return:
        '''
        # z是两组样本中抽取数据共同的均值
        conbimed = np.append(x_data, y_data)
        z = np.mean(conbimed)
        x_mean = np.mean(x_data)
        y_mean = np.mean(y_data)

        x_size = len(x_data)
        y_size = len(y_data)
        for t in range(self.times):
            x_samples = np.random.choice(x_data, x_size)  #choice方法可以模拟有放回的抽取，size的值等于len(data)
            y_samples = np.random.choice(y_data, y_size)

            self.x_total_samples = np.append(self.x_total_samples, x_samples)
            self.y_total_samples = np.append(self.y_total_samples, y_samples)
            # x样本组抽取的样本数据预处理,下面两行语句中，后者比前者高效是矢量计算
            # x_samples = list(map(lambda i: i - x_mean + z, x_samples))
            x_samples = x_samples - x_mean + z
            x_samples_mean = np.mean(x_samples)
            # y样本组抽取的样本数据预处理，下面两行语句中，后者比前者高效是矢量计算
            # y_samples  = list(map(lambda i: i - y_mean + z, y_samples))
            y_samples = y_samples - y_mean + z
            y_samples_mean = np.mean(y_samples)
            x_samples_std = np.std(x_samples)
            y_samples_std = np.std(y_samples)

            t_val = Bootstrap_test.calc_t_val(x_samples_mean, y_samples_mean, x_size, y_size, x_samples_std, y_samples_std, self.side)
            self.t_val_list.append(t_val)
        # print(self.t_val_list)
        p_val = self.calc_p_val(t_obs, self.t_val_list)

        #定义双尾检验时 p为 α/2；单尾检验时 p为α
        if self.side == "one-side":
            p_ref = 1 - self.conf
        elif self.side == "two-side":
            p_ref = (1 - self.conf)/2
        if p_val < p_ref:
            res = "Significant difference"
        else:
            res = "Insignificant difference"
        #求取x和y样本组经过self.time次抽取后形成的总体的均值和标准差
        x_total_mean = np.mean(self.x_total_samples)
        x_total_std = np.std(self.x_total_samples)
        y_total_mean = np.mean(self.y_total_samples)
        y_total_std = np.std(self.y_total_samples)

        return res, p_val, p_ref, x_total_mean, y_total_mean, x_total_std, y_total_std

    @staticmethod
    def calc_t_val(x_mean, y_mean, x_size, y_size, x_std, y_std, side):
        '''计算 bootstrap抽样样本的 t值，此处两组样本的t值使用的是不具备方差齐性时的 t'检验公式
        :param x_mean: 经过预处理后样本组一的均值
        :param y_mean: 经过预处理后样本组二的均值
        :param x_size: 样本组一的样本量
        :param y_size: 样本组二的样本量
        :param x_std: 经过预处理后样本组一的标准差
        :param y_std: 经过预处理后样本组二的标准差
        '''
        if side == "one-side":
            numerator = x_mean - y_mean
        elif side == "two-side":
            numerator = np.abs(x_mean - y_mean)
        else:
            print("invalid input")
            exit()
        denominator = np.sqrt(np.square(x_std)/x_size + np.square(y_std)/y_size)
        t_val = numerator/denominator
        return t_val

    def calc_p_val(self, t_obs, t_boot_list):
        '''通过bootstrap抽样所得的数据求得的t值，与observed数据测得的t值比较，求出p值
            如果置信区间取95%，当H0正确时，拒绝零假设所犯的第一类错误α的值为0.05；
        当t_boot<t_obs,说明拒绝零假设时犯第一类错误的概率α>0.05，求取抽样的
        所有次数中犯第一类错误的次数占总次数的比值，求取的实际上就是bootstrap
        抽样后的α值。
        :param t_obs: 是preprocess方法中计算出来的两个样本组原始值的 t检验值
        :param t_boot_list: 是一个列表，是经过 self.time次抽样后得到的self.time个 t'值组成的列表
        '''
        count = 0.0
        for t in range(self.times):
            if t_obs > t_boot_list[t]:  #t_boot<t_obs表示拒绝H0时会犯第一类错误
                count += 1.0
        p_val = count/self.times  #self.time次中有 count次拒绝H0后发生第一类错误的概率
        return p_val

    def main(self, x_data, y_data):
        '''将之前所有方法进行整合形成一套完整的bootstrap检验流程
        :param x_data: 参数类型是list，是样本组一的数据值组成的列表
        :param y_data: 参数类型是list，是样本组二的数据值组成的列表
        :return res是最终的检验结果 "Significant difference"或是"Insignificant difference"；
                p_val是发生第一类错误的概率； p_ref是α值，即 1-self.conf；
                x_mean是样本组x的均值； y_mean是样本组y的均值；
                x_std是样本组x的标准差；y_std是样本组y的标准差。
        '''
        x_data = np.array(x_data)
        y_data = np.array(y_data)
        t_obs = self.preprocess(x_data, y_data)
        res, p_val, p_ref, x_mean, y_mean, x_std, y_std = self.bootstrap(x_data, y_data, t_obs)
        return res, p_val, p_ref, x_mean, y_mean, x_std, y_std


if __name__ == "__main__":
    warnings.filterwarnings("ignore")