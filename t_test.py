# -*- coding:utf-8 -*-
import numpy as np
from scipy import stat


class T_test:
    def __init__(self, miu1, sigma1, miu2, sigma2, size1, size2, alpha=0.05):
        self.mean1 = miu1
        self.std1 = sigma1
        self.mean2 = miu2
        self.std2 = sigma2
        self.size1 = size1
        self.size2 = size2
        self.df = self.size1 + self.size2 - 2
        self.alpha = alpha

    def homogeneity_of_variance_test(self):
        variance1 = np.square(self.std1)
        variance2 = np.square(self.std2)
        if variance1 > variance2:
            f_val = variance1 / variance2
        else:
            f_val = variance2 / variance1
        return f_val

    def t_test(self):
        if self.size1 < 60 or self.size2 < 60:
            #本次boostrap模拟的数据样本都大于60
            # f_val = self.homogeneity_of_variance_test()
            # f_standard = stats.f.cdf()
            pass
        else:
            s_combine = ((self.size1 - 1) * np.square(self.std1) + (self.size2 - 1) * np.square(self.std2))/self.df
            s = np.sqrt(s_combine * (1/self.size1 + 1/self.size2))
            t_val = np.abs(self.mean1 - self.mean2)/s
            p_val = stats.t.cdf(t_val, self.df)
        if p_val < self.alpha:
            res = "Significant difference"
        else:
            res = "Insignificant difference"
        return [res, t_val, p_val]