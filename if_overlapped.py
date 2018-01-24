'''
Description:
    用于检验两个序列是否有重叠。
===========================
@author:hbs
@date:2018-1-20
@version:1.1
'''
# -*- coding:utf-8 -*-


def if_overlapped(interval1, interval2):
    '''判断同一个蛋白上两个domain区间是否重叠
    :param interval1: 传入第一个区间
    :param interval2: 传入第二个区间
    '''
    overlapped = False
    # 相交的情况1:       p1 --------------------.min([interval1[1],interval2[1]])
    # max([[interval1[0],interval2[0]]) .------------ p2
    if max([interval1[0], interval2[0]]) <= min([interval1[1], interval2[1]]):
        overlapped = True

    # 相交的情况2:   不满足 p1 -------         或               p1 --------的情况
    #                              p2 ------     --------- p2
    elif not (interval1[0] > interval2[1] or interval1[1] < interval2[0]):
        overlapped = True
    return overlapped