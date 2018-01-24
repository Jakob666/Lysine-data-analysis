#-*- coding:utf-8 -*-
'''
Description:
    记录annovar产生的注释数据（包括蛋白质突变保守型打分）的文件目录：
        ① /data1/data/zengyanru/LysineTCGA/organ_cancerabbr/tcga
        ② /data1/data/zengyanru/LysineTCGA/organ_cancerabbr/icgc
        ③ /data1/data/zengyanru/LysineTCGA/organ_cancerabbr/cosmic
    存放蛋白质p值注释的文件目录：
        E:/protein_p_val (目前没有上传到服务器）
        文件目录中按照不同癌症及不同的修饰进行分类，每个文件中的蛋白都有相应的p值、对应基因名等信息
    用于作为索引的文件目录（但是由于存在历史遗留问题，一定少用）：
        /data1/data/mnliu/ModifyOut/lysine_ptm
        文件目录中按照不同的修饰进行分类，里面有各种修饰蛋白的 UniprotID、ENG和EST基因

    用上面的几个文件目录生成文件的内容：
        ① signficant protein--> in_domain  out_domain两个文件
        ② insignificant protein--> in domain  out_domain两个文件   （先做出sig的，再做insig的）
    每个 in_domain out_domain文件中存在如下列：
        protein_id  modification  seq  site  from  to  patient  con_score  source(icgc or tcga or cosmic)
===================================================================
@author: hbs
@date: 2018-1-18（由于原始文件有问题，暂时停滞）
@version: 1.0
'''