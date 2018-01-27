# -*- coding:utf-8 -*-
import os
from domain_bootstrap import DomainBootstrapTest
import multiprocessing as mp


if __name__ == "__main__":

    def bootstrap_test(in_domain_file, out_domain_file, cancer, mod_type, output_doc):
        d = DomainBootstrapTest(in_domain_file, out_domain_file,mod_type)
        content = d.main()
        if content.find("higher") != -1:
            print(cancer + ":" + mod_type)
        with open(output_doc, "w") as f:
            f.write(content)
            f.close()
        return None

    pool = mp.Pool(7)

    cancer_kind = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP',
                   'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM',
                   'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']
    # cancer_kind = ["PCPG"]

    calc_dir = "/data1/hbs/all_cancer_analysis/all_cancer_modification_in_domain_or_not/seven_modify_test_result/"
    result_dir = "/data1/hbs/all_cancer_analysis/all_cancer_modification_in_domain_or_not/seven_modify_test_result/"

    for c in cancer_kind:
        in_domain_file = calc_dir + c + "/" + "in_domain_calc.txt"
        out_domain_file = calc_dir + c + "/" + "out_domain_calc.txt"
        mod_type = "all combined 7 modification"
        output_doc = calc_dir + c + "/" + "significant_test.txt"
        # bootstrap_test(in_domain_file, out_domain_file, c, mod_type, output_doc)
        pool.apply_async(bootstrap_test, (in_domain_file, out_domain_file, c, mod_type, output_doc))
    pool.close()
    pool.join()

    #测试代码
    # bootstrap_test("C:/Users/hbs/Desktop/domain_test/in_domain_calc.txt", "C:/Users/hbs/Desktop/domain_test/out_domain_calc.txt","acc","all", "")