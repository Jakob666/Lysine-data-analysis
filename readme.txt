数据分析流程：
=========================== 将突变位点原始数据进行整理 ====================================
    ① PTMSiteSelection.py程序首先生成如下形式的文件，
-------------------------------------------------------------------------------------------------------
*protein_id	Canonical	*cancer_type	*K position	position in motif	*left flank	*right flank	mutated position	*from	*to	sample_num	*direct mutation
P51532-1	P51532	BLCA	1027	2	7	7	1030	G	C	TCGA-DK-AA77	no
P51532-1	P51532	BLCA	1027	2	7	7	1033	K	N	TCGA-CU-A3KJ	no
P51532-1	P51532	BLCA	1027	2	7	7	1030	G	C	TCGA-DK-AA77	no
P51532-1	P51532	BLCA	1027	2	7	7	1033	K	N	TCGA-CU-A3KJ	no
O94832	O94832	BLCA	758	1	7	7	757	G	V	TCGA-G2-A3VY	no
O94832	O94832	BLCA	758	1	7	7	754	R	Q	TCGA-FD-A3SP	no
O94832	O94832	BLCA	758	1	7	7	757	G	V	TCGA-G2-A3VY	no
O94832	O94832	BLCA	758	1	7	7	754	R	Q	TCGA-FD-A3SP	no
Q00973-1	Q00973	BLCA	95	1	7	7	89	L	F	TCGA-DK-A6AW	no
--------------------------------------------------------------------------------------------------------
这是ICGC、TCGA和COSMIC文件整理得到的，重要的信息列用*标记，其中 
protein_id是对应蛋白名称（细致到亚形）；
cancer_type是要针对某一种癌症进行分析；
K position是赖氨酸的位置
left flank是中心位点左侧序列长度
right flank是中心位点右侧序列长度
from、to是从什么aa突变到什么aa
direct mutation表示是否为中心突变
存放位置 /data1/hbs/all_cancer_analysis/all_cancer_sig_prot_motif_modify/  和  /data1/hbs/all_cancer_analysis/all_cancer_insig_prot_motif_modify/

    ② 将Total.elm（fastr序列文件在total_fastr目录下），将所有的fastr转化为fasta格式如下，
--------------------------------------------------------------------------------------------------------
>PLMD-2@O00116@Ubiquitination@Homo sapiens@102
MAEAAAAAGGTGLGAGASYGSAADRDRDPDPDRAGRRLRVLSGHLLGRPREALSTNECKARRAASAATAAPTATPAAQESGTIPKKRQEVMKWNGWGYNDSKFIFNKKGQIELTGKRYPLSGMGLPTFKEWIQNTLGVNVEHKTTSKASLNPSDTPPSVVNEDFLHDLKETNISYSQEADDRVFRAHGHCLHEIFLLREGMFERIPDIVLWPTCHDDVVKIVNLACKYNLCIIPIGGGTSVSYGLMCPADETRTIISLDTSQMNRILWVDENNLTAHVEAGITGQELERQLKESGYCTGHEPDSLEFSTVGGWVSTRASGMKKNIYGNIEDLVVHIKMVTPRGIIEKSCQGPRMSTGPDIHHFIMGSEGTLGVITEATIKIRPVPEYQKYGSVAFPNFEQGVACLREIAKQRCAPASIRLMDNKQFQFGHALKPQVSSIFTSFLDGLKKFYITKFKGFDPNQLSVATLLFEGDREKVLQHEKQVYDIAAKFGGLAAGEDNGQRGYLLTYVIAYIRDLALEYYVLGESFETSAPWDRVVDLCRNVKERITRECKEKGVQFAPFSTCRVTQTYDAGACIYFYFAFNYRGISDPLTVFEQTEAAAREEILANGGSLSHHHGVGKLRKQWLKESISDVGFGMLKSVKEYVDPNNIFGNRNLL
>PLMD-2@O00116@Acetylation@Homo sapiens@107
MAEAAAAAGGTGLGAGASYGSAADRDRDPDPDRAGRRLRVLSGHLLGRPREALSTNECKARRAASAATAAPTATPAAQESGTIPKKRQEVMKWNGWGYNDSKFIFNKKGQIELTGKRYPLSGMGLPTFKEWIQNTLGVNVEHKTTSKASLNPSDTPPSVVNEDFLHDLKETNISYSQEADDRVFRAHGHCLHEIFLLREGMFERIPDIVLWPTCHDDVVKIVNLACKYNLCIIPIGGGTSVSYGLMCPADETRTIISLDTSQMNRILWVDENNLTAHVEAGITGQELERQLKESGYCTGHEPDSLEFSTVGGWVSTRASGMKKNIYGNIEDLVVHIKMVTPRGIIEKSCQGPRMSTGPDIHHFIMGSEGTLGVITEATIKIRPVPEYQKYGSVAFPNFEQGVACLREIAKQRCAPASIRLMDNKQFQFGHALKPQVSSIFTSFLDGLKKFYITKFKGFDPNQLSVATLLFEGDREKVLQHEKQVYDIAAKFGGLAAGEDNGQRGYLLTYVIAYIRDLALEYYVLGESFETSAPWDRVVDLCRNVKERITRECKEKGVQFAPFSTCRVTQTYDAGACIYFYFAFNYRGISDPLTVFEQTEAAAREEILANGGSLSHHHGVGKLRKQWLKESISDVGFGMLKSVKEYVDPNNIFGNRNLL
--------------------------------------------------------------------------------------------------------
文件存放位置 /data1/hbs/total_fastr

    ③ interproscan处理
    将全部的fastr文件转为fasta文件（存放在total_fastr目录下），将所有的序列放入interproscan中预测蛋白上的
domain，生成到 Total_interproscan.tsv文件中（该文件也在total_fastr目录下）。
interproscan处理后的文件Total_interproscan.tsv文件中内容形式为，
--------------------------------------------------------------------------------------------------------
*title										*seq_len 			*start	*end
PLMD-3268@P04908@Succinylation@Homo	eb9a1e9f1c207300757bcd20505f7ab5	130	SMART	SM00414		3	123	9.1E-81	T	11-01-2018	IPR002119	Histone H2A	GO:0000786|GO:0003677|GO:0005634
PLMD-3268@P04908@Succinylation@Homo	eb9a1e9f1c207300757bcd20505f7ab5	130	Pfam	PF00125	Core histone H2A/H2B/H3/H4	5	89	3.6E-17	T	11-01-2018	IPR007125	Histone H2A/H2B/H3	GO:0003677
PLMD-3268@P04908@Succinylation@Homo	eb9a1e9f1c207300757bcd20505f7ab5	130	Pfam	PF16211	C-terminus of histone H2A	92	126	1.5E-20	T	11-01-2018	IPR032454	Histone H2A, C-terminal domain	
PLMD-11001@P50502@Succinylation@Homo	fb2fc9e524790584f6b64ff7fef621b2	369	SMART	SM00727		319	358	3.2E-10	T	11-01-2018	IPR006636	Heat shock chaperonin-binding	
PLMD-11001@P50502@Succinylation@Homo	fb2fc9e524790584f6b64ff7fef621b2	369	SMART	SM00028		114	147	30.0	T	11-01-2018	IPR019734	Tetratricopeptide repeat	GO:0005515
PLMD-11001@P50502@Succinylation@Homo	fb2fc9e524790584f6b64ff7fef621b2	369	SMART	SM00028		148	181	8.9E-4	T	11-01-2018	IPR019734	Tetratricopeptide repeat	GO:0005515
PLMD-11001@P50502@Succinylation@Homo	fb2fc9e524790584f6b64ff7fef621b2	369	SMART	SM00028		182	215	64.0	T	11-01-2018	IPR019734	Tetratricopeptide repeat	GO:0005515
PLMD-7767@P30101@Succinylation@Homo	1f72603c28a716eea7cb211291fdc8cd	505	Pfam	PF13848	Thioredoxin-like domain	160	355	5.3E-28	T	11-01-2018
PLMD-7767@P30101@Succinylation@Homo	1f72603c28a716eea7cb211291fdc8cd	505	Pfam	PF00085	Thioredoxin	27	130	1.2E-33	T	11-01-2018	IPR013766	Thioredoxin domain	GO:0045454
PLMD-7767@P30101@Succinylation@Homo	1f72603c28a716eea7cb211291fdc8cd	505	Pfam	PF00085	Thioredoxin	378	482	2.5E-31	T	11-01-2018	IPR013766	Thioredoxin domain	GO:0045454
PLMD-46924@D6RBL5@Succinylation@Homo	f08e10a58bddc7cc68cced4e31d3506d	260	SMART	SM00335		203	255	1.2E-22	T	11-01-2018	IPR018502	Annexin repeat	GO:0005509|GO:0005544
PLMD-46924@D6RBL5@Succinylation@Homo	f08e10a58bddc7cc68cced4e31d3506d	260	SMART	SM00335		44	96	1.3E-23	T	11-01-2018	IPR018502	Annexin repeat	GO:0005509|GO:0005544
PLMD-46924@D6RBL5@Succinylation@Homo	f08e10a58bddc7cc68cced4e31d3506d	260	SMART	SM00335		128	180	4.7E-15	T	11-01-2018	IPR018502	Annexin repeat	GO:0005509|GO:0005544
PLMD-46924@D6RBL5@Succinylation@Homo	f08e10a58bddc7cc68cced4e31d3506d	260	Pfam	PF00191	Annexin	191	255	3.4E-23	T	11-01-2018	IPR018502	Annexin repeat	GO:0005509|GO:0005544
--------------------------------------------------------------------------------------------------------
重要的field已经使用*标记，其中
title：要通过Uniprot Accession与上一步生成的motif中的条目进行匹配；
seq_len：序列的长度
start、end：domain的起始位点和终止位点，这两个位点在最终的显著性差异检验的时候要将overlap的区域进行合并；
存放位置 /data1/hbs/all_cancer_analysis/all_domain_classify_by_modify/

================================ 数据的基本处理工作结束 ===========================================

====================== 对domain区和非domain区域的突变数目进行显著性检验（bootstrap） ========================
    ⑤ counting.py程序处理①中产生的文件统计各个修饰motif的样本数目，使用如下规则去除样本重复和motif重复
       当同一个样本因为motif和序列不同归类在 sig=yes和sig≠yes中时
            motif优先级 sig=yes -> sig≠ yes
       当相同的样本因为选取的蛋白序列不同分到 central和surround 两类时，优先保留 central的
            修饰位置优先级 central -> surround -> out
···最终形成
	insignificant/significant protein
        |-----------> 某种cancer -->某修饰 ------->in motif mutation----->central mutation
        |                       	   |                       |----->surround mutation
        |                       	   |------>out motif mutation
        |
        |-----------> 某种cancer -->某修饰 ------->in motif mutation----->central mutatu=ion
                               		  |                        |----->surround mutation
                               		  |------->out motif mutation
文件的内容形式为：
---------------------------------------------------------------------------------------------------------
protein_id	Canonical	K position	left flank	right flank	mutated position	from	to	sample_num	sample_count
C9JN71		C9JN71		   352		    7		    7			347		  A	 T	TCGA-PK-A5HB	     1
O95602		O95602		   222		    7		    7			229		  A	 T	TCGA-OR-A5J4	     1
P61006-1	P61006		   3		    7		    7			5		  Y	 C	TCGA-OR-A5KU	     1
Q12778		Q12778		   248		    7		    7			241		  P	 Q	TCGA-OR-A5JA	     1
Q13085-1	Q13085		  2127		    7		    7			2132		  D	 Y	TCGA-OR-A5JA	     1
Q5SW79-1	Q5SW79		   440		    7		    7			439		  G	 W	TCGA-OR-A5KB	     1
Q5TA45-1	Q5TA45		   535		    7		    7			538		  L	 M	TCGA-PK-A5HB	     1
Q7Z460-1	Q7Z460		   570		    7		    7			572		  S	 C	TCGA-OR-A5J4	     1
Q86VS8		Q86VS8		   324		    7		    7			320		  V	 L	TCGA-OR-A5LJ	     1 
Q8N0Z3		Q8N0Z3		   260	 	    7		    7			261		  R	 G	TCGA-OR-A5JV	     1
---------------------------------------------------------------------------------------------------------
存放位置 /data1/hbs/all_cancer_analysis/all_cancer_modification_calculate/sig_prot 和 /data1/hbs/all_cancer_analysis/all_cancer_modification_calculate/insig_prot


    ⑥ ifInDomain.py程序将所有 in motif和 out motif突变位点定位到 in_domain和 out_domain中，这个检测的是某一癌症、某种修饰的显著性蛋白
       在domain 和 非domain区域 每千个氨基酸序列的突变数 是否有显著差异。相应的文件生成在目录下，
       insignificant/significant protein
        |-----------> 某种癌症 ---> 某种修饰------->in domain mutation     
        |                       	   |                       
        |                       	   |------>out domain mutation      (但是此次不检测非显著蛋白，所以只生成了一个文件目录存放显著性蛋白的数据）
        |
        |-----------> 某种癌症 ---> 某种修饰------->in domain mutation
                               		   |                        
                               		   |------->out domain mutation
生成的文件内容形式为
-----------------------------------------------------------------------------------------------------------
Uniprot Accession	protein length	in_domain	non_domain length	sample_count
C9JN71			     531.0	    no			508.0		    1.0
O95602			    1720.0	    no			928.0		    1.0
Q12778			     655.0	    no			564.0		    1.0
Q13085-1		    2346.0	    no			1792.0		    1.0
Q15911-1		    3703.0	    no			3703.0		    1.0
-----------------------------------------------------------------------------------------------------------
文件存放位置   /data1/hbs/all_cancer_analysis/all_cancer_modification_in_domain_or_not/calc_result
检测结果存放于 /data1/hbs/all_cancer_analysis/all_cancer_modification_in_domain_or_not/test_result

    ⑦ ifInDomainV2.py程序是将所有 in motif突变位点定位到 in domain 和 out domain中，不同的是此次不区分不同的修饰，一个癌症产生一个文件。

	|-----------> 某种癌症 ---> 7种修饰------->in domain mutation     
        |                       	   |                       
        |                       	   |------>out domain mutation      (但是此次不检测非显著蛋白，所以只生成了一个文件目录存放显著性蛋白的数据）
        |
        |-----------> 某种癌症 ---> 7种修饰------->in domain mutation
                               		   |                        
                               		   |------->out domain mutation
文件形式与⑥中相同，
文件存放位置   /data1/hbs/all_cancer_analysis/all_cancer_modification_in_domain_or_not/seven_modify_test_result
检测结果存放于 /data1/hbs/all_cancer_analysis/all_cancer_modification_in_domain_or_not/seven_modify_test_result
==================================== bootstrap检验结束 ==================================
 

==================================== 保守性分析开始 ======================================
分析一：
  ⑧ 
    分别将突变分为 显著蛋白 in domain、显著蛋白 out domain、非显著蛋白 in domain和非显著蛋白 out domain
四种，分别存放在以下四个文件目录中，这个过程使用的是/data1/hbs/all_cancer_analysis/all_cancer_mutation_site_conservation/ifInDomain.py
	/data1/hbs/all_cancer_analysis/all_cancer_mutation_site_conservation/all_cancer_in_domain_mut_site/
	/data1/hbs/all_cancer_analysis/all_cancer_mutation_site_conservation/all_cancer_out_domain_mut_site/
    

   ⑨ 之后的工作移步至 222.200.186.117，端口号为27的服务器上，
将蛋白质突变位点匹配到相应的基因组坐标上（annovar获得），一系列文件存放于
    /data/zengyanru/LysineTCGA/find_mutsite_wholegenome/find_extra_score 文件目录下
生成的文件形式为，
--------------------------------------------------------------------------------------------------------------
Cancer	protein	mut_in_out	ENST	Isoform	    ptm_type	ptm_site	ptm_motif	mutation_site	 Chr	 Start	 phyloP100way_vertebrate	phyloP100way_vertebrate_rankscore	phastCons100way_vertebrate	phastCons100way_vertebrate_rankscore
ACC	P68431	  in	ENST00000621411	  nan	  Succinylation	  123	      KRVTIMPKDIQLA	   D124H	  6	26031691	7.62				0.821					1.0					0.715
ACC	P00338	  in	ENST00000379412	P00338-1  Succinylation	  118	      QRNVNIFKFIIPN	   V114M	  11	18400932	7.968				0.874					1.0					0.715
ACC	P51649	  in	ENST00000357578	P51649-1  Succinylation	  365	      KAFAEAMKKNLRVG	   V370A	  6	24522861	5.889				0.693					1.0					0.715
ACC	P0CG47	  in	ENST00000395837	  nan	  Succinylation	  124	      QRLIFAGKQLEDG	   D128G	  17	16382290	8.273				0.897					1.0					0.715										  （mutation site的位置）
-----------------------------------------------------------------------------------------------------------------
再将⑧中得到的四类突变点与其对应的保守性分值相匹配，得到尽可能多的突变点的保守性得分。

    ⑩ 将一个癌症的 sig prot in domain突变 和 sig prot out domain突变
		    insig prot in domain突变 和 insig prot out domain突变 四种位点的保守性打分绘制CDF图
       这四种位点绘制在一个图中作为参照。

分析二：
    由于domain区域本身保守性就比非 domain区域的高，所以改变方案，只比较显著性蛋白赖氨酸相关突变（in motif）和 赖氨酸无关突变（out motif）的 in domain和 out domain保守性分数

    所以需要整理的是赖氨酸无关突变的 in domain和 out domain保守分数
    之前做过的 in motif 和 out motif 无重复的文件统一放在 /data1/hbs/all_cancer_analysis/all_cancer_modification_calculate/sig_prot目录中，从这里提取即可，
    将最后的结果存放在
	/data1/hbs/all_cancer_analysis/all_cancer_mutation_site_conservation/sig_prot_lysine_irrelated_mut_site/
	/data1/hbs/all_cancer_analysis/all_cancer_mutation_site_conservation/sig_prot_lysine_related_mut_site/
    这两个文件目录中。两个目录的分级相同
	related/irralated ----> in domain ----> cancer_name -----> mod_type_in_domain.txt
			   |
			   |--> out domain ----> cancer_name -----> mode_type_out_domain.txt

===================================== 保守性分析工作完成 ======================================

========================================= 挑选癌症突变频次高的蛋白及其特征分析==================================
    ? 根据文件 all_sig_pro1.txt中的数据，挑选出在某个癌症中出现两次及以上的蛋白，该文件的形式为
-----------------------------------------------------------------------------------------------
Q0VDD8-1	4	[['LUAD_sumo_all_tes_result', '0.05'], ['LUAD_ubi_all_tes_result', '0.04885897435897436'], ['SKCM_sumo_all_tes_result', '0.0'], ['STAD_ubi_all_tes_result', '0.042928571428571434']]	DNAH14	F	F
Q9NPG3-1	3	[['CHOL_suc_all_tes_result', '0.0'], ['ESCA_sumo_all_tes_result', '0.0'], ['KIRC_suc_all_tes_result', '0.0']]	UBN1	F	F
P38606-1	1	[['BLCA_gly_all_tes_result', '0.0']]	ATP6V1A	F	F
------------------------------------------------------------------------------------------------
比如文件的第一行，Q0VDD8-1 这个蛋白在LUAD这个癌症中出现两次，将这种蛋白视为高频蛋白；而像Q9NPG3-1和P38606-1即使出现在癌症中，但是同一个癌症中没有出现多次，则视为低频蛋白。
选出这些高频蛋白的突变点在 in motif和 out motif中出现次数（当同一个位点因为选取的序列不同，同时出现在in motif 和 out motif中时，保留 in motif的记录）

    ? 高频蛋白计算 motif区域 和 非motif区域的长度  












