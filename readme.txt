数据分析流程：
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
    ② 将Total.elm（fastr序列文件在total_fastr目录下），将所有的fastr转化为fasta格式如下，
--------------------------------------------------------------------------------------------------------
>PLMD-2@O00116@Ubiquitination@Homo sapiens@102
MAEAAAAAGGTGLGAGASYGSAADRDRDPDPDRAGRRLRVLSGHLLGRPREALSTNECKARRAASAATAAPTATPAAQESGTIPKKRQEVMKWNGWGYNDSKFIFNKKGQIELTGKRYPLSGMGLPTFKEWIQNTLGVNVEHKTTSKASLNPSDTPPSVVNEDFLHDLKETNISYSQEADDRVFRAHGHCLHEIFLLREGMFERIPDIVLWPTCHDDVVKIVNLACKYNLCIIPIGGGTSVSYGLMCPADETRTIISLDTSQMNRILWVDENNLTAHVEAGITGQELERQLKESGYCTGHEPDSLEFSTVGGWVSTRASGMKKNIYGNIEDLVVHIKMVTPRGIIEKSCQGPRMSTGPDIHHFIMGSEGTLGVITEATIKIRPVPEYQKYGSVAFPNFEQGVACLREIAKQRCAPASIRLMDNKQFQFGHALKPQVSSIFTSFLDGLKKFYITKFKGFDPNQLSVATLLFEGDREKVLQHEKQVYDIAAKFGGLAAGEDNGQRGYLLTYVIAYIRDLALEYYVLGESFETSAPWDRVVDLCRNVKERITRECKEKGVQFAPFSTCRVTQTYDAGACIYFYFAFNYRGISDPLTVFEQTEAAAREEILANGGSLSHHHGVGKLRKQWLKESISDVGFGMLKSVKEYVDPNNIFGNRNLL
>PLMD-2@O00116@Acetylation@Homo sapiens@107
MAEAAAAAGGTGLGAGASYGSAADRDRDPDPDRAGRRLRVLSGHLLGRPREALSTNECKARRAASAATAAPTATPAAQESGTIPKKRQEVMKWNGWGYNDSKFIFNKKGQIELTGKRYPLSGMGLPTFKEWIQNTLGVNVEHKTTSKASLNPSDTPPSVVNEDFLHDLKETNISYSQEADDRVFRAHGHCLHEIFLLREGMFERIPDIVLWPTCHDDVVKIVNLACKYNLCIIPIGGGTSVSYGLMCPADETRTIISLDTSQMNRILWVDENNLTAHVEAGITGQELERQLKESGYCTGHEPDSLEFSTVGGWVSTRASGMKKNIYGNIEDLVVHIKMVTPRGIIEKSCQGPRMSTGPDIHHFIMGSEGTLGVITEATIKIRPVPEYQKYGSVAFPNFEQGVACLREIAKQRCAPASIRLMDNKQFQFGHALKPQVSSIFTSFLDGLKKFYITKFKGFDPNQLSVATLLFEGDREKVLQHEKQVYDIAAKFGGLAAGEDNGQRGYLLTYVIAYIRDLALEYYVLGESFETSAPWDRVVDLCRNVKERITRECKEKGVQFAPFSTCRVTQTYDAGACIYFYFAFNYRGISDPLTVFEQTEAAAREEILANGGSLSHHHGVGKLRKQWLKESISDVGFGMLKSVKEYVDPNNIFGNRNLL
--------------------------------------------------------------------------------------------------------

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

    ⑤ counting.py程序处理①中产生的文件统计各个修饰motif的样本数目，使用如下规则去除样本重复和motif重复
       当同一个样本因为motif和序列不同归类在 sig=yes和sig≠yes中时
            motif优先级 sig=yes -> sig≠ yes
       当相同的样本因为选取的蛋白序列不同分到 central和surround 两类时，优先保留 central的
            修饰位置优先级 central -> surround -> out


    ⑥ ifInDomain.py程序将所有突变位点定位到 in_domain和 out_domain中
