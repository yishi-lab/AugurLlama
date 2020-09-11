# AugurLlama: a tool for proteomics based Nanobody discovery
This tool is designed to assist the automated proteomic discovery of diverse nanobody repertoires. 

The tool was developed by Zhe Sang, Yufei Xiang and Yi Shi. The rule/algorithm in this tool was summarized frmo manual inspectino of MS data and exploration of features of nanobody sequence database. The tool is still in the very early stage and not well optimized. More integrated and sophisticated versioon is currently under development.


The input files are:

1) PSM table(excel format)  
2) annotated spectra(txt)
3) MS raw file
4) Sequence database
where 1) and 2) are exported from Proteome Discovery.

Workflow of the tool:

 1) In-silico digestion of the database
 2) Annotation of nanobody sequence(idenfication of FR annd CDR regions)
 3) Reading PSMs and removal of low-quality CDR3 peptides
 4) Quantification of CDR3(the fingerprint)
 5) Affinity group classification(0:low, 1:mediocre, 2: high)

Requirements:

1) Windows operating system
2) pymsfilereader(https://github.com/frallain/pymsfilereader)
3) Biopython,Scipy,Numpy,Pandas

Contact: Zhe Sang, sang.zhe(at)pitt.edu
