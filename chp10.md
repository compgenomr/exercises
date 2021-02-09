## Exercises

### Differential methylation
The main objective of this exercise is getting differential methylated cytosines between two groups of samples: IDH-mut (AML patients with IDH mutations) vs. NBM (normal bone marrow samples).

1. Download methylation call files from GEO. These files are readable by methlKit using default `methRead` arguments. [Difficulty: **Beginner**]

samples     Link      
-------     ------  
IDH1_rep1    [link](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM919990&format=file&file=GSM919990%5FIDH%2Dmut%5F1%5FmyCpG%2Etxt%2Egz) 
IDH1_rep2    [link](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM919991&format=file&file=GSM919991%5FIDH%5Fmut%5F2%5FmyCpG%2Etxt%2Egz) 
NBM_rep1     [link](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM919982&format=file&file=GSM919982%5FNBM%5F1%5FmyCpG%2Etxt%2Egz)
NBM_rep2     [link](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM919984&format=file&file=GSM919984%5FNBM%5F2%5FRep1%5FmyCpG%2Etxt%2Egz)

Example code for reading a file:
```{r exaReadMethExercise, eval=FALSE}
library(methylKit)
m=methRead("~/Downloads/GSM919982_NBM_1_myCpG.txt.gz",
           sample.id = "idh",assembly="hg18")
```

2. Find differentially methylated cytosines. Use chr1 and chr2 only if you need to save time. You can subset it after you download the files either in R or Unix. The files are for hg18 assembly of human genome. [Difficulty: **Beginner**]
3. Describe the general differential methylation trend, what is the main effect for most CpGs? [Difficulty: **Intermediate**]
4. Annotate differentially methylated cytosines (DMCs) as promoter/intron/exon? [Difficulty: **Beginner**]
5. Which genes are the nearest to DMCs? [Difficulty: **Intermediate**]
6. Can you do gene set analysis either in R or via web-based tools? [Difficulty: **Advanced**]



### Methylome segmentation
The main objective of this exercise is to learn how to do methylome segmentation and the downstream analysis for annotation and data integration.

1. Download the human embryonic stem-cell (H1 Cell Line) methylation bigWig files from the [Roadmap Epigenomics website](http://egg2.wustl.edu/roadmap/web_portal/processed_data.html#MethylData). It may take a while to understand how the website is structured and which bigWig file to use. That is part of the exercise. The files you will download are for hg19 assembly unless stated otherwise. [Difficulty: **Beginner**]
2. Do segmentation on hESC methylome. You can only use chr1 if using the whole genome takes too much time. [Difficulty: **Intermediate**]
3. Annotate segments and the kinds of gene-based features each segment class overlaps with (promoter/exon/intron). [Difficulty: **Beginner**]
4. For each segment type, annotate the segments with chromHMM annotations from the Roadmap Epigenome database available [here](https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html#core_15state). The specific file you should use is [here](https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E003_15_coreMarks_mnemonics.bed.gz). This is a bed file with chromHMM annotations. chromHMM annotations are parts of the genome identified by a hidden-Markov-model-based machine learning algorithm. The segments correspond to active promoters, enhancers, active transcription, insulators, etc. The chromHMM model uses histone modification ChIP-seq and potentially other ChIP-seq data sets to annotate the genome.[Difficulty: **Advanced**]
