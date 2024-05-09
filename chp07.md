## Exercises and solutions for Chapter 7

For this set of exercises, we will use the `chip_1_1.fq.bz2` and `chip_2_1.fq.bz2` files from the `QuasR` package. You can reach the folder that contains the files as follows:
```{r seqProcessEx,eval=FALSE}
folder=(system.file(package="QuasR", "extdata"))
dir(folder) # will show the contents of the folder
```
1. Plot the base quality distributions of the ChIP-seq samples `Rqc` package.
**HINT**: You need to provide a regular expression pattern for extracting the right files from the folder. `"^chip"` matches the files beginning with "chip". [Difficulty: **Beginner/Intermediate**]


**solution:**
```{r,echo=FALSE,eval=FALSE}
library(Rqc)
qcRes=rqc(path = folder, pattern = "^chip", openBrowser=FALSE)
rqcCycleQualityBoxPlot(qcRes)

```

2. Now we will trim the reads based on the quality scores. Let's trim 2-4 bases on the 3' end depending on the quality scores. You can use the `QuasR::preprocessReads()` function for this purpose. [Difficulty: **Beginner/Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
# Since preprocessReads() doesn't trim by score, you must look at the QC report above (see previous exercise) to determine number of bases to trim. Looking at the previous plot, the last base (#36) is often low-quality (<20) for both samples, while the last 3 bases are marginal in chip_2_1.fq. So, the truncateEndBases argument was set to 3.

# sample files
infiles <- system.file(package="QuasR", "extdata",
                       c("chip_1_1.fq.bz2","chip_2_1.fq.bz2"))
outfiles <- paste(pattern=c("chip_1_1.","chip_2_1."),"trim3",".fastq",sep="")
unlink(outfiles) # this removes old files, which allows next command to work
preprocessReads(infiles, outfiles, truncateEndBases = 3)

# Alternative Approach: Instead of trimming all reads to the same extent with preprocessReads(), the ShortRead::trimTails() function can trim *individual* reads based on their quality scores.
library(ShortRead)
# Unfortunately, ShortRead:readFastq can't read in .bz2 files!
# readFastq(dirPath = "extdata/", pattern = "chip_._1.fq.bz2")

# but preprocessReads() can read in .fq.bz2 files and output reads into .fastq format
outfiles <- paste(pattern=c("chip_1_1.","chip_2_1."),"fastq",sep="")
unlink(outfiles) # this removes old files, which allows next command to work
preprocessReads(infiles, outfiles) # we want the original reads, so no trimming
# now read each .fastq file into a ShortReadQ object
fq1 = readFastq(outfiles[1]) # here's the first one
fq1 # 2597 reads originally
# use ShortRead::trimTails to trim each read separately; for help, run ?TrimTails
fq1_1 <- trimTails(fq1, 1, "4") # k = 1, a = "4"; trim at first base with quality < 20 
                                # ASCII of "4" is 52, minus 33 equals quality of 19
fq1_1 # 2361 reads after this trimming
# in contrast, ShortRead::trimTailw uses a window to trim; here trim where 2/5 bases <20
fq1_2w <- trimTailw(fq1, 2, "4", 2) # here, window is set to trim where 2/5 bases <20
fq1_2w # 2151 reads after this trimming
 
```

3. Align the trimmed and untrimmed reads using `QuasR` and plot alignment statistics, did the trimming improve alignments? [Difficulty: **Intermediate/Advanced**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
# Mapping/aligning reads to the genome
# IF "extdata" folder is missing from current working directory, uncomment next line to copy example data
# file.copy(system.file(package="QuasR", "extdata"), ".", recursive=TRUE)

# path to genome file in fasta format
genomeFile <- "extdata/hg19sub.fa"

# alignment with UNTRIMMED reads
# text file with list of fastq file paths and sample names (for grouping)
sampleFile <- "extdata/samples_chip_single.txt"

# check names of files in this list:
read.delim(sampleFile)
# create alignments 
proj_trim <- qAlign(sampleFile, genomeFile)
proj_trim
# save a quality control report to the "qc_report_untrim.pdf" file using the qQCReport function
qQCReport(proj_trim, "qc_report_untrim.pdf")

# alignment with TRIMMED reads
# duplicate "samples_chip_single.txt" file, rename as "samples_trimmed.txt", move to folder with trimmed files
# and edit FileNames with file paths to trimmed files and sample names (for grouping)
sampleFile <- "samples_trimmed.txt" 
# check that FileNames in this list are for trimmed files:
read.delim(sampleFile)
# create alignments 
proj_trim <- qAlign(sampleFile, genomeFile)
# save a quality control report to the "qc_report_trim.pdf" file using the qQCReport function
qQCReport(proj_trim, "qc_report_trim.pdf")
 
```

*Comparing the qc_report for the sequences trimmed by 3 bases and the qc_report for the untrimmed sequences, this trimming didn’t really improve the alignments significantly. For example, with both the untrimmed and trimmed reads, only 88-90% mapped to the genome (p. 4 of qc_report).*

