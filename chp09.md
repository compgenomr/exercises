## Exercises and solutions for Chapter 9

### Quality control

1. Apply the fragment size estimation procedure to all ChIP and Input available datasets. [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
# Some objects must be prepared before starting this exercise. First, get the paths to the data files of interest:
# get path to all chip-seq datasets and subsets by file format
data_path <- system.file('extdata/chip-seq',package='compGenomRData') 
chip_files <- list.files(data_path, full.names=TRUE)
bam_files <- list.files(data_path, full.names=TRUE, pattern='bam$')
bw_files <- list.files(data_path, full.names=TRUE, pattern='bw$')

# Next, get human genome data, focusing on chr21, since all files are limited to this chromosome.
library(GenomeInfoDb) # load the chromosome info package
hg_chrs <- getChromInfoFromUCSC('hg38') # get human chromosome lengths
hg_chrs <- subset(hg_chrs, grepl('chr21$',chrom)) # get length of chromosome 21
seqlengths = with(hg_chrs, setNames(size, chrom)) # convert to named vector

# Now follow the process described in subsection 9.5.4, Plus and minus strand cross-correlation. The maximum cross-correlation value when shifting the strands should correspond to the average DNA fragment length present in the library.* 
library(GenomicAlignments)
# use apply function to perform commands on each bam file
reads_list <- lapply(bam_files, function(x){ # for each bam file, 
  reads <-  readGAlignments(x)  # load the reads 
  reads <- granges(reads)  # and convert to a GRanges object
  reads <- resize(reads, width=1, fix='start') # then set each range to start position
  reads <- keepSeqlevels(reads, 'chr21', pruning.mode='coarse') # remove extra levels
})
# calculate the coverage profile for plus and minus strand for each dataset
wsize <- 1:400 # define shift range (for later use)
jaccard = function(x,y)sum((x & y)) / sum((x | y)) # calculates jaccard similarity
cc_list <- lapply(reads_list, function(x){ # for each GRanges object (dataset), 
  reads <- split(x, strand(x)) # split each object based on strand
  cov <- lapply(reads, function(y){
    coverage(y, width = seqlengths)[[1]] > 0 # convert coverage vector to boolean vector
    }) 
  cov <- lapply(cov, as.vector)
  cc  <- shiftApply(SHIFT = wsize, # shift the + vector by 1 - 400 nucleotides and after each shift
                    X     = cov[['+']], 
                    Y     = cov[['-']], 
                    FUN   = jaccard) # calculate the similarity between strands 
  })

# convert the results into data frames 
cc_list <- lapply(cc_list, function(cc){
  data.frame(fragment_size = wsize, cross_correlation = cc)
})
 
```

2. Visualize the resulting distributions. [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
library(ggplot2)
# make list of experiment names to label plots
expt_names <- sub('.chr21.bam','', basename(bam_files)) # remove path and common suffix
expt_names = sub('GM12878_hg38_','', expt_names) # remove common prefix
# add experiment name to each dataframe in cc_list
names(cc_list) <- expt_names

lapply(cc_list, function(cc){
  ggplot(data = cc, aes(fragment_size, cross_correlation)) +
    geom_point() +
    geom_vline(xintercept = which.max(cc$cross_correlation), 
               size=2, color='red', linetype=2) +
    theme_bw() +
    theme(
        axis.text = element_text(size=10, face='bold'),
        axis.title = element_text(size=14,face="bold"),
        plot.title = element_text(hjust = 0.5)) +
    labs(x = 'Shift in base pairs',
         y = ('Jaccard similarity')
         )
}) 
```

3. How does the Input sample distribution differ from the ChIP samples? [Difficulty: **Beginner**]

**solution:**

*The Jaccard similarity plots for the Input samples (7-11) are more likely to be diffuse with no defined peaks; the major exception is Input_r5 (#11). In contrast, the Jaccard similarity plots for the ChIP samples generally have a discrete peak, indicative of the fragment size; the major exception is H3K4me1, which is quite diffuse. (#4).*


4. Write a function which converts the bam files into bigWig files. [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
library(rtracklayer)
bam2bw <- function(bam_file){
  reads <-  readGAlignments(bam_file) # read in bam file
  reads = granges(reads) # convert to GRanges object
  reads = resize(reads, width=200, fix='start') # extend each read toward 3' end
  seqnames <- unique(seqnames(reads)) # get names (once) of each chr in object  
  reads <- keepSeqlevels(reads, seqnames, pruning.mode='coarse') # remove unused levels
  cov = coverage(reads, width = seqlengths) # convert reads into signal profile
  bw_file <- sub('.bam','.bigWig', bam_file) # change filename extension to .bigWig
  work_dir <- getwd()
  output_file <- file.path(work_dir, basename(bw_file))
  export.bw(cov, output_file) # export as bigWig file
} 
```

5. Apply the function to all files, and visualize them in the genome browser.
Observe the signal profiles. What can you notice, about the similarity of the samples? [Difficulty: **Beginner**]
  
**solution:**
```{r,echo=FALSE,eval=FALSE}
# convert bam reads to bigWig coverage profile
lapply(bam_files, bam2bw)
# visualize profiles in separate tracks
library(Gviz)
axis   = GenomeAxisTrack(
    range = GRanges('chr21', IRanges(1, width=seqlengths))
)
# convert each signal into genomic ranges and define each signal track
# read in .bigWig files with coverage profiles
bw_files <- sub('.bam','.bigWig', bam_files) 
work_dir <- getwd()
bw_files <- file.path(work_dir, basename(bw_files)) # set path to working dir
bw_list <- lapply(bw_files, import.bw)
# name each bigWig object in the list with its experiment name
expt_names <- sub('.chr21.bam','', basename(bam_files)) # remove path and common suffix
expt_names = sub('GM12878_hg38_','', expt_names) # remove common prefix
names(bw_list) <- expt_names
# visualize each profile in a genome browser
# dtrack_list <- lapply(bw_list, function(cov){
#   gcov <- as(cov, 'GRanges')
#   dtrack <- DataTrack(gcov, name = names(cov), type='l')
#   return(dtrack)
# })
# alternative; b/c when using lapply to produce dtrack_list, plotTracks threw error
i=1
for(i in 1:length(bw_list)) {
  gcov <- as(bw_list[[i]], 'GRanges')
  dtrack <- DataTrack(gcov, name = names(bw_list[i]), type='l')
  plotTracks(trackList = list(axis, dtrack), 
             sizes            = c(.1,1), 
             background.title = "black")
  } 
```

*All of the Input samples (#7-11) have a low level of reads (peaks rarely exceed 20 counts) that are distributed rather uniformly; again, Input_r5 is the exception. The ChIP samples have more peaks with high levels of reads (often hundreds per peak); SMC3_r1 is an exception.*

6. Use `GViz` to visualize the profiles for CTCF, SMC3 and ZNF143. [Difficulty: **Beginner/Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
# Visualization of profiles are generated by code in previous exercise.
 
```

7. Calculate the cross correlation for both CTCF replicates, and the input samples. How does the profile look for the control samples? [Difficulty: **Intermediate**]
  
**solution:**
```{r,echo=FALSE,eval=FALSE}
# following the example in 9.6.2
# use tileGenome() to return list of GRanges of given width, spanning whole chromosome
tilling_window = tileGenome(seqlengths, tilewidth=1000)
tilling_window = unlist(tilling_window) # convert list to one GRanges object
so = summarizeOverlaps(tilling_window, bam_files) # count reads in each window
counts = assays(so)[[1]] # extract counts from SummarizedExperiment
cpm = t(t(counts)*(1000000/colSums(counts))) # calculate cpm from counts matrix
cpm = cpm[rowSums(cpm) > 0,] # remove tiles with no reads
colnames(cpm) = sub('.chr21.bam','', colnames(cpm)) # shorten column names
colnames(cpm) = sub('GM12878_hg38_','', colnames(cpm))
# calculate the pearson correlation coefficient between CTCF and input samples
cor(cpm[,c(1,2,7:11)], method='pearson') # subset for specified columns
 
```

*The correlation between the two CTCF replicates is very high (r = 0.94) indicating that the enriched fragments were consistent. As expected, the CTCF replicates do not exhibit strong correlation with any input samples, which do not contain enriched fragments. The input samples exhibit weak to moderate cross correlation (r = 0.39-0.64), consistent with similar composition but weak signals.*

8. Calculate the cross correlation coefficients for all samples and 
visualize them as a heatmap. [Difficulty: **Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
correlation_matrix <- cor(cpm, method='pearson') # use entire matrix/all samples

library(ComplexHeatmap)
library(circlize)
heatmap_col = circlize::colorRamp2( # define color palette where 0 is white
    breaks = c(-1,0,1),
    colors = c('blue','white','red')
)

# plot the heatmap using the Heatmap function
Heatmap(
    matrix = correlation_matrix, 
    col    = heatmap_col
)
 
```

#### Peak calling

1. Use `normR` to call peaks for all SMC3, CTCF, and ZNF143 samples. [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
# rather than arbitrarily select a single Input file (with fewer reads), 
# let's merge all 5 Input files into one
library(Rsamtools)
mergeBam(files = bam_files[7:11],
         overwrite = TRUE,
         destination = file.path(work_dir, "GM12878_hg38_Inputs_merged.chr21.bam"),
         indexDestination = TRUE)

# using code in 9.6.2 as template:
library(normr)
# peak calling using CTCF (1,2), SMC3 (12,13), and ZNF143 (14,15) samples with merged Input
peaks_list <- lapply(bam_files[c(1,2,12:15)], function(chip_file){
  ctcf_fit = enrichR(treatment = chip_file,
                     control   = file.path(work_dir, "GM12878_hg38_Inputs_merged.chr21.bam"),
                     genome    = "hg38",
                     verbose   = FALSE)
  })
# show summary for all samples
lapply(peaks_list, summary)
 
```

2. Calculate the percentage of reads in peaks for the CTCF experiment. [Difficulty: **Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
# get regions with peaks for CTCF replicates (1,2)
CTCF_peaks <- lapply(peaks_list[1:2], function(fit){
  ctcf_peaks = getRanges(fit) # extract all ranges for an experiment
  ctcf_peaks$qvalue = getQvalues(fit) # annotate ranges with adjusted p values
  ctcf_peaks$enrichment = getEnrichment(fit) # annotate ranges with calculated enrichment
  ctcf_peaks = subset(ctcf_peaks, !is.na(component)) # selects ranges corresponding to enriched class
  ctcf_peaks = subset(ctcf_peaks, qvalue < 0.01) # filter by stringent q value threshold
  ctcf_peaks = ctcf_peaks[order(ctcf_peaks$qvalue)] # sort peaks based on q values
  ctcf_peaks = GenomicRanges::reduce(ctcf_peaks) # collapse neighboring regions
})
# clean up levels to avoid problems later
CTCF_peaks <- lapply(CTCF_peaks, function(sample){
  keepSeqlevels(sample, 'chr21', pruning.mode='coarse') 
  })

# find which reads overlap peak regions for each replicate
# CTCF_r1 sample
CTCF_r1_reads <- resize(reads_list[[1]], 200) # resize each read assuming 200 bp fragments
# calculate percentage of reads overlapping a peak
100*sum(countOverlaps(CTCF_r1_reads, CTCF_peaks[[1]]))/length(CTCF_r1_reads)

# CTCF_r2 sample
CTCF_r2_reads <- resize(reads_list[[2]], 200) # resize each read assuming 200 bp fragments
# calculate percentage of reads overlapping a peak
100*sum(countOverlaps(CTCF_r2_reads, CTCF_peaks[[2]]))/length(CTCF_r2_reads)

# combined for both replicates
100*(sum(countOverlaps(CTCF_r1_reads, CTCF_peaks[[1]])) + sum(countOverlaps(CTCF_r2_reads, CTCF_peaks[[2]])))/(length(CTCF_r1_reads) + length(CTCF_r2_reads))

*Roughly one-third of the reads in the CTCF_r1 (26%) and CTCF_r2 (35%) map to CTCF peaks.*
 
```

3. Download the blacklisted regions corresponding to the hg38 human genome, and calculate
the percentage of CTCF peaks falling in such regions. [Difficulty: **Advanced**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

4. Unify the biological replicates by taking an intersection of peaks.
How many peaks are specific to each biological replicate, and how many peaks overlap. [Difficulty: **Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
CTCF_peaks_unified <- CTCF_peaks[[1]][countOverlaps(CTCF_peaks[[1]], CTCF_peaks[[2]]) > 0]
length(CTCF_peaks_unified)
length(CTCF_peaks[[1]]) - length(CTCF_peaks_unified)
length(CTCF_peaks[[2]]) - length(CTCF_peaks_unified)
 
```

*Altogether 429 peaks overlap between the two CTCF replicates, while 47 peaks are specific to r1 and 286 are specific to r2.*

5. Plot a scatter plot of signal strengths for biological replicates. Do intersecting peaks have equal signal strength in both samples? [Difficulty: **Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
# for the 429 peaks shared by the CTCF replicates, count # reads per peak in each replicate
reads_CTCFpeak <- tibble(CTCF_r1 = countOverlaps(CTCF_peaks_unified, CTCF_r1_reads),
                         CTCF_r2 = countOverlaps(CTCF_peaks_unified, CTCF_r2_reads))
# normalize to account for difference in total number of reads in each sample
reads_CTCFpeak <- reads_CTCFpeak %>% 
  mutate(cpm_r1 = 10^6*CTCF_r1/length(CTCF_r1_reads),
         cpm_r2 = 10^6*CTCF_r2/length(CTCF_r2_reads))
# plot the normalized count for each CTCF peak
reads_CTCFpeak %>% 
  ggplot(aes(x=cpm_r1, y=cpm_r2)) +
  geom_point() +
  geom_abline(slope = 1, color = "red")
 
```

*After normalizing the counts, the signal strengths of the intersecting peaks are generally comparable in the two samples, as indicated by the close fit with the red line for y=x. However, there does seems to be some excess of normalized read counts in r2 vs. r1.*

6. Quantify the combinatorial binding of all three proteins. Find
the number of places which are bound by all three proteins, by 
a combination of two proteins, and exclusively by one protein.
Annotate the different regions based on their genomic location. [Difficulty: **Advanced**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

7. Correlate the normR enrichment score for CTCF with peak presence/absence
(create boxplots of enrichment for peaks which contain and do not contain CTCF motifs). [Difficulty: **Advanced**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

8. Explore the co-localization of CTCF and ZNF143. Where are the co-bound
regions located? Which sequence motifs do they contain? Download the ChIA-pet
data for the GM12878 cell line, and look at the 3D interaction between different
classes of binding sites. [Difficulty: **Advanced**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

#### Motif discovery

1. Repeat the motif discovery analysis on peaks from the ZNF143 transcription factor.
How many motifs do you observe? How do the motifs look (visualize the motif logs)? [Difficulty: **Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
ZNF143_peaks <- lapply(peaks_list[5:6], function(fit){
  peaks = getRanges(fit) # extract all ranges for an experiment
  peaks$qvalue = getQvalues(fit) # annotate ranges with adjusted p values
  peaks$enrichment = getEnrichment(fit) # annotate ranges with calculated enrichment
  peaks = subset(peaks, !is.na(component)) # selects ranges corresponding to enriched class
  peaks = subset(peaks, qvalue < 0.01) # filter by stringent q value threshold
  peaks = peaks[order(peaks$qvalue)] # sort peaks based on q values
  peaks = GenomicRanges::reduce(peaks) # collapse neighboring regions
})
# clean up levels to avoid problems later
ZNF143_peaks <- lapply(ZNF143_peaks, function(sample){
  keepSeqlevels(sample, 'chr21', pruning.mode='coarse') 
  })
# focus each peak at its center
ZNF143_peaks_resized <- lapply(ZNF143_peaks, function(peak){
  resize(peak, width = 50, fix='center')
  })
library(rGADEM)
# load the human genome sequence
library(BSgenome.Hsapiens.UCSC.hg38)
# run GADEM() for motif discovery 
novel_motifs <- GADEM(unlist(as(ZNF143_peaks_resized, "GRangesList")), # use peaks of both replicates
                    verbose=1, # print results to screen
                    genome=Hsapiens) # need genome to get base composition
consensus(novel_motifs) # view the consensus sequence of each motif
nOccurrences(novel_motifs) # count occurrences of each motif
getPWM(novel_motifs) # get position weight matrix for each motif
plot(novel_motifs) # visualize each motif
 
```

*Three motifs are observed.*

2. Scan the ZNF143 peaks with the top motifs found in the previous exercise. 
Where are the motifs located? [Difficulty: **Advanced**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

3. Scan the CTCF peaks with the top motifs identified in the **ZNF143** peaks.
Where are the motifs located? What can you conclude from the previous exercises?
[Difficulty: **Advanced**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```





