## Exercises and solutions for Chapter 8

### Exploring the count tables

Here, import an example count table and do some exploration of the expression data. 

```{r exSetup1, eval=FALSE}
counts_file <- system.file("extdata/rna-seq/SRP029880.raw_counts.tsv",
                           package = "compGenomRData")
coldata_file <- system.file("extdata/rna-seq/SRP029880.colData.tsv", 
                            package = "compGenomRData")
```

1. Normalize the counts using the TPM approach. [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
# first compare total counts for each sample
counts <- as.matrix(read.table(counts_file, header = T, sep = '\t'))
colSums(counts) / 1E9 # display billions of reads for each sample 
# calculate TPM (transcripts per million)
#find gene length normalized values 
geneLengths <- as.vector(subset(counts, select = c(width)))
rpk <- apply(subset(counts, select = c(-width)), 
             2, 
             function(x) {x/(geneLengths/1000)}
            )
#normalize by the sample size using rpk values
tpm <- apply(rpk, 
             2, 
             function(x) {x / sum(as.numeric(x)) * 10^6}
             )
colSums(tpm) # confirm that normalized values sum to 10^6
head(tpm) # glimpse data

```

2. Plot a heatmap of the top 500 most variable genes. Compare with the heatmap obtained using the 100 most variable genes. [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#compute the variance of each gene across samples
V <- apply(X = tpm, 
           MARGIN = 1, 
           FUN = var)

# sort the results by variance in decreasing order 
selectedGenes <- names(V[order(V, decreasing = T)])
# make heatmaps clustering genes and samples
library(pheatmap)
# add annotations for samples 
colData <- read.table(coldata_file, header = T, sep = '\t', 
                      stringsAsFactors = TRUE)
# heatmap with the top 100 genes 
pheatmap(tpm[selectedGenes[1:100], ], scale = 'row', show_rownames = FALSE,
         annotation_col = colData)
# heatmap with the top 500 genes 
pheatmap(tpm[selectedGenes[1:500], ], scale = 'row', show_rownames = FALSE)
 
```

*Clustering with the TPM values of the top 100 genes or top 500 genes produces a tree that separately clusters the control and case samples. However, there is some difference in which samples cluster together between the two heatmaps. For example, Case 1 is most closely related to case 4 in the 100 gene tree, but Case 1 is most closely related to Case 3 in the 500 gene tree.*


3. Re-do the heatmaps setting the `scale` argument to `none`, and `column`. Compare the results with `scale = 'row'`. [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
pheatmap(tpm[selectedGenes[1:500], ], show_rownames = FALSE)
pheatmap(tpm[selectedGenes[1:500], ], scale = 'column', show_rownames = FALSE)
 
```

*Scaling by rows allows you to compare relative expression levels of each gene across the samples. Without any scaling, only extreme differences in absolute expression levels are noticeable; in particular, the variation among lowly-level genes (TPM < 50000) is not noticeable.  When scaling by columns, the heatmap indicates that nearly all of these genes with similar absolute values show average levels (log = 0) of expression among the nearly 20,000 genes overall.*


4. Draw a correlation plot for the samples depicting the sample differences as 'ellipses', drawing only the upper end of the matrix, and order samples by hierarchical clustering results based on `average` linkage clustering method. [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
correlationMatrix <- cor(tpm)
# visualize
library(corrplot)
corrplot(correlationMatrix, 
         method = "ellipse", type = "upper",
         order = 'hclust', 
         hclust.method = "average",
         # addrect = 2, 
         addCoef.col = 'black', 
         number.cex = 0.7) 
 
```

5. How else could the count matrix be subsetted to obtain quick and accurate clusters? Try selecting the top 100 genes that have the highest total expression in all samples and re-draw the cluster heatmaps and PCA plots. [Difficulty: **Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
# Heatmaps
# compute the sum of each gene across samples
TE <- apply(X = tpm, 
           MARGIN = 1, 
           FUN = sum)
# sort the results by total of expression values in decreasing order 
selectedGenes_TE <- names(TE[order(TE, decreasing = T)])
# heatmap with the top 100 genes 
pheatmap(tpm[selectedGenes_TE[1:100], ], scale = 'row', show_rownames = FALSE)

# PCA
library(ggfortify)
# transpose the matrix so genes become columns
M <- t(tpm[selectedGenes,])
# transform the counts to log2 scale 
M <- log2(M + 1) # + 1 to avoid log(0)
# compute PCA 
pcaResults <- prcomp(M)

# plot PCA results making use of ggplot2's autoplot function
autoplot(pcaResults, data = colData, colour = 'group')
 
```

6. Add an additional column to the annotation data.frame object to annotate the samples and use the updated annotation data.frame to plot the heatmaps. (Hint: Assign different batch values to CASE and CTRL samples). Make a PCA plot and color samples by the added variable (e.g. batch). [Difficulty: Intermediate]

**solution:**
```{r,echo=FALSE,eval=FALSE}
colData$batch <- c(rep(letters[1:2],times = 5)) # simulate batch effect
pheatmap(tpm[selectedGenes_TE[1:100],], scale = 'row', 
         show_rownames = FALSE, 
         annotation_col = colData)
 
```

7. Try making the heatmaps using all the genes in the count table, rather than sub-selecting. [Difficulty: **Advanced**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

8. Use the [`Rtsne` package](https://cran.r-project.org/web/packages/Rtsne/Rtsne.pdf) to draw a t-SNE plot of the expression values. Color the points by sample group. Compare the results with the PCA plots. [Difficulty: **Advanced**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```


### Differential expression analysis

Firstly, carry out a differential expression analysis starting from raw counts.
Use the following datasets:

```
counts_file <- system.file("extdata/rna-seq/SRP029880.raw_counts.tsv", 
                            package = "compGenomRData")
coldata_file <- system.file("extdata/rna-seq/SRP029880.colData.tsv", 
                            package = "compGenomRData")
```

- Import the read counts and colData tables.
- Set up a DESeqDataSet object.
- Filter out genes with low counts.
- Run DESeq2 contrasting the `CASE` sample with `CONTROL` samples. 

Now, you are ready to do the following exercises: 

1. Make a volcano plot using the differential expression analysis results. (Hint: x-axis denotes the log2FoldChange and the y-axis represents the -log10(pvalue)). [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
# First create needed objects (see list above)
# 1. read count data: remove the 'width' column from previous matrix
counts <- as.matrix(read.table(counts_file, header = T, sep = '\t'))
countData <- as.matrix(subset(counts, select = c(-width)))
# 2. get table with experimental setup 
colData <- read.table(coldata_file, header = T, sep = '\t', 
                      stringsAsFactors = TRUE)
# 3. define the design formula, indicating variable of interest in colData
designFormula <- "~ group"
# set up DESeqDataSet object
library(DESeq2)
# Now create a DESeq dataset object from the three elements above 
dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = colData, 
                              design = as.formula(designFormula))
# For each gene, we count the total number of reads for that gene in all samples 
# and only keep those that have at least 2 reads
dds <- dds[rowSums(DESeq2::counts(dds)) > 1, ]
# run analysis on this filtered data set
dds <- DESeq(dds)
# compute the contrast (difference) for the 'group' variable 
# where 'CTRL' samples are used as the control group. 
DEresults = results(dds, contrast = c("group", 'CASE', 'CTRL')) 
# sort rows by increasing p-value
DEresults <- DEresults[order(DEresults$pvalue),]

# show the top 10 genes to briefly check results
head(DEresults, n = 10) 

# make volcano plot
ggplot(data = as.data.frame(DEresults), 
       aes(x = log2FoldChange, y = -log10(pvalue))) + 
  geom_point()

```

2. Use DESeq2::plotDispEsts to make a dispersion plot and find out the meaning of this plot. (Hint: Type ?DESeq2::plotDispEsts) [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
plotDispEsts(dds, ymin = 1e-03, # removes empty plot space due to outlier
             # finalcol = NULL # uncomment to see all original values in black
             )
 
```

3. Explore `lfcThreshold` argument of the `DESeq2::results` function. What is its default value? What does it mean to change the default value to, for instance, `1`? [Difficulty: **Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
?DESeq2::results 

```

*Default value is 0. If the default value is changed to 1, then this performs a hypothesis test that the absolute values of log2 fold changes are less than or equal to 1 (assuming that altHypothesis is at the default setting of greaterAbs).*


4. What is independent filtering? What happens if we don't use it? Google `independent filtering statquest` and watch the online video about independent filtering. [Difficulty: **Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
 
```

*Independent filtering is a method to reduce the number of tested genes to improve the power of statistical testing, the ability to detect true positives. The method removes genes with low expression values, which are more likely to yield false positives due to higher dispersion (see above). By reducing the number of false positives, the proportion of true positives is increased.*


5. Re-do the differential expression analysis using the `edgeR` package. Find out how much DESeq2 and edgeR agree on the list of differentially expressed genes. [Difficulty: **Advanced**] 

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

6. Use the `compcodeR` package to run the differential expression analysis using at least three different tools and compare and contrast the results following the `compcodeR` vignette. [Difficulty: **Advanced**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```


### Functional enrichment analysis

1. Re-run gProfileR, this time using pathway annotations such as KEGG, REACTOME, and protein complex databases such as CORUM, in addition to the GO terms. Sort the resulting tables by columns `precision` and/or `recall`. How do the top GO terms change when sorted for `precision`, `recall`, or `p.value`? [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
library(gprofiler2) # updated version of gProfileR 
library(knitr) # to make nice tables

# filter for genes with significant changes in expression
#remove genes with NA values 
DE <- DEresults[!is.na(DEresults$padj),]
#select genes with adjusted p-values below 0.05
DE <- DE[DE$padj < 0.05,]
#select genes with absolute log2 fold change above 1 (two-fold change)
DE <- DE[abs(DE$log2FoldChange) > 1,]

#get the list of genes of interest
genesOfInterest <- rownames(DE)

#calculate enriched terms amongst this gene list for multiple data sources
multiResults <- gost(query = genesOfInterest, # new function with gprofiler2
                       organism = 'hsapiens', 
                       sources = c('GO', "KEGG", "REACTOME", "CORUM"))
kable(multiResults$result) # nicely formats data in table

#sort GO terms for precision, recall or p.value
goResults <- gost(query = genesOfInterest, # new function with gprofiler2
                       organism = 'hsapiens', 
                       sources = 'GO')

# sorted by increasing p_value 
# based on sizes of query (gene list), term (gene set), and intersection (overlap)
kable(head(goResults$result[order(goResults$result$p_value),])) # nicely formatted table of top lists

# sorted by decreasing precision, intersection_size/query_size
kable(head(goResults$result[order(goResults$result$precision),])) # nicely formatted table of top lists

# sorted by decreasing recall, intersection_size/term_size
kable(head(goResults$result[order(goResults$result$recall, decreasing = T),])) # nicely formatted table of top lists
 
```

*With this large gene list (query_size), the results sorted by decreasing precision is similar to the results sorted by increasing p-value. In both cases, enriched terms with very low p-values and gene sets with many entries (large term sizes) are at the top of the list. These gene sets are higher level, such as "protein binding" or "cytoplasm". In contrast, when the results are sorted by decreasing recall, gene sets with fewer entries (small term sizes) are at the top of the list. The latter sort can be useful for focusing on more specific biological features, such as "blood coagulation, fibrin clot formation" and "complement activation, alternative pathway".*


2. Repeat the gene set enrichment analysis by trying different options for the `compare` argument of the `GAGE:gage`
function. How do the results differ? [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#use the normalized counts to carry out a GSEA
normalizedCounts <- DESeq2::counts(dds, normalized = TRUE)

library(gage)
# use the top term from the GO results to create 1st gene set 
# but first restrict goResults to terms that have < 100 genes 
go <- goResults$result[goResults$result$term_size < 100,] 
geneSet1 <- gconvert(go$term_id[1], organism = "hsapiens") 

# randomly select 25 genes from the counts table for 2nd gene set
geneSet2 <- sample(rownames(normalizedCounts), 25)

geneSets <- list('top_GO_term' = geneSet1$name,
                 'random_set' = geneSet2)

# perform GSEA with different options for compare argument
# 'paired' compares each reference with its matching sample (like paired t-test) 
# appropriate here since each reference and sample is from same person
print("compare = 'paired'")
gseaResults <- gage(exprs = log2(normalizedCounts+1), # avoids log(0)
                    ref = match(rownames(colData[colData$group == 'CTRL',]), 
                                colnames(normalizedCounts)), 
                    samp = match(rownames(colData[colData$group == 'CASE',]), 
                                 colnames(normalizedCounts)),
                    gsets = geneSets, 
                    compare = 'paired')
gseaResults$greater
# 'unpaired' compares all possible pairs of ref and samples 
print("compare = 'unpaired'")
gseaResults <- gage(exprs = log2(normalizedCounts+1), # avoids log(0)
                    ref = match(rownames(colData[colData$group == 'CTRL',]), 
                                colnames(normalizedCounts)), 
                    samp = match(rownames(colData[colData$group == 'CASE',]), 
                                 colnames(normalizedCounts)),
                    gsets = geneSets, 
                    compare = 'unpaired')
gseaResults$greater
# '1ongroup' compares each sample to average of all references
print("compare = '1ongroup'")
gseaResults <- gage(exprs = log2(normalizedCounts+1), # avoids log(0)
                    ref = match(rownames(colData[colData$group == 'CTRL',]), 
                                colnames(normalizedCounts)), 
                    samp = match(rownames(colData[colData$group == 'CASE',]), 
                                 colnames(normalizedCounts)),
                    gsets = geneSets, 
                    compare = '1ongroup')
gseaResults$greater
# 'as.group' compares the two groups (like two-sample t-test) 
print("compare = 'as.group'")
gseaResults <- gage(exprs = log2(normalizedCounts+1), # avoids log(0)
                    ref = match(rownames(colData[colData$group == 'CTRL',]), 
                                colnames(normalizedCounts)), 
                    samp = match(rownames(colData[colData$group == 'CASE',]), 
                                 colnames(normalizedCounts)),
                    gsets = geneSets, 
                    compare = 'as.group')
gseaResults$greater
 
```

*The paired comparison matches each case with its control sample, which accounts for the variation between individuals. The unpaired comparison matches each case with all control samples and produces a similar p-value. The 1ongroup comparison matches each case to the average of all controls.  Finally, the as.group comparison doesn't account for individual variation, and this is reflected in the higher p-value compared to the other comparisons.*


3. Make a scatter plot of GO term sizes and obtained p-values by setting the `gProfiler::gprofiler` argument `significant = FALSE`. Is there a correlation of term sizes and p-values? (Hint: Take -log10 of p-values). If so, how can this bias be mitigated? [Difficulty: **Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
goResults_all <- gost(query = genesOfInterest, # new function with gprofiler2
                  organism = 'hsapiens', 
                  significant = FALSE,
                  sources = 'GO')
ggplot(goResults_all$result, aes(x= sqrt(term_size), y= -log10(p_value))) + 
  geom_point()
cor(x = sqrt(goResults_all$result$term_size), # square-root transformation explained below
    y = -log10(goResults_all$result$p_value))
 
```

*The scatter plot (and Pearson correlation coefficient) indicates that there is a correlation between term sizes and p-values. After trying several mathematical transformations of the explanatory variable, this was most evident with a square root transformation of term_size (and a negative log transformation of p_value). To mitigate for this bias, -log10(p-value) could be divided by the square root of the term size.*


4. Do a gene-set enrichment analysis using gene sets from top 10 GO terms. [Difficulty: **Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
# retrieve the top 10 GO term identifiers
go.term_ids <- go$term_id[1:10] 
# construct gene sets by multiply applying gconvert() for each GO term_id
geneSets_10 <- mapply(function(x) {
    gconvert(x, organism = "hsapiens")
  }, go.term_ids
)
geneSets_10
# geneSets_10 <- go.term_ids %>% 
library(dplyr) # used for pull() and %>% pipe
library(purrr) # used for map, an alternative to mapply
geneSets_10 <- go.term_ids %>% map(function(x) {
    gconvert(x, organism = "hsapiens") %>%
        pull("name")
    }
)
names(geneSets_10) <- go.term_ids # name each list element with its GO term_id
gseaResults_10 <- gage(exprs = log2(normalizedCounts+1), # avoids log(0)
                       ref = match(rownames(colData[colData$group == 'CTRL',]), 
                                   colnames(normalizedCounts)), 
                       samp = match(rownames(colData[colData$group == 'CASE',]), 
                                    colnames(normalizedCounts)),
                       gsets = geneSets_10, 
                       compare = 'paired')
gseaResults_10$greater
 
```

5. What are the other available R packages that can carry out gene set enrichment analysis for RNA-seq datasets? [Difficulty: **Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
 
```

*A general search for "gene set enrichment analysis" at bioconductor.org returned an overwhelming number of entries.  At the bottom of the bioconductor.org page on [gage](http://bioconductor.org/packages/release/bioc/html/gage.html), the biocViews Details includes a link to the (GeneSetEnrichment)(http://bioconductor.org/packages/release/BiocViews.html#___GeneSetEnrichment) category, listing a more manageable number (139) of packages. An internet search for "gene set enrichment analysis with R" identifies additional packages, including gprofiler2.* 


6.  Use the topGO package (https://bioconductor.org/packages/release/bioc/html/topGO.html) to re-do the GO term analysis. Compare and contrast the results with what has been obtained using the `gProfileR` package. Which tool is faster, `gProfileR` or topGO? Why? [Difficulty: **Advanced**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

7. Given a gene set annotated for human, how can it be utilized to work on _C. elegans_ data? (Hint: See `biomaRt::getLDS`). [Difficulty: **Advanced**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

8. Import curated pathway gene sets with Entrez identifiers from the [MSIGDB database](http://software.broadinstitute.org/gsea/msigdb/collections.jsp) and re-do the GSEA for all curated gene sets. [Difficulty: **Advanced**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```


### Removing unwanted variation from the expression data

For the exercises below, use the datasets at: 
```
counts_file <- system.file('extdata/rna-seq/SRP049988.raw_counts.tsv', 
                           package = 'compGenomRData')
colData_file <- system.file('extdata/rna-seq/SRP049988.colData.tsv', 
                           package = 'compGenomRData')
```

1. Run RUVSeq using multiple values of `k` from 1 to 10 and compare and contrast the PCA plots obtained from the normalized counts of each RUVSeq run. [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

2. Re-run RUVSeq using the `RUVr()` function. Compare PCA plots from `RUVs`, `RUVg` and `RUVr` using the same `k` values and find out which one performs the best. [Difficulty: **Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

3. Do the necessary diagnostic plots using the differential expression results from the EHF count table. [Difficulty: **Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

4. Use the `sva` package to discover sources of unwanted variation and re-do the differential expression analysis using variables from the output of `sva` and compare the results with `DESeq2` results using `RUVSeq` corrected normalization counts. [Difficulty: **Advanced**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```


