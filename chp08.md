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
#coming soon
 
```

2. Plot a heatmap of the top 500 most variable genes. Compare with the heatmap obtained using the 100 most variable genes. [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

3. Re-do the heatmaps setting the `scale` argument to `none`, and `column`. Compare the results with `scale = 'row'`. [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

4. Draw a correlation plot for the samples depicting the sample differences as 'ellipses', drawing only the upper end of the matrix, and order samples by hierarchical clustering results based on `average` linkage clustering method. [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

5. How else could the count matrix be subsetted to obtain quick and accurate clusters? Try selecting the top 100 genes that have the highest total expression in all samples and re-draw the cluster heatmaps and PCA plots. [Difficulty: **Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

6. Add an additional column to the annotation data.frame object to annotate the samples and use the updated annotation data.frame to plot the heatmaps. (Hint: Assign different batch values to CASE and CTRL samples). Make a PCA plot and color samples by the added variable (e.g. batch). [Difficulty: Intermediate]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
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
#coming soon
 
```

2. Use DESeq2::plotDispEsts to make a dispersion plot and find out the meaning of this plot. (Hint: Type ?DESeq2::plotDispEsts) [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

3. Explore `lfcThreshold` argument of the `DESeq2::results` function. What is its default value? What does it mean to change the default value to, for instance, `1`? [Difficulty: **Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

4. What is independent filtering? What happens if we don't use it? Google `independent filtering statquest` and watch the online video about independent filtering. [Difficulty: **Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

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
#coming soon
 
```

2. Repeat the gene set enrichment analysis by trying different options for the `compare` argument of the `GAGE:gage`
function. How do the results differ? [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

3. Make a scatter plot of GO term sizes and obtained p-values by setting the `gProfiler::gprofiler` argument `significant = FALSE`. Is there a correlation of term sizes and p-values? (Hint: Take -log10 of p-values). If so, how can this bias be mitigated? [Difficulty: **Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

4. Do a gene-set enrichment analysis using gene sets from top 10 GO terms. [Difficulty: **Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

5. What are the other available R packages that can carry out gene set enrichment analysis for RNA-seq datasets? [Difficulty: **Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

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


