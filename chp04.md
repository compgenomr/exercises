## Exercises and solutions for Chapter 4

For this set of exercises we will be using the expression data shown below:
```{r dataLoadClu,eval=FALSE}
expFile=system.file("extdata",
                    "leukemiaExpressionSubset.rds",
                    package="compGenomRData")
mat=readRDS(expFile)

```

### Clustering

1. We want to observe the effect of data transformation in this exercise. Scale the expression matrix with the `scale()` function. In addition, try taking the logarithm of the data with the `log2()` function prior to scaling. Make box plots of the unscaled and scaled data sets using the `boxplot()` function. [Difficulty: **Beginner/Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
# transform data matrix
scaled_mat <- scale(mat) # by default, center=TRUE, scale=TRUE
logscaled_mat <- scale(log2(mat))
# make boxplots
boxplot(mat)
boxplot(scaled_mat)
boxplot(logscaled_mat) 
```

2. For the same problem above using the unscaled data and different data transformation strategies, use the `ward.d` distance in hierarchical clustering and plot multiple heatmaps. You can try to use the `pheatmap` library or any other library that can plot a heatmap with a dendrogram. Which data-scaling strategy provides more homogeneous clusters with respect to disease types? [Difficulty: **Beginner/Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
library(pheatmap) # See https://towardsdatascience.com/pheatmap-draws-pretty-heatmaps-483dab9a3cc for a helpful tutorial 

# The leukemia type for each sample is indicated by the first 3 letters of each ID, which we'll store to use in pheatmap
annotation_col <-  data.frame(LeukemiaType =substr(colnames(mat),1,3)) # store first 3 letters from each ID in dataframe
rownames(annotation_col)=colnames(mat) # add rownames to this dataframe to crossreference each ID to its type in next step

# generate heatmap with unscaled data
pheatmap(mat = mat,
  show_rownames = FALSE, show_colnames = FALSE, # names are too long to print
  annotation_col = annotation_col, # to show leukemia type for each sample
  clustering_method = "ward.D2", # distance metric is euclidean by default
  main = "Unscaled heatmap") # adds title
# generate heatmap with scaled data
pheatmap(mat = scaled_mat, show_rownames = FALSE, show_colnames = FALSE, annotation_col = annotation_col, clustering_method = "ward.D2", main = "Scaled heatmap")
# generate heatmap with log-scaled data
pheatmap(mat = logscaled_mat, show_rownames = FALSE, show_colnames = FALSE, annotation_col = annotation_col, clustering_method = "ward.D2", main = "Log-scaled heatmap")
```

The clusters seem comparable with all three strategies: In all three heatmaps, one case doesn't cluster with its type, but otherwise the remaining cases do cluster as expected.


3. For the transformed and untransformed data sets used in the exercise above, use the silhouette for deciding number of clusters using hierarchical clustering. [Difficulty: **Intermediate/Advanced**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

4. Now, use the Gap Statistic for deciding the number of clusters in hierarchical clustering. Is it the same number of clusters identified by two methods? Is it similar to the number of clusters obtained using the k-means algorithm in the chapter. [Difficulty: **Intermediate/Advanced**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

### Dimension reduction
We will be using the leukemia expression data set again. You can use it as shown in the clustering exercises.

1. Do PCA on the expression matrix using the `princomp()` function and then use the `screeplot()` function to visualize the explained variation by eigenvectors. How many top components explain 95% of the variation? [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

2. Our next tasks are to remove eigenvectors and reconstruct the matrix using SVD, then calculate the reconstruction error as the difference between original and reconstructed matrix. Remove a few eigenvectors, reconstruct the matrix and calculate the reconstruction error. Reconstruction error can be euclidean distance between original and reconstructed matrices. HINT: You have to use the `svd()` function and equalize eigenvalue to $0$ for the component you want to remove. [Difficulty: **Intermediate/Advanced**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

3. Produce a 10-component ICA from the expression data set. Remove each component and measure the reconstruction error without that component. Rank the components by decreasing reconstruction-error. [Difficulty: **Advanced**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

4. In this exercise we use the `Rtsne()` function on the leukemia expression data set. Try to increase and decrease perplexity t-sne, and describe the observed changes in 2D plots. [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```
