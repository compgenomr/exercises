## Exercises and solutions for Chapter 11

### Matrix factorization methods

1. Find features associated with iCluster and MFA factors, and visualize the feature weights. [Difficulty: **Beginner**]

**solution:**

iClusterPlus stores feature LASSO weights in the `beta` element of the result. Recall we saved the results of running `iClusterPlus` in the `r.icluster` variable. The data type is a list, with one element per data type, i.e. `r.icluster$beta[[1]]` will house the coefficients for the first data type, etc.

FactoMineR stores MFA coefficients in exactly the same way as NMF, in the `mfa.w <- r.mfa$quanti.var$coord` variable. Visualization can follow the same logic as the NMF example in the text. 

2. Normalizing the data matrices by their $\lambda_1$'s as in MFA supposes we wish to assign each data type the same importance in the down-stream analysis. This leads to a natural generalization whereby the different data types may be differently weighted. Provide an implementation of weighed-MFA where the different data types may be assigned individual weights. [Difficulty: **Intermediate**]

**solution:**

This can be achieved by using R's `prcomp` and foregoing the `FactoMineR` package. E.g.

```{r}
+big_x <- rbind(x1 / factor1, x2 / factor2, x3 / factor3)
+result <- prcomp(big_x) 
```

where `factor_x` can be an arbitrary scaling factor chosen for each input modality.

3. In order to use NMF algorithms on data which can be negative, we need to split each feature into two new features, one positive and one negative. Implement the following function, and see that the included test does not fail: [Difficulty: **Intermediate/Advanced**]

```{r,moNMFExerciseColumnSplitting,eval=FALSE, echo=TRUE}
# Implement this function
split_neg_columns <- function(x) {
    # your code here
}

# a test that shows the function above works
test_split_neg_columns <- function() {
    input <- as.data.frame(cbind(c(1,2,1),c(0,1,-2)))
    output <- as.data.frame(cbind(c(1,2,1), c(0,0,0), c(0,1,0), c(0,0,2)))
    stopifnot(all(output == split_neg_columns(input)))
}

# run the test to verify your solution
test_split_neg_columns()
```


**solution:**
```{r,echo=FALSE,eval=FALSE}
# For NMF to be usable on data types which allow negativity, we'll use a trick,
# where we split each column into two non-negative columns, one for the negative
# numbers and one for the positive numbers. Here's a function that does that:
split_neg_columns <- function(x) {
    new_cols <- list()
    for(i in seq_len(dim(x)[2])) {
        new_cols[[paste0(colnames(x)[i],'+')]] <- sapply(X = x[,i], 
                                                         function(x) max(0,x))
        new_cols[[paste0(colnames(x)[i],'-')]] <- sapply(X = -x[,i],
                                                         function(x) max(0,x))
    }
    new_cols
    return(do.call(cbind, new_cols))
}

# and here's a test that shows the function above works
test_split_neg_columns <- function() {
    test.input <- as.data.frame(cbind(c(1,2,1),c(0,1,-2)))
    expected.output <- as.data.frame(cbind(c(1,2,1),
                                           c(0,0,0),
                                           c(0,1,0),
                                           c(0,0,2)))
    stopifnot(all(as.matrix(expected.output) == as.matrix(split_neg_columns(
      test.input))))
}
test_split_neg_columns()

# Here's a function that undoes the previous transformation,
# so it takes each pair of columns and merges them to one, where the values
# from the first column are taken as positives, and the values from the second
# column as negatives.
merge_neg_columns <- function(x) {
    new_cols <- list()
    for(i in seq_len(dim(x)[2]/2)) {
        pos_col <- x[,2*i-1]
        neg_col <- x[,2*i]
        merged_col <- pos_col
        merged_col[neg_col>0] <- -neg_col[neg_col>0]
        new_cols[[i]] <- merged_col
    }
    return(do.call(cbind, new_cols))
}

# And here's a test for the merging function.
test_merge_neg_columns <- function() {
    input1 <- as.data.frame(cbind(c(1,2,1),c(0,1,-2)))
    input2 <- split_neg_columns(input1)
    output <- merge_neg_columns(input2)
    stopifnot(all(output == input1))
}
test_merge_neg_columns()
 
```


4. The iCluster+ algorithm has some parameters which may be tuned for maximum performance. The `iClusterPlus` package has a method, `iClusterPlus::tune.iClusterPlus`, which does this automatically based on the Bayesian Information Criterion (BIC). Run this method on the data from the examples above and find the optimal lambda and alpha values. [Difficulty: **Beginner/Intermediate**]

**solution:**



```{r,echo=FALSE,eval=FALSE}
r.icluster.tuned <- iClusterPlus::tune.iClusterPlus(
    cpus=4,
    t(x1), # Providing each omics type
    t(x2),
    t(x3),
    type=c("gaussian", "binomial", "multinomial"), # Providing the distributions
    K=2, # provide the number of factors to learn
    alpha=c(1,1,1), # as well as other model parameters
    n.lambda=35)
 
```

### Clustering using latent factors

1. Why is one-hot clustering more suitable for NMF than iCluster? [Difficulty: **Intermediate**]

**solution:**

This is best approached by looking at the NMF and iCluster latent factor values, in figures 11.9 and 11.12. There, it is apparent that NMF latent factors tend to be disentangled, i.e. only one of them tends to be active for each data point. This is exactly the motivation for one-hot clustering. Contrast this with the latent factor values for iCluster.

2. Which clustering algorithm produces better results when combined with NMF, K-means, or one-hot clustering? Why do you think that is? [Difficulty: **Intermediate/Advanced**]

**solution:**

For the examples in the text above, it is clear that one-hot clustering provides a more intuitive and easier to explain clustering result in combination with NMF, because this disentangled NMF attempts to assign one factor to each sample, which also simplifies explaining what drives each cluster (the corresponding NMF factor). However, in general, one might have specific downstream analysis tasks for which K-means or other clustering algorithms might be more appropriate.

### Biological interpretation of latent factors

1. Another covariate in the metadata of these tumors is their _CpG island methylator Phenotype_ (CIMP). This is a phenotype carried by a group of colorectal cancers that display hypermethylation of promoter CpG island sites, resulting in the inactivation of some tumor suppressors. This is also assayed using an external test. Do any of the multi-omics methods surveyed find a latent variable that is associated with the tumor's CIMP phenotype? [Difficulty: **Beginner/Intermediate**]


**solution:**

The solution is similar to what was done in the text for MSI status:

```{r,moNMFCIMP,echo=FALSE, eval=FALSE}
a <- data.frame(age=covariates$age, gender=as.factor(covariates$gender), msi=covariates$msi, cimp=as.factor(covariates$cimp))
b <- nmf.h
colnames(b) <- c('factor1', 'factor2')
cov_factor <- cbind(a,b)
ggplot2::ggplot(cov_factor, ggplot2::aes(x=cimp, y=factor1, group=cimp)) + ggplot2::geom_boxplot() + ggplot2::ggtitle("NMF factor 1 and CIMP status")
ggplot2::ggplot(cov_factor, ggplot2::aes(x=cimp, y=factor2, group=cimp)) + ggplot2::geom_boxplot() + ggplot2::ggtitle("NMF factor 2 and CIMP status")
```

At this point, examine the two figures. Is it apparent which factor is associated with CIMP status?

2. Does MFA give a disentangled representation? Does `iCluster` give disentangled representations? Why do you think that is? [Difficulty: **Advanced**]

**solution:**

MFA and iCluster tend not to give disentangled representations. MFA is based on principal component analysis, which finds variance-minimizing representations, wherein most features get assigned some nonzero weight most of the time. iCluster, with its expectation maximization approach, also does the same.

3. Figures \@ref(fig:moNMFClinicalCovariates) and \@ref(fig:moNMFClinicalCovariates2) show that MSI/MSS tumors have different values for NMF factors 1 and 2. Which NMF factor is associated with microsatellite instability? [Difficulty: **Beginner**]

**solution:**

Figure \@ref(fig:moNMFClinicalCovariates) shows that NMF factor 1 is high for samples which are microsatellite instable, and figure \@ref(fig:moNMFClinicalCovariates2) shows that NMF factor 2 is high for samples which are microsatellite stable. Hence, NMF factor 1 may be said to be associated with microsatellite instability.

4. Microsatellite instability (MSI) is associated with hyper-mutated tumors. As seen in Figure \@ref(fig:momutationsHeatmap), one of the subtypes has tumors with significantly more mutations than the other. Which subtype is that? Which NMF factor is associated with that subtype? And which NMF factor is associated with MSI? [Difficulty: **Advanced**]

**solution:**

From figure \@ref(fig:momutationsHeatmap) one can see that CMS1 tumors harbor many more mutations than do CMS3 tumors. Figure 11.13 shows that CMS1 is predicted by NMF factor 1, the same factor associated with Microsatellite instability (see previous exercise).
