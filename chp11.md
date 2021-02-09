## Exercises

### Matrix factorization methods

1. Find features associated with iCluster and MFA factors, and visualize the feature weights. [Difficulty: **Beginner**]

2. Normalizing the data matrices by their $\lambda_1$'s as in MFA supposes we wish to assign each data type the same importance in the down-stream analysis. This leads to a natural generalization whereby the different data types may be differently weighted. Provide an implementation of weighed-MFA where the different data types may be assigned individual weights. [Difficulty: **Intermediate**]

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

4. The iCluster+ algorithm has some parameters which may be tuned for maximum performance. The `iClusterPlus` package has a method, `iClusterPlus::tune.iClusterPlus`, which does this automatically based on the Bayesian Information Criterion (BIC). Run this method on the data from the examples above and find the optimal lambda and alpha values. [Difficulty: **Beginner/Intermediate**]

### Clustering using latent factors

1. Why is one-hot clustering more suitable for NMF than iCluster? [Difficulty: **Intermediate**]

2. Which clustering algorithm produces better results when combined with NMF, K-means, or one-hot clustering? Why do you think that is? [Difficulty: **Intermediate/Advanced**]

### Biological interpretation of latent factors

1. Another covariate in the metadata of these tumors is their _CpG island methylator Phenotype_ (CIMP). This is a phenotype carried by a group of colorectal cancers that display hypermethylation of promoter CpG island sites, resulting in the inactivation of some tumor suppressors. This is also assayed using an external test. Do any of the multi-omics methods surveyed find a latent variable that is associated with the tumor's CIMP phenotype? [Difficulty: **Beginner/Intermediate**]

```{r,moNMFCIMP,echo=FALSE, eval=FALSE}
a <- data.frame(age=covariates$age, gender=as.factor(covariates$gender), msi=covariates$msi, cimp=as.factor(covariates$cimp))
b <- nmf.h
colnames(b) <- c('factor1', 'factor2')
cov_factor <- cbind(a,b)
ggplot2::ggplot(cov_factor, ggplot2::aes(x=cimp, y=factor1, group=cimp)) + ggplot2::geom_boxplot() + ggplot2::ggtitle("NMF factor 1 and CIMP status")
ggplot2::ggplot(cov_factor, ggplot2::aes(x=cimp, y=factor2, group=cimp)) + ggplot2::geom_boxplot() + ggplot2::ggtitle("NMF factor 2 and CIMP status")
```

2. Does MFA give a disentangled representation? Does `iCluster` give disentangled representations? Why do you think that is? [Difficulty: **Advanced**]

3. Figures \@ref(fig:moNMFClinicalCovariates) and \@ref(fig:moNMFClinicalCovariates2) show that MSI/MSS tumors have different values for NMF factors 1 and 2. Which NMF factor is associated with microsatellite instability? [Difficulty: **Beginner**]

4. Microsatellite instability (MSI) is associated with hyper-mutated tumors. As seen in Figure \@ref(fig:momutationsHeatmap), one of the subtypes has tumors with significantly more mutations than the other. Which subtype is that? Which NMF factor is associated with that subtype? And which NMF factor is associated with MSI? [Difficulty: **Advanced**]

