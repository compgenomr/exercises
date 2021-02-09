## Exercises and solutions for chapter 3

### How to summarize collection of data points: The idea behind statistical distributions

1. Calculate the means and variances 
of the rows of the following simulated data set, and plot the distributions
of means and variances using `hist()` and `boxplot()` functions. [Difficulty: **Beginner/Intermediate**]  
```{r getDataChp3Ex,eval=FALSE}
set.seed(100)

#sample data matrix from normal distribution
gset=rnorm(600,mean=200,sd=70)
data=matrix(gset,ncol=6)
```


2. Using the data generated above, calculate the standard deviation of the
distribution of the means using the `sd()` function. Compare that to the expected
standard error obtained from the central limit theorem keeping in mind the
population parameters were  $\sigma=70$ and $n=6$. How does the estimate from the random samples change if we simulate more data with
`data=matrix(rnorm(6000,mean=200,sd=70),ncol=6)`? [Difficulty: **Beginner/Intermediate**] 

3. Simulate 30 random variables using the `rpois()` function. Do this 1000 times and calculate the mean of each sample. Plot the sampling distributions of the means
using a histogram. Get the 2.5th and 97.5th percentiles of the
distribution. [Difficulty: **Beginner/Intermediate**] 
4. Use the `t.test()` function to calculate confidence intervals
of the mean on the first random sample `pois1` simulated from the `rpois()` function below. [Difficulty: **Intermediate**] 
```{r exRpoisChp3,eval=FALSE}
#HINT
set.seed(100)

#sample 30 values from poisson dist with lamda paramater =30
pois1=rpois(30,lambda=5)

```
5. Use the bootstrap confidence interval for the mean on `pois1`. [Difficulty: **Intermediate/Advanced**] 
6. Compare the theoretical confidence interval of the mean from the `t.test` and the bootstrap confidence interval. Are they similar? [Difficulty: **Intermediate/Advanced**]  
7. Try to re-create the following figure, which demonstrates the CLT concept.[Difficulty: **Advanced**] 
```{r,echo=FALSE,message=FALSE,warning=FALSE}
set.seed(101)
require(mosaic)
par(mfcol=c(4,3))
par(mar=c(5.1-2,4.1-1,4.1,2.1-2))
d=c(rnorm(1000,mean=10,sd=8))
hist(d,main="",
     col="black",border="white",breaks=20,xlab="",ylab=""
     )
abline(v=mean(d),col="red")
mtext(expression(paste(mu,"=10")),cex=0.6)
mtext("normal",cex=0.8,line=1)

bimod10=rowMeans(do(1000)*rnorm(5,mean=10,sd=8))
bimod30=rowMeans(do(1000)*rnorm(15,mean=10,sd=8))
bimod100=rowMeans(do(1000)*rnorm(50,mean=10,sd=8))
hist(bimod10,xlim=c(0,20),main="",xlab="",ylab="",breaks=20,col="gray",
     border="gray")
mtext("n=10",side=2,cex=0.8,line=2)
hist(bimod30,xlim=c(0,20),main="",xlab="",ylab="",breaks=20,col="gray",
     border="gray")
mtext("n=30",side=2,cex=0.8,line=2)
hist(bimod100,xlim=c(0,20),main="",xlab="",ylab="",breaks=20,col="gray",
     border="gray")
mtext("n=100",side=2,cex=0.8,line=2)

d=rexp(1000)
hist(d,main="",
     col="black",border="white",breaks=20,xlab="",ylab=""
     )
abline(v=mean(d),col="red")
mtext(expression(paste(mu,"=1")),cex=0.6)
mtext("exponential",cex=0.8,line=1)
mtext("Distributions of different populations",line=2)

exp10 =rowMeans(do(2000)*rexp(10))
exp30 =rowMeans(do(2000)*rexp(30))
exp100=rowMeans(do(2000)*rexp(100))
hist(exp10,xlim=c(0,2),main="",xlab="",ylab="",breaks=20,col="gray",
     border="gray")
mtext("Sampling distribution of sample means",line=2)
hist(exp30,xlim=c(0,2),main="",xlab="",ylab="",breaks=20,col="gray",
     border="gray")
hist(exp100,xlim=c(0,2),main="",xlab="",ylab="",breaks=20,col="gray",
     border="gray")

d=runif(1000)
hist(d,main="",
     col="black",border="white",breaks=20,xlab="",ylab=""
     )
abline(v=mean(d),col="red")
mtext(expression(paste(mu,"=0.5")),cex=0.6)

mtext("uniform",cex=0.8,line=1)
unif10 =rowMeans(do(1000)*runif(10))
unif30 =rowMeans(do(1000)*runif(30))
unif100=rowMeans(do(1000)*runif(100))
hist(unif10,xlim=c(0,1),main="",xlab="",ylab="",breaks=20,col="gray",
     border="gray")
hist(unif30,xlim=c(0,1),main="",xlab="",ylab="",breaks=20,col="gray",
     border="gray")
hist(unif100,xlim=c(0,1),main="",xlab="",ylab="",breaks=20,col="gray",
     border="gray")
```

### How to test for differences in samples
1. Test the difference of means of the following simulated genes
using the randomization, `t-test()`, and `wilcox.test()` functions.
Plot the distributions using histograms and boxplots. [Difficulty: **Intermediate/Advanced**] 
```{r exRnorm1chp3,eval=FALSE}
set.seed(101)
gene1=rnorm(30,mean=4,sd=3)
gene2=rnorm(30,mean=3,sd=3)

```

2. Test the difference of the means of the following simulated genes
using the randomization, `t-test()` and `wilcox.test()` functions.
Plot the distributions using histograms and boxplots. [Difficulty: **Intermediate/Advanced**] 
```{r exRnorm2chp3,eval=FALSE}
set.seed(100)
gene1=rnorm(30,mean=4,sd=2)
gene2=rnorm(30,mean=2,sd=2)

```

3. We need an extra data set for this exercise. Read the gene expression data set as follows:
`gexpFile=system.file("extdata","geneExpMat.rds",package="compGenomRData") data=readRDS(gexpFile)`. The data has 100 differentially expressed genes. The first 3 columns are the test samples, and the last 3 are the control samples. Do 
a t-test for each gene (each row is a gene), and record the p-values.
Then, do a moderated t-test, as shown in section "Moderated t-tests" in this chapter, and record 
the p-values. Make a p-value histogram and compare two approaches in terms of the number of significant tests with the $0.05$ threshold.
On the p-values use FDR (BH), Bonferroni and q-value adjustment methods.
Calculate how many adjusted p-values are below 0.05 for each approach.
[Difficulty: **Intermediate/Advanced**] 

### Relationship between variables: Linear models and correlation

Below we are going to simulate X and Y values that are needed for the 
rest of the exercise.
```{r exLM1chp3,eval=FALSE}
# set random number seed, so that the random numbers from the text
# is the same when you run the code.
set.seed(32)

# get 50 X values between 1 and 100
x = runif(50,1,100)

# set b0,b1 and variance (sigma)
b0 = 10
b1 = 2
sigma = 20
# simulate error terms from normal distribution
eps = rnorm(50,0,sigma)
# get y values from the linear equation and addition of error terms
y = b0 + b1*x+ eps
```


1. Run the code then fit a line to predict Y based on X. [Difficulty:**Intermediate**] 
2. Plot the scatter plot and the fitted line. [Difficulty:**Intermediate**] 
3. Calculate correlation and R^2. [Difficulty:**Intermediate**] 
4. Run the `summary()` function and 
try to extract P-values for the model from the object
returned by `summary`. See `?summary.lm`. [Difficulty:**Intermediate/Advanced**] 
5. Plot the residuals vs. the fitted values plot, by calling the `plot()` 
function with `which=1` as the second argument. First argument
is the model returned by `lm()`. [Difficulty:**Advanced**] 


6. For the next exercises, read the data set histone modification data set. Use the following to get the path to the file:
```
hmodFile=system.file("extdata",
                    "HistoneModeVSgeneExp.rds",
                     package="compGenomRData")`
```
There are 3 columns in the dataset. These are measured levels of H3K4me3,
H3K27me3 and gene expression per gene. Once you read in the data, plot the scatter plot for H3K4me3 vs. expression. [Difficulty:**Beginner**] 

7. Plot the scatter plot for H3K27me3 vs. expression. [Difficulty:**Beginner**] 

8. Fit the model for prediction of expression data using: 1) Only H3K4me3 as explanatory variable, 2) Only H3K27me3 as explanatory variable, and 3) Using both H3K4me3 and H3K27me3 as explanatory variables. Inspect the `summary()` function output in each case, which terms are significant. [Difficulty:**Beginner/Intermediate**] 

10. Is using H3K4me3 and H3K27me3 better than the model with only H3K4me3? [Difficulty:**Intermediate**] 

11. Plot H3k4me3 vs. H3k27me3. Inspect the points that do not
follow a linear trend. Are they clustered at certain segments 
of the plot? Bonus: Is there any biological or technical interpretation
for those points? [Difficulty:**Intermediate/Advanced**] 