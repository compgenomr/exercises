## Exercises and solutions for Chapter 9

### Quality control

1. Apply the fragment size estimation procedure to all ChIP and Input available datasets. [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

2. Visualize the resulting distributions. [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

3. How does the Input sample distribution differ from the ChIP samples? [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

4. Write a function which converts the bam files into bigWig files. [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

5. Apply the function to all files, and visualize them in the genome browser.
Observe the signal profiles. What can you notice, about the similarity of the samples? [Difficulty: **Beginner**]
  
**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

6. Use `GViz` to visualize the profiles for CTCF, SMC3 and ZNF143. [Difficulty: **Beginner/Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

7. Calculate the cross correlation for both CTCF replicates, and
the input samples. How does the profile look for the control samples? [Difficulty: **Intermediate**]
  
**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

8. Calculate the cross correlation coefficients for all samples and 
visualize them as a heatmap. [Difficulty: **Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

#### Peak calling

1. Use `normR` to call peaks for all SMC3, CTCF, and ZNF143 samples. [Difficulty: **Beginner**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

2. Calculate the percentage of reads in peaks for the CTCF experiment. [Difficulty: **Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
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
#coming soon
 
```

5. Plot a scatter plot of signal strengths for biological replicates. Do intersecting
peaks have equal signal strength in both samples? [Difficulty: **Intermediate**]

**solution:**
```{r,echo=FALSE,eval=FALSE}
#coming soon
 
```

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
#coming soon
 
```

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





