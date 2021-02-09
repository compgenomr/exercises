## Exercises

### Quality control

1. Apply the fragment size estimation procedure to all ChIP and Input available datasets. [Difficulty: **Beginner**]

2. Visualize the resulting distributions. [Difficulty: **Beginner**]

3. How does the Input sample distribution differ from the ChIP samples? [Difficulty: **Beginner**]

4. Write a function which converts the bam files into bigWig files. [Difficulty: **Beginner**]

5. Apply the function to all files, and visualize them in the genome browser.
Observe the signal profiles. What can you notice, about the similarity of the samples? [Difficulty: **Beginner**]
  
6. Use `GViz` to visualize the profiles for CTCF, SMC3 and ZNF143. [Difficulty: **Beginner/Intermediate**]

7. Calculate the cross correlation for both CTCF replicates, and
the input samples. How does the profile look for the control samples? [Difficulty: **Intermediate**]
  
8. Calculate the cross correlation coefficients for all samples and 
visualize them as a heatmap. [Difficulty: **Intermediate**]

#### Peak calling

1. Use `normR` to call peaks for all SMC3, CTCF, and ZNF143 samples. [Difficulty: **Beginner**]

2. Calculate the percentage of reads in peaks for the CTCF experiment. [Difficulty: **Intermediate**]

3. Download the blacklisted regions corresponding to the hg38 human genome, and calculate
the percentage of CTCF peaks falling in such regions. [Difficulty: **Advanced**]

4. Unify the biological replicates by taking an intersection of peaks.
How many peaks are specific to each biological replicate, and how many peaks overlap. [Difficulty: **Intermediate**]

5. Plot a scatter plot of signal strengths for biological replicates. Do intersecting
peaks have equal signal strength in both samples? [Difficulty: **Intermediate**]

6. Quantify the combinatorial binding of all three proteins. Find
the number of places which are bound by all three proteins, by 
a combination of two proteins, and exclusively by one protein.
Annotate the different regions based on their genomic location. [Difficulty: **Advanced**]

7. Correlate the normR enrichment score for CTCF with peak presence/absence
(create boxplots of enrichment for peaks which contain and do not contain CTCF motifs). [Difficulty: **Advanced**]

8. Explore the co-localization of CTCF and ZNF143. Where are the co-bound
regions located? Which sequence motifs do they contain? Download the ChIA-pet
data for the GM12878 cell line, and look at the 3D interaction between different
classes of binding sites. [Difficulty: **Advanced**]

#### Motif discovery

1. Repeat the motif discovery analysis on peaks from the ZNF143 transcription factor.
How many motifs do you observe? How do the motifs look (visualize the motif logs)? [Difficulty: **Intermediate**]

2. Scan the ZNF143 peaks with the top motifs found in the previous exercise. 
Where are the motifs located? [Difficulty: **Advanced**]

3. Scan the CTCF peaks with the top motifs identified in the **ZNF143** peaks.
Where are the motifs located? What can you conclude from the previous exercises?
[Difficulty: **Advanced**]

```{r include=FALSE, eval=TRUE, echo=FALSE}
rm(list=ls())
gc()
```




