# chromVAR 

## Step 1: `creat_count_object_for_chromVAR.R`

- create SummarizedExperiment object using `counts per peak per cell
  matrix` and `reads depth per cell` 
- filterred sample with mininal depth `min_depth=1500` and
  `min_in_peaks=0.15` by using `filterSamplesPlot` function 
- filtered peaks using `filterPeaks` 
- added GC bias by `addGCBias` using genome `BSgenome.Hsapiens.UCSC.hg19`
- lastly added cell_type assignment for each cell

## Step 2: `run_chromVAR_over_jaspar.R` and select variable motifs 

- use default motif database by `getJapserMotifs`
- Then `computeDeviations` function to calculate Z score matrix
  (motifs by cell)
- variability was calculated by `computeVariability` function and
variable motifs are selected by variability larger than 1.2. Then the
mean Z scores for each motif on each celltype are calculated and
scaled, presenting in heatmap.


## Step 3: select celltype-specific variable motifs 

- Student's T tests were performed for each motif to find out motifs
  whose Z scores are significant different between sub-cell-types
  (alpha1 and alpha2 or beta1 and beta2). The raw p-value from T tests
  are adjusted by BH procedure to produce FDR. Those motifs with FDR
  less than 0.01 and absolute delta Z values larger than 0.5 are
  selected. (`calc_varability.R`) 

- These sub-cell-type specific motifs are selected and calculated the
  average Z score across the sub-celltypes.
  
## Step 4: get motif by pesudostate matrix

- The previous derivated motifs by cell Z score matrix was used to get
  the motifs by pesudostat matrix by coverting the cell id to the
  corresponding pesudostates. 
  

# Differential promotor openning 
## Data preparation (`summarize data at transcript level`)
- get all the TSS from gencode and extend to -500bp to +500bp regions defined as promoter
- find peaks overlap with any promoter and annotate that peak by that promoter's gene 
- check the binary peak-cell matrix and get the cell-gene matrix whose
  value is whether the cell has a reads in peaks overlap with gene's
  promoter region. 

## Differential opening promoters between cell subtypes (alpha1
  vs. alpha2 and beta 1 vs beta 2) 

- (`Binary_comp`) for all the cells labeled alpha cells, we performed
one-sided Fisher's exact text for each promoter to test whether this 
gene is more open in one subtype than the other. We adjusted p value
from the test using Bonferroni correction. The promoters with adjusted
p value less than 0.01 were selected as significant differential
openning. The odds ratios from the Fisher's exact test were used to
rank transcripts in GSEA analysis.

- The unique genes from the previous transcrit lists are input into
  EnrichR via R interphase to perform GO Biological Processes (BP)
  enrichment analysis using the 2018 version  db. Only BP terms
  contain less than 100 and adjusted p value less than 0.1 are
  reported. 


  
  
