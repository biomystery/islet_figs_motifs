* chromVAR 

Step 1: `creat_count_object_for_chromVAR.R`

- create SummarizedExperiment object using `counts per peak per cell matrix` and `reads depth per cell`
- filterred sample with mininal depth `min_depth=1500` and `min_in_peaks=0.15` by using `filterSamplesPlot` function
- filtered peaks using `filterPeaks` 
- added GC bias by `addGCBias` using genome `BSgenome.Hsapiens.UCSC.hg19`
- lastly added cell_type assignment for each cell

Step 2: `run_chromVAR_over_jaspar.R` and select variable motifs 

- use default motif database by `getJapserMotifs`
- Then `computeDeviations` function to calculate Z score matrix
- variability was calculated by `computeVariability` function and variable motifs are selected by variability larger than 1.2. Then the mean Z scores for each motif on each celltype are calculated and scaled, presenting in heatmap.  


Step 3: select celltype-specific variable motifs 

- Student's T tests were performed for each motif to find out motifs whose Z scores are significant different between sub-cell-types (alpha1 and alpha2 or beta1 and beta2). The raw p-value from T tests are adjusted by BH procedure to produce FDR. Those motifs with FDR less than 0.01 and absolute delta Z values larger than 0.5 are selected. (`calc_varability.R`)

- These sub-cell-type specific motifs are selected and calculated the average Z score across the sub-celltypes.   

* Differential promotor openning 

* 