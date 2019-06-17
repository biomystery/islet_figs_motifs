
# Run ChromVar 

## Steps: 

0. <filter_reads_in_final_cells.R>
   * Input: 1. the raw reads file;  2. the final cell and cell_type file.
   * Output: 1. a sparse matrix for peaks counts.; 2. a sumarizeReads
   
1. <count_reads_in_peaks.R>
   * input: 1. peaks (255k) 2. filtered reads (128M reads) 3. cell-celltype dict (15k cells)
   * output: 1. peaks x cells - sparse counts

2. <cout_reads_depth_per_cell.R>
   * input: 1.  filtered reads (o128M reads) 3. cell-celltype dict (15k cells)
   * output: 1.# reads per cells
   
3. <creat_count_object_for_chromVAR.R>
   * input: 1. peaks x cells - sparse counts 2. cell,read_depth 3. cell_type assignement
   * output: 1. summarizedExperiement(SE) obj for chromVAR
   
4. <run_chromVAR_over_jaspar.R>
	*  input: 1. summarizedExperiement(SE) obj for chromVAR 2. Jaspar matrix 
	*  output: 1. motif x cell (z score) 2. plot: ranked

   
