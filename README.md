# GNS_pm_relative

## Obtaining Relative Proper Motions from GNS1 and GNS2 Using GNS2 as the Reference Frame

Proper motions calculated from GNS1 and GNS2 data.

0. **`gns1_B_lists.py`**  
   Combines the H and Ks lists for the new reduction.
   
1. **`GNS1_GNS2_alignment.py`**
Performs an initial alignment (GNS1 to GNS2) using `aa`, transforms the xy coordinates of GNS1, and generates an alignment list. This script also creates first- and second-order grids for alignment (used in the next script) using IDL code.  

1. **`aa_GNS_GNS2.py`**  
   Performs an initial alignment (GNS1 to GNS2) using `aa`, transforms the xy coordinates of GNS1, and generates an alignment list. This script also creates first- and second-order grids for alignment (used in the next script) using IDL code.  
   > **Note**: Ban lists are located in TEAPOT -> `/home/data/results/HAWK-I/data/list_2015/field12/chip1/cat_Ban_12_1.txt`.  
   > **Note 2**: `aa_fg_GNS_GNS2.py` filters for foreground stars only. This is useful for testing the alignment using only foreground stars, allowing comparison with Gaia, which includes only foreground stars in the relative PM frame.

2. **`IDL_BSdeg2_align_gns_1_2.pro`**  
   Creates lists of the aligned star xy positions using bootstrapping, with degrees from 1 to 3.  
   > **Note**: This script uses modified `compare_lists_BS.pro`, which accounts for repeated values. This allows for lists containing repeated points, facilitating bootstrap statistics.

___

2.1 **`IDL_JK_align_gns_1_2.pro`**  
   Computes jackknife alignment using the grid reference stars and stores the resulting lists.

___

3. **`IDL_aligns_gns_1_2.pro`**  
   Aligns GNS1 with GNS2 using all common stars and iterates with polynomial transformations.

3.1 **`IDL_grid_aligns_gns1_2.pro`**  
   Similar to step #3, but uses only the reference stars from the grid generated in step #1.

4. **`pm_analysis_relative.py`**  
   Generates overlapping tiles and searches for clusters within each tile, creating corresponding lists.

5. **`candidates_plotting_gns.py`**  
   Plots all identified clusters in the cluster frame.

## Quality Check

6. **`pm_uncerts_relative.py`**  
   Plots position, alignment, and proper motion uncertainties.
