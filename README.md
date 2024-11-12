# GNS_pm_relative

## Obtaining Relative Proper Motions from GNS1 and GNS2 Using GNS2 as the Reference Frame

Proper motions calculated from GNS1 and GNS2 data.

0. **`gns1_B_lists.py`**  
   Combines the H and Ks lists for the new reduction.
   
1. **`GNS1_GNS2_alignment.py`**
   - Perferm the alignment between GNS1 and GNS2 with an afine transformation frisr and the with a polynomial, using a python version of polywarp and compute the proper motions
   - Comapares the pm with those of Gaia and Hosek.