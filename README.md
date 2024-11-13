# GNS_pm_relative

## Obtaining Relative Proper Motions from GNS1 and GNS2 Using GNS2 as the Reference Frame

This process calculates proper motions using data from GNS1 and GNS2.

0. **`gns1_B_lists.py`**  
   - Combines the H and Ks lists for GNS2.

1. **`GNS1_GNS2_alignment.py`**
   - Performs alignment between GNS1 and GNS2 with an affine transformation first, followed by a polynomial transformation.
   - Uses a Python version of `polywarp` to compute proper motions.

2. **`gns_pm_comparison.py`**   
   - Compares the proper motions with those from Gaia and Hosek.