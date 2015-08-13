#UniOrder - Ordering by Minimum Entropy Approach.

# Introduction #

Inference of linear order of elements using Minimum Entropy Approach and Best-Fit Extension.

# Details #

Python\_UniOrder - script to infer linear order of elements by analysis of two dimensional matrices of pairwise distances. Python\_UniOrder tries to find the 2D matrix which has minimal total sum of differences between adjacent cells. The script calculates so called 'delta' for each pair of adjacent cells by subtracting one pairwise score from the adjacent one. Then it calculates the sum of absolute values of all deltas and chooses that matrix which has a lowest value. In other words, Python\_UniOrder is searching for the matrix with the lowest entropy among a set of available matrices.

All values in a two dimensional matrix of pairwise distances are taken into account even for pairs of very distant items. In other words, contribution of pairwise distances between any pair of elements are equal to find the best optimal order with lowest entropy.

Finding of the optimal linear order for N elements has (N!)/2 complexity (factorial problem). For ten items it is (10!)/2 = 1,814,400 (almost two million of different orders is available for 10 elements). So, we can not check all available matrices for a set of 12 elements or higher in a real time (using single CPU of 2 to 4 GHz).

To minimize a number of matrices to analyze, Python\_UniOrder script uses an initial frame work order of elements and tries to insert all other elements one by one into the frame map calculating 'delta' for each iteration. New map with lowest delta is selected after each iteration. Run time using this approach does not exceed N x N x N time.

If initial frame work order is not available, Python\_UniOrder can check ALL POSSIBLE COMBINATIONS of orders for up to 10 chosen elements. Then it adds items one by one from the 'markers to map list'. In this case run time for the script is (N!)/2 + (M x M x M) where N is a number of frame work markers, M is a number of all markers to map.

Note: UniOrder was derived from Python\_MadMapper\_V248\_XDELTA\_119.py script (see atgc-map project)