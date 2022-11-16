# HZLS - Hager–Zhang line search
A MATLAB implementation of the Hager–Zhang line-search technique [1, 2], which uses the approximate Wolfe conditions.

This MATLAB version of the Hager–Zhang bracketing has been implemented by following the existing Julia code available at:
https://julianlsolvers.github.io/LineSearches.jl/stable/index.html

If there are any problems or bugs, feel free to email me at msutti (at) ncts.tw

## I) Version History
This code was originally included in the RMGLS package, released in June 2020. We decided to make this new version available independently of the RMGLS code.
* Ver 2, November 2022.
* Ver 1, June 2020. See [here](https://yareta.unige.ch/home/detail/833fad78-f35c-497c-b282-c62cae3172c9). 

## II) Main scripts
* `Driver_HZ2005_Figure_41` reproduces Figure 4.1 from [1].
* `Driver_HZLS_quartic` illustrates how to use steepest descent with the Hager–Zhang line search to minimize a simple quadratic function.
* `Driver_SD_WW_vs_HZLS` compares steepest descent with classical weak Wolfe conditions against steepest descent using the Hager–Zhang line search technique. This script can be used to reproduce the figures from Chapter 5.3 of [3].

## III) License
Code written by me is GPL licensed.

## IV) References
- [1] W. W. Hager and H. Zhang, "A new conjugate gradient method with guaranteed descent and an efficient line search", SIAM J. Optim., 16 (2005), pp. 170–192.
- [2] W. W. Hager and H. Zhang, "Algorithm 851: CG DESCENT, a Conjugate Gradient Method with Guaranteed Descent", ACM Trans. Math. Softw., 32 (2006), pp. 113–137.
- [3] M. Sutti, "Riemannian Algorithms on the Stiefel and the Fixed-Rank Manifold", Ph.D. thesis, University of Geneva, December 2020.
