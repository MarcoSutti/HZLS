# HZLS - Hager--Zhang line search
MATLAB implementation of the Hager-Zhang line-search technique.

Line search by Hager and Zhang, which uses the approximate Wolfe conditions.
This MATLAB version of the Hager-Zhang bracketing has been implemented by
following the existing Julia code available at:
https://julianlsolvers.github.io/LineSearches.jl/stable/index.html

## Drivers
* `Driver_HZ2005_Figure_41` reproduces Figure 4.1 in [1].
* `Driver_HZLS_quartic` illustrates how to use steepest descent with the Hager--Zhang line search to minimize a simple quadratic function.
* `Driver_SD_WW_vs_HZLS` compares steepest descent with classical weak Wolfe conditions against steepest descent using the Hager--Zhang line search technique.

## References:
- [1] Hager and Zhang, "A new conjugate gradient method with guaranteed descent and an efficient line search", SIAM J. Optim., 16 (2005), pp. 170–192.
- [2] Hager and Zhang, "Algorithm 851: CG DESCENT, a Conjugate Gradient Method with Guaranteed Descent", ACM Trans. Math. Softw., 32 (2006), pp. 113–137.
