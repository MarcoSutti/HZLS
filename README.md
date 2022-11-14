# HZLS
MATLAB implementation of the Hager-Zhang line-search technique.

Line search by Hager and Zhang, which uses the approximate Wolfe conditions.
This MATLAB version of the Hager-Zhang bracketing has been implemented by
following the existing Julia code available at:
https://julianlsolvers.github.io/LineSearches.jl/stable/index.html

References: [1] Hager and Zhang, A new conjugate gradient method with
                 guaranteed descent and an efficient line search, SIAM
                 J. Optim., 16 (2005), pp. 170–192.
            [2] Hager and Zhang, Algorithm 851: CG DESCENT, a Conjugate
                 Gradient Method with Guaranteed Descent, ACM Trans. Math.
                 Softw., 32 (2006), pp. 113–137.