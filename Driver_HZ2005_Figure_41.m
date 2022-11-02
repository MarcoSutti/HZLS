%==========================================================================
% Driver for reproducing Figure 4.1 in [1].
% References: [1] Hager and Zhang, A new conjugate gradient method with
%                 guaranteed descent and an efficient line search, SIAM
%                 Journal on Optimization, 16 (2005), pp. 170???192.
% Created:     29.08.2020
% Last change: 29.08.2020
%==========================================================================
close all; clear; clc;

options_plot;
%--------------------------------------------------------------------------

a = 1 - 2.5e-8;
b = 1 + 2.5e-8;


h = (b - a)/10000;


x=a:h:b;
F = 1-2.*x+x.^2;

plot( x, F, 'k-', 'LineWidth', 2 )