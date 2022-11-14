%==========================================================================
% Driver for reproducing Figure 4.1 in [1].
% References: [1] Hager and Zhang, A new conjugate gradient method with
%                 guaranteed descent and an efficient line search, SIAM
%                 Journal on Optimization, 16 (2005), pp. 170???192.
% Created:     2020.08.29
% Last change: 2020.08.29
%==========================================================================
close all; clear; clc;

% Add folder and its subfolders to MATLAB path for the current session:
addpath(genpath('../HZLS'))

options_plot;
%--------------------------------------------------------------------------

a = 1 - 2.5e-8;
b = 1 + 2.5e-8;


h = (b - a)/10000;


x=a:h:b;
F = 1-2.*x+x.^2;

plot( x, F, '-', 'Color', red, 'LineWidth', 2.5 )