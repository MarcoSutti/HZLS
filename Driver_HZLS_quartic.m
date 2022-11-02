%==========================================================================
% Driver for steepest descent with Hager-Zhang line search.

% Created:     2021.02.17
% Last change: 2022.11.02

%   Nov 2, 2022:
%       Deep cleanup.
%   Feb 17, 2021:
%       Created.
%==========================================================================
close all; clear; clc;

% Add folder and its subfolders to MATLAB path for the current session:
addpath(genpath('../HZLS_code_scorporato'))

%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------
% Fix stream of random number generator:
rng(1)
%--------------------------------------------------------------------------
% Options for steepest descent:
options_sd.maxiter = 9;     % number of steepest descent iterations
options_sd.verbosity = 2;    % verbosity of the output of steepest descent
                             % verbosity levels: 0, 1, 2.
%--------------------------------------------------------------------------
% Load Hager-Zhang linesearch options
% This loads the struct 'options_hzls' with the default parameters for
% using Hager-Zhang line search. If you wish to change the default
% parameters, you need to do so in 'options_hagerzhang_ls.m'.
options_hagerzhang_ls;
%--------------------------------------------------------------------------
% Define function phi(t) and its derivative for the Hager-Zhang line search
options_hzls.eval_phi = 'get_phi';
options_hzls.eval_phidphi = 'get_phidphi';
%--------------------------------------------------------------------------
% DEFINE PROBLEM STRUCTURE
%--------------------------------------------------------------------------
% Definition of cost function and gradient
problem.cost = @(x) x.^4 - 5 * x.^2 + 20;
problem.egrad = @(x) 4*x.^3 - 10 * x;
%--------------------------------------------------------------------------
% Initial guess
x0 = 2;
%--------------------------------------------------------------------------
% Exact value of the minimizer (for reference):
xex = sqrt(5/2);
% Get the reference value of Fh
Fh_ref = 55/4;
%--------------------------------------------------------------------------
% End of data
%--------------------------------------------------------------------------
pars.Fh_ref = Fh_ref;
pars.xex = xex;
problem.xex = xex;

% MS, added 25.01.2021:
F0 = problem.cost(x0);
pars.F0 = F0;
%--------------------------------------------------------------------------
% LAUNCH STEEPEST DESCENT WITH APPROXIMATE WOLFE CONDITIONS
%--------------------------------------------------------------------------
[ ~, info_HZLS ] = steepest_descent_hzls( x0, problem, options_sd, options_hzls );
%--------------------------------------------------------------------------
% Plot results:
%--------------------------------------------------------------------------
plot_sd_hzls( info_HZLS, Fh_ref )

% info_HZLS.opt_variable
