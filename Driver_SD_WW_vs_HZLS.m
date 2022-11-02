%==========================================================================
% Driver for comparing two versions of steepest descent which use two
% different line-search techniques:
%    1. Line search with weak Wolfe conditions: this version uses the line
%       search routine "linesch_ww" provided by Michael Overton in his
%       HANSO 2.0 code package;
%    2. Line search by Hager and Zhang, which uses the approximate Wolfe
%       conditions. The Hager-Zhang line search has been implemented by
%       following the Julia code on:
%       https://julianlsolvers.github.io/LineSearches.jl/stable/index.html

% Use this script as it is provided to generate Fig. 5 in:
%       Riemannian multigrid line search for low-rank problems,
%       M. Sutti and B. Vandereycken, SIAM J. Sci. Comput., 43(3), A1803–A1831, 2021.
%       https://epubs.siam.org/doi/10.1137/20M1337430

% References: [1] Hager and Zhang, A new conjugate gradient method with
%                 guaranteed descent and an efficient line search, SIAM
%                 Journal on Optimization, 16 (2005), pp. 170–192.
%             [2] Hager and Zhang, Algorithm 851: CG DESCENT, a Conjugate
%                 Gradient Method with Guaranteed Descent, ACM Trans. Math.
%                 Softw., 32 (2006), pp. 113–137.

% Created:     2020.05.01
% Last change: 2022.11.02

%   Aug 29, 2020:
%       Added Goldstein-Price function.
%   Aug 19, 2020:
%       Added Rosenbrock's function.
%   May 27, 2020:
%       Added 'problem.eval_phi' and 'problem.eval_phidphi'.
%   May 23, 2020:
%       Code cleanup.
%   May 20, 2020:
%       The struct 'options_hzls' is now assigned as a field of the struct
%       'problem'.
%==========================================================================
close all; clear; clc;

% Add folder and its subfolders to MATLAB path for the current session:
addpath(genpath('../HZLS'))

%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------
% Fix stream of random number generator:
rng(1)
%--------------------------------------------------------------------------
% Choose function type: 'easy_quadratic', 'rosenbrock_function'
%                       'goldstein_price'
pars.fgname = 'easy_quadratic';
% We use a modified Rosenbrock function in which the optimal cost is
% nonzero; this is in the same order of ideas as the remark at p. 189-190
% in reference [1].

z0 = 1;
%--------------------------------------------------------------------------
% Only for quadratic function
% Choose variable type: 'vector', 'matrix'
pars.var_type = 'matrix';
%--------------------------------------------------------------------------
% Only for Rosenbrock function
% Choose Rosenbrock function dimension: dim = 2 for 2D, dim = 3 for 3D:
pars.dim = 2;
%--------------------------------------------------------------------------
n = 100;       % size of the problem
cond_A = 10;   % condition number of matrix A
%--------------------------------------------------------------------------
% Options for steepest descent:
options_sd.maxiter = 200;     % number of steepest descent iterations
options_sd.verbosity = 1;    % verbosity of the output of steepest descent
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
% End of data
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% DEFINE PROBLEM STRUCTURE
%--------------------------------------------------------------------------
% Definition of cost function and gradient
if strcmp(pars.fgname,'easy_quadratic')
    
    % Build the matrix A
    [Q,~] = qr(rand(n));
    S = diag(logspace(0,log10(1/cond_A),n));
    A = Q*S*Q';
    
    % Prescribe the exact solution
    if strcmp(pars.var_type,'vector')
        xex = rand(n,1);
    elseif strcmp(pars.var_type,'matrix')
        xex = rand(n,n);
    end
    
    % Build the matrix B from the exact solution
    B = A*xex;
    %----------------------------------------------------------------------
    % Store xex, A and B in a struct that is passed to steepest_descent_hzls
    pars.xex = xex;
    pars.A = A;
    pars.B = B;
    %----------------------------------------------------------------------
    % Get the reference value of Fh
    [ Fh_ref, ~ ] = easy_quadratic( xex, pars );
    %----------------------------------------------------------------------
    % Initial guess
    if strcmp(pars.var_type,'vector')
        x0 = rand(n,1);
    elseif strcmp(pars.var_type,'matrix')
        x0 = rand(n,n);
    end
    
    % MS, added 25.01.2021:
    [ F0, ~ ] = easy_quadratic( x0, pars );
    
    % fileName = ['initial_random_variables.mat'];
    % save( fileName, 'A', 'xex', 'x0' )
    %----------------------------------------------------------------------
    if strcmp(pars.var_type,'vector')
        problem.cost = @(x)((0.5*A*x - B)'*x);
        problem.egrad = @(x)(A*x - B);
    elseif strcmp(pars.var_type,'matrix')
        problem.cost = @(X) 0.5*trace(X' * A * X) - trace(X' * B);
        problem.egrad = @(X) A*X - B;
    end
    problem.xex = xex;
    
elseif strcmp(pars.fgname,'rosenbrock_function')
    if pars.dim == 2
        % Rosenbrock function, dimension dim = 2
        % We use a modified Rosenbrock function in which the optimal cost is nonzero:
        problem.cost = @(x) 100*(x(2)-x(1)^2)^2+(1-x(1))^2 + z0;
        problem.egrad = @(x) [ 2*(200*x(1)^3-200*x(1)*x(2)+x(1)-1);
            200*(x(2)-x(1)^2) ];
    elseif pars.dim == 3
        % Rosenbrock function, dimension dim = 3
        problem.cost = @(x) 100*(x(2)-x(1)^2)^2+(1-x(1))^2 + 100*(x(3)-x(2)^2)^2+(1-x(2))^2;
        problem.egrad = @(x) [ 2*(200*x(1)^3-200*x(1)*x(2)+x(1)-1);
            -200*x(1)^2 + 400*x(2)^3 + x(2)*(202 - 400*x(3))-2;
            200*(x(3)-x(2)^2) ];
    end
    problem.rgrad = @(x) problem.M.egrad2rgrad( x, problem.egrad(x) );
    
    % Exact solution:
    xex = ones(pars.dim,1);
    problem.xex = xex;
    pars.xex = xex;
   
    pars.z0 = z0;
    
    %----------------------------------------------------------------------
    % Get the reference value of Fh
    [ Fh_ref, ~ ] = rosenbrock_function( xex, pars );
    %----------------------------------------------------------------------
    
    % Generate initial guess close to the exact solution:
    x0 = xex + rand(pars.dim,1)/2.25e5;
    
    % MS, added 25.01.2021:
    [ F0, ~ ] = rosenbrock_function( x0, pars );
    
elseif strcmp(pars.fgname,'goldstein_price')
    problem.cost = @(x) (1 + (x(1)+x(2)+1)^2 * (19-14*x(1)+3*x(1)^2-14*x(2)+6*x(1)*x(2)+3*x(2)^2) ) ...
        * (30 + (2*x(1)-3*x(2))^2*(18-32*x(1)+12*x(1)^2+48*x(2)-36*x(1)*x(2)+27*x(2)^2) );
    problem.egrad = @(x) [ 24*(8*x(1)^3 - 4*x(1)^2*(9*x(2) + 4) + 6*x(1)*(9*x(2)^2 + 8*x(2) + 1) ...
        - 9*x(2)*(3*x(2)^2 + 4*x(2) + 1))*((3*x(1)^2 + 2*x(1)*(3*x(2) - 7) ...
        + 3*x(2)^2 - 14*x(2) + 19)*(x(1) + x(2) + 1)^2 + 1) + 12*(x(1)^3 + x(1)^2*(3*x(2) - 2) ...
        + x(1)*(3*x(2)^2 - 4*x(2) - 1) + x(2)^3 - 2*x(2)^2 - x(2) + 2)*((12*x(1)^2 - 4*x(1)*(9*x(2) + 8) ...
        + 3*(9*x(2)^2 + 16*x(2) + 6))*(2*x(1) - 3*x(2))^2 + 30);
        12*(x(1)^3 + x(1)^2*(3*x(2) - 2) + x(1)*(3*x(2)^2 - 4*x(2) - 1) ...
        + x(2)^3 - 2*x(2)^2 - x(2) + 2)*((12*x(1)^2 - 4*x(1)*(9*x(2) + 8) ...
        + 3*(9*x(2)^2 + 16*x(2) + 6))*(2*x(1) - 3*x(2))^2 + 30) - 36*(8*x(1)^3 ...
        - 4*x(1)^2*(9*x(2) + 4) + 6*x(1)*(9*x(2)^2 + 8*x(2) + 1) ...
        - 9*x(2)*(3*x(2)^2 + 4*x(2) + 1))*((3*x(1)^2 + 2*x(1)*(3*x(2) - 7) ...
        + 3*x(2)^2 - 14*x(2) + 19)*(x(1) + x(2) + 1)^2 + 1) ];
    % problem.rgrad = @(x) problem.M.egrad2rgrad( x, problem.egrad(x) );
    
    % Exact solution:
    xex = [ 0; -1 ];
    problem.xex = xex;
    pars.xex = xex;
    
    %----------------------------------------------------------------------
    % Get the reference value of Fh
    % [ Fh_ref, ~ ] = goldstein_price( xex );
    Fh_ref = 3;
    %----------------------------------------------------------------------
    
    % Generate initial guess close to the exact solution:
    x0 = xex + rand(2,1);
    
     % MS, added 25.01.2021:
    [ F0, ~ ] = goldstein_price( x0, pars );
end

% MS, added 25.01.2021:
pars.F0 = F0;

%--------------------------------------------------------------------------
% LAUNCH STEEPEST DESCENT WITH WEAK WOLFE CONDITIONS
%--------------------------------------------------------------------------
[ ~, info_WWLS ] = steepest_descent_ww( x0, options_sd, pars );
%--------------------------------------------------------------------------
% LAUNCH STEEPEST DESCENT WITH APPROXIMATE WOLFE CONDITIONS
%--------------------------------------------------------------------------
[ ~, info_HZLS ] = steepest_descent_hzls( x0, problem, options_sd, options_hzls );
%--------------------------------------------------------------------------
% Plot results:
plot_sd_ww_vs_hzls( info_WWLS, info_HZLS, Fh_ref, pars );
