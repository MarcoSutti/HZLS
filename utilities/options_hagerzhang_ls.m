%   +-----------------------------------------------------------------+
%   | Options for Hager-Zhang line-search.                            |
%   +-----------------------------------------------------------------+
%   | This script collects all the parameters and options used in the |
%   | Hager-Zhang [HZ] line search procedure.                         |
%   | Created:     16.10.2019                                         |
%   | Last change: 28.04.2020                                         |
%   +-----------------------------------------------------------------+

% warning off;

% c_1 Wolfe sufficient decrease condition
options_hzls.delta = 0.1;   % (Nocedal & Wright recommends 0.01?)
                            % The value from HZ paper is 0.1.

% c_2 Wolfe curvature condition (Recommend 0.1 for GradientDescent)
options_hzls.sigma = 0.9;   % (Nocedal & Wright recommends 0.1 for GradientDescent)
                            % The value from HZ paper is 0.9.
% MS: the value 0.1 yields a more regular convergence behavior of the
%     gradient norm.

options_hzls.epsilon = 1e-6;    % The value from HZ2005, p. 186.
options_hzls.alphamax = 100;
options_hzls.rho = 5;           % Hager and Zhang propose rho = 5
options_hzls.gamma = 0.66;      % gamma is the decay factor for the bracketing interval (2/3)
options_hzls.linesearchmax = 100;
options_hzls.display = false;

%--------------------------------------------------------------------------
% Parameters for initialguess routine
%--------------------------------------------------------------------------
options_hzls.psi0 = 0.01; 
options_hzls.psi1 = 0.2;
% psi2 is the coefficient used to define the initial guess for c when it is
% not the very first iteration. HZ, and also the authors of the Julia code,
% use psi2 = 2.0.
options_hzls.psi2 = 2.0;
options_hzls.psi3 = 0.1;
options_hzls.quadstep = false;