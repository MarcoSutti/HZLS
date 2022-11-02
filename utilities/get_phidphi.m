function [ phi_c, dphi_c ] = get_phidphi( problem, x, s, c )

% function [ phi_c, dphi_c, stat_hzls ] = get_phidphi( problem, x, s, c, stat_hzls )
% Purpose: Returns phi_c = f(x+c*s) and its derivative dphi_c = ∇f(x+c*s).s
% Created:     2020.05.27
% Last change: 2022.11.02

%   May 27, 2020:
%       This new version replaces the previous one, which made use of
%       'function_phi' and 'function_phidphi'.

% Retraction
x_new = x + c * s;

% Compute phi_c
phi_c = problem.cost(x_new);

% Increase counter for the number of function evaluations
% stat_hzls.nf = stat_hzls.nf + 1;

% Compute dphi_c
grad_x_new = problem.egrad(x_new);
dphi_c = dot( s(:), grad_x_new(:) );

end