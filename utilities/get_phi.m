function [ phi_c ] = get_phi( problem, x, s, c )

% function [ phi_c, stat_hzls ] = get_phi( problem, x, s, c )
% Purpose: Returns phi_c = f(x+c*s).
% Created:     27.05.2020
% Last change: 27.05.2020

%   May 27, 2020:
%       This new version replaces the previous one, which made use of
%       'function_phi' and 'function_phidphi'.

% Retraction
x_new = x + c * s;

% Compute phi(c)
phi_c = problem.cost(x_new);

% Increase counter for the number of function evaluations
% stat_hzls.nf = stat_hzls.nf + 1;
%
end