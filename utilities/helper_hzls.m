function [ c, stat_hzls ] = helper_hzls( problem, x, s, options_hzls )

% function [ c, stat_hzls ] = helper_hzls( problem, x, s, options_hzls )
% Purpose: Wraps the 'hagerzhang_linesearch' routine.

% Created:     2019.10.16
% Last change: 2022.11.02

%   May 20, 2020:
%       The struct 'options_hzls' is now passed as an input parameter to
%       'hagerzhang_linesearch'.

% Perform Hager-Zhang linesearch:
[ alphas, ~, ~, stat_hzls ] = hagerzhang_linesearch( problem, x, s, options_hzls );

% The stepsize 'c' we look for is the last entry in the alphas vector:
c = alphas(end);

if options_hzls.display
    fprintf( 'c = %.4e\n', c );
end

end