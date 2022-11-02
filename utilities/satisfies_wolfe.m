function [ outputArg ] = satisfies_wolfe( c, phi_c, dphi_c, phi_0, dphi_0, phi_lim, options_hzls )

% function [ outputArg ] = satisfies_wolfe( c, phi_c, dphi_c, phi_0, dphi_0, phi_lim, options_hzls )
% Purpose: Check Wolfe & approximate Wolfe.
% Created:     27.09.2019
% Last change: 06.12.2019

% Equation (22): original Wolfe conditions (T1)
wolfe1 = ( options_hzls.delta * dphi_0 >= (phi_c - phi_0) / c ) && ( dphi_c >= options_hzls.sigma * dphi_0 );

% Equation (23): Approximate Wolfe conditions (T2)
wolfe2 = ( ( (2 * options_hzls.delta - 1) * dphi_0 >= dphi_c ) && ( dphi_c >= options_hzls.sigma * dphi_0 ) ) && ( phi_c <= phi_lim );

outputArg = ( wolfe1 || wolfe2 );

end