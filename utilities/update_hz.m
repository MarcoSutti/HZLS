function [ ia, ib, alphas, values, slopes ] = update_hz( alphas, values, slopes, ia, ib, ic, phi_lim, x, s, options_hzls )

% function [ ia, ib, alphas, values, slopes, stat_hzls ] = update_hz( alphas, values, slopes, ia, ib, ic, phi_lim, display, problem, x, s, stat_hzls )
% Purpose: Implements [ \bar{a}, \bar{b} ] = update(a,b,c); see Hager and
%          Zhang, stages U0-U3, p. 123.
%          Given a third point, pick the best two that retain the bracket
%          around the minimum (as defined by HZ, eq. 29) (opposite slope
%          condition). b will be the upper bound, and a the lower bound.
% References: [1] Hager and Zhang, Algorithm 851: CG DESCENT, a Conjugate
%                 Gradient Method with Guaranteed Descent, ACM Trans. Math.
%                 Softw., 32 (2006), pp. 113–137.
%             [2] https://julianlsolvers.github.io/LineSearches.jl/stable/index.html
% Created:     27.09.2019
% Last change: 2022.10.24

%   Oct 24, 2022:
%       Adaptation to the code for Riem. CG.
%   May 23, 2020:
%       Code cleanup.

a = alphas(ia);
b = alphas(ib);
assert( values(ia) <= phi_lim );
assert( slopes(ib) >= 0 );
assert( b > a );
c = alphas(ic);
phi_c = values(ic);
dphi_c = slopes(ic);
if options_hzls.display
    fprintf('update: ia = %d, a = %.3e, ib = %d, b = %.3e, c = %.3e, phi_c = %.3e, dphi_c = %.3e.\n', ia, a, ib, b, c, phi_c, dphi_c);
end
%--------------------------------------------------------------------------
% U0
%--------------------------------------------------------------------------
if c < a || c > b
    if options_hzls.display
        fprintf('U0: c is out of the bracketing interval: a = %.3e, ib = %d, b = %.3e, c = %.3e.\n', a, b, c);
    end
    % keep the same a and b
    return;
end
%--------------------------------------------------------------------------
% U1
%--------------------------------------------------------------------------
if dphi_c >= 0
    ib = ic;   % replace b with a closer point
    return;
end
%--------------------------------------------------------------------------
% U2
% NB: When we are in U2 and U3, the opposite slope condition is naturally
%     NOT satisfied.
% We know dphi_c < 0. However, phi may not be monotonic between a and c, so
% check that the value is also smaller than phi_0. (It's more dangerous to
% replace a than b, since we're leaving the secure environment of alpha=0;
% that's why we didn't check this above.)
%--------------------------------------------------------------------------
if phi_c <= phi_lim
    ia = ic;   % replace a
    return;
else
    %----------------------------------------------------------------------
    % U3
    %----------------------------------------------------------------------
    % phi_c is bigger than phi_0, which implies that the minimum
    % lies between a and c. Find it via bisection.
    
    % Debugging
    assert( dphi_c < 0 );
    assert( phi_c > phi_lim );
    assert( alphas(ia) < alphas(ic) );
    assert( values(ic) > values(ia) );
    
    [ ia, ib, alphas, values, slopes ] = bisect(alphas, values, slopes, ...
        ia, ic, phi_lim, x, s, options_hzls );
end

end