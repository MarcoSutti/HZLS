function [ ia, ib, alphas, values, slopes ] = bisect( alphas, values, ...
    slopes, ia, ib, phi_lim, x, s, options_hzls )

% function [ ia, ib, alphas, values, slopes ] = bisect( alphas, values, ...
%     slopes, ia, ib, phi_lim, display, x, s, options_hzls )
% Purpose: Implements stage U3 of the update routine (with theta=0.5).
%          This is the loop encoded in the lines U3a-U3c, p. 123, of:
%             [1] Hager and Zhang, Algorithm 851: CG DESCENT, a Conjugate
%                 Gradient Method with Guaranteed Descent, ACM Trans. Math.
%                 Softw., 32 (2006), pp. 113–137.
% Created:     27.09.2019
% Last change: 2022.10.24

%   Oct 24, 2022:
%       Adaptation to the code for Riem. CG.
%   May 27, 2020:
%       Replaced function_phidphi with feval(problem.eval_phidphi, ...).
%   May 23, 2020:
%       Code cleanup.

a = alphas(ia);
b = alphas(ib);

% Debugging (HZ, conditions shown following U3)
% assert( slopes(ia) < 0 );
% assert( values(ia) <= phi_lim );
% assert( slopes(ib) < 0 );
% assert( values(ib) > phi_lim );
% assert( b > a );

% This is the loop encoded in the lines U3a-U3c of HZ pseudocode.
while b - a > eps(b)
    
    if display
        fprintf('U3a, bisect: a = %.3e, b = %.3e, (b - a) = %3e.\n', a, b, b-a );
    end

    %----------------------------------------------------------------------
    % U3a: set d = ( \bar{a} + \bar{b} )/2
    %----------------------------------------------------------------------
    d = (a + b) / 2;

    % [ phi_d, dphi_d, stat_hzls ] = feval(problem.eval_phidphi, problem, x, Dg, d, stat_hzls );
    % 2022.10.20
    % [ phi_d, dphi_d ] = eval_phidphi(Da, Db, s, d);
    [ phi_d, dphi_d ] = feval( options_hzls.eval_phidphi, x, s, d );

    if display
        fprintf('        phi_d = %.4e,   dphi_d = %.4e.\n', phi_d, dphi_d );
    end

    assert( isfinite(phi_d) && isfinite(dphi_d) );
    
    alphas = [alphas; d];
    values = [values; phi_d];
    slopes = [slopes; dphi_d];
    
    id = length(alphas);
    if dphi_d >= 0   % || abs(dphi_d) < eps   % MS, 15.10.2019: I relaxed this condition
        ib = id;     % replace b, return
        return;
    end
    % dphi < 0, otherwise we would not reach this line of the code
    %----------------------------------------------------------------------
    % U3b: If phi_d <= phi_lim, then \bar{a} = d, and go to U3a
    %----------------------------------------------------------------------
    if phi_d <= phi_lim
        a = d;   % replace a, but keep bisecting until dphi_b > 0
        ia = id;
    else
        %------------------------------------------------------------------
        % U3c: If phi_d > phi_lim, then \bar{b} = d
        %------------------------------------------------------------------
        b = d;
        ib = id;
    end
end

% If bisect terminates without enforcing dphi_d >= 0, throw a warning
% if dphi_d < 0
%     warning('U3 ended without enforcing the condition dphi_d >= 0 !!!!')
% end

end