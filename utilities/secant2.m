function [ iswolfe, ia, ib, alphas, values, slopes ] = secant2( alphas, values, slopes, ia, ib, phi_lim, x, s, problem, options_hzls )

% function [ iswolfe, ia, ib, alphas, values, slopes, stat_hzls ] = secant2( alphas, values, slopes, ia, ib, phi_lim, x, Dg, options_hzls )
% Purpose: Implementation of [\bar{a},\bar{b}] = secant^{2}(a,b).
% References: [1] Hager and Zhang, Algorithm 851: CG DESCENT, a Conjugate
%                 Gradient Method with Guaranteed Descent, ACM Trans. Math.
%                 Softw., 32 (2006), pp. 113–137.
%             [2] https://julianlsolvers.github.io/LineSearches.jl/stable/index.html
% Created:     2019.09.27
% Last change: 2022.10.24

%   Oct 24, 2022:
%       Adaptation to the code for Riem. CG.
%   May 27, 2020:
%       Replaced function_phidphi with feval(problem.eval_phidphi, ...).
%   May 23, 2020:
%       Code cleanup.

% MS, 27.11.2019: Introduced these global variables to use quadstep.
% global info_global;
% global column_idx;
% global current_lev;

phi_0 = values(1);
dphi_0 = slopes(1);
a = alphas(ia);
b = alphas(ib);
dphi_a = slopes(ia);
dphi_b = slopes(ib);

if not( dphi_a < 0 && dphi_b >= 0 )
    warning('The opposite slope condition is not satisfied: dphi(a) = %.4e; dphi(b) = %.4e.', dphi_a, dphi_b );
end
%--------------------------------------------------------------------------
% S1
%--------------------------------------------------------------------------
c = secant_formula(a, b, dphi_a, dphi_b);
if options_hzls.display
    fprintf('S1: a = %.3e, b = %.3e, c = %.3e.\n', a, b, c);
end
assert( isfinite(c) );

% 2022.10.20
% [ phi_c, dphi_c ] = eval_phidphi( x, Dg, c );
[ phi_c, dphi_c ] = feval( options_hzls.eval_phidphi, problem, x, s, c );

assert( isfinite(phi_c) && isfinite(dphi_c) );

alphas = [alphas; c];
values = [values; phi_c];
slopes = [slopes; dphi_c];

ic = length(alphas);

if satisfies_wolfe( c, phi_c, dphi_c, phi_0, dphi_0, phi_lim, options_hzls )
    if options_hzls.display
        fprintf('S1: first c satisfied Wolfe conditions.\n');
    end
        
    % MS, 27.11.2019: Save the accepted c for future use.
    % info_global(column_idx) = c;
    
    iswolfe = true;
    ia = ic;
    ib = ic;
    return;   % return true, ic, ic
end

if options_hzls.display
    fprintf('S1: The opposite slope condition is enforced: dphi(a) = %.4e; dphi(b) = %.4e.\n', slopes(ia), slopes(ib) );
    fprintf('S1: Values of phi: phi(a) = %.4e; phi(b) = %.4e.\n', values(ia), values(ib) );
    fprintf('S1: update_hz is called by secant2\n')
end

[ iA, iB, alphas, values, slopes ] = update_hz(alphas, values, slopes, ...
    ia, ib, ic, phi_lim, x, s, options_hzls );
if options_hzls.display
    fprintf('S1: iA = %d, iB = %d, ic = %d.\n', iA, iB, ic);
end
a = alphas(iA);
b = alphas(iB);
%--------------------------------------------------------------------------
% S2
%--------------------------------------------------------------------------
if iB == ic
    % we updated b, make sure we also update a
    c = secant(alphas, slopes, ib, iB);
%--------------------------------------------------------------------------
% S3
%--------------------------------------------------------------------------
elseif iA == ic
    % we updated a, do it for b too
    c = secant(alphas, slopes, ia, iA);
end
%--------------------------------------------------------------------------
% S4
%--------------------------------------------------------------------------
if (iA == ic || iB == ic) && ( ( a <= c) && ( c <= b ) )
    if options_hzls.display
        fprintf('S4: second c = %.3e.\n', c);
    end

    % 2022.10.24
    [ phi_c, dphi_c ] = feval( options_hzls.eval_phidphi, problem, x, s, c );

    assert( isfinite(phi_c) && isfinite(dphi_c) );

    alphas = [alphas; c];
    values = [values; phi_c];
    slopes = [slopes; dphi_c];
    
    ic = length(alphas);

    if satisfies_wolfe( c, phi_c, dphi_c, phi_0, dphi_0, phi_lim, options_hzls )
        if options_hzls.display
            fprintf('S4: second c satisfied Wolfe conditions.\n')
        end
        
        % MS, 27.11.2019: Save the accepted c for future use.
        % info_global(current_lev,column_idx) = c;
        
        iswolfe = true;
        ia = ic;
        ib = ic;
        return;   % return true, ic, ic
    end
    
    if options_hzls.display
        fprintf('S4: update_hz is called by secant2\n');
    end
    [ iA, iB, alphas, values, slopes ] = update_hz(alphas, values, slopes, ...
        iA, iB, ic, phi_lim, x, s, options_hzls );
end
if options_hzls.display
    fprintf('secant2 output: a = %.4e, b = %.4e.\n', alphas(iA), alphas(iB) );
    fprintf('                dphi(a) = %.4e; dphi(b) = %.4e.\n', slopes(iA), slopes(iB) );
end
iswolfe = false;
ia = iA;
ib = iB;
return;   % return false, iA, iB
end

%--------------------------------------------------------------------------
% OTHER FUNCTIONS
%--------------------------------------------------------------------------
function [ c ] = secant( alphas, slopes, ia, ib )

% Created:     27.09.2019
% Last change: 27.09.2019

c = secant_formula( alphas(ia), alphas(ib), slopes(ia), slopes(ib) );

end

function [ c ] = secant_formula( a, b, dphi_a, dphi_b )

% Purpose: Secant formula (p. 123 of HZ).
% Created:     27.09.2019
% Last change: 27.09.2019

c = (a * dphi_b - b * dphi_a) / (dphi_b - dphi_a);

end