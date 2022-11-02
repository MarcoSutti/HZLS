function [ state, stat_hzls ] = initialguess( problem, x, Dg, options_hzls )

% function [ state, stat_hzls ] = initialguess( x, s, options_hzls )
% Purpose: Implements the routine [c] = initial(k) to generate the starting
%          guess c used by the bracket routine, as described in HZ paper,
%          p. 124.
% References: [1] Hager and Zhang, Algorithm 851: CG DESCENT, a Conjugate
%                 Gradient Method with Guaranteed Descent, ACM Trans. Math.
%                 Softw., 32 (2006), pp. 113–137.
%             [2] https://julianlsolvers.github.io/LineSearches.jl/stable/index.html

% If alpha0 is NaN, then procedure I0 is called at the first iteration,
% otherwise, we select according to procedure I1-2, with starting value alpha0.

% Created:     2019.10.11
% Last change: 2022.10.24

%   Oct 24, 2022:
%       Adaptation to the code for Riem. CG.
%   May 27, 2020:
%       Replaced function_phi and function_phidphi with
%       feval(problem.eval_phi, ... ) and feval(problem.eval_phidphi, ...).
%   May 23, 2020:
%       Code cleanup.

% % MS, 20.11.2019: Introduced these global variables to use quadstep.
% global info_global;
% global column_idx;
% global current_lev;x.Dc

% Statistics for line-search (MS, 14.10.2019)
% If we are here, it means that we just started a brand new steepest
% descent iteration, so it is the moment to reset the counter for the
% number of function evaluations, stat_hzls.nf
stat_hzls.nf = 0;

% if isnan(info_global(current_lev,column_idx))
% If we're at the first iteration (info_global(current_lev,column_idx)
% is NaN), then we pick the initial step size according to HZ, I0

% MS, 26.11.2019: Debugging:
% fprintf('Debugging message: Entering hzI0.\n' );
state.alpha = hzI0( problem, x, options_hzls, stat_hzls );

%**********************************************************************
% MS, 26.11.2019: Check that dphi_0 is strictly negative
% [ ~, dphi_0, stat_hzls ] = feval(problem.eval_phidphi, problem, x, s, state.alpha, stat_hzls );

% 2022.10.24
[ ~, dphi_0 ] = feval( options_hzls.eval_phidphi, problem, x, Dg, state.alpha );

% If dphi_0 is not strictly negative, then we do this kind of hack to
% force the code to return an alpha such that dphi_0 < 0.
if dphi_0 >= 0
    while dphi_0 >= 0
        state.alpha = 0.5 * state.alpha;
        % [ ~, dphi_0, stat_hzls ] = feval(problem.eval_phidphi, problem, x, Dg, state.alpha, stat_hzls );

        % 2022.10.20
        [ ~, dphi_0 ] = feval( options_hzls.eval_phidphi, problem, x, Dg, state.alpha );
        % [ ~, dphi_0 ] = eval_phidphi( x, Dg, state.alpha );

        if state.alpha == 0
%             dphi_0
            state.alpha = eps;
            return;
        end
    end
end
% We need to have dphi_0 < 0 before proceeding further:
assert( dphi_0 < 0 );
%**********************************************************************

% Keep track that the first iteration has been done:
%     info_global(current_lev,1) = 1;

% else
%     % If we are here, it means that this is not the first iteration of the
%     % optimization method, so we can reuse the alpha0_global from the previous iteration
%     state.alpha = info_global(current_lev,column_idx);
%     [ phi_0, dphi_0 ] = eval_phidphi(Da, Db, Dg, state.alpha);
%
%     % Pick the initial step size according to HZ #I1-2
%     [ state.alpha, stat_hzls ] = hzI12( state.alpha, problem, Da, Db, Dg, phi_0, dphi_0, options_hzls, stat_hzls );
% end

end

%--------------------------------------------------------------------------
% OTHER FUNCTIONS
%--------------------------------------------------------------------------
% I1
%--------------------------------------------------------------------------
% Generate initial guess for step size (HZ, stage I0)
function [ alpha ] = hzI0( problem, x, options_hzls, stat_hzls )

if options_hzls.display
    fprintf( 'I0\n' );
end

f_x = feval( options_hzls.eval_phi, problem, x, 0, 0 );
% f_x = eval_LS( x, 0, 0);

% Increase counter for the number of function evaluations
stat_hzls.nf = stat_hzls.nf + 1;

alpha = options_hzls.psi0 * abs(f_x);

alpha = min(alpha, options_hzls.alphamax);

end


%--------------------------------------------------------------------------
% Pick the initial step size (HZ I1-I2)
%--------------------------------------------------------------------------
function [ alphatest, stat_hzls ] = hzI12( alpha, x, Dg, phi_0, dphi_0, options_hzls, stat_hzls )

if options_hzls.display
    fprintf( 'I1-2\n' );
end

% Prevent values of 'x_new' that are likely to make phi(x_new) infinite
iterfinitemax = ceil(-log2(eps));

alphatest = options_hzls.psi1 * alpha;   % MS: This is psi1*alpha_{k-1} in HZ paper.
alphatest = min(alphatest, options_hzls.alphamax);

% Reminder: df is problem.cost and s is given by -problem.grad(xk), xk = x
% MS: this is denoted as phi(psi1*alpha_{k-1}) in HZ paper.
% [ phitest, stat_hzls ] = feval(problem.eval_phi, problem, x, s, alphatest, stat_hzls );

% 2022.10.24
[ phitest, ~ ] = feval(options_hzls.eval_phi, problem, x, Dg, alphatest );

[ alphatest, phitest ] = get_finite( alphatest, phitest, options_hzls.psi3, iterfinitemax, phi_0, x, Dg );

quadstep_success = false;
%--------------------------------------------------------------------------
% I1
%--------------------------------------------------------------------------
if options_hzls.display
    fprintf( 'I1\n' );
end
if options_hzls.quadstep   % If quadstep is true...

    a = ((phitest - phi_0) / alphatest - dphi_0) / alphatest;  % quadratic fit

    if options_hzls.display
        fprintf( 'quadfit: alphatest =  %.4e, phi_0 = %.4e, dphi_0 = %.4e, phitest = %.4e, quadcoef = %.4e.\n', ...
            alphatest, phi_0, dphi_0, phitest, a );
    end

    % ... and if phi(psi1*alpha_{k-1}) ≤ phi(0), and the quadratic
    % interpolant q() is strongly convex (i.e., a > 0), then ...
    if isfinite(a) && a > 0 && phitest <= phi_0
        % ... choose minimum of quadratic interpolant q()
        alphatest2 = -dphi_0 / 2 / a;
        alphatest2 = min(alphatest2, options_hzls.alphamax);
        % [ phitest2, stat_hzls ] = feval(problem.eval_phi, problem, x, s, alphatest2, stat_hzls );
        % 2022.10.24
        phitest2 = feval( options_hzls.eval_phidphi, problem, x, Dg, alphatest2 );
        if isfinite(phitest2)
            alphatest = alphatest2;
            phitest = phitest2;
            if options_hzls.display
                fprintf('alpha guess (quadratic): %.4e.\n', alphatest );
            end
        end
    end
end
%--------------------------------------------------------------------------
% I2
%--------------------------------------------------------------------------
if ( ~options_hzls.quadstep || ~quadstep_success ) && ( phitest <= phi_0 )
    % If no quadstep or it fails, expand the interval.
    % While the phitest <= phi_0 condition was not in the paper, it gives
    % a significant boost to the speed. The rationale behind it is that
    % since the slope at alpha = 0 is negative, if phitest > phi_0 then a
    % local minimum must be between alpha = 0 and alpha = alphatest, so
    % alpha_test is good enough to return.
    alphatest = options_hzls.psi2 * alpha;   % MS: This is c = psi2*alpha_{k-1} in HZ paper.
    alphatest = min(alphatest, options_hzls.alphamax);

    % [ phitest, stat_hzls ] = feval(problem.eval_phi, problem, x, s, alphatest, stat_hzls );

    % 2022.10.24
    phitest = feval( options_hzls.eval_phi, problem, x, Dg, alphatest );

    [ alphatest, ~ ] = get_finite( alphatest, phitest, options_hzls.psi3, iterfinitemax, phi_0, x, Dg );
    if options_hzls.display
        fprintf('alpha guess (expand): %.4e.\n', alphatest );
    end
end
% return alphatest
end   % end of function hzI12


%--------------------------------------------------------------------------
% Function get_finite
%--------------------------------------------------------------------------
function [ alpha, phi ] = get_finite( alpha, phi, psi3, itermax, phi_0, x, Dg )

% Returns alpha, phi.

iter = 1;
while ~isfinite(phi)
    if iter >= itermax
        phi = phi_0;
        return;   % returns phi_0
    end
    alpha = psi3 * alpha;

    % 2022.10.24
    phi = feval( options_hzls.eval_phi, problem, x, Dg, alpha );

    iter = iter + 1;
end

end