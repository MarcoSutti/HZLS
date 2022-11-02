function [ xk, info ] = steepest_descent_ww( xk, options_sd, pars )

% function [ xk, info ] = steepest_descent_ww( xk, options_sd, pars )
% Purpose: Performs Euclidean steepest descent.
%          This version uses the line search function 'linesch_ww' provided
%          in the HANSO 2.0 code by Michael Overton.

% Created:     2020.05.01
% Last change: 2022.11.02

%   Aug 29, 2020:
%       Added Goldstein-Price function.
%   Aug 19, 2020:
%       Added Rosenbrock's function.
%   May 3, 2020:
%       Added field nfeval in the info struct.

% NB: The length is options_sd.maxiter+1 because we take into account also
%     the values corresponding to the 0th iteration.
info.cost = zeros( 1, options_sd.maxiter+1 );
info.gradnorm = zeros( 1, options_sd.maxiter+1 );
info.opt_variable = zeros( 1, options_sd.maxiter+1 );

% Added 03.05.2020:
info.nfeval = zeros( 1, options_sd.maxiter );

norm_xex = norm( pars.xex, 'fro' );

if options_sd.verbosity >= 1
    fprintf('+-----------------------------------------------------------------+\n' );
    fprintf('|                          Steepest Descent                       |\n');
    fprintf('|                with original weak Wolfe conditions              |\n');
    fprintf('+-----------------------------------------------------------------+\n');
end

if strcmp(pars.fgname,'easy_quadratic')
    [ fk, gk ] = easy_quadratic( xk, pars );
elseif strcmp(pars.fgname,'rosenbrock_function')
    [ fk, gk ] = rosenbrock_function( xk, pars );
elseif strcmp(pars.fgname,'goldstein_price')
    [ fk, gk ] = goldstein_price( xk, pars );
end

info.cost(1) = fk;
info.gradnorm(1) = norm( gk, 'fro' );
info.opt_variable(1) = norm( pars.xex - xk, 'fro' )/norm_xex;

% Display iteration information.
if options_sd.verbosity >= 2
    fprintf(' iter\t               cost val\t    grad. norm\t  ||x-x*||/||x*||\n');
    fprintf('%5d\t%+.16e\t%.8e\t   %.8e\n', 0, info.cost(1), info.gradnorm(1), info.opt_variable(1) );
end

tic

for iter=2:options_sd.maxiter+1
    
    if strcmp(pars.fgname,'easy_quadratic')
        [ fk, gk ] = easy_quadratic( xk, pars );
    elseif strcmp(pars.fgname,'rosenbrock_function')
        [ fk, gk ] = rosenbrock_function( xk, pars );
    elseif strcmp(pars.fgname,'goldstein_price')
        [ fk, gk ] = goldstein_price( xk, pars );
    end
    
    % Set steepest descent direction to minus gk
    dk = -gk;
    
    % Use the line search from HANSO 2.0 by Michael Overton
    % [alpha, xk, fk, gk, fail, beta, gradbeta, fevalrec, nfeval] = linesch_ww(xk, fk, gk, dk, pars, 1e-4, 0.5, -inf, 1);
    [~, xk, fk, gk, ~, ~, ~, ~, nfeval] = linesch_ww(xk, fk, gk, dk, pars, 1e-4, 0.5, -inf, 1);
    
    info.nfeval(iter-1) = nfeval;
    
    info.cost(iter) = fk;
    info.gradnorm(iter) = norm( gk, 'fro' );
    info.opt_variable(iter) = norm( pars.xex - xk, 'fro' )/norm_xex;
    
    % Display iteration information
    if options_sd.verbosity >= 2
        fprintf('%5d\t%+.16e\t%.8e\t   %.8e\n', iter-1, info.cost(iter), info.gradnorm(iter), info.opt_variable(iter) );
    end
    
end

time_end = toc;

if iter == options_sd.maxiter+1
    output.msg = 'maximum number of mg iterations reached';
end
%--------------------------------------------------------------------------

if options_sd.verbosity >= 1
    fprintf('Termination message: %s.\n', output.msg);
    fprintf('Total time: %.2f s.\n', time_end );
    fprintf('Last grad. norm: %.4e.\n', info.gradnorm(iter) );
    fprintf('Total number of function evaluations: %.d.\n', sum(info.nfeval) );
end

end