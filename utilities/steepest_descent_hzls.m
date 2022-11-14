function [ xk, info ] = steepest_descent_hzls( xk, problem, options_sd, options_hzls )

% function [ xk, info ] = steepest_descent_hzls( xk, problem, options_sd, options_hzls )
% Purpose: Performs Euclidean steepest descent using the Hager-Zhang line-
%          search procedure (HZLS).

% Created:     2020.05.01
% Last change: 2022.11.02

%   May 3, 2020:
%       Added field nfeval in the info struct.

% NB: The length is options_sd.maxiter+1 because we take into account also
%     the values corresponding to the 0th iteration.
info.cost = zeros( 1, options_sd.maxiter+1 );
info.gradnorm = zeros( 1, options_sd.maxiter+1 );
info.opt_variable = zeros( 1, options_sd.maxiter+1 );

% Added 03.05.2020:
info.nfeval = zeros( 1, options_sd.maxiter );

norm_xex = norm( problem.xex, 'fro' );

if options_sd.verbosity >= 1
    fprintf('+-----------------------------------------------------------------+\n' );
    fprintf('|                          Steepest Descent                       |\n');
    fprintf('|                    with Hager-Zhang line search                 |\n');
    fprintf('+-----------------------------------------------------------------+\n');
end

fk = problem.cost( xk );
gk = problem.egrad( xk );

info.cost(1) = fk;
info.gradnorm(1) = norm( gk, 'fro' );
info.opt_variable(1) = norm( problem.xex - xk, 'fro' )/norm_xex;

% Display iteration information.
if options_sd.verbosity >= 2
    fprintf(' iter\t               cost val\t    grad. norm\t  ||x-x*||/||x*||\n');
    fprintf('%5d\t%+.16e\t%.8e\t   %.8e\n', 0, info.cost(1), info.gradnorm(1), info.opt_variable(1) );
end

tic

for iter=2:options_sd.maxiter+1
    
    gk = problem.egrad( xk );
    
    % Set steepest descent direction to minus gk
    dk = -gk;
    
    % Perform Hager-Zhang line search
    [ tmin, stat_hzls ] = helper_hzls( problem, xk, dk, options_hzls );
    
    info.nfeval(iter-1) = stat_hzls.nf;
    
    % Update to new xk
    xk = xk + tmin * dk;
    
    info.cost(iter) = problem.cost( xk );
    info.gradnorm(iter) = norm( gk, 'fro' );
    info.opt_variable(iter) = norm( problem.xex - xk, 'fro' )/norm_xex;
    
    % Display iteration information
    if options_sd.verbosity >= 2
        fprintf('%5d\t%+.16e\t%.8e\t   %.8e\n', iter-1, info.cost(iter), info.gradnorm(iter), info.opt_variable(iter) );
    end
    
end

time_end = toc;

if iter == options_sd.maxiter+1
    output.msg = 'maximum number of sd iterations reached';
end
%--------------------------------------------------------------------------

if options_sd.verbosity >= 1
    fprintf('Termination message: %s.\n', output.msg);
    fprintf('Total time: %.2f s.\n', time_end );
    fprintf('Last grad. norm: %.4e.\n', info.gradnorm(iter) );
    fprintf('Total number of function evaluations: %.d.\n', sum(info.nfeval) );
end

end