function [ f, g ] = easy_quadratic( x, pars )

% function [ f, g ] = easy_quadratic( x, pars )
% Created:     01.05.2020
% Last change: 01.05.2020

if strcmp(pars.var_type,'vector')
    f = (0.5*pars.A*x - pars.B)'*x;
    g = pars.A*x - pars.B;
elseif strcmp(pars.var_type,'matrix')
    f = 0.5*trace(x'*pars.A*x) - trace(x'*pars.B);
    g = pars.A*x - pars.B;
end

end