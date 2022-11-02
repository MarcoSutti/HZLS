function [ f, g ] = rosenbrock_function( x, pars )

% function [ f, g ] = rosenbrock_function( x, pars )
% Created:     19.08.2020
% Last change: 28.08.2020

% if pars.dim == 2
% Rosenbrock function, dimension dim = 2
f = 100*(x(2)-x(1)^2)^2 + (1-x(1))^2 + pars.z0;
g =  [ 2*(200*x(1)^3-200*x(1)*x(2)+x(1)-1);
    200*(x(2)-x(1)^2) ];
% elseif pars.dim == 3
%     % Rosenbrock function, dimension dim = 3
%     f = @(x) 100*(x(2)-x(1)^2)^2+(1-x(1))^2 + 100*(x(3)-x(2)^2)^2+(1-x(2))^2;
%     g = @(x) [ 2*(200*x(1)^3-200*x(1)*x(2)+x(1)-1);
%         -200*x(1)^2 + 400*x(2)^3 + x(2)*(202 - 400*x(3))-2;
%         200*(x(3)-x(2)^2) ];
% end

end