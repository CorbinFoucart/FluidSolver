% Poisson Source Function 2D
% Corbin Foucart

function z = poisson_source_2d(x,y,D)


% constant source function
%  z = 1;
 
% z = 100*y^3;

%z = -(2*(y-1).*y + 2*(x-1).*x);

% < 1, 0 >
%z = -D*(2.0*(y-1.0).*y + 2.0*(x-1.0).*x) + (1.0-2.0*x).*y.*(1.0-y); 

% < -x, y >
%z = x .*(x - y) .* y - 2.0 * D * (-x + x.^2.0 + (-1.0 + y).*y);

% Crouzeix <u,v>
z = 2.0.*D.*(x - x.^2.0 + y - y.^2.0);

end