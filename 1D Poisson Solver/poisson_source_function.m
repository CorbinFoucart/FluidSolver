% Poisson Source Function
% Corbin Foucart
% ------------------------------------------------------------------------
% Takes in a sorted vector of source points and returns the corresponding
% results of the 1D source function. 

function s_points = poisson_source_function(x_points)


% Constant source
% constant_value = 10;
% s_points = zeros(length(x_points), 1) + constant_value;

% Linear source
% s_points = 1000 * x_points;

% Quadratic source
 s_points = 50 * x_points.^2;

% Other Sample Function
% s_points = 50 * abs(sin(2*pi*x_points)) .* x_points.^2;

end