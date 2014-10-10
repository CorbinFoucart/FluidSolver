% Voronoi Generator 1D Function
% Corbin Foucart
% ------------------------------------------------------------------------
% Function takes in the interval data [a,b] as well as number of randomly
% placed (uniform distribution) mesh points N, including boundary points.
%
% Function returns the randomly generated mesh points, which include the
% endpoints of the interval as well as the Vornoi boundary points that
% separate the mesh points. 

function [points, v_boundaries] = voronoi_gen_1D_function(a, b, N)

% ----------- creation of points ------------ %
N = N - 2; % correction, since boundary points should be included in N
% linear method
% points = linspace(a,b,N);

% random method, sorted points
inner_points = sort(a + (b - a).*rand(N, 1));
points = [a, inner_points(1:end)', b]';
% check for repeating points, recompute if a repeat
while unique(points) ~= points 
    points = a + (b-a).*rand(N,1); 
end

%------- calculation of boundaries -------- %
v_boundaries = zeros(length(points) - 1, 1);
    for i = 1:length(points) - 1
        v_boundaries(i) = 1/2 * (points(i) + points(i + 1));
    end
end





