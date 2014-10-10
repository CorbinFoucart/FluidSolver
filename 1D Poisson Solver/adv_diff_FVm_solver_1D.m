% 1D advection, diffusion equation
% Corbin Foucart
% -------------------------------------------------------------

clear all; close all; clc; 

% Input
% Interval
a = 0;
b = 10;
N = 40; % Number of mesh nodes, including endpoints

% Constants
D = 2;
v = 3;

% Boundary Conditions
c0 = 0;
c1 = 5;

% ------------------------ Mesh Generation ----------------------%
[x, v_boundaries] = voronoi_gen_1D_function(a, b, N);

% ----------------- Source Function, Integration ----------------%
% For each control volume, compute the source function at N_quad more inner
% points between Voronoi boundaries and numerically integrate the source
% function using a trapezoidal quadrature scheme

int_boundaries = [a, v_boundaries(1:end)', b]';

% Quadrature resolution
N_quad = 10;
approx_integral = zeros(N, 1);
for i = 2:length(int_boundaries) - 1
    source_interval = linspace(int_boundaries(i),...
        int_boundaries(i+1), N_quad);
    s_x = poisson_source_function(source_interval);
    approx_integral(i) = trapz(source_interval, s_x);
end

% ----------------- Calculation of Solution ----------------%
% Solve the system Ax = b

A = zeros(N);
% Matrix without boundary terms
for j = 2 : N - 1
    for k = 2 : N - 1
        if j == k
            A(j, k-1) = -D/(x(j)-x(j-1)) - v/2;
            A(j, k) = D/(x(j+1)-x(j)) + D/(x(j)-x(j-1));
            A(j, k+1) = -D/(x(j+1)-x(j)) + v/2;
        end
    end
end

% Boundary Terms
% Boundary Terms
A(1, 1) = c0;
A(N, N) = c1;

% Solution
c = A\approx_integral;

% ----------------------- Plotting ------------------------ %
Y_points = zeros(length(x), 1);

% Visual source function
vis_source_points = linspace(a, b, 100)';
vis_source_fn = poisson_source_function(vis_source_points);


plot(x, Y_points, 'ko')
hold on
for i = 1:length(v_boundaries)
    plot(line([v_boundaries(i), v_boundaries(i)], ...
        [-0.05*max(vis_source_fn), 0.05*max(vis_source_fn)]))
end
plot(vis_source_points, vis_source_fn, 'b-')
axis([a, b, -1, 1.1*max([max(vis_source_fn), max(c)])])
plot(x, c, 'k-')
