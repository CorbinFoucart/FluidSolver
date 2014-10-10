% Voronoi Generator 1D
% Corbin Foucart
% ----------------------
% Script takes in a interval [a,b], and either randomly or regularly 
% momputes a specified number of N points within the interval. 
% The script then calculates and plots the Voronoi ranges.

clear all; close all; clc;

% input
a = 0;
b = 100;
N = 8;

% ----------- creation of points ------------ %
interval = [a,b];

% linear method
% points = linspace(a,b,N);

% random method, sorted points
inner_points = sort(a + (b - a).*rand(N, 1));
points = [a, inner_points(1:end)', b]';
% check for repeating points
while unique(points) ~= points 
    points = a + (b-a).*rand(N,1);
end

%------- creation of Vornoi intervals -------- %
boundaries = zeros(length(points) - 1, 1);
for i = 1:length(points) - 1
    boundaries(i) = 1/2 * (points(i) + points(i + 1));
end


% ----------------- plotting ----------------- %
Y_points = zeros(length(points), 1);

plot(points, Y_points, 'ko')
hold on
for i = 1:length(boundaries)
    plot(line([boundaries(i), boundaries(i)], [-1, 1]))
end
axis([a, b, -5, 5]);





