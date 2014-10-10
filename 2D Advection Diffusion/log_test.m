% Test code for stuff, scratch work
% Corbin Foucart

clear all; 
close all; clc;

% Finding test code
%{
points = [1,2,3;
          3,5,7;
          5,8,3;
          0,1,2;
          10,12,13]
      
A = (points == 3)      
B = (points == 5)
C = A + B
trues = sum(C,2)

ind = [find(trues >= 2)]'

%}
% Barycenter test code
%{
tri_points = [0.5,0;
              0.5, 0.5;
              0, 0.5];      
r = [0.625, 0.25];
[bary, i] = adif2d_barycentric_coords(tri_points, r)

figure()
plot(tri_points(:,1), tri_points(:,2), 'ko')
hold all
plot(r(1), r(2), 'rx')

pts = [0,0; 0.5,0; 0,0.5]
test = [0.375, 0.25]

[lamb tru] = adif2d_barycentric_coords(pts,test)
%}

% Inpoly test code
%{
n = 6;            
            
bound_points = ...
    [0.5000,         0;
     1.0000,         0;
     1.0000,    0.5000;
     1.0000,    1.0000;
     0.5000,    1.0000;
     0.5000,    0.5000]   
 
 inner_pt = ...
    [0.7500    0.2500;
     1.0000    0.2500]

list = [[1:n]',[[2:n]';1]]
[in, on] = inpolygon(inner_pt, bound_points, list)

truth = in+on

isempty(find(truth == 0))
%}
% Return triangles test code
%{
vfilename =  'voronoi_data_test.mat';
load(vfilename, 'v_lines', 'neighbor_mat', 'edge_matrix', 'points', 'conn_list');

outer_points = [5 6 7 8 9]'
      
sum_mat = zeros(size(conn_list))
for i = 1:length(outer_points)  
    sum_mat = sum_mat + (conn_list == outer_points(i))
end

% return indices in the matrix that sum to 3, because we want all the
% triangles that contain exactly 3 points in our point list.
ind = find((sum(sum_mat,2) == 3))
%}          
% Getting points of a delaunay trinagle test code
%{
vfilename =  'voronoi_data_test.mat';
load(vfilename, 'v_lines', 'neighbor_mat', 'edge_matrix', 'points', 'conn_list');

% triangle nums
tris = [1 2 6]
% extract point numbers conn_list
pt_nums = conn_list(tris, :)
% get actual points for each row
acpoints = points(pt_nums(1,:),:)
%}
% Return edges test code
%{
vfilename =  'voronoi_data_test.mat';
load(vfilename, 'v_lines', 'neighbor_mat', 'edge_matrix', 'points', 'conn_list');

region_points = [1 2 4 5]
edges = adif2d_edge_returner(region_points, edge_matrix)
%}


%Initial norm data save;
%{
l2_data = []
linf_data= []
h1_data = []
filename = 'norm_data_sqrgrid'
save(filename, 'h1_data', 'linf_data', 'l2_data')
%}


% Norm Plotting Scratch Work
%
load('norm_data_sqrgrid', 'l2_data','linf_data', 'h1_data')
h = 1./sqrt(l2_data(:,1));
%figure()
%loglog(h, l2_data(:,2))
%title('mesh size h vs l2, log')
%xlabel('mesh size h')
%ylabel('L^2 norm')
pl2 = polyfit(log(h), log(l2_data(:,2)),1)

%figure()
%loglog(h,linf_data(:,2))
%title('mesh size h against L^\infty')
%xlabel('mesh size h')
%ylabel('L^\infty norm')
plinf = polyfit(log(h), log(linf_data(:,2)), 1)

%figure()
%loglog(h,abs(h1_data(:,2)))
%title('mesh size h against h1 norm')
%xlabel('mesh size h')
%ylabel('h1 norm')
ph1 = polyfit(log(h), log(h1_data(:,2)), 1)
%




%{
% order points in distance from the first one
seg_points = [-2,-1; 2,1; 1,0.5; 8,4; -2,-1; 6,3; 4,2]
seg_points(:,3) = sqrt((seg_points(:,1) - seg_points(1,1)).^2 +...
    (seg_points(:,2) - seg_points(1,2)).^2)
seg_points = sortrows(seg_points,3)
%}
%{
% Point selection
vfilename =  'voronoi_data_test.mat';
load(vfilename, 'v_lines', 'neighbor_mat', 'edge_matrix', 'points',...
    'conn_list','CC','dtri', 'adj_triangles');

A = [1 2]
B = conn_list(A,:)

points(B(1,:),:)
points(B(2,:),:)

%}

% Experimenting with cells
%{
C{1} = {magic(3), 17}
C{2} = {magic(3), 17}

C{1}{1}
%}
%{
% Read in Data, initialize column vectors and matrices. 
vfilename =  'voronoi_data_10pt.mat';
load(vfilename);

% Get mesh crossing data
cross_dat = adif2d_edge_segmentation_fun(v_lines, neighbor_mat,...
            edge_matrix, points,conn_list, adj_triangles);
    
    %}

% Eigenvalue plotting
% Clearing Code
%{
D = [];
filename = 'eigSave1'
save(filename, 'D')
%}