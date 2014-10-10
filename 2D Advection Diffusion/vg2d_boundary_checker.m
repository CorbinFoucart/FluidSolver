% Voronoi 2D mesh generator, boundary checker
% Corbin Foucart
%-------------------------------------------------------------------------
% Takes in the list of triangulation points and checks which ones are 
% on the boundary. Makes sense to do this now rather than initially because
% of the necessary addition of steiner points.

function rn_dtri_pts = vg2d_boundary_checker(a,b,c,d,dtri_points)

% Mark which points are boundary points, append to dtri_points
boundary_tracker = zeros(length(dtri_points(:,1)), 1);
for i = 1:length(dtri_points(:,1))
    if dtri_points(i,1) == a || dtri_points(i,1) == b...
         || dtri_points(i,2) == c || dtri_points(i,2) == d
        boundary_tracker(i) = 1;
    end
end

% 3rd column whether boundary point, 0,1
rn_dtri_pts = [dtri_points, boundary_tracker]; 
end