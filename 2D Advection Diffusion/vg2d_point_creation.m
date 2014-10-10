% Voronoi 2D mesh generator, random creation of points
% Corbin Foucart
%-------------------------------------------------------------------------
% Input : a,b,c,d -> [a,b]x[c,d] domain boundaries
%         N_bp, N_ip, N_cp -> num. of boundary points, inner points and
%         corner points. 
%         0 < sf < 0.5 keeps points certain ratio of the boundary away from
%         the edge
% 

function return_pts = vg2d_point_creation(a,b,c,d,N_bp,...
    N_ip, N_cp, sf)

% ------------------------ Creation of Points ------------------------ %

% Test Case 1
%{
return_pts = [0,0;
          0.5,0;
          1,0;
          0,0.5;
          0.5,0.5;
          1,0.5;
          0,1;
          0.5,1;
          1,1];
%}
% Test Case 2
%{
return_pts = [0,0; 0.25,0; 0.5,0; 0.75,0; 1,0;
              0,0.25; 0.25,0.25; 0.5,0.25; 0.75,0.25; 1,0.25;
              0,0.5; 0.25,0.5; 0.5,0.5; 0.75,0.5; 1,0.5;
              0,0.75; 0.25,0.75; 0.5,0.75; 0.75,0.75; 1,0.75;
              0,1; 0.25,1; 0.5,1; 0.75,1; 1,1]
%}

% Randomly generated Case
total_points = N_bp + N_ip + N_cp;

% initialize point matrix
points = zeros(total_points, 2);
% boundary points
bound_points = zeros(N_bp + 4, 2);
bound_points(1:4,:) = [[a,c];[b,c];[b,d];[a,d]];

% Place points randomly on edges
for i = 5:length(bound_points)
    side = randi([1,4]);
    placement = rand(1,1);
    switch side
        case 1
            bound_points(i,:) = [a, c + (d-c)*placement];
        case 2
            bound_points(i,:) = [a + (b-a)*placement, d];
        case 3
            bound_points(i,:) = [b, c + (d-c)*placement];
        case 4
            bound_points(i,:) = [a + (b-a)*placement, c];
    end
end

% Sort the points counter clockwise
%{
sort_bound_points = [[a,c];
                     sortrows(bound_points(bound_points(:,2) == c, :),1); ...
                     [b,c];
                     sortrows(bound_points(bound_points(:,1) == b, :),2); ...
                     [b,d];
                    flipud(sortrows(bound_points(bound_points(:,2) == d, :),1)); ...
                     [a,d];
                    flipud(sortrows(bound_points(bound_points(:,1) == a, :),2))]
                
bound_points = sort_bound_points;
%}

% random inner points
 
inner_x_points = a + sf*(b - a) + (1-2*sf)*(b - a).*rand(N_ip, 1);
inner_y_points = c + sf*(d - c) + (1-2*sf)*(d - c).*rand(N_ip,1);

% check for repeating points, only need to check 1 dimension
while unique(inner_x_points) ~= inner_x_points 
    inner_x_points = a + sf*(b - a) + (1-2*sf)*(b - a).*rand(N_ip, 1);
end

% Assignment
% Assign boundary points last so when steiner points added, they will be
% grouped with other boundary points. 
points(1:N_ip,1:2) = [inner_x_points, inner_y_points];
points(N_ip + 1:end, 1:2) = bound_points;

% redundancy check
points = unique(points, 'rows');

% ----------------------- Delaunay Triangulation ------------------------%

tri = delaunay(points);
dtri = delaunayTriangulation(points);
dtri_points = dtri.Points;              % row is vertex num
conn_list = dtri.ConnectivityList;      % row is triangle num
CC = circumcenter(dtri);                % row is triangle num
[mcon ncon] = size(conn_list);

% -------------------- Enforce boundary conformation ------------------- %
    while(true)
        % prob_cc: rows are triangle_num, columns are coordinates
        prob_cc = zeros(length(CC(:,1)),2);
        repair_pts = zeros(length(CC(:,1)),2); % must have 0s removed
        for i = 1:length(CC(:,1))
            cct = CC(i,:);
            % check if boundary triangles need steiner points, add them. 
            if cct(1) < a 
                prob_cc(i,:) = CC(i,:);
                prob_pts =  sortrows([conn_list(i,:)', dtri_points(conn_list(i,:),:)], 2);
                if prob_pts(1,2) == a && prob_pts(2,2) == a % boundary point check
                   x1 = prob_pts(1,2:3)'; x2 = prob_pts(2,2:3)'; x3 = prob_pts(3,2:3)';
                   theta_out = acos(dot((x2-x3),(x1-x3))/((norm(x2-x3))*norm(x1-x3)));
                   if theta_out > pi/2 % triangle is obtuse
                       repair_pts(i,:) = [a, 1/2*(x1(2)+x2(2))];
                   end
                end
            end
            if cct(1) > b
                prob_cc(i,:) = CC(i,:);
                prob_pts =  flipud(sortrows([conn_list(i,:)', dtri_points(conn_list(i,:),:)], 2));
                if prob_pts(1,2) == b && prob_pts(2,2) == b % boundary point check
                   x1 = prob_pts(1,2:3)'; x2 = prob_pts(2,2:3)'; x3 = prob_pts(3,2:3)';
                   theta_out = acos(dot((x2-x3),(x1-x3))./((norm(x2-x3))*norm(x1-x3)));
                   if theta_out > pi/2 
                       repair_pts(i,:) = [b, 1/2*(x1(2)+x2(2))];
                   end
                end
            end
            if cct(2) < c
                prob_cc(i,:) = CC(i,:);
                prob_pts =  sortrows([conn_list(i,:)', dtri_points(conn_list(i,:),:)], 3);
                if prob_pts(1,3) == c && prob_pts(2,3) == c % boundary point check
                   x1 = prob_pts(1,2:3)'; x2 = prob_pts(2,2:3)'; x3 = prob_pts(3,2:3)';
                   theta_out = acos(dot((x2-x3),(x1-x3))./((norm(x2-x3))*norm(x1-x3)));
                   if theta_out > pi/2
                       repair_pts(i,:) = [1/2*(x1(1)+x2(1)), c];
                   end
                end
            end
            if cct(2) > d
                prob_cc(i,:) = CC(i,:);
                prob_pts =  flipud(sortrows([conn_list(i,:)', dtri_points(conn_list(i,:),:)], 3));
                if prob_pts(1,3) == d && prob_pts(2,3) == d % boundary point check
                   x1 = prob_pts(1,2:3)'; x2 = prob_pts(2,2:3)'; x3 = prob_pts(3,2:3)';
                   theta_out = acos(dot((x2-x3),(x1-x3))./((norm(x2-x3))*norm(x1-x3)));
                   if theta_out > pi/2 
                       repair_pts(i,:) = [1/2*(x1(1)+x2(1)), d];
                   end
                end
            end
        end

        % check to see if no more points have been repaired - implies that no
        % more steiner points must be added. 
        if repair_pts == zeros(length(CC(:,1)),2)
            break
        end

        % add repair points, remove redundancies
        points = unique([points; repair_pts],'rows');

        % re triangulate 
        tri = delaunay(points);
        dtri = delaunayTriangulation(points);
        dtri_points = dtri.Points;              % row is vertex num
        conn_list = dtri.ConnectivityList;      % row is triangle num
        CC = circumcenter(dtri);                % row is triangle num
        [mcon ncon] = size(conn_list);
    end

return_pts = points;

end
