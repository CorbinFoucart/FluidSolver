% Barycentric Coordinate Calculator
% Corbin Foucart
% -------------------------------------------------------------------------
% Takes in a 3x3 matrix, and a 2D point r = [x,y] in a row vector
%                                           
%
% R = [ x1, y1]  
%     [ x2, y2]
%     [ x3, y3]
%
%  Normalizing the w components ensure that we have affine coordinates
%  Barycentric coordinates returned as a column vector 

function [rn_lamb rn_in] = adif2d_barycent(R,r)

% Use matlab's convex hull to orient points counterclockwise
vi = convhull(R);
R = [R(vi,1), R(vi,2)]; % Re order R according to convex hull
%R = unique(R, 'rows'); % Remove doubled entry

% Columns
R = R';
r = r';

A = [R(1,1) - R(1,3),  R(1,2) - R(1,3);
     R(2,1) - R(2,3),  R(2,2) - R(2,3)];
b =  [r(1) - R(1,3); r(2)- R(2,3)];

lamb = A\b;
rn_lamb = [lamb(1); lamb(2); 1 - lamb(1) - lamb(2)];

% Tolerance 
for i = 1:3
    if rn_lamb(i) < 0
        if abs(rn_lamb(i)) < 0.000001
            rn_lamb(i) = 0;
        end
    end
end

% check if point is within triangle or on boundary
% no lambdas less than 1 -> point is inside or on boundary
if sum(find(rn_lamb < 0)) == 0
    rn_in = 1;
else
    rn_in = 0;
end


end