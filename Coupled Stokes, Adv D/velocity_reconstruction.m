% Velocity Reconstruction function
% Corbin Foucart

% Reconstruction of velocity
% note: above, the boolean isRT0mode determined how the right hand side of
% the matrix equation Ax = b was computed. This implies two different right
% hand sides, which implies a different [usol, vsol]. Regardless, we
% perform a reconstruction of the velocity solution to pass to the
% Convection diffusion solver. This reconstruction operator is consistent
% for both velocity solutions. 
%
% The rationale behind this is that when the stokes solver is coupled to
% the convection diffusion solver, the reconstruction changes the accuracy
% of the scheme by a negligible amount, but hugely improves stability of
% the finite volumes scheme. 

% Input: RT is a cell containing all the geometric information necessary to
% compute the reconstruction. See Crouzeix_fn.m for specifics, but basici
% structure is RT{i,1}{1,k}
%
% i: triangle number
% k: 1:3 for edges of the triangle.
% RT{i,1}{1,k} contains an array:
%   [aTx aTy b unormx unormy]
%
% usol, vsol are the velocity solutions from the stokes solver (i.e. the
% x vector in the Ax = b problem).
%
% dim is used only as a counting parameter used to know the number of
% triangles in our mesh.

function [RVx, RVy] = velocity_reconstruction(RT, usol, vsol, dim) 

RV_Ax = zeros(dim,1);
RV_Ay = zeros(dim,1);
RT_b  = zeros(dim,1);
for i = 1:dim
    for j = 1:3
        localv = [usol(RT{i,1}{1,j}(6)), vsol(RT{i,1}{1,j}(6))];
        % second row of cell relevant velocities
        RT{i,1}{2,j} = localv;
        % 3rd row dot product between normals and velocities
        RT{i,1}{3,j} = dot(RT{i,1}{1,j}(4:5), localv);
        % 4th row dot products multiplied with the respective factors
        % scalar * array
        RT{i,1}{4,j} = RT{i,1}{3,j}*RT{i,1}{1,j}(1:3);
    end
    % next row in cell the values [aTx aTy b] for each triangle
    RT{i,2} = RT{i,1}{4,1} + RT{i,1}{4,2} + RT{i,1}{4,3}; 

    % Save values to more easily accessible vectors
    RV_Ax(i,1) =  RT{i,2}(1);
    RV_Ay(i,1) =  RT{i,2}(2);
    RT_b(i,1)  =  RT{i,2}(3);  % necessary?
end

% return
RVx = RV_Ax;
RVy = RV_Ay;

end