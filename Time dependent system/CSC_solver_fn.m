% Coupled Stokes, Convection Diffusion Solver, time independent
% Corbin Foucart
% -------------------------------------------------------------------------
%
function [rTsol, rpsol, rusol, rvsol, rNu, rTriA, rXmps, rYmps, rXbcs, rYbcs] ...
    = CSC_solver_fn(vn, isRT0mode, Ra, T0, dt)

% load mesh data
load(vn);

TguessMode = 0;

    % mesh points
    x = points(:,1);
    y = points(:,2);

    % triangulations for plotting
    triplot = delaunay(x,y);
    

% --------------------------- first iteration -------------------------- %

% Read in 'inital condition' from last time step

%{

  figure()
    trisurf(triplot,x,y,Tsol)
    title('Temperature Profile')
    colormap jet
    axis square
    alpha(.7)
    xlabel('x')
    ylabel('y')
    zlabel('T')
%}

% generate corresponding pressure and velocity solutions to that inital
% condition, as well as a first iteration of temperature
[psol, Rtx, Rty, S, Ssol, Srhs, usol, vsol, xbcs, ybcs, xmps, ymps, s2e, triA] ...
    = crouzeix_fn(Ra, isRT0mode, T0, vn);

[Tsol, A, b] = CD_solver_fn(Rtx, Rty, vn, T0, dt);


% intialization of matrix cells (for tracking residuals)
% note: we don't need to save the vector b in the convection diffusion
% equation Ax = b, because it is the zero vector, so the residual is simply
% res = Ax. For the stokes equations, we must save the RHS.
%
% note also that the second main slot of the cell, res_mat{2}, always
% contains the most recent information from the solvers in the iteration.
res_mats{2,1} = {A, Tsol, b};
res_mats{2,2} = {S, Ssol, Srhs};

% ------------------------- coupled system ------------------------------ %
% no longer use first iteration guesses
iterationNumber = 1;

% intialize residual data save
res_data = [];
 
%{
    figure()
    trisurf(triplot,x,y,Tsol)
    title('Temperature Profile')
    colormap jet
    axis square
    alpha(.7)
    xlabel('x')
    ylabel('y')
    zlabel('T')
%}


while true

    % display iteration number
    iterationNumber
    

    % first step: move the most recent solve data into the first slot of
    % res_mat{1} to allow solvers to populate res_mat{2} with new solution data
    res_mats{1,1} = res_mats{2,1};
    res_mats{1,2} = res_mats{2,2};

    % step forward in solution iteration
    [psol, Rtx, Rty, S, Ssol, Srhs, usol, vsol, xbcs, ybcs, xmps, ymps, s2e, triA]...
        = crouzeix_fn(Ra, isRT0mode, Tsol, vn);
    [Tsol, A, b] = CD_solver_fn(Rtx, Rty, vn, T0, dt);
%{    
    figure()
    trisurf(triplot,x,y,Tsol)
    title('Temperature Profile')
    colormap jet
    axis square
    alpha(.7)
    xlabel('x')
    ylabel('y')
    zlabel('T')
    %}
    
    % save new solution data
    res_mats{2,1} = {A, Tsol, b};
    res_mats{2,2} = {S, Ssol, Srhs};

    % computation of residual data
    % number of entries
    N_stokes = length(Srhs);
    N_CD = length(Tsol);

    % residuals (basically computing Ax - b)
    res_stokes = 1/N_stokes * sum(abs(res_mats{1,2}{1,1}*res_mats{1,2}{1,2} ...
        - res_mats{2,2}{1,3}))

    res_CD = 1/N_CD * sum(abs(res_mats{2,1}{1,1}*res_mats{1,1}{1,2}...
        - res_mats{2,1}{1,3}))

    v_mag = mean(sqrt(usol.^2 + vsol.^2))

    if iterationNumber > 1 && res_stokes < 1e-8 && res_CD < 1e-8
        break
    end
    
    
    iterationNumber = iterationNumber + 1;

end


% --------------------- Nusselt number Calc ----------------------------- %
gTmp = zeros(length(edge_matrix(:,1)) , 2);
% for each triangle
Nu = 0;
for i = 1:length(CC(:,1))
    
    % get vertex numbers of triangle
    vs = conn_list(i, :);
    
    % get edge numbers of triangle
    es = s2e(i,:);    
    
    % get vertex info
    for j = 1:3
        % r1 = [x1 y1 T1], r2 = [x2 y2 T2]... etc
        r(j,:) = [points(vs(j), 1), points(vs(j), 2), Tsol(vs(j))];
    end
    
    av_T = sum(r(:,3))/3;
    
    % compute normal vector
    nrm = cross(r(2,:) - r(1,:), r(3,:) - r(1,:));
    
    % compute gradient from normal vector 
    grad = -1*[nrm(1)/nrm(3), nrm(2)/nrm(3)];
    
    for j = 1:3
        gTmp(es(j), :) = grad;
    end

    % phi = 1 - y
    gPhi = [0, -1];
    
    Nu = Nu + dot(gPhi,(-grad + av_T*[Rtx(i), Rty(i)]))*triA(i);
end

Nu = abs(Nu);


% ------------------------- returning data ----------------------------- %
rTsol = Tsol;
rpsol = psol;
rusol = usol;
rvsol = vsol;

rNu = Nu;
rTriA = triA;
rXmps = xmps;
rYmps = ymps;
rXbcs = xbcs;
rYbcs = ybcs;
end












