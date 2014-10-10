% Coupled Stokes, Convection Diffusion Solver, time independent
% Corbin Foucart
% -------------------------------------------------------------------------
%
clear all; 
close all;
clc;
format long

% ----------------------- user input / options -------------------------- %
% choice of mesh to read in from data file
%vn =  'zVrm10.mat';
vn =  'zVg50.mat';
mnn = vn(4:5);
load(vn);

% options
isRT0mode = 1;
Ra = 1e4;
RaStr = int2str(Ra);

% picture save filenames
tempSave = strcat('Ra', RaStr, 'mesh', mnn,'temp','.jpg');
velSave = strcat('Ra', RaStr, 'mesh', mnn,'vel','.jpg');


% -------------- first iteration (guess of v = 0 everywhere)------------- %
first_it = 1;
Rtx = 1;
Rty = 1;

% Generate Temperature 'guess' *Must change BCs here as well*
TguessMode = 1;


[Tsol, A, b] = CD_solver_fn(Rtx, Rty, first_it, vn, TguessMode);
[psol, Rtx, Rty, S, Ssol, Srhs, usol, vsol, xbcs, ybcs, xmps, ymps, s2e, triA] ...
    = crouzeix_fn(Ra, isRT0mode, Tsol, vn);

% turn guessing off, for future iteration
TguessMode = 0;

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
first_it = 0;
iterationNumber = 1;

% intialize residual data save
res_data = [];

% see steps
%{
% mesh points
x = points(:,1);
y = points(:,2);

% triangulations for plotting
triplot = delaunay(x,y);LE6_vol_26_iss_2_024105_1.pdf
ptri = delaunay(xbcs, ybcs);
%}

while true

    % display iteration number
    iterationNumber
   
    %close all;

    % see steps
    %{
    figure('Renderer','zbuffer');
    trisurf(triplot,x,y,Tsol)
    title('Temperature Profile')
    colormap jet
    axis square
    set(gca,'NextPlot','replaceChildren');
    F(iterationNumber) = getframe;
    %}
    
    

    % first step: move the most recent solve data into the first slot of
    % res_mat{1} to allow solvers to populate res_mat{2} with new solution data
    res_mats{1,1} = res_mats{2,1};
    res_mats{1,2} = res_mats{2,2};

    % step forward in solution iteration
    [Tsol, A, b] = CD_solver_fn(Rtx, Rty, first_it, vn, TguessMode);
    [psol, Rtx, Rty, S, Ssol, Srhs, usol, vsol, xbcs, ybcs, xmps, ymps, s2e, triA]...
        = crouzeix_fn(Ra, isRT0mode, Tsol, vn);

    % save new solution data
    res_mats{2,1} = {A, Tsol, b};
    res_mats{2,2} = {S, Ssol, Srhs};

    % computation of residual data
    % number of entries
    N_stokes = length(Srhs);
    N_CD = length(Tsol);

    % residuals
    res_stokes = 1/N_stokes * sum(abs(res_mats{1,2}{1,1}*res_mats{1,2}{1,2} ...
        - res_mats{2,2}{1,3}))

    res_CD = 1/N_CD * sum(abs(res_mats{2,1}{1,1}*res_mats{1,1}{1,2}...
        - res_mats{2,1}{1,3}))

    v_mag = mean(sqrt(usol.^2 + vsol.^2))

    if iterationNumber > 1 && res_stokes < 1e-9 && res_CD < 1e-9
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

Nu = abs(Nu)

% -------------------- visualization, plots ----------------------------- %

% mesh points
x = points(:,1);
y = points(:,2);

% triangulations for plotting
triplot = delaunay(x,y);
ptri = delaunay(xbcs, ybcs);

% Temperature Profile
figure()
trisurf(triplot,x,y,Tsol)
title('Temperature Profile')
colormap jet
axis square
alpha(.7)
xlabel('x')
ylabel('y')
zlabel('T')

h = gcf;
saveas(h,tempSave);

% Velocity Field
figure()
quiver(xmps, ymps, usol, vsol);
title('Velocity Profile')
axis square
alpha(.6)
xlabel('x')
ylabel('y')

h = gcf
saveas(h,velSave)

% Pressure Profile
figure()
trisurf(ptri, xbcs, ybcs, psol)
title('Pressure Profile')
colormap winter
axis square
alpha(.6)
xlabel('x')
ylabel('y')
zlabel('P')












