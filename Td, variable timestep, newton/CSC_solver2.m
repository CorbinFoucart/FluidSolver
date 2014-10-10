% Coupled Stokes, Convection Diffusion Solver, time independent
% Corbin Foucart
% -------------------------------------------------------------------------
%{
clear all; 
close all;
clc;



% ----------------------- user input / options -------------------------- %
% choice of mesh to read in from data file
%vn =  'zVrm10.mat';
vn =  'zVg20.mat';
mnn = vn(4:5);
load(vn);

cd 'RT1 - data files'
load('Ra50000mesh20temp.mat')
cd ..

% options
isRT0mode = 0;
Ra = 5e5;
RaStr = int2str(Ra);

% picture save filenames
tempSave = strcat('Ra', RaStr, 'mesh', mnn,'temp','.jpg');
velSave = strcat('Ra', RaStr, 'mesh', mnn,'vel','.jpg');


% -------------- first iteration (guess of v = 0 everywhere)------------- %
first_it = 1;
Rtx = 1;
Rty = 1;

% Generate Temperature 'guess' *Must change BCs here as well*
TguessMode = 0;


[psol, Rtx, Rty, S, Ssol, Srhs, usol, vsol, xbcs, ybcs, xmps, ymps] ...
    = crouzeix_fn(Ra, isRT0mode, Tsol, vn);
[Tsol, A, b] = CD_solver_fn(Rtx, Rty, first_it, vn, TguessMode);

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
    [psol, Rtx, Rty, S, Ssol, Srhs, usol, vsol, xbcs, ybcs, xmps, ymps]...
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


% -------------------- visualization, plots ----------------------------- %
%}
% mesh points
x = points(:,1);
y = points(:,2);

% triangulations for plotting
triplot = delaunay(x,y);
ptri = delaunay(xbcs, ybcs);

% Gridfitting for countour plots
contNodesX = 0:0.1:1;
contNodesY = 0:0.1:1;

% magnitude of vectors
%contZ = sqrtm(usol*usol' + vsol*vsol');
%{
contZ = usol + vsol;

[contx, conty, contz] = gridfit(xmps, ymps, contZ,...
    contNodesX, contNodesY);
%}
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

%{
% Contour Plot of Velocity Field
figure()
contour(contx, conty, contz);
title('Velocity Contour Profile')
axis square
xlabel('x')
ylabel('y')
%}












