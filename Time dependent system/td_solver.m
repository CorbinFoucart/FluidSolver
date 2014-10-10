% time dependent solver
% corbin foucart


% allows for iteration of previous CSC_solver files in time.
 clear all;
 close all; clc;

% ----------------------- user input / options -------------------------- %
% choice of mesh to read in from data file
%vn =  'zVrm10.mat';
vn =  'zVg30.mat';
load(vn);

dt = 0.001; % time step [s]

isRT0mode = 1;
Ra = 1;
N = 10;

% --------------------------- Computation ------------------------------- %
first_it = 1;



% solution record
Sol_keep = zeros(length(points(:,1)), N+1);
solu_keep = zeros(length(adj_triangles(:,1)), N+1);
solv_keep = zeros(length(adj_triangles(:,1)), N+1);
for T = 1:N
    
    if first_it
    % Generate temperature guess, or initial condition
    [Tsol, A, b] = CD_solver_guess_fn(vn);
    % save data for plotting
    Sol_keep(:,1) = Tsol;
    
    first_it = 0;
    time_iteration = 0
    end

    if ~first_it
        % count iteration
        time_iteration = time_iteration + 1
        % iterate within timestep for accuracy
        [Tsol, psol, usol, vsol, Nu, TriA, xmps, ymps, xbcs, ybcs] ...
            = CSC_solver_fn(vn, isRT0mode, Ra, Tsol, dt);
        % save solution data for plotting
        Sol_keep(:, T+1) = Tsol;   
        solu_keep(:, T+1) = usol;
        solv_keep(:, T+1) = vsol;
    end



end



% -------------------- visualization, plots ----------------------------- %

% mesh points
x = points(:,1);
y = points(:,2);

% triangulations for plotting
triplot = delaunay(x,y);
%ptri = delaunay(xbcs, ybcs);


% Preallocate movie structure.
mov(1:N) = struct('cdata', [],...
                        'colormap', []);

% Create movie.
Z = Sol_keep(:,1); trisurf(triplot,x,y,Z); 
axis tight
set(gca,'nextplot','replacechildren');
for k = 1:N + 1 
   trisurf(triplot,x, y, Sol_keep(:,k))
   view([-80 36])
   mov(k) = getframe(gcf);
end

% Create AVI file.
movie2avi(mov, 'myTemp.avi', 'compression', 'None', 'fps', 10);


%{
% Temperature Profile
figure()
trisurf(triplot,x,y,Sol_keep(:,i))
title('Temperature Profile')
colormap jet
axis square
alpha(.7)
xlabel('x')
ylabel('y')
zlabel('T')
%}

% scale velocities
maxU = max(solu_keep(:));
maxV = max(solv_keep(:));
maxV = max(maxU, maxV);
largest_v = num2str(maxV);

solu_keep = solu_keep/(10*maxV);
solv_keep = solv_keep/(10*maxV);

% Preallocate movie structure.
movV(1:N) = struct('cdata', [],...
                        'colormap', []);

% Create movie.
quiver(xmps, ymps, solu_keep(:,1), solv_keep(:,1),'Autoscale', 'off');
axis tight
set(gca,'nextplot','replacechildren');
for k = 1:N + 1 
   quiver(xmps, ymps, solu_keep(:,k), solv_keep(:,k),'Autoscale', 'off');
   title(strcat('Velocity Profile MaxV: ',largest_v))
    axis square
    axis([0 1 0 1])
    alpha(.6)
    xlabel('x')
    ylabel('y')
   movV(k) = getframe(gcf);
end

% Create AVI file.
movie2avi(movV, 'myVel.avi', 'compression', 'None', 'fps', 10);


%{
for i = 1:N+1
    
% Velocity Field
figure()
quiver(xmps, ymps, solu_keep(:,i), solv_keep(:,i),'Autoscale', 'off');
title(strcat('Velocity Profile MaxV: ',largest_v))
axis square
axis([0 1 0 1])
alpha(.6)
xlabel('x')
ylabel('y')

end
%}


 %{
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
%}

