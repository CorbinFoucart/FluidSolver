clear all;
close all;
clc;

%video plotter#
 load('solkeep.mat')
 
 
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
   view([-100 36])
   mov(k) = getframe(gcf);
end

% Create AVI file.
movie2avi(mov, 'myTemp.avi', 'compression', 'None', 'fps', 10);
 
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