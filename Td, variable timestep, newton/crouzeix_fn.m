% crouzeix.m in 2D function

% output format
% format longEng
% format long

function [Pfsol, rRTx, rRTy, rS, rSsol, rSrhs, rUsol, rVsol,...
    rXbcs, rYbcs, rXmps, rYmps, s2e, triA] ...
    = crouzeix_fn(Ra, isRT0mode, Tsol, vfilename)

%'voronoi_sqrgrid_data_5pts.mat' 'voronoi_data_5pt' 
load(vfilename);

% machine eps
cmpEps = 1.0e-12;

nx = points(:,1);
ny = points(:,2);
nodenum = length(nx);

% mesh generation
trirep = delaunayTriangulation(nx, ny);
e = trirep.edges();
edgenum = length(e);
trinum = trirep.size(1);

simplex2edges = zeros(trinum, 3);
hi = zeros(trinum, 1);

rhs = zeros(2 * edgenum + trinum, 1);

% defining the boundary
isBoundaryEdge = zeros(edgenum, 1);
for i=1:edgenum
    xa = nx(e(i, 1));
    ya = ny(e(i, 1));

    xb = nx(e(i, 2));
    yb = ny(e(i, 2));

    if ((abs(xa) < cmpEps) && (abs(xb) < cmpEps)) ...
            || ((abs(ya) < cmpEps) && (abs(yb) < cmpEps)) ...
            || ((abs(xa - 1.0) < cmpEps) && (abs(xb - 1.0) < cmpEps)) ...
            | ((abs(ya - 1.0) < cmpEps) && (abs(yb - 1.0) < cmpEps))
        isBoundaryEdge(i) = 1;
    end
end

% Read in quadrature formulas, weights
[qmat, ws] = quad_formulas();
pdelibqmat = zeros(54, 2);
for i=1:54
    pdelibqmat(i, 1) = -qmat(i, 1) + qmat(i, 2) - qmat(i, 3);
    pdelibqmat(i, 2) = -qmat(i, 1) - qmat(i, 2) + qmat(i, 3);
end

% construct adjacency matrix returning for each triangle
% all adjacent edges
for i=1:edgenum
    tris = edgeAttachments(trirep, e(i, 1), e(i, 2));
    trisvec=cell2mat(tris);
    for j=1:length(trisvec)
        k=1;
        t = trisvec(j);
        while(k > 0 && k <= hi(t))
            if simplex2edges(t, k+1) == i
              k = -1;
            else
              k = k+1;
            end
        end
        if k > 0
            hi(t) = hi(t) + 1;
            simplex2edges(t, hi(t)) = i;
        end
    end
end

simplex2nodes = zeros(trinum, 3);
hi = zeros(trinum, 1);

% construct adjacency matrix returning for each triangle
% all adjacent nodes
for i=1:nodenum
    tris = vertexAttachments(trirep, i);
    trisvec=cell2mat(tris);
    for j=1:length(trisvec)
        k=1;
        t = trisvec(j);
        while(k > 0 && k <= hi(t))
            if simplex2nodes(t, k+1) == i
              k = -1;
            else
              k = k+1;
            end
        end
        if k > 0
            hi(t) = hi(t) + 1;
            simplex2nodes(t, hi(t)) = i;
        end
    end
end

% compute triangle areas by Heron's formula
triareas = zeros(trinum, 1);
for i=1:trinum
    n1 = simplex2nodes(i, 1);
    n2 = simplex2nodes(i, 2);
    n3 = simplex2nodes(i, 3);

    l12 = sqrt((nx(n1) - nx(n2)).^2 + (ny(n1) - ny(n2)).^2);
    l13 = sqrt((nx(n1) - nx(n3)).^2 + (ny(n1) - ny(n3)).^2);
    l23 = sqrt((nx(n2) - nx(n3)).^2 + (ny(n2) - ny(n3)).^2);

    s = (l12 + l13 + l23) / 2;

    triareas(i) = sqrt(s * (s - l12) * (s - l13) * (s - l23));
end

xbcs = zeros(trinum, 1);
ybcs = zeros(trinum, 1);

% for RT0 reconstruction and pressure visualization
for i=1:trinum
    n1 = simplex2nodes(i, 1);
    n2 = simplex2nodes(i, 2);
    n3 = simplex2nodes(i, 3);

    xbcs(i) = 1/3 * (nx(n1) + nx(n2) + nx(n3));
    ybcs(i) = 1/3 * (ny(n1) + ny(n2) + ny(n3));
end

%simplex2edges
%simplex2nodes
%triareas

%triplot(simplex2nodes, nx, ny);

% allocating the stiffness matrix
is = [];
js = [];
ss = [];
mat = sparse(is, js, ss, 2 * edgenum + trinum, 2 * edgenum + trinum);

phi_as = zeros(3, 1);
phi_bs = zeros(3, 1);
phi_cs = zeros(3, 1);

phi_mat = zeros(3);
phi_rhs = zeros(3, 1);
phi_sol = zeros(3, 1);


for i=1:trinum
    e1 = simplex2edges(i, 1);
    e2 = simplex2edges(i, 2);
    e3 = simplex2edges(i, 3);
    
    emp1x = (nx(e(e1, 1)) + nx(e(e1, 2))) / 2;
    emp1y = (ny(e(e1, 1)) + ny(e(e1, 2))) / 2;
    
    emp2x = (nx(e(e2, 1)) + nx(e(e2, 2))) / 2;
    emp2y = (ny(e(e2, 1)) + ny(e(e2, 2))) / 2;
    
    emp3x = (nx(e(e3, 1)) + nx(e(e3, 2))) / 2;
    emp3y = (ny(e(e3, 1)) + ny(e(e3, 2))) / 2;
    
    phi_mat(1, 1) = 1;
    phi_mat(2, 1) = 1;
    phi_mat(3, 1) = 1;
    phi_mat(1, 2) = emp1x;
    phi_mat(2, 2) = emp2x;
    phi_mat(3, 2) = emp3x;
    
    phi_mat(1, 3) = emp1y;
    phi_mat(2, 3) = emp2y;
    phi_mat(3, 3) = emp3y;
    
    phi_sol = phi_mat \ [1 0 0]';
    phi_as(1) = phi_sol(1);
    phi_bs(1) = phi_sol(2);
    phi_cs(1) = phi_sol(3);
    
    phi_sol = phi_mat \ [0 1 0]';
    phi_as(2) = phi_sol(1);
    phi_bs(2) = phi_sol(2);
    phi_cs(2) = phi_sol(3);
    
    phi_sol = phi_mat \ [0 0 1]';
    phi_as(3) = phi_sol(1);
    phi_bs(3) = phi_sol(2);
    phi_cs(3) = phi_sol(3);
    
    %phi_as
    %phi_bs
    %phi_cs
    
    le = [e1 e2 e3]';
    
    % triangle corner data, for barycentric coords
    corner_nodes = simplex2nodes(i,:);
    corner_x = nx(corner_nodes);
    corner_y = ny(corner_nodes);
    
    % assembly of the Laplacian operators
    for ii=1:3
        for jj=1:3
            if(~isBoundaryEdge(le(ii)) && ~isBoundaryEdge(le(jj)))
              mat(le(jj), le(ii)) = mat(le(jj), le(ii)) + triareas(i) * (phi_bs(ii) * phi_bs(jj) + phi_cs(ii) * phi_cs(jj));
              mat(edgenum+le(jj), edgenum+le(ii)) = mat(edgenum+le(jj), edgenum+le(ii)) + triareas(i) * (phi_bs(ii) * phi_bs(jj) + phi_cs(ii) * phi_cs(jj));
            end
        end
    end
    
    % assembly of the divergence and the pressure gradient operators
    for ii=1:3
        % (-u_x, q)
        if (~isBoundaryEdge(le(ii))) && (i ~= 1)
          mat(2 * edgenum + i, le(ii)) = mat(2 * edgenum + i, le(ii)) + triareas(i) * (-phi_bs(ii));
          mat(le(ii), 2 * edgenum + i) = mat(le(ii), 2 * edgenum + i) + triareas(i) * (-phi_bs(ii));
        end
        
        % (-v_y, q)
        if (~isBoundaryEdge(le(ii))) && (i ~= 1)
          mat(2 * edgenum + i, edgenum + le(ii)) = mat(2 * edgenum + i, edgenum + le(ii)) + triareas(i) * (-phi_cs(ii));
          mat(edgenum + le(ii), 2 * edgenum + i) = mat(edgenum + le(ii), 2 * edgenum + i) + triareas(i) * (-phi_cs(ii));
        end
    end
    
    % assembly of the right hand side
    [qN, dummy] = size(qmat);
    
    if ~isRT0mode  
      for qi=1:qN
          xc = barycentricToCartesian(trirep, [i]', [qmat(qi, 1) qmat(qi, 2) qmat(qi, 3)]);
          xl = xc(1);
          yl = xc(2);
          
          % linear interpolation over the triangle to find temperature at
          % present quadrature point.
          Tl = temp_interp(Tsol, corner_nodes, corner_x, corner_y, xl, yl);
          
          for ii=1:3
            if ~isBoundaryEdge(le(ii))
                
              % note: 0 used to be f1(xl,yl), now zero do to bouyancy,
              % Ra * T_l used to be f2(xl,yl)
              rhs(le(ii)) = rhs(le(ii)) + triareas(i) * ws(qi) ...
                  * 0 * (phi_as(ii) + xl * phi_bs(ii) + yl * phi_cs(ii));
              rhs(edgenum + le(ii)) = rhs(edgenum + le(ii)) + triareas(i) * ws(qi) ...
                  * Ra * Tl * (phi_as(ii) + xl * phi_bs(ii) + yl * phi_cs(ii));
              
                    % For velocity reconstruction
                    xa = nx(e(le(ii), 1));
                    ya = ny(e(le(ii), 1));
                    
                    xb = nx(e(le(ii), 2));
                    yb = ny(e(le(ii), 2));
                    
                    % get edge midpoint
                    empx = 1/2 * (xa + xb);
                    empy = 1/2 * (ya + yb);
                    
                    % compute triangle edge length sigma
                    sigma = sqrt((xb-xa)^2 + (yb-ya)^2);
                    
                    % compute unit normal to triangle edge sigma:
                    nsigmax = -(yb - ya) / sigma;
                    nsigmay = (xb - xa) / sigma;
                    
                    % get right orientation of edge normal
                    if (empx - xbcs(i)) * nsigmax ...
                            + (empy - ybcs(i)) * nsigmay < 0
                        nsigmax = -nsigmax;
                        nsigmay = -nsigmay;
                    end
                                    
                    % aTx, aTy, b parameter Raviart
                    aTx = 1/triareas(i) * sigma * (empx - xbcs(i));
                    aTy = 1/triareas(i) * sigma * (empy - ybcs(i));
                    b = sigma/triareas(i);
                    
                    % Save RT data in cell for later use
                    % Row contains [aTx aTy b unormx unormy]
                    RT{i,1}{ii} = [aTx, aTy, b, nsigmax, nsigmay, le(ii)];    
            else
                    % not boundary edge case
                    RT{i,1}{ii} = [0, 0, 0, 0, 0, le(ii)];
            end
          end
      end
    end % if ~isRT0mode
    
    if isRT0mode
        % RT0 reconstruction of the form v = a_T + c_T ( x - x_T)
        for qi=1:qN
            xc = barycentricToCartesian(trirep, [i]', [qmat(qi, 1) qmat(qi, 2) qmat(qi, 3)]);
            xl = xc(1);
            yl = xc(2);
            
          % linear interpolation over the triangle to find temperature at
          % present quadrature point.
          Tl = temp_interp(Tsol, corner_nodes, corner_x, corner_y, xl, yl);
            
            for ii=1:3
                if ~isBoundaryEdge(le(ii))
                    xa = nx(e(le(ii), 1));
                    ya = ny(e(le(ii), 1));
                    
                    xb = nx(e(le(ii), 2));
                    yb = ny(e(le(ii), 2));
                    
                    % get edge midpoint
                    empx = 1/2 * (xa + xb);
                    empy = 1/2 * (ya + yb);
                    
                    % compute triangle edge length sigma
                    sigma = sqrt((xb-xa)^2 + (yb-ya)^2);
                    
                    % b parameter Raviart
                    b = sigma/triareas(i);
                    
                    % compute unit normal to triangle edge sigma:
                    nsigmax = -(yb - ya) / sigma;
                    nsigmay = (xb - xa) / sigma;
                    
                    % get right orientation of edge normal
                    if (empx - xbcs(i)) * nsigmax ...
                            + (empy - ybcs(i)) * nsigmay < 0
                        nsigmax = -nsigmax;
                        nsigmay = -nsigmay;
                    end
                    
                    % nsigmax
                    % nsigmay
                    
                    % horizontal degree of freedom
                    aTx = 1/triareas(i) * sigma * nsigmax * (empx - xbcs(i));
                    aTy = 1/triareas(i) * sigma * nsigmax * (empy - ybcs(i));
                    
                    % cT is 1/2 of the divergence of hat function (0, 1)
                    cT = phi_bs(ii) / 2;
                    
                    % note: 0 used to be f1(xl,yl), now zero do to bouyancy,
                    % Ra * T_l used to be f2(xl,yl)
                    rhs(le(ii)) = rhs(le(ii)) + triareas(i) * ws(qi) ...
                        * (0 * (aTx + cT * (xl - xbcs(i))) ...
                    + Ra * Tl * (aTy + cT * (yl - ybcs(i))));
                
                    % vertical degree of freedom
                    aTx = 1/triareas(i) * sigma * nsigmay * (empx - xbcs(i));
                    aTy = 1/triareas(i) * sigma * nsigmay * (empy - ybcs(i));
                    
                    % cT is 1/2 of the divergence of hat function (0, 1)
                    cT = phi_cs(ii) / 2;
                    rhs(edgenum + le(ii)) = rhs(edgenum + le(ii)) + triareas(i) * ws(qi) ...
                        * (0 * (aTx + cT * (xl - xbcs(i))) ...
                    + Ra * Tl * (aTy + cT * (yl - ybcs(i))));
                
                    % Construct Raviart Thomas Array for each edge
                    % Row contains [aTx aTy b unormx unormy]
                    
                    aTx = 1/triareas(i) * sigma * (empx - xbcs(i));
                    aTy = 1/triareas(i) * sigma * (empy - ybcs(i));
                    
                    RT{i,1}{ii} = [aTx, aTy, b, nsigmax, nsigmay, le(ii)];
                else
                    RT{i,1}{ii} = [0, 0, 0, 0, 0, le(ii)];
                end
            end
        end
    end % if isRT0mode
end

% assembly of the boundary behavior
for i=1:edgenum
    if(isBoundaryEdge(i))
        mat(i, i) = 1;
        mat(edgenum + i, edgenum + i) = 1;
    end
end
mat(2 * edgenum + 1, 2 * edgenum + 1) = -1;

% linear solve
sol = mat \ rhs;

% flow visualization
xmps = zeros(edgenum, 1);
ymps = zeros(edgenum, 1);
usol = zeros(edgenum, 1);
vsol = zeros(edgenum, 1);
psol = zeros(trinum,  1);

for i=1:edgenum
    xmps(i) = 1/2 * (nx(e(i, 1)) + nx(e(i,2)));
    ymps(i) = 1/2 * (ny(e(i, 1)) + ny(e(i,2)));
    usol(i) = sol(i);
    vsol(i) = sol(edgenum + i);
end

% Reconstruction of velocity computation
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
%
% In this case, reconstruction is handled by an external function.
[RV_Ax, RV_Ay] = velocity_reconstruction(RT, usol, vsol, length(triareas));

% change pressure level
pavg = 0;
mesomega = 0;
for i=1:trinum
    mesomega = mesomega + triareas(i);
    pavg = pavg + triareas(i) * sol(2 * edgenum + i);
end
pavg = pavg / mesomega;

for i=1:trinum
    psol(i) = sol(2 * edgenum + i) - pavg;
end

Pfsol = psol;


% Return reconstructed velocity elements, or normal velocity components
rRTx = RV_Ax;
rRTy = RV_Ay;

% Return stokes matrix and solution vector, rhs of equation
rS = mat;
rSsol = sol;
rSrhs = rhs;

% return velocity field solution
rUsol = usol;
rVsol = vsol;

% return barycenters
rXbcs = xbcs;
rYbcs = ybcs;

% return midpoints
rXmps = xmps;
rYmps = ymps;

s2e = simplex2edges;
triA = triareas;

% norm stuff (still needed?)
%{
l2normExact = sqrt(l2normExact)
l2norm = sqrt(l2norm)
l2errnorm = sqrt(l2errnorm)

h1normExact = sqrt(h1normExact)
h1norm = sqrt(h1norm)
h1errnorm = sqrt(h1errnorm)

divnorm = sqrt(divnorm)

pnormExact = sqrt(pnormExact)
pnorm = sqrt(pnorm)
perrnorm = sqrt(perrnorm)
%}

end







