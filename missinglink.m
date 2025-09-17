% Square lattice with one missing central link (same method as your script)
clear; clc;

% ---------------- Parameters (match your style) ----------------
lp0 = 0.005;                 % half edge length (m)
h   = 2*lp0;                 % square grid spacing in x and y
dx  = h; dy = h;

L = 2.0;                     % domain length (x)
W = 2.0;                     % domain width (y)

nx = round(L/dx) + 1;
ny = round(W/dy) + 1;

sigma_target = 6e4;          % sheet conductivity in J = -sigma * grad(phi)  (S/m)

% ---------------- Generate nodes & connectivity (square grid) ----------------
nodes    = zeros(nx*ny,3);   % [id, x, y]
node_map = zeros(ny,nx);
id = 1;
for j = 1:ny
  for i = 1:nx
    x = (i-1)*dx;
    y = (j-1)*dy;
    nodes(id,:)   = [id, x, y];
    node_map(j,i) = id;
    id = id + 1;
  end
end
coords = nodes(:,2:3);
N = size(coords,1);

% 4-neighbor edges (+x, +y only; undirected handled in assembly)
edges = [];
for j = 1:ny
  for i = 1:nx
    n = node_map(j,i);
    if i < nx
      edges(end+1,:) = [n, node_map(j,i+1)]; %#ok<AGROW>
    end
    if j < ny
      edges(end+1,:) = [n, node_map(j+1,i)]; %#ok<AGROW>
    end
  end
end

% ---------------- Remove a single central link ----------------
cx = L/2; cy = W/2;
midpts = 0.5*(coords(edges(:,1),:) + coords(edges(:,2),:));
[~, kmin] = min(vecnorm(midpts - [cx,cy], 2, 2));
central_edge = edges(kmin,:);           % store for plotting
edges(kmin,:) = [];                     % remove that link

% ---------------- Assemble Laplacian (Kirchhoff) ----------------
M = sparse(N,N);
b = zeros(N,1);
K = 1;                                  % unit edge weight (we scale J by sigma later)

for e = 1:size(edges,1)
  i = edges(e,1); j = edges(e,2);
  M(i,i) = M(i,i) + K; 
  M(j,j) = M(j,j) + K;
  M(i,j) = M(i,j) - K; 
  M(j,i) = M(j,i) - K;
end

% ---------------- Dirichlet BC: left -5, right +5 ----------------
left_nodes  = find( abs(coords(:,1)-0) < 1e-12 );
right_nodes = find( abs(coords(:,1)-L) < 1e-12 );

bc_nodes = [left_nodes; right_nodes];
bc_vals  = [-5*ones(size(left_nodes)); 5*ones(size(right_nodes))];

for k = 1:numel(bc_nodes)
  n = bc_nodes(k);
  M(n,:) = 0; M(n,n) = 1;
  b(n)   = bc_vals(k);
end

% ---------------- Solve potentials ----------------
V = M \ b;

% ---------------- Interpolate phi onto grid (for plotting V) ----------------
xg = linspace(0,L,400);
yg = linspace(0,W,400);
[X,Y] = meshgrid(xg,yg);
Z = griddata(coords(:,1), coords(:,2), V, X, Y, 'natural');

% ---------------- Current density J = -sigma * grad(phi) ----------------
% Triangulation on all nodes (same idea as your script)

DT = delaunayTriangulation(coords(:, 1), coords(:, 2));  % constrained by edges
tri = DT.ConnectivityList;
P   = DT.Points; 

% ---------------- Remove triangles containing the missing central edge ----------------
mask_bad_tri = any(tri == central_edge(1) | tri == central_edge(2), 2);
tri_clean = tri(~mask_bad_tri, :);

% Extract triangle points for remaining triangles
t1 = tri_clean(:,1); t2 = tri_clean(:,2); t3 = tri_clean(:,3);
x1 = P(t1,1);  y1 = P(t1,2);
x2 = P(t2,1);  y2 = P(t2,2);
x3 = P(t3,1);  y3 = P(t3,2);

V1 = V(t1); V2 = V(t2); V3 = V(t3);

twiceA = (x2 - x1).*(y3 - y1) - (x3 - x1).*(y2 - y1);
twiceA(abs(twiceA) < 1e-15) = NaN;   % guard degenerate tris

% Gradients of phi in each triangle (piecewise linear)
dphidx = ( V1.*(y2 - y3) + V2.*(y3 - y1) + V3.*(y1 - y2) ) ./ twiceA;
dphidy = ( V1.*(x3 - x2) + V2.*(x1 - x3) + V3.*(x2 - x1) ) ./ twiceA;

% Current density per triangle
Jx_tri = -sigma_target .* dphidx;
Jy_tri = -sigma_target .* dphidy;

% Triangle centroids
xc = (x1 + x2 + x3)/3;
yc = (y1 + y2 + y3)/3;

% Flatten triangles to nodes

Jx_node = accumarray([t1;t2;t3], [Jx_tri;Jx_tri;Jx_tri], [N,1], @mean);
Jy_node = accumarray([t1;t2;t3], [Jy_tri;Jy_tri;Jy_tri], [N,1], @mean);

% Interpolate node-averaged J onto grid using linear

FJx = scatteredInterpolant(coords(:,1), coords(:,2), Jx_node, 'linear','none');
FJy = scatteredInterpolant(coords(:,1), coords(:,2), Jy_node, 'linear','none');

%{
FJx = scatteredInterpolant(xc, yc, Jx_tri, 'nearest','none');
FJy = scatteredInterpolant(xc, yc, Jy_tri, 'nearest','none');
%}

Jx = FJx(X,Y);
Jy = FJy(X,Y);




% ---------------- Plots ----------------

% 1) Heatmap of Jx (show missing link)
figure('Color','w');
pcolor(X, Y, Jx); shading interp; axis equal tight; box on
colormap(cool); cb=colorbar; cb.Label.String='J_x (A/m)';
xlabel('x (m)'); ylabel('y (m)');
title('Current density J_x (A/m)');
hold on;
cu = coords(central_edge(1),:); 
cv = coords(central_edge(2),:);
plot([cu(1) cv(1)], [cu(2) cv(2)], 'k-', 'LineWidth', 3);

% 2) Jx along midline (y = W/2) and quarterline (y = W/4)
y_mid = W/2;
y_quarter = W/4;

JxI = griddedInterpolant({yg,xg}, Jx, 'linear','none'); % array order = (y,x)
Jx_mid     = JxI( y_mid    * ones(size(xg)), xg );
Jx_quarter = JxI( y_quarter* ones(size(xg)), xg );

figure('Color','w');
plot(xg, Jx_mid, 'LineWidth',2, 'DisplayName', sprintf('y = %.2f (midline)', y_mid));
hold on;
plot(xg, Jx_quarter, 'LineWidth',2, 'DisplayName', sprintf('y = %.2f (quarterline)', y_quarter));
grid on; box on; axis tight;
xlabel('x (m)'); ylabel('J_x (A/m)');
title('Current density J_x along midline and quarterline');
legend('Location','best'); hold off;

% (Optional) Heatmap of potential with contours (useful sanity check)
figure('Color','w');
pcolor(X, Y, Z); shading interp; hold on;
contour(X, Y, Z, 20, 'k', 'LineWidth', 1);
axis equal tight; box on;
colormap(cool); cb=colorbar; cb.Label.String='V (Volts)';
xlabel('x (m)'); ylabel('y (m)');
title('Potential V with contours');
