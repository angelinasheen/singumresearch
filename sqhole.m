
% Parameters 
lp0 = 0.005;               % half edge length (m)
dx  = 2*lp0;               % horizontal spacing
dy  = sqrt(3)*lp0;         % vertical spacing

nx = round(2/dx) + 1;
ny = round(2/dy) + 1;
L = 2.0;
W  = 2.0;

% Physical sheet conductivity 
sigma_target = 6e4;        

% Generate nodes and connectivity 
nodes    = zeros(nx*ny,3);
node_map = zeros(ny,nx);
id = 1;
for j = 1:ny
  for i = 1:nx
    x = (i-1)*dx + mod(j,2)*lp0;
    y = (j-1)*dy;
    nodes(id,:)    = [id, x, y];
    node_map(j,i)  = id;
    id = id + 1;
  end
end
coords = nodes(:,2:3);

edges = [];
for j = 1:ny
  for i = 1:nx
    n = node_map(j,i);
    if i < nx
      edges(end+1,:) = [n, node_map(j,i+1)];
    end
    if j < ny
      if mod(j,2)==1
        if i<nx, edges(end+1,:) = [n, node_map(j+1,i+1)]; end
        edges(end+1,:) = [n, node_map(j+1,i)];
      else
        if i>1, edges(end+1,:) = [n, node_map(j+1,i-1)]; end
        edges(end+1,:) = [n, node_map(j+1,i)];
      end
    end
  end
end

N = size(coords,1);
M = sparse(N,N);
b = zeros(N,1);
K = 1;                     

% Square hole 0.4×0.4 m 
hole_w = 0.4; 
hole_h = 0.4;
cx = L/2; 
cy = W/2;
hole_nodes = find( coords(:,1)>=cx-hole_w/2 & coords(:,1)<=cx+hole_w/2 & ...
                   coords(:,2)>=cy-hole_h/2 & coords(:,2)<=cy+hole_h/2 );
valid = true(N,1);
valid(hole_nodes) = false;

% Assemble M only for valid nodes 
for e = 1:size(edges,1)
  i = edges(e,1); j = edges(e,2);
  if valid(i) && valid(j)
    M(i,i) = M(i,i) + K; 
    M(j,j) = M(j,j) + K;
    M(i,j) = M(i,j) - K; 
    M(j,i) = M(j,i) - K;
 

  end
end

% Dirichlet BC on left/right boundaries 
left_nodes  = find(valid & coords(:,1)<1e-6);
right_nodes = find(valid & abs(coords(:,1)-L)<1e-6);
bc_nodes = [left_nodes; right_nodes];
bc_vals  = [-5*ones(size(left_nodes)); 5*ones(size(right_nodes))];
for k = 1:numel(bc_nodes)
  n = bc_nodes(k);
  M(n,:) = 0; M(n,n) = 1;
  b(n)   = bc_vals(k);
end

% Solve potentials 
V = M\b;

% Interpolate phi onto grid 
xg = linspace(0,L,400);
yg = linspace(0,W,400);
[X,Y] = meshgrid(xg,yg);
Z = griddata(coords(valid,1), coords(valid,2), V(valid), X, Y, 'natural');


% Mask hole region
hole_mask = (X>=cx-hole_w/2 & X<=cx+hole_w/2 & Y>=cy-hole_h/2 & Y<=cy+hole_h/2);
Z(hole_mask) = NaN;

% ---------------- Current density J = -K * grad(V) ----------------
% Build a triangulation on the valid nodes
coords_v = coords(valid,:);          % Nx2
V_v      = V(valid);                 % Nx1
DT = delaunayTriangulation(coords_v(:,1), coords_v(:,2));
tri = DT.ConnectivityList;           % Tx3, triangles by vertex ids (into coords_v)
% Extract triangle vertices (vectors of length T = #triangles)
t1 = tri(:,1); t2 = tri(:,2); t3 = tri(:,3);
x1 = coords_v(t1,1);  y1 = coords_v(t1,2);
x2 = coords_v(t2,1);  y2 = coords_v(t2,2);
x3 = coords_v(t3,1);  y3 = coords_v(t3,2);

% Potentials at triangle vertices
V1 = V_v(t1);  V2 = V_v(t2);  V3 = V_v(t3);

% Gradient of phi is constant inside each linear triangle
% twiceA = signed 2*area; safe for mixed orientation
twiceA = (x2 - x1).*(y3 - y1) - (x3 - x1).*(y2 - y1);
twiceA(abs(twiceA) < 1e-15) = NaN;  % guard degenerate tris

dphidx = ( V1.*(y2 - y3) + V2.*(y3 - y1) + V3.*(y1 - y2) ) ./ twiceA;
dphidy = ( V1.*(x3 - x2) + V2.*(x1 - x3) + V3.*(x2 - x1) ) ./ twiceA;

% Current density per triangle: J = -sigma * grad(phi)
Jx_tri = -sigma_target .* dphidx;
Jy_tri = -sigma_target .* dphidy;

% Triangle centroids (for interpolation to a continuous field)
xc = (x1 + x2 + x3)/3;
yc = (y1 + y2 + y3)/3;


% Continuous field via scatteredInterpolant (no ringing)
FJx = scatteredInterpolant(xc, yc, Jx_tri, 'natural','none');
FJy = scatteredInterpolant(xc, yc, Jy_tri, 'natural','none');

% Sample on your plot grid / lines
Jx = FJx(X, Y);
Jy = FJy(X, Y);

% Mask the hole for plotting
Jx(hole_mask) = NaN;  Jy(hole_mask) = NaN;



% ---------------- Jx along bottom quarter line (y = W/4) ----------------
y_quarter = W/4;
JxI = griddedInterpolant({yg,xg}, Jx, 'linear','none');  % array order = (y,x)
Jx_quarter = JxI( y_quarter * ones(size(xg)), xg);
y_mid = W/2;
Jx_mid = JxI(y_mid*ones(size(xg)), xg);
inHole_mid = (xg>=cx-hole_w/2 & xg<=cx+hole_w/2) & ...
             (abs(y_mid-cy) <= hole_h/2);
Jx_mid(inHole_mid) = NaN;


figure('Color','w');
plot(xg, Jx_mid, 'LineWidth', 2, 'DisplayName', 'y = 1 (centerline)');
hold on;
plot(xg, Jx_quarter, 'LineWidth', 2, 'DisplayName', 'y = 0.5 (quarterline)');
grid on; box on; axis tight;
xlabel('x (m)'); ylabel('J_x (A/m)');
legend('Location','best')

% ---------------- Heatmap of Jx ----------------
figure('Color','w');
pcolor(X, Y, Jx); shading interp; axis equal tight; box on
colormap(cool); cb=colorbar; cb.Label.String='A/m';
xlabel('x (m)'); ylabel('y (m)');
xticks(0:0.5:L); yticks(0:0.5:W);
title('Current density J_x (A/m)');

% ---------------- Heatmap of V ----------------
figure('Color','w');
pcolor(X, Y, Z); 
shading interp; 
hold on;
contour(X, Y, Z, 20, 'k', 'LineWidth', 1);  % contour lines
axis equal tight; 
box on
colormap(cool); cb=colorbar; cb.Label.String='A/m';
xlabel('x (m)'); ylabel('y (m)');
xticks(0:0.5:L); yticks(0:0.5:W);
title('Current density J_x (A/m)');




