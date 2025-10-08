% Square lattice with central unit zero conductivity
clear; clc;

% ---------------- Parameters (match your style) ----------------
lp0 = 0.005;                 % half edge length (m)
h   = 2*lp0;                 % square grid spacing in x and y
dx  = h; dy = h;

L = 2.0;                     % domain length (x)
W = 2.0;                     % domain width (y)

nx = round(L/dx) + 1;
ny = round(W/dy) + 1;
sigma_target = 6e4; % sheet conductivity in J = -sigma * grad(phi) (S/m)
sigma_x_out = 6e4;   % conductivity used by FEM in x outside the box
sigma_y_out = 6e4;   % conductivity used by FEM in y outside the box

use_inner_box = true;         % set to false if FEM has no special patch
xL=0.9925; xR=1.0075; yB=0.9925; yT=1.0075;
sigma_x_in = 2e4;             % FEM x-conductivity inside the box
sigma_y_in = 2e4;             % FEM y-conductivity inside the box

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

E  = size(edges,1);
Ge = zeros(E,1);

% ---------------- Assign conductivities, central unit zero ----------------
for e = 1:E
    i  = edges(e,1); 
    j  = edges(e,2);
    p0 = coords(i,:); 
    p1 = coords(j,:);
    pm = 0.5*(p0+p1);      % edge midpoint for region lookup
    dx_e = p1(1)-p0(1);
    dy_e = p1(2)-p0(2);

    % Default outside conductivity
    sx = sigma_x_out; sy = sigma_y_out;

    % Inner central box: assign very small conductivities
    if use_inner_box && pm(1)>=xL && pm(1)<=xR && pm(2)>=yB && pm(2)<=yT
        sx = sigma_x_in;
        sy = sigma_y_in;
    end

    if abs(pm(1)-1.0) < 0.02 && abs(pm(2)-1.0) < 0.02
        fprintf('Edge %5d: midpoint=(%.4f, %.4f), sx=%.2e, sy=%.2e\n', ...
            e, pm(1), pm(2), sx, sy);
    end

    % Orientation: horizontal → x, vertical → y
    if abs(dy_e) < 1e-15
        Ge(e) = sx;
    elseif abs(dx_e) < 1e-15
        Ge(e) = sy;
    else
        error('Non axis-aligned edge encountered.');
    end
end
% ---------------- Assemble Laplacian (Kirchhoff) with anisotropic Ge -------
M = sparse(N,N);
b = zeros(N,1);

for e = 1:size(edges,1)
  i = edges(e,1); j = edges(e,2);
  g = Ge(e);
  if g == 0, continue; end
  M(i,i) = M(i,i) + g; 
  M(j,j) = M(j,j) + g;
  M(i,j) = M(i,j) - g; 
  M(j,i) = M(j,i) - g;
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
DT = delaunayTriangulation(coords(:, 1), coords(:, 2));  
tri = DT.ConnectivityList;
P   = DT.Points; 

% Extract triangle points for remaining triangles
t1 = tri(:,1); t2 = tri(:,2); t3 = tri(:,3);
x1 = P(t1,1);  y1 = P(t1,2);
x2 = P(t2,1);  y2 = P(t2,2);
x3 = P(t3,1);  y3 = P(t3,2);

V1 = V(t1); V2 = V(t2); V3 = V(t3);

twiceA = (x2 - x1).*(y3 - y1) - (x3 - x1).*(y2 - y1);
twiceA(abs(twiceA) < 1e-15) = NaN;

dphidx = ( V1.*(y2 - y3) + V2.*(y3 - y1) + V3.*(y1 - y2) ) ./ twiceA;
dphidy = ( V1.*(x3 - x2) + V2.*(x1 - x3) + V3.*(x2 - x1) ) ./ twiceA;

valid = isfinite(dphidx) & isfinite(dphidy);
t1v = t1(valid); t2v = t2(valid); t3v = t3(valid);

xc = (x1(valid) + x2(valid) + x3(valid))/3;
yc = (y1(valid) + y2(valid) + y3(valid))/3;

if exist('sigma_x_out','var') && exist('sigma_y_out','var')
    sigma_x_tri = sigma_x_out*ones(sum(valid),1);
    sigma_y_tri = sigma_y_out*ones(sum(valid),1);
    if exist('use_inner_box','var') && use_inner_box
        inside = (xc>=xL & xc<=xR & yc>=yB & yc<=yT);
        sigma_x_tri(inside) = sigma_x_in;
        sigma_y_tri(inside) = sigma_y_in;
    end
    Jx_tri = -sigma_x_tri .* dphidx(valid);
    Jy_tri = -sigma_y_tri .* dphidy(valid);
else
    Jx_tri = -sigma_target .* dphidx(valid);
    Jy_tri = -sigma_target .* dphidy(valid);
end

%{
Jx_node = accumarray([t1;t2;t3], [Jx_tri;Jx_tri;Jx_tri], [N,1], @mean);
Jy_node = accumarray([t1;t2;t3], [Jy_tri;Jy_tri;Jy_tri], [N,1], @mean);

FJx = scatteredInterpolant(coords(:,1), coords(:,2), Jx_node, 'linear','none');
FJy = scatteredInterpolant(coords(:,1), coords(:,2), Jy_node, 'linear','none');
%}
FJx = scatteredInterpolant(xc, yc, Jx_tri, 'natural','none');
FJy = scatteredInterpolant(xc, yc, Jy_tri, 'natural','none');

Jx = FJx(X,Y);
Jy = FJy(X,Y);

% ---------------- Plots ----------------
figure('Color','w');
pcolor(X, Y, Jx); shading interp; axis equal tight; box on
colormap(cool); cb=colorbar; cb.Label.String='J_x (A/m)';
xlabel('x (m)'); ylabel('y (m)');
title('Current density J_x (A/m)');

y_mid = W/2;
y_quarter = W/4;

JxI = griddedInterpolant({yg,xg}, Jx, 'linear','none');
Jx_mid     = JxI( y_mid    * ones(size(xg)), xg );
Jx_quarter = JxI( y_quarter* ones(size(xg)), xg );

data_case1 = readmatrix('current_density_Jx_center_fig4.csv');
data_case2 = readmatrix('current_density_Jx_quarter_fig4.csv');

x1 = data_case1(:,1); Jx1 = data_case1(:,3);
x2 = data_case2(:,1); Jx2 = data_case2(:,3);

figure('Color','w');
plot(xg, Jx_mid, 'b-', 'LineWidth',2, 'DisplayName', sprintf('y = %.2f (midline)', y_mid));
hold on;
plot(xg, Jx_quarter, 'y--', 'LineWidth',2, 'DisplayName', sprintf('y = %.2f (quarterline)', y_quarter));

plot(x1, Jx1, 'r-', 'LineWidth', 1.5, 'DisplayName','FEM Centerline');
plot(x2, Jx2, 'g--', 'LineWidth', 1.5, 'DisplayName','FEM Quarterline');

grid on; box on; axis tight;
xlabel('x (m)'); ylabel('J_x (A/m)');
title('Overlay: Square Mesh vs FEM current density J_x');
legend('Location','best');
%{
grid on; box on; axis tight;
xlabel('x (m)'); ylabel('J_x (A/m)');
title('Current density J_x along midline and quarterline');
legend('Location','best'); hold off;
%}
figure('Color','w');
pcolor(X, Y, Z); shading interp; hold on;
contour(X, Y, Z, 20, 'k', 'LineWidth', 1);
axis equal tight; box on;
colormap(cool); cb=colorbar; cb.Label.String='V (Volts)';
xlabel('x (m)'); ylabel('y (m)');
title('Potential V with contours');
